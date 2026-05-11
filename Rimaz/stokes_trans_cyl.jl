# Transient Stokes flow around a cylinder
# Based on the Navier-Stokes tutorial (ns_tut_cyl.jl) with the convective term dropped.
# The governing equation is:
#   M * du/dt = K * u
# where K is the Stokes operator (viscous diffusion + pressure-velocity coupling)
# and M is the velocity mass matrix.

using Ferrite, SparseArrays, BlockArrays, LinearAlgebra, WriteVTK

using DiffEqBase
using OrdinaryDiffEqRosenbrock: Rodas5P

ν = 1.0 / 1000.0 # kinematic viscosity

using FerriteGmsh
using FerriteGmsh: Gmsh
Gmsh.initialize()
gmsh.option.set_number("General.Verbosity", 2)
dim = 2

# Build the channel geometry with a circular obstacle
rect_tag = gmsh.model.occ.add_rectangle(0, 0, 0, 1.1, 0.41)
circle_tag = gmsh.model.occ.add_circle(0.2, 0.2, 0, 0.05)
circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
gmsh.model.occ.cut([(dim, rect_tag)], [(dim, circle_surf_tag)])

gmsh.model.occ.synchronize()

bottomtag = gmsh.model.model.add_physical_group(dim - 1, [6], -1, "bottom")
lefttag   = gmsh.model.model.add_physical_group(dim - 1, [7], -1, "left")
righttag  = gmsh.model.model.add_physical_group(dim - 1, [8], -1, "right")
toptag    = gmsh.model.model.add_physical_group(dim - 1, [9], -1, "top")
holetag   = gmsh.model.model.add_physical_group(dim - 1, [5], -1, "hole")

gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)

gmsh.model.mesh.generate(dim)
grid = togrid()
Gmsh.finalize()

# Taylor-Hood Q2/Q1 elements (same as NS tutorial)
ip_v = Lagrange{RefQuadrilateral, 2}()^dim
qr   = QuadratureRule{RefQuadrilateral}(4)
cellvalues_v = CellValues(qr, ip_v)

ip_p = Lagrange{RefQuadrilateral, 1}()
cellvalues_p = CellValues(qr, ip_p)

dh = DofHandler(grid)
add!(dh, :v, ip_v)
add!(dh, :p, ip_p)
close!(dh)

# -----------------------------------------------------------------------
# Boundary conditions (identical to NS tutorial)
# -----------------------------------------------------------------------
ch = ConstraintHandler(dh)

noslip_facet_names = ["top", "bottom", "hole"]
∂Ω_noslip = union(getfacetset.((grid,), noslip_facet_names)...)
noslip_bc = Dirichlet(:v, ∂Ω_noslip, (x, t) -> Vec((0.0, 0.0)), [1, 2])
add!(ch, noslip_bc)

∂Ω_inflow = getfacetset(grid, "left")

vᵢₙ(t) = min(t * 1.5, 1.5) # ramp inflow velocity up to 1.5

parabolic_inflow_profile(x, t) = Vec((4 * vᵢₙ(t) * x[2] * (0.41 - x[2]) / 0.41^2, 0.0))
inflow_bc = Dirichlet(:v, ∂Ω_inflow, parabolic_inflow_profile, [1, 2])
add!(ch, inflow_bc)

∂Ω_free = getfacetset(grid, "right") # outflow: left free (no BC)

close!(ch)
update!(ch, 0.0)

# -----------------------------------------------------------------------
# Assemble the velocity mass matrix M
# Only the velocity block is non-zero (pressure has no time derivative).
# -----------------------------------------------------------------------
function assemble_mass_matrix(cellvalues_v::CellValues, cellvalues_p::CellValues, M::SparseMatrixCSC, dh::DofHandler)
    n_basefuncs_v = getnbasefunctions(cellvalues_v)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    n_basefuncs   = n_basefuncs_v + n_basefuncs_p
    v▄, p▄ = 1, 2
    Mₑ = BlockedArray(zeros(n_basefuncs, n_basefuncs),
                      [n_basefuncs_v, n_basefuncs_p],
                      [n_basefuncs_v, n_basefuncs_p])

    mass_assembler = start_assemble(M)
    for cell in CellIterator(dh)
        fill!(Mₑ, 0)
        Ferrite.reinit!(cellvalues_v, cell)

        for q_point in 1:getnquadpoints(cellvalues_v)
            dΩ = getdetJdV(cellvalues_v, q_point)
            # Only the velocity-velocity block: M_vv = ∫ φᵢ · φⱼ dΩ
            for i in 1:n_basefuncs_v
                φᵢ = shape_value(cellvalues_v, q_point, i)
                for j in 1:n_basefuncs_v
                    φⱼ = shape_value(cellvalues_v, q_point, j)
                    Mₑ[BlockIndex((v▄, v▄), (i, j))] += φᵢ ⋅ φⱼ * dΩ
                end
            end
        end
        assemble!(mass_assembler, celldofs(cell), Mₑ)
    end
    return M
end

# -----------------------------------------------------------------------
# Assemble the Stokes stiffness matrix K
# K_vv = -ν ∫ ∇φᵢ : ∇φⱼ dΩ   (viscous diffusion, negative → on RHS)
# K_vp = +∫ (∇·φᵢ) ψⱼ dΩ      (pressure gradient)
# K_pv = +∫ ψᵢ (∇·φⱼ) dΩ      (incompressibility constraint)
# -----------------------------------------------------------------------
function assemble_stokes_matrix(cellvalues_v::CellValues, cellvalues_p::CellValues, ν, K::SparseMatrixCSC, dh::DofHandler)
    n_basefuncs_v = getnbasefunctions(cellvalues_v)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    n_basefuncs   = n_basefuncs_v + n_basefuncs_p
    v▄, p▄ = 1, 2
    Kₑ = BlockedArray(zeros(n_basefuncs, n_basefuncs),
                      [n_basefuncs_v, n_basefuncs_p],
                      [n_basefuncs_v, n_basefuncs_p])

    stiffness_assembler = start_assemble(K)
    for cell in CellIterator(dh)
        fill!(Kₑ, 0)
        Ferrite.reinit!(cellvalues_v, cell)
        Ferrite.reinit!(cellvalues_p, cell)

        for q_point in 1:getnquadpoints(cellvalues_v)
            dΩ = getdetJdV(cellvalues_v, q_point)

            # Viscous block: -ν ∫ ∇φᵢ : ∇φⱼ dΩ
            for i in 1:n_basefuncs_v
                ∇φᵢ = shape_gradient(cellvalues_v, q_point, i)
                for j in 1:n_basefuncs_v
                    ∇φⱼ = shape_gradient(cellvalues_v, q_point, j)
                    Kₑ[BlockIndex((v▄, v▄), (i, j))] -= ν * ∇φᵢ ⊡ ∇φⱼ * dΩ
                end
            end

            # Pressure-velocity coupling blocks
            for j in 1:n_basefuncs_p
                ψ = shape_value(cellvalues_p, q_point, j)
                for i in 1:n_basefuncs_v
                    divφ = shape_divergence(cellvalues_v, q_point, i)
                    Kₑ[BlockIndex((v▄, p▄), (i, j))] += (divφ * ψ) * dΩ
                    Kₑ[BlockIndex((p▄, v▄), (j, i))] += (ψ * divφ) * dΩ
                end
            end
        end
        assemble!(stiffness_assembler, celldofs(cell), Kₑ)
    end
    return K
end

# -----------------------------------------------------------------------
# Time integration parameters
# -----------------------------------------------------------------------
T      = 20.0
Δt₀    = 0.001
Δt_save = 0.1

# Allocate and assemble the system matrices
M = allocate_matrix(dh)
M = assemble_mass_matrix(cellvalues_v, cellvalues_p, M, dh)

K = allocate_matrix(dh)
K = assemble_stokes_matrix(cellvalues_v, cellvalues_p, ν, K, dh)

u₀ = zeros(ndofs(dh))
apply!(u₀, ch)

jac_sparsity = sparse(K)

# Apply BCs to mass matrix (zero rows/cols of constrained dofs)
apply!(M, ch)

# -----------------------------------------------------------------------
# ODE parameter struct (same pattern as NS tutorial)
# -----------------------------------------------------------------------
struct RHSparams
    K::SparseMatrixCSC
    ch::ConstraintHandler
    dh::DofHandler
    cellvalues_v::CellValues
    u::Vector
end

p = RHSparams(K, ch, dh, cellvalues_v, copy(u₀))

# -----------------------------------------------------------------------
# Step limiter: enforce Dirichlet BCs after each accepted step
# (same as NS tutorial)
# -----------------------------------------------------------------------
function ferrite_limiter!(u, _, p, t)
    update!(p.ch, t)
    return apply!(u, p.ch)
end

# -----------------------------------------------------------------------
# Transient Stokes RHS  (NS tutorial minus the convective term)
#
#   M * du/dt = K * u
#
# The convective term  -∫ (v·∇v) · φⱼ dΩ  is simply omitted.
# -----------------------------------------------------------------------
function stokes!(du, u_uc, p::RHSparams, t)
    (; K, ch, u) = p

    u .= u_uc
    update!(ch, t)
    apply!(u, ch)

    # Linear Stokes operator only — no convective term
    mul!(du, K, u)
    return
end

# -----------------------------------------------------------------------
# Jacobian of the transient Stokes RHS
# Since the problem is linear the Jacobian is just K (no nonlinear part).
# -----------------------------------------------------------------------
function stokes_jac!(J, u_uc, p::RHSparams, t)
    (; K, ch, u) = p

    u .= u_uc
    update!(ch, t)
    apply!(u, ch)

    # Jacobian equals the Stokes operator K (no nonlinear contribution)
    nonzeros(J) .= nonzeros(K)

    return apply!(J, ch)
end

# -----------------------------------------------------------------------
# Build and solve the ODE problem
# -----------------------------------------------------------------------
struct FreeDofErrorNorm
    ch::ConstraintHandler
end
(fe_norm::FreeDofErrorNorm)(u::Union{AbstractFloat, Complex}, t) = DiffEqBase.ODE_DEFAULT_NORM(u, t)
(fe_norm::FreeDofErrorNorm)(u::AbstractArray, t) = DiffEqBase.ODE_DEFAULT_NORM(u[fe_norm.ch.free_dofs], t)

rhs = ODEFunction(stokes!, mass_matrix = M; jac = stokes_jac!, jac_prototype = jac_sparsity)
problem = ODEProblem(rhs, u₀, (0.0, T), p)

timestepper = Rodas5P(autodiff = false, step_limiter! = ferrite_limiter!)

integrator = init(
    problem, timestepper; initializealg = NoInit(), dt = Δt₀,
    adaptive = true, abstol = 1.0e-4, reltol = 1.0e-5,
    progress = true, progress_steps = 1,
    verbose = true, internalnorm = FreeDofErrorNorm(ch), d_discontinuities = [1.0]
)

# -----------------------------------------------------------------------
# Time-stepping loop with VTK output
# -----------------------------------------------------------------------
pvd = paraview_collection("stokes_trans_cyl_20")
for (step, (u, t)) in enumerate(intervals(integrator))
    VTKGridFile("stokes_trans_cyl_20-$step", dh) do vtk
        write_solution(vtk, dh, u)
        pvd[t] = vtk
    end
end
vtk_save(pvd)

println("Done!")
