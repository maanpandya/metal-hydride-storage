using Ferrite, SparseArrays, BlockArrays, LinearAlgebra, WriteVTK
using DiffEqBase
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqRosenbrock: Rodas5P
using Plots

using FerriteGmsh
using FerriteGmsh: Gmsh

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
const nu = 1.0 / 1000.0      # kinematic viscosity [m^2/s]
const rho = 1.0              # density [kg/m^3]
const p0_inlet = 1.0         # prescribed total pressure at inlet [Pa]

# Total-pressure penalty parameters
const alpha = 100.0
const t_ramp = 1.0
alpha_of_t(t) = alpha * min(t / t_ramp, 1.0)

# Channel geometry (used to skip inlet corners)
const H_channel = 0.41
const tol_corner = 0.0011

# ---------------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------------
dim = 2
Gmsh.initialize()
gmsh.option.set_number("General.Verbosity", 2)

rect_tag = gmsh.model.occ.add_rectangle(0, 0, 0, 1.1, 0.41)
circle_tag = gmsh.model.occ.add_circle(0.2, 0.2, 0, 0.05)
circle_curve_tag = gmsh.model.occ.add_curve_loop([circle_tag])
circle_surf_tag = gmsh.model.occ.add_plane_surface([circle_curve_tag])
gmsh.model.occ.cut([(dim, rect_tag)], [(dim, circle_surf_tag)])
gmsh.model.occ.synchronize()

gmsh.model.model.add_physical_group(dim - 1, [6], -1, "bottom")
gmsh.model.model.add_physical_group(dim - 1, [7], -1, "left")
gmsh.model.model.add_physical_group(dim - 1, [8], -1, "right")
gmsh.model.model.add_physical_group(dim - 1, [9], -1, "top")
gmsh.model.model.add_physical_group(dim - 1, [5], -1, "hole")

gmsh.option.setNumber("Mesh.Algorithm", 11)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 20)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)

gmsh.model.mesh.generate(dim)
grid = togrid()
Gmsh.finalize()

println("Mesh: $(getncells(grid)) cells, $(getnnodes(grid)) nodes")

# ---------------------------------------------------------------------------
# FE spaces
# ---------------------------------------------------------------------------
ip_v = Lagrange{RefQuadrilateral, 2}()^dim
ip_p = Lagrange{RefQuadrilateral, 1}()
ip_g = Lagrange{RefQuadrilateral, 1}()

qr = QuadratureRule{RefQuadrilateral}(4)
cellvalues_v = CellValues(qr, ip_v, ip_g)
cellvalues_p = CellValues(qr, ip_p, ip_g)

qr_facet = FacetQuadratureRule{RefQuadrilateral}(4)
facetvalues_v = FacetValues(qr_facet, ip_v, ip_g)
facetvalues_p = FacetValues(qr_facet, ip_p, ip_g)

dh = DofHandler(grid)
add!(dh, :v, ip_v)
add!(dh, :p, ip_p)
close!(dh)

println("Total DOFs: $(ndofs(dh))")

# ---------------------------------------------------------------------------
# Boundary conditions
# ---------------------------------------------------------------------------
ch = ConstraintHandler(dh)

noslip_facet_names = ["top", "bottom", "hole"]
noslip_facets = union(getfacetset.((grid,), noslip_facet_names)...)
add!(ch, Dirichlet(:v, noslip_facets, (x, t) -> Vec((0.0, 0.0)), [1, 2]))

# Inlet has no Dirichlet velocity; it is driven by total-pressure penalty.
inlet_facets = getfacetset(grid, "left")

# Outlet pressure gauge to remove pressure null-space.
add!(ch, Dirichlet(:p, getfacetset(grid, "right"), (x, t) -> 0.0))

close!(ch)
update!(ch, 0.0)

# ---------------------------------------------------------------------------
# Matrix assembly
# ---------------------------------------------------------------------------
function assemble_mass_matrix(cellvalues_v::CellValues, cellvalues_p::CellValues, M::SparseMatrixCSC, dh::DofHandler)
    n_v = getnbasefunctions(cellvalues_v)
    n_p = getnbasefunctions(cellvalues_p)
    n = n_v + n_p
    vblk, pblk = 1, 2
    M_e = BlockedArray(zeros(n, n), [n_v, n_p], [n_v, n_p])

    assembler = start_assemble(M)
    for cell in CellIterator(dh)
        fill!(M_e, 0)
        Ferrite.reinit!(cellvalues_v, cell)

        for qp in 1:getnquadpoints(cellvalues_v)
            dOmega = getdetJdV(cellvalues_v, qp)
            for i in 1:n_v
                phi_i = shape_value(cellvalues_v, qp, i)
                for j in 1:n_v
                    phi_j = shape_value(cellvalues_v, qp, j)
                    M_e[BlockIndex((vblk, vblk), (i, j))] += phi_i ⋅ phi_j * dOmega
                end
            end
        end

        assemble!(assembler, celldofs(cell), M_e)
    end

    return M
end

function assemble_stokes_matrix(cellvalues_v::CellValues, cellvalues_p::CellValues, nu, K::SparseMatrixCSC, dh::DofHandler)
    n_v = getnbasefunctions(cellvalues_v)
    n_p = getnbasefunctions(cellvalues_p)
    n = n_v + n_p
    vblk, pblk = 1, 2
    K_e = BlockedArray(zeros(n, n), [n_v, n_p], [n_v, n_p])

    assembler = start_assemble(K)
    for cell in CellIterator(dh)
        fill!(K_e, 0)
        Ferrite.reinit!(cellvalues_v, cell)
        Ferrite.reinit!(cellvalues_p, cell)

        for qp in 1:getnquadpoints(cellvalues_v)
            dOmega = getdetJdV(cellvalues_v, qp)

            for i in 1:n_v
                grad_phi_i = shape_gradient(cellvalues_v, qp, i)
                for j in 1:n_v
                    grad_phi_j = shape_gradient(cellvalues_v, qp, j)
                    K_e[BlockIndex((vblk, vblk), (i, j))] -= nu * grad_phi_i ⊡ grad_phi_j * dOmega
                end
            end

            for j in 1:n_p
                psi_j = shape_value(cellvalues_p, qp, j)
                for i in 1:n_v
                    div_phi_i = shape_divergence(cellvalues_v, qp, i)
                    K_e[BlockIndex((vblk, pblk), (i, j))] += div_phi_i * psi_j * dOmega
                    K_e[BlockIndex((pblk, vblk), (j, i))] += psi_j * div_phi_i * dOmega
                end
            end
        end

        assemble!(assembler, celldofs(cell), K_e)
    end

    return K
end

# ---------------------------------------------------------------------------
# Time-stepping parameters
# ---------------------------------------------------------------------------
const T_end = 5.0

M = allocate_matrix(dh)
assemble_mass_matrix(cellvalues_v, cellvalues_p, M, dh)

K = allocate_matrix(dh)
assemble_stokes_matrix(cellvalues_v, cellvalues_p, nu, K, dh)

u0 = zeros(ndofs(dh))
apply!(u0, ch)

apply!(M, ch)
jac_sparsity = sparse(K)

# ---------------------------------------------------------------------------
# Inlet total-pressure model
# r(ux,p) = ux^2 - 2(p0 - p)/rho
# ---------------------------------------------------------------------------
@inline total_p_residual(ux, p_static) = ux^2 - 2.0 * (p0_inlet - p_static) / rho
@inline dtotal_p_residual_dux(ux) = 2.0 * ux
@inline dtotal_p_residual_dp() = 2.0 / rho

# ---------------------------------------------------------------------------
# RHS/Jacobian parameter container
# ---------------------------------------------------------------------------
struct RHSparams
    K::SparseMatrixCSC
    ch::ConstraintHandler
    dh::DofHandler
    cellvalues_v::CellValues
    facetvalues_v::FacetValues
    facetvalues_p::FacetValues
    inlet_facets
    u::Vector
end

params = RHSparams(K, ch, dh, cellvalues_v, facetvalues_v, facetvalues_p, inlet_facets, copy(u0))

function ferrite_limiter!(u, _, p, t)
    update!(p.ch, t)
    return apply!(u, p.ch)
end

# ---------------------------------------------------------------------------
# Nonlinear convection term
# ---------------------------------------------------------------------------
function navierstokes_rhs_element!(dv_e, v_e, cellvalues_v)
    n_basefuncs = getnbasefunctions(cellvalues_v)
    for qp in 1:getnquadpoints(cellvalues_v)
        dOmega = getdetJdV(cellvalues_v, qp)
        grad_v = function_gradient(cellvalues_v, qp, v_e)
        v = function_value(cellvalues_v, qp, v_e)
        for j in 1:n_basefuncs
            phi_j = shape_value(cellvalues_v, qp, j)
            dv_e[j] -= v ⋅ grad_v' ⋅ phi_j * dOmega
        end
    end
    return
end

function navierstokes_jac_element!(J_e, v_e, cellvalues_v)
    n_basefuncs = getnbasefunctions(cellvalues_v)
    for qp in 1:getnquadpoints(cellvalues_v)
        dOmega = getdetJdV(cellvalues_v, qp)
        grad_v = function_gradient(cellvalues_v, qp, v_e)
        v = function_value(cellvalues_v, qp, v_e)

        for j in 1:n_basefuncs
            phi_j = shape_value(cellvalues_v, qp, j)
            for i in 1:n_basefuncs
                phi_i = shape_value(cellvalues_v, qp, i)
                grad_phi_i = shape_gradient(cellvalues_v, qp, i)
                J_e[j, i] -= (phi_i ⋅ grad_v' + v ⋅ grad_phi_i') ⋅ phi_j * dOmega
            end
        end
    end
    return
end

# ---------------------------------------------------------------------------
# Inlet total-pressure residual and Jacobian additions
# ---------------------------------------------------------------------------
function add_inlet_total_pressure_rhs!(du, u, p::RHSparams, t)
    (; dh, facetvalues_v, facetvalues_p, inlet_facets) = p

    v_range = dof_range(dh, :v)
    p_range = dof_range(dh, :p)

    for facet in FacetIterator(dh, inlet_facets)
        Ferrite.reinit!(facetvalues_v, facet)
        Ferrite.reinit!(facetvalues_p, facet)

        celld = celldofs(facet)
        v_dofs = @view celld[v_range]
        p_dofs = @view celld[p_range]
        v_e = u[v_dofs]
        p_e = u[p_dofs]

        nphi_v = getnbasefunctions(facetvalues_v)
        Re = zeros(length(v_dofs))

        for qp in 1:getnquadpoints(facetvalues_v)
            x = spatial_coordinate(facetvalues_v, qp, getcoordinates(facet))
            if x[2] <= tol_corner || x[2] >= H_channel - tol_corner
                continue
            end

            dGamma = getdetJdV(facetvalues_v, qp)
            ux = function_value(facetvalues_v, qp, v_e)[1]
            p_static = function_value(facetvalues_p, qp, p_e)
            res = total_p_residual(ux, p_static)

            for i in 1:nphi_v
                phi_x_i = shape_value(facetvalues_v, qp, i)[1]
                Re[i] -= alpha_of_t(t) * res * phi_x_i * dGamma
            end
        end

        assemble!(du, v_dofs, Re)
    end

    return
end

function add_inlet_total_pressure_jac!(J, u, p::RHSparams, t)
    (; dh, facetvalues_v, facetvalues_p, inlet_facets) = p

    v_range = dof_range(dh, :v)
    p_range = dof_range(dh, :p)

    assembler = start_assemble(J; fillzero = false)

    for facet in FacetIterator(dh, inlet_facets)
        Ferrite.reinit!(facetvalues_v, facet)
        Ferrite.reinit!(facetvalues_p, facet)

        celld = celldofs(facet)
        v_dofs = @view celld[v_range]
        p_dofs = @view celld[p_range]
        v_e = u[v_dofs]

        nphi_v = getnbasefunctions(facetvalues_v)
        nphi_p = getnbasefunctions(facetvalues_p)

        Jvv = zeros(length(v_dofs), length(v_dofs))

        for qp in 1:getnquadpoints(facetvalues_v)
            x = spatial_coordinate(facetvalues_v, qp, getcoordinates(facet))
            if x[2] <= tol_corner || x[2] >= H_channel - tol_corner
                continue
            end

            dGamma = getdetJdV(facetvalues_v, qp)
            ux = function_value(facetvalues_v, qp, v_e)[1]
            dres_dux = dtotal_p_residual_dux(ux)
            dres_dp = dtotal_p_residual_dp()

            for i in 1:nphi_v
                phi_x_i = shape_value(facetvalues_v, qp, i)[1]

                for j in 1:nphi_v
                    phi_x_j = shape_value(facetvalues_v, qp, j)[1]
                    Jvv[i, j] -= alpha_of_t(t) * dres_dux * phi_x_i * phi_x_j * dGamma
                end

                global_i = v_dofs[i]
                for j in 1:nphi_p
                    psi_j = shape_value(facetvalues_p, qp, j)
                    global_j = p_dofs[j]
                    J[global_i, global_j] -= alpha_of_t(t) * dres_dp * phi_x_i * psi_j * dGamma
                end
            end
        end

        assemble!(assembler, v_dofs, Jvv)
    end

    return
end

# ---------------------------------------------------------------------------
# Full transient Navier-Stokes RHS/Jacobian
# ---------------------------------------------------------------------------
function navierstokes!(du, u_uc, p::RHSparams, t)
    (; K, ch, dh, cellvalues_v, u) = p

    u .= u_uc
    update!(ch, t)
    apply!(u, ch)

    # Linear Stokes contribution
    mul!(du, K, u)

    # Nonlinear volume convection
    v_range = dof_range(dh, :v)
    n_basefuncs = getnbasefunctions(cellvalues_v)
    v_e = zeros(n_basefuncs)
    du_e = zeros(n_basefuncs)

    for cell in CellIterator(dh)
        Ferrite.reinit!(cellvalues_v, cell)
        v_celldofs = @view celldofs(cell)[v_range]

        v_e .= @views u[v_celldofs]
        fill!(du_e, 0.0)
        navierstokes_rhs_element!(du_e, v_e, cellvalues_v)
        assemble!(du, v_celldofs, du_e)
    end

    # Nonlinear inlet total-pressure contribution
    add_inlet_total_pressure_rhs!(du, u, p, t)

    return
end

function navierstokes_jac!(J, u_uc, p::RHSparams, t)
    (; K, ch, dh, cellvalues_v, u) = p

    u .= u_uc
    update!(ch, t)
    apply!(u, ch)

    # Start from linear Stokes operator
    nonzeros(J) .= nonzeros(K)

    assembler = start_assemble(J; fillzero = false)

    # Convection Jacobian
    n_basefuncs = getnbasefunctions(cellvalues_v)
    J_e = zeros(n_basefuncs, n_basefuncs)
    v_e = zeros(n_basefuncs)
    v_range = dof_range(dh, :v)

    for cell in CellIterator(dh)
        Ferrite.reinit!(cellvalues_v, cell)
        v_celldofs = @view celldofs(cell)[v_range]

        v_e .= @views u[v_celldofs]
        fill!(J_e, 0.0)
        navierstokes_jac_element!(J_e, v_e, cellvalues_v)
        assemble!(assembler, v_celldofs, J_e)
    end

    # Inlet total-pressure Jacobian
    add_inlet_total_pressure_jac!(J, u, p, t)

    return apply!(J, ch)
end

rhs = ODEFunction(
    navierstokes!,
    mass_matrix = M;
    jac = navierstokes_jac!,
    jac_prototype = jac_sparsity,
)
problem = ODEProblem(rhs, u0, (0.0, T_end), params)

struct FreeDofErrorNorm
    ch::ConstraintHandler
end
(fe_norm::FreeDofErrorNorm)(u::Union{AbstractFloat, Complex}, t) = DiffEqBase.ODE_DEFAULT_NORM(u, t)
(fe_norm::FreeDofErrorNorm)(u::AbstractArray, t) = DiffEqBase.ODE_DEFAULT_NORM(u[fe_norm.ch.free_dofs], t)

timestepper = Rodas5P(autodiff = false, step_limiter! = ferrite_limiter!)

# Let the solver choose adaptive time steps automatically.
solve_time = @elapsed sol = solve(
    problem,
    timestepper;
    adaptive = true,
    abstol = 1.0e-4,
    reltol = 1.0e-5,
    progress = true,
    progress_steps = 1,
    verbose = true,
    internalnorm = FreeDofErrorNorm(ch),
    d_discontinuities = [t_ramp],
    save_everystep = true,
)

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
import Base.Filesystem: mkpath
mkpath("sol/sol_ns_p0_inlet_auto_dt")

println("Total solver wall-clock time: $(round(solve_time; digits = 3)) seconds")
println("Adaptive solver chose $(length(sol.t)) saved time points.")

println("Accepted adaptive time steps:")
if length(sol.t) >= 2
    for i in 2:length(sol.t)
        dt_i = sol.t[i] - sol.t[i - 1]
        println("Step $(i - 1): t = $(sol.t[i]), dt = $(dt_i)")
    end
else
    println("Only one time point was saved; no dt history is available.")
end

if length(sol.t) >= 2
    dt_history = diff(sol.t)
    t_history = sol.t[2:end]

    dt_plot = plot(
        t_history,
        dt_history;
        xlabel = "time",
        ylabel = "dt",
        title = "Adaptive time-step evolution (Transient Navier-Stokes + p0 inlet)",
        lw = 2,
        marker = :circle,
        legend = false,
    )
    savefig(dt_plot, "sol/sol_ns_p0_inlet_auto_dt/navier-stokes-p0-inlet-auto-dt-history.png")
    println("Saved time-step plot to 'sol/sol_ns_p0_inlet_auto_dt/navier-stokes-p0-inlet-auto-dt-history_20.png'.")
end

pvd = paraview_collection("sol/sol_ns_p0_inlet_auto_dt/navier-stokes-p0-inlet-auto-dt_6")
for (step, (u, t)) in enumerate(zip(sol.u, sol.t))
    VTKGridFile("sol/sol_ns_p0_inlet_auto_dt/navier-stokes-p0-inlet-auto-dt-$(lpad(step, 4, '0'))", dh) do vtk
        write_solution(vtk, dh, u)
        pvd[t] = vtk
    end
end
vtk_save(pvd)

println("Done! Open 'sol/sol_ns_p0_inlet_auto_dt/navier-stokes-p0-inlet-auto-dt.pvd' in ParaView.")
