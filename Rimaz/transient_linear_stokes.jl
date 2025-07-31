using Ferrite             # Finite element framework
using SparseArrays         # Sparse matrix operations
using LinearAlgebra        # Basic linear algebra
using BlockArrays          # Blocked matrix operations
using DifferentialEquations    # ODE/DAE solving framework
using OrdinaryDiffEq          # Advanced ODE solvers
using DiffEqCallbacks        # Callback functions for DifferentialEquations
using Statistics             # Statistical functions
using Plots                  # Plotting and visualization
using UnPack               # For @unpack macro
using WriteVTK             # For VTK output and paraview_collection


# Import specific functions from Ferrite
using Ferrite: update!, apply!
using OrdinaryDiffEq: TimeChoiceIterator



# Domain Geometry
H = 0.25                        # Channel height [m]
L = 4 * H                       # Channel length [m] (aspect ratio 4:1)

# Mesh Parameters
nelem_base = 8                  # Reduced for faster testing
nels = (4 * nelem_base, nelem_base)  # Elements in (x,y) directions

# Physical Properties
μ = 1e-3                        # Dynamic viscosity [Pa·s]
ρ = 1                      #

# Inlet Velocity Parameters
v_max = 20.0                     # Maximum inlet velocity [m/s]
ramp_time = 3                   # Time to reach full velocity [s]

# Finite Element Parameters
dim = 2                         # Spatial dimension
degree = 1                      # Base polynomial degree

# Solver Tolerances
reltol = 1e-4                   # Relative tolerance
abstol = 1e-5                   # Absolute tolerance


# Generate mesh and setup finite elements
left_corner = Vec((0.0, 0.0))     # Bottom-left corner
right_corner = Vec((L, H))        # Top-right corner

# Generate structured quadrilateral mesh
grid = generate_grid(Quadrilateral, nels, left_corner, right_corner)
addvertexset!(grid, "corner", (x) -> x[1] ≈ 0.0 && x[2] ≈ 0.0)
# Define interpolation spaces (Taylor-Hood elements: Q₂-Q₁)
ipu = Lagrange{RefQuadrilateral, degree+1}()^dim     # Q₂ × Q₂ (vector)
ipp = Lagrange{RefQuadrilateral, degree}()           # Q₁ (scalar)

# Create DOF handler
dh = DofHandler(grid)
add!(dh, :u, ipu)    # Add velocity field
add!(dh, :p, ipp)    # Add pressure field
close!(dh)

# Quadrature rules for integration
qr = QuadratureRule{RefQuadrilateral}(2*degree+1)    # Higher order for accuracy
ipg = Lagrange{RefQuadrilateral, 1}()                # Linear geometric mapping

# Cell values for volume integration
cvu = CellValues(qr, ipu, ipg)    # Velocity cell values
cvp = CellValues(qr, ipp, ipg)    # Pressure cell values
qr_facet = FacetQuadratureRule{RefQuadrilateral}(2)
fvp = FacetValues(qr_facet, ipp, ipg) # required for pressure constraint 

function setup_mean_constraint(dh, fvp)
    assembler = Ferrite.COOAssembler()
    # All external boundaries
    set = union(
            getfacetset(dh.grid, "left"),
            getfacetset(dh.grid, "right"),
            getfacetset(dh.grid, "bottom"),
            getfacetset(dh.grid, "top"),
    )
    # Allocate buffers
    range_p = dof_range(dh, :p)
    element_dofs = zeros(Int, ndofs_per_cell(dh))
    element_dofs_p = view(element_dofs, range_p)
    element_coords = zeros(Vec{2}, 4) # assuming 2D mesh with quadrilaterals only 
    Ce = zeros(1, length(range_p)) # Local constraint matrix (only 1 row)
    # Loop over all the boundaries
    for (ci, fi) in set
        Ce .= 0
        getcoordinates!(element_coords, dh.grid, ci)
        Ferrite.reinit!(fvp, element_coords, fi)
        celldofs!(element_dofs, dh, ci)
        for qp in 1:getnquadpoints(fvp)
            dΓ = getdetJdV(fvp, qp)
            for i in 1:getnbasefunctions(fvp)
                Ce[1, i] += shape_value(fvp, qp, i) * dΓ
            end
        end
        # Assemble to row 1
        assemble!(assembler, [1], element_dofs_p, Ce)
    end
    C, _ = finish_assemble(assembler)
    # Create an AffineConstraint from the C-matrix
    _, J, V = findnz(C)
    _, constrained_dof_idx = findmax(abs2, V)
    constrained_dof = J[constrained_dof_idx]
    V ./= V[constrained_dof_idx]
    mean_value_constraint = AffineConstraint(
        constrained_dof,
        Pair{Int,Float64}[J[i] => -V[i] for i in 1:length(J) if J[i] != constrained_dof],
        0.0,
    )

    return mean_value_constraint
end

println("Transient Stokes problem setup with Taylor-Hood elements")

# boundary conditions
# Create constraint handler
ch = ConstraintHandler(dh)

# 1. INLET BOUNDARY CONDITION (LEFT) - Time-dependent parabolic profile

# Time-dependent inlet velocity function with smooth ramp-up

# Parabolic velocity profile function
function parabolic_inflow_profile(x, t, v_max, H, ramp_time)

    y = x[2]
    
    # Get time-dependent velocity magnitude
    v_t = min(t * v_max / ramp_time, v_max)  # Linear ramp to v_max
    
    # Parabolic profile: u_x = v_t * 4*y*(H-y)/H^2
    # This gives maximum velocity v_t at y = H/2 and zero at y = 0, H
    u_x = v_t * 4.0 * y * (H - y) / (H^2)
    u_y = 0.0  # No y-component at inlet
    
    return Vec((u_x, u_y))
end
# vᵢₙ(t) = min(t * 5, 5) #inflow velocity

# parabolic_inflow_profile(x, t) = Vec((10 * vᵢₙ(t) * x[2] * (0.25 - x[2]) / 0.25^2, 0.0))

inlet_bc = Dirichlet(:u, getfacetset(grid, "left"), 
                    (x, t) -> parabolic_inflow_profile(x, t, v_max,H, ramp_time), [1,2])
add!(ch, inlet_bc)

# 2. NO-SLIP WALL CONDITIONS
# Right wall (outlet/wall)
right_wall_bc = Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> Vec{2}((0.0, 0.0)))
add!(ch, right_wall_bc)

# Top wall  
top_wall_bc = Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> Vec{2}((0.0, 0.0)))
add!(ch, top_wall_bc)

# Bottom wall
bottom_wall_bc = Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> Vec{2}((0.0, 0.0)))
add!(ch, bottom_wall_bc)

# 3. PIN PRESSURE (remove nullspace for closed domain)
# Pin pressure at first pressure node (typically a corner)
# pressure_pin = Dirichlet(:p, Set([1]), (x, t) -> 0.0)
# add!(ch, pressure_pin)
mean_value_constraint = setup_mean_constraint(dh, fvp)
add!(ch, mean_value_constraint)

# Close constraint handler
close!(ch)

# Apply constraints to matrices (at t=0 initially)
update!(ch, 0.0)

println("Constraint handler setup complete")
# Assemble mass matrix

function assemble_mass_matrix!(M::SparseMatrixCSC, dh::DofHandler, cvu::CellValues, cvp::CellValues)
    """Assemble the mass matrix for the transient Stokes system."""
    
    # Get dimensions for block structure
    n_basefuncs_v = getnbasefunctions(cvu)
    n_basefuncs_p = getnbasefunctions(cvp)
    n_basefuncs = n_basefuncs_v + n_basefuncs_p
    
    # Block indices: 1 = velocity, 2 = pressure
    v▄, p▄ = 1, 2
    
    # Local element mass matrix (blocked)
    Mₑ = BlockedArray(zeros(n_basefuncs, n_basefuncs), 
                     [n_basefuncs_v, n_basefuncs_p], 
                     [n_basefuncs_v, n_basefuncs_p])

    # Initialize assembler
    mass_assembler = start_assemble(M)
    
    # Element loop
    for cell in CellIterator(dh)
        fill!(Mₑ, 0)
        Ferrite.reinit!(cvu, cell)

        # Quadrature loop
        for q_point in 1:getnquadpoints(cvu)
            dΩ = getdetJdV(cvu, q_point)
            
            # Velocity-velocity mass matrix: ∫ Nᵢ · Nⱼ dΩ
            for i in 1:n_basefuncs_v
                φᵢ = shape_value(cvu, q_point, i)
                for j in 1:n_basefuncs_v
                    φⱼ = shape_value(cvu, q_point, j)
                    # Dot product for vector shape functions
                    Mₑ[BlockIndex((v▄, v▄), (i, j))] += φᵢ ⋅ φⱼ * dΩ
                end
            end
        end
        
        # Assemble local matrix to global matrix
        assemble!(mass_assembler, celldofs(cell), Mₑ)
    end

    return M
end

function assemble_stokes_matrix!(K, dh, cvu, cvp, viscosity)
    """
    Assemble the stiffness matrix for the Stokes system.
    
    The stiffness matrix K has the block structure:
    K = [K_uu  K_up]
        [K_pu   0  ]
    
    where:
    - K_uu: viscous term (∫ μ ∇u : ∇v dΩ)
    - K_up: pressure gradient term (∫ -p ∇·v dΩ)
    - K_pu: continuity constraint (∫ -q ∇·u dΩ)
    """
    
    assembler = start_assemble(K)
    ke = zeros(ndofs_per_cell(dh), ndofs_per_cell(dh))
    
    # Get DOF ranges for velocity and pressure
    range_u = dof_range(dh, :u)
    ndofs_u = length(range_u)
    range_p = dof_range(dh, :p)
    ndofs_p = length(range_p)
    
    # Pre-allocate arrays for shape function values
    ϕᵤ = Vector{Vec{2,Float64}}(undef, ndofs_u)        # Velocity shape functions
    ∇ϕᵤ = Vector{Tensor{2,2,Float64,4}}(undef, ndofs_u) # Velocity gradients (2×2 tensor)
    divϕᵤ = Vector{Float64}(undef, ndofs_u)             # Velocity divergences
    ϕₚ = Vector{Float64}(undef, ndofs_p)                # Pressure shape functions
    
    # Element loop
    for cell in CellIterator(dh)
        Ferrite.reinit!(cvu, cell)
        Ferrite.reinit!(cvp, cell)
        ke .= 0
        
        # Quadrature loop
        for qp in 1:getnquadpoints(cvu)
            dΩ = getdetJdV(cvu, qp)
            
            # Evaluate shape functions at quadrature point
            for i in 1:ndofs_u
                ϕᵤ[i] = shape_value(cvu, qp, i)
                ∇ϕᵤ[i] = shape_gradient(cvu, qp, i)
                divϕᵤ[i] = shape_divergence(cvu, qp, i)
            end
            for i in 1:ndofs_p
                ϕₚ[i] = shape_value(cvp, qp, i)
            end
            
            # K_uu: Viscous term ∫ μ ∇u : ∇v dΩ
            for (i, I) in pairs(range_u), (j, J) in pairs(range_u)
                ke[I, J] += viscosity * (∇ϕᵤ[i] ⊡ ∇ϕᵤ[j]) * dΩ
            end
            
            # K_up: Pressure gradient term ∫ -p ∇·v dΩ
            for (i, I) in pairs(range_u), (j, J) in pairs(range_p)
                ke[I, J] += (-divϕᵤ[i] * ϕₚ[j]) * dΩ
            end
            
            # K_pu: Continuity constraint ∫ -q ∇·u dΩ
            for (i, I) in pairs(range_p), (j, J) in pairs(range_u)
                ke[I, J] += (-divϕᵤ[j] * ϕₚ[i]) * dΩ
            end
        end
        
        # Assemble local matrix to global matrix
        assemble!(assembler, celldofs(cell), ke)
    end
    
    return K
end

T = 15.0
Δt₀ = 0.001
Δt_save = 0.1

M = allocate_matrix(dh)  # Mass matrix
M = assemble_mass_matrix!(M, dh, cvu, cvp)

# K = allocate_matrix(dh)  # Stiffness matrix
# K = assemble_stokes_matrix!(K, dh, cvu, cvp, μ)
coupling = [true true; true false] # no coupling between pressure test/trial functions
K = allocate_matrix(dh, ch; coupling = coupling)
u₀ = zeros(ndofs(dh))
apply!(u₀, ch);
println("stiffness matrix...")

jac_sparsity = sparse(K);

apply!(M, ch)

struct RHSparams
    K::SparseMatrixCSC
    ch::ConstraintHandler
    dh::DofHandler
    cellvalues_v::CellValues
    u::Vector
end
p = RHSparams(K, ch, dh, cvu, copy(u₀))

function ferrite_limiter!(u, _, p, t)
    update!(p.ch, t)
    return apply!(u, p.ch)
end


function stokes_rhs!(du, u_uc, p::RHSparams, t)
    
    @unpack K, ch, dh, cellvalues_v, u = p

    u .= u_uc  # Update solution with current values
    update!(ch, t)  # Update constraints for current time
    apply!(u, ch)  # Apply constraints to solution

    # mul!(du, K, u)   # du = -K*u
    mul!(du, K, u, -1.0, 0.0)

    return
end


function stokes_jacobian!(J, u_uc, p, t)

    @unpack K, ch, dh, cellvalues_v, u = p

    u .= u_uc  # Update solution with current values
    update!(ch, t)  # Update constraints for current time
    apply!(u, ch)  # Apply constraints to solution

    J .= -K
    apply!(J, ch)
    return
end

f_ode = ODEFunction(stokes_rhs!, mass_matrix = M; jac = stokes_jacobian!, jac_prototype=sparse(K))

problem = ODEProblem(f_ode, u₀, (0.0, T), p)

struct FreeDofErrorNorm
    ch::ConstraintHandler
end
(fe_norm::FreeDofErrorNorm)(u::Union{AbstractFloat, Complex}, t) = DiffEqBase.ODE_DEFAULT_NORM(u, t)
(fe_norm::FreeDofErrorNorm)(u::AbstractArray, t) = DiffEqBase.ODE_DEFAULT_NORM(u[fe_norm.ch.free_dofs], t)

timestepper = Rodas5P(autodiff = false, step_limiter! = ferrite_limiter!);

integrator = init(
    problem, timestepper; initializealg = NoInit(), dt = Δt₀,
    adaptive = true, abstol = 1.0e1, reltol = 1.0e1,
    progress = true, progress_steps = 1,
    verbose = true, internalnorm = FreeDofErrorNorm(ch), d_discontinuities = [1.0]
);
pvd = paraview_collection("ns_transient_r")
for (step, (u, t)) in enumerate(intervals(integrator))
    VTKGridFile("ns_transient_r-$step", dh) do vtk
        write_solution(vtk, dh, u)
        pvd[t] = vtk
    end
end
vtk_save(pvd);
