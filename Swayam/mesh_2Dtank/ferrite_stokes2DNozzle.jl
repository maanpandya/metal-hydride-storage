# ============================================================================
# 2D Tank-Nozzle Stokes Flow Simulation using Ferrite.jl
# ============================================================================
# Solves Stokes flow in the 2D tank-nozzle geometry created by gmsh_2DtankNozzle.jl
# Uses the mesh with proper boundary condition labels: "walls" and "inlet"

using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

function setup_tank_grid()
    # Check if the tank mesh file exists
    mesh_file = "tank_2D_with_BCs.msh"
    if !isfile(mesh_file)
        error("Tank mesh file $(mesh_file) not found. Please run gmsh_2DtankNozzle.jl first to generate the mesh.")
    end
    
    println("Loading 2D tank-nozzle mesh from $(mesh_file)...")
    
    # Read the tank mesh file as a Ferrite Grid
    grid = FerriteGmsh.togrid(mesh_file)
    
    println("✅ Tank grid loaded successfully!")
    println("   - Number of cells: $(getncells(grid))")
    println("   - Number of nodes: $(getnnodes(grid))")
    
    # Print available boundary sets (facesets in Ferrite.jl)
    println("\n📋 Available boundary sets:")
    try
        # Try to access facesets - the exact field name may vary by Ferrite version
        if hasfield(typeof(grid), :facesets)
            for (name, set) in grid.facesets
                println("   - $(name): $(length(set)) faces")
            end
        elseif hasfield(typeof(grid), :facetsets)
            for (name, set) in grid.facetsets
                println("   - $(name): $(length(set)) faces")
            end
        else
            println("   - Boundary set information not directly accessible")
            println("   - Grid loaded successfully, proceeding with simulation")
        end
    catch e
        println("   - Could not access boundary sets: $(e)")
        println("   - Grid loaded successfully, proceeding with simulation")
    end
    
    return grid
end

function setup_fevalues(ipu, ipp, ipg)
    # Note: Tank mesh uses quadrilaterals, so we need RefQuadrilateral
    qr = QuadratureRule{RefQuadrilateral}(2)
    cvu = CellValues(qr, ipu, ipg)
    cvp = CellValues(qr, ipp, ipg)
    qr_facet = FacetQuadratureRule{RefQuadrilateral}(2)
    fvp = FacetValues(qr_facet, ipp, ipg)
    return cvu, cvp, fvp
end

function setup_dofs(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)  # Velocity field
    add!(dh, :p, ipp)  # Pressure field
    close!(dh)
    return dh
end

function setup_mean_constraint(dh, fvp)
    assembler = Ferrite.COOAssembler()
    
    # For the tank geometry, use all external boundaries (walls + inlet) for pressure constraint
    set = union(
        getfacetset(dh.grid, "walls"),
        getfacetset(dh.grid, "inlet"),
    )
    
    # Allocate buffers
    range_p = dof_range(dh, :p)
    element_dofs = zeros(Int, ndofs_per_cell(dh))
    element_dofs_p = view(element_dofs, range_p)
    element_coords = zeros(Vec{2}, 4)  # 4 nodes for quadrilateral elements
    Ce = zeros(1, length(range_p)) # Local constraint matrix (only 1 row)
    
    # Loop over all the boundaries
    for (ci, fi) in set
        Ce .= 0
        getcoordinates!(element_coords, dh.grid, ci)
        reinit!(fvp, element_coords, fi)
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
        Pair{Int, Float64}[J[i] => -V[i] for i in 1:length(J) if J[i] != constrained_dof],
        0.0,
    )
    return mean_value_constraint
end

function setup_constraints(dh, fvp)
    ch = ConstraintHandler(dh)
    
    # Get boundary sets from the tank mesh
    walls_set = getfacetset(dh.grid, "walls")
    inlet_set = getfacetset(dh.grid, "inlet")
    
    println("\n🔧 Setting up boundary conditions:")
    println("   - Walls: $(length(walls_set)) faces with no-slip condition (u = 0)")
    println("   - Inlet: $(length(inlet_set)) faces with prescribed velocity")
    
    # No-slip boundary condition on walls (u = 0)
    dbc_walls = Dirichlet(:u, walls_set, (x, t) -> [0.0, 0.0], [1, 2])
    add!(ch, dbc_walls)
    
    # Inlet boundary condition - prescribe a simple inflow profile
    # For now, let's use a uniform inflow velocity
    inlet_velocity = 1.0  # m/s inflow in x-direction
    dbc_inlet = Dirichlet(:u, inlet_set, (x, t) -> [inlet_velocity, 0.0], [1, 2])
    add!(ch, dbc_inlet)
    
    # Add mean pressure constraint to fix pressure level
    mean_value_constraint = setup_mean_constraint(dh, fvp)
    add!(ch, mean_value_constraint)
    
    # Finalize
    close!(ch)
    update!(ch, 0)
    return ch
end

function assemble_system!(K, f, dh, cvu, cvp)
    assembler = start_assemble(K, f)
    ke = zeros(ndofs_per_cell(dh), ndofs_per_cell(dh))
    fe = zeros(ndofs_per_cell(dh))
    range_u = dof_range(dh, :u)
    ndofs_u = length(range_u)
    range_p = dof_range(dh, :p)
    ndofs_p = length(range_p)
    ϕᵤ = Vector{Vec{2, Float64}}(undef, ndofs_u)
    ∇ϕᵤ = Vector{Tensor{2, 2, Float64, 4}}(undef, ndofs_u)
    divϕᵤ = Vector{Float64}(undef, ndofs_u)
    ϕₚ = Vector{Float64}(undef, ndofs_p)
    
    # Material properties
    μ = 1.0  # Dynamic viscosity (Pa·s)
    
    for cell in CellIterator(dh)
        reinit!(cvu, cell)
        reinit!(cvp, cell)
        ke .= 0
        fe .= 0
        for qp in 1:getnquadpoints(cvu)
            dΩ = getdetJdV(cvu, qp)
            for i in 1:ndofs_u
                ϕᵤ[i] = shape_value(cvu, qp, i)
                ∇ϕᵤ[i] = shape_gradient(cvu, qp, i)
                divϕᵤ[i] = shape_divergence(cvu, qp, i)
            end
            for i in 1:ndofs_p
                ϕₚ[i] = shape_value(cvp, qp, i)
            end
            
            # Stokes equations assembly
            # u-u: viscous term μ∇u⋅∇v
            for (i, I) in pairs(range_u), (j, J) in pairs(range_u)
                ke[I, J] += μ * (∇ϕᵤ[i] ⊡ ∇ϕᵤ[j]) * dΩ
            end
            
            # u-p: pressure term -p∇⋅v
            for (i, I) in pairs(range_u), (j, J) in pairs(range_p)
                ke[I, J] += (-divϕᵤ[i] * ϕₚ[j]) * dΩ
            end
            
            # p-u: continuity term -q∇⋅u
            for (i, I) in pairs(range_p), (j, J) in pairs(range_u)
                ke[I, J] += (-divϕᵤ[j] * ϕₚ[i]) * dΩ
            end
            
            # Right-hand side: body force (gravity or other external forces)
            # For this example, let's add a small body force to drive flow
            for (i, I) in pairs(range_u)
                x = spatial_coordinate(cvu, qp, getcoordinates(cell))
                # Small gravitational force in y-direction
                body_force = Vec{2}((0.0, -0.1))
                fe[I] += (ϕᵤ[i] ⋅ body_force) * dΩ
            end
        end
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    return K, f
end

function main()
    println("🚀 Starting 2D Tank-Nozzle Stokes Flow Simulation")
    println("=" ^ 60)
    
    # Load the tank grid
    grid = setup_tank_grid()
    
    # Interpolations - use quadrilateral elements to match the tank mesh
    ipu = Lagrange{RefQuadrilateral, 2}()^2  # Quadratic velocity (biquadratic)
    ipp = Lagrange{RefQuadrilateral, 1}()    # Linear pressure (bilinear)
    
    println("\n🔧 Setting up finite element spaces:")
    println("   - Velocity: Quadratic (Q2) elements")
    println("   - Pressure: Linear (Q1) elements")
    
    # DOF handler
    dh = setup_dofs(grid, ipu, ipp)
    println("   - Total DOFs: $(ndofs(dh))")
    
    # Calculate DOFs for each field manually
    u_range = dof_range(dh, :u)
    p_range = dof_range(dh, :p)
    println("   - Velocity DOFs: $(length(u_range))")
    println("   - Pressure DOFs: $(length(p_range))")
    
    # FE values
    ipg = Lagrange{RefQuadrilateral, 1}()  # Linear geometric interpolation
    cvu, cvp, fvp = setup_fevalues(ipu, ipp, ipg)
    
    # Boundary conditions
    ch = setup_constraints(dh, fvp)
    
    # Global tangent matrix and rhs
    println("\n🔧 Assembling system matrix...")
    coupling = [true true; true false] # no coupling between pressure test/trial functions
    K = allocate_matrix(dh, ch; coupling = coupling)
    f = zeros(ndofs(dh))
    
    # Assemble system
    assemble_system!(K, f, dh, cvu, cvp)
    println("   ✅ System assembled")
    
    # Apply boundary conditions and solve
    println("\n🔧 Solving linear system...")
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)
    println("   ✅ Solution computed")
    
    # Export the solution
    println("\n💾 Exporting solution...")
    VTKGridFile("tank-stokes-flow", grid) do vtk
        write_solution(vtk, dh, u)
    end
    println("   ✅ Solution exported to tank-stokes-flow.vtu")
    
    # Print solution statistics
    u_vals = u[dof_range(dh, :u)]
    p_vals = u[dof_range(dh, :p)]
    
    println("\n📊 Solution Statistics:")
    println("   - Max velocity magnitude: $(round(maximum(norm, reshape(u_vals, 2, :)), digits=6))")
    println("   - Max pressure: $(round(maximum(p_vals), digits=6))")
    println("   - Min pressure: $(round(minimum(p_vals), digits=6))")
    
    println("\n✅ Tank flow simulation completed successfully!")
    println("   Open 'tank-stokes-flow.vtu' in ParaView to visualize the results.")
    
    return grid, dh, u
end

# Run the simulation
main()
