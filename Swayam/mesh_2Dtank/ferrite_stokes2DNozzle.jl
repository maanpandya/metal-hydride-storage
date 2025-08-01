# 2D Tank-Nozzle Stokes Flow Simulation using Ferrite.jl

using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays

const INPUT_MESH_FILE_NAME = "meshes_generated/2D_meshes/quadtrilateral_meshes/tank_2D_quad_fine.msh" # add .msh extension
const OUTPUT_SOLUTION_FILE_NAME = "tank_2D_quad_fine" # dont add .vtu extension

function setup_tank_grid()
    mesh_file = INPUT_MESH_FILE_NAME
    if !isfile(mesh_file)
        error("Tank mesh file $(mesh_file) not found. Please run gmsh_2DtankNozzle.jl first to generate the mesh.")
    end
    
    grid = FerriteGmsh.togrid(mesh_file)
    return grid
end

function setup_fevalues(ipu, ipp, ipg)
    qr = QuadratureRule{RefQuadrilateral}(2)
    cvu = CellValues(qr, ipu, ipg)
    cvp = CellValues(qr, ipp, ipg)
    qr_facet = FacetQuadratureRule{RefQuadrilateral}(2)
    fvp = FacetValues(qr_facet, ipp, ipg)
    return cvu, cvp, fvp
end

function setup_dofs(grid, ipu, ipp)
    dh = DofHandler(grid)
    add!(dh, :u, ipu)
    add!(dh, :p, ipp)
    close!(dh)
    return dh
end

function setup_mean_constraint(dh, fvp)
    assembler = Ferrite.COOAssembler()
    
    set = union(
        getfacetset(dh.grid, "walls"),
        getfacetset(dh.grid, "inlet"),
    )
    
    range_p = dof_range(dh, :p)
    element_dofs = zeros(Int, ndofs_per_cell(dh))
    element_dofs_p = view(element_dofs, range_p)
    element_coords = zeros(Vec{2}, 4)
    Ce = zeros(1, length(range_p))
    
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
        assemble!(assembler, [1], element_dofs_p, Ce)
    end
    C, _ = finish_assemble(assembler)
    
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
    
    walls_set = getfacetset(dh.grid, "walls")
    inlet_set = getfacetset(dh.grid, "inlet")
    
    # No-slip boundary condition on walls
    dbc_walls = Dirichlet(:u, walls_set, (x, t) -> [0.0, 0.0], [1, 2])
    add!(ch, dbc_walls)
    
    # Inlet boundary condition
    inlet_velocity = 1.0
    dbc_inlet = Dirichlet(:u, inlet_set, (x, t) -> [inlet_velocity, 0.0], [1, 2])
    add!(ch, dbc_inlet)
    
    # Add mean pressure constraint
    mean_value_constraint = setup_mean_constraint(dh, fvp)
    add!(ch, mean_value_constraint)
    
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
    
    μ = 1.0  # Dynamic viscosity
    
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
            
            # u-u: viscous term
            for (i, I) in pairs(range_u), (j, J) in pairs(range_u)
                ke[I, J] += μ * (∇ϕᵤ[i] ⊡ ∇ϕᵤ[j]) * dΩ
            end
            
            # u-p: pressure term
            for (i, I) in pairs(range_u), (j, J) in pairs(range_p)
                ke[I, J] += (-divϕᵤ[i] * ϕₚ[j]) * dΩ
            end
            
            # p-u: continuity term
            for (i, I) in pairs(range_p), (j, J) in pairs(range_u)
                ke[I, J] += (-divϕᵤ[j] * ϕₚ[i]) * dΩ
            end
            
            # Right-hand side: body force
            for (i, I) in pairs(range_u)
                x = spatial_coordinate(cvu, qp, getcoordinates(cell))
                body_force = Vec{2}((0.0, -0.1))
                fe[I] += (ϕᵤ[i] ⋅ body_force) * dΩ
            end
        end
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    return K, f
end

function main()
    # Load the tank grid
    grid = setup_tank_grid()
    
    # Interpolations
    ipu = Lagrange{RefQuadrilateral, 2}()^2
    ipp = Lagrange{RefQuadrilateral, 1}()
    
    # DOF handler
    dh = setup_dofs(grid, ipu, ipp)
    
    # FE values
    ipg = Lagrange{RefQuadrilateral, 1}()
    cvu, cvp, fvp = setup_fevalues(ipu, ipp, ipg)
    
    # Boundary conditions
    ch = setup_constraints(dh, fvp)
    
    # Global tangent matrix and rhs
    coupling = [true true; true false]
    K = allocate_matrix(dh, ch; coupling = coupling)
    f = zeros(ndofs(dh))
    
    # Assemble system
    assemble_system!(K, f, dh, cvu, cvp)
    
    # Apply boundary conditions and solve
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)
    
    # Export the solution
    VTKGridFile(OUTPUT_SOLUTION_FILE_NAME, grid) do vtk
        write_solution(vtk, dh, u)
    end
    
    return grid, dh, u
end

# Run the simulation
main()
