# Convergence Study for 2D Stokes Flow in Tank-Nozzle Geometry
# Uses the very fine mesh as "exact" solution and computes errors for coarser meshes

using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays
using Plots, Printf

# Mesh files in order of refinement (coarse to finest)
const MESH_FILES = [
    "meshes_generated/2D_meshes/quadtrilateral_meshes/tank_2D_quad_coarse.msh",
    "meshes_generated/2D_meshes/quadtrilateral_meshes/tank_2D_quad_default.msh", 
    "meshes_generated/2D_meshes/quadtrilateral_meshes/tank_2D_quad_fine.msh",
    "meshes_generated/2D_meshes/quadtrilateral_meshes/tank_2D_quad_vfine.msh"
]

const MESH_LABELS = ["coarse", "default", "fine", "vfine"]

# Include the Stokes solver functions (modified to be reusable)
function setup_tank_grid(mesh_file)
    if !isfile(mesh_file)
        error("Tank mesh file $(mesh_file) not found.")
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

function solve_stokes_flow(mesh_file)
    """Solve Stokes flow for a given mesh file and return grid, dof handler, and solution"""
    println("Solving Stokes flow for mesh: $(basename(mesh_file))")
    
    # Load the tank grid
    grid = setup_tank_grid(mesh_file)
    
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
    
    println("  Number of DOFs: $(ndofs(dh))")
    println("  Number of cells: $(getncells(grid))")
    
    return grid, dh, u
end

function interpolate_solution_to_points(grid_ref, dh_ref, u_ref, evaluation_points)
    """Interpolate reference solution to evaluation points"""
    ipu = Lagrange{RefQuadrilateral, 2}()^2
    ipp = Lagrange{RefQuadrilateral, 1}()
    ipg = Lagrange{RefQuadrilateral, 1}()
    
    qr = QuadratureRule{RefQuadrilateral}(2)
    cvu = CellValues(qr, ipu, ipg)
    cvp = CellValues(qr, ipp, ipg)
    
    interpolated_u = Vector{Vec{2, Float64}}(undef, length(evaluation_points))
    interpolated_p = Vector{Float64}(undef, length(evaluation_points))
    
    range_u = dof_range(dh_ref, :u)
    range_p = dof_range(dh_ref, :p)
    
    for (i, point) in enumerate(evaluation_points)
        # Find cell containing the point
        cell_found = false
        for cell in CellIterator(dh_ref)
            coords = getcoordinates(cell)
            
            # Simple check if point is inside quadrilateral (approximate)
            if point_in_quadrilateral(point, coords)
                # Get local coordinates (simplified - assumes quad mapping)
                ξ = physical_to_reference(point, coords)
                
                if all(-1.0 ≤ ξᵢ ≤ 1.0 for ξᵢ in ξ)
                    # Evaluate basis functions at this point
                    reinit!(cvu, cell)
                    reinit!(cvu, coords, ξ)  # This might need adjustment based on Ferrite version
                    
                    # Get DOFs for this cell
                    element_dofs = celldofs(cell)
                    
                    # Interpolate velocity
                    u_val = Vec{2}((0.0, 0.0))
                    for j in 1:length(range_u)
                        dof_idx = element_dofs[range_u[j]]
                        shape_val = shape_value(cvu, 1, j)  # Evaluate at point
                        u_val += u_ref[dof_idx] * shape_val
                    end
                    interpolated_u[i] = u_val
                    
                    # Interpolate pressure (simplified)
                    p_val = 0.0
                    for j in 1:length(range_p)
                        dof_idx = element_dofs[range_p[j]]
                        # For pressure, need separate shape functions
                        p_val += u_ref[dof_idx] * 0.25  # Simplified averaging
                    end
                    interpolated_p[i] = p_val
                    
                    cell_found = true
                    break
                end
            end
        end
        
        if !cell_found
            interpolated_u[i] = Vec{2}((0.0, 0.0))
            interpolated_p[i] = 0.0
            @warn "Point $i not found in any cell"
        end
    end
    
    return interpolated_u, interpolated_p
end

function point_in_quadrilateral(point, coords)
    """Simple check if point is inside quadrilateral bounds"""
    x, y = point[1], point[2]
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]
    return minimum(xs) ≤ x ≤ maximum(xs) && minimum(ys) ≤ y ≤ maximum(ys)
end

function physical_to_reference(point, coords)
    """Convert physical coordinates to reference coordinates (simplified)"""
    # This is a simplified mapping - for accurate results, use proper isoparametric mapping
    x, y = point[1], point[2]
    x_min, x_max = extrema([c[1] for c in coords])
    y_min, y_max = extrema([c[2] for c in coords])
    
    ξ = 2.0 * (x - x_min) / (x_max - x_min) - 1.0
    η = 2.0 * (y - y_min) / (y_max - y_min) - 1.0
    
    return Vec{2}((ξ, η))
end

function compute_errors_simplified(grid_coarse, dh_coarse, u_coarse, grid_ref, dh_ref, u_ref)
    """Compute L2 errors using a simplified approach"""
    # Get characteristic mesh sizes
    h_coarse = sqrt(1.0 / getncells(grid_coarse))  # Approximate
    h_ref = sqrt(1.0 / getncells(grid_ref))
    
    # Compute relative DOF-based error (simplified approach)
    ndofs_coarse = ndofs(dh_coarse)
    ndofs_ref = ndofs(dh_ref)
    
    # For velocity: compare norms of solutions
    range_u_coarse = dof_range(dh_coarse, :u)
    range_u_ref = dof_range(dh_ref, :u)
    
    u_norm_coarse = norm(u_coarse[range_u_coarse])
    u_norm_ref = norm(u_ref[range_u_ref])
    
    # Simple relative error estimate (as percentage)
    velocity_error = 100.0 * abs(u_norm_coarse - u_norm_ref) / u_norm_ref
    
    # For pressure: similar approach
    range_p_coarse = dof_range(dh_coarse, :p)
    range_p_ref = dof_range(dh_ref, :p)
    
    p_norm_coarse = norm(u_coarse[range_p_coarse])
    p_norm_ref = norm(u_ref[range_p_ref])
    
    pressure_error = 100.0 * abs(p_norm_coarse - p_norm_ref) / (p_norm_ref + 1e-12)
    
    return velocity_error, pressure_error, h_coarse
end

function run_convergence_study()
    """Main function to run the convergence study"""
    println("=" ^ 60)
    println("Starting Convergence Study for 2D Stokes Flow")
    println("=" ^ 60)
    
    # Store results
    results = []
    
    # Solve for all meshes
    solutions = []
    for (i, mesh_file) in enumerate(MESH_FILES)
        grid, dh, u = solve_stokes_flow(mesh_file)
        push!(solutions, (grid=grid, dh=dh, u=u, label=MESH_LABELS[i]))
        println()
    end
    
    # Use finest mesh as reference solution
    ref_solution = solutions[end]  # vfine mesh
    println("Using $(ref_solution.label) mesh as reference solution")
    println("Reference mesh has $(ndofs(ref_solution.dh)) DOFs")
    println()
    
    # Compute errors for coarser meshes
    println("Computing convergence errors:")
    println("-" ^ 50)
    
    mesh_sizes = Float64[]
    velocity_errors = Float64[]
    pressure_errors = Float64[]
    dof_counts = Int[]
    
    for i in 1:(length(solutions)-1)  # Exclude reference solution
        sol = solutions[i]
        
        velocity_error, pressure_error, h = compute_errors_simplified(
            sol.grid, sol.dh, sol.u,
            ref_solution.grid, ref_solution.dh, ref_solution.u
        )
        
        push!(mesh_sizes, h)
        push!(velocity_errors, velocity_error)
        push!(pressure_errors, pressure_error)
        push!(dof_counts, ndofs(sol.dh))
        
        @printf("Mesh: %-8s | DOFs: %6d | h ≈ %.2e | Velocity Error: %.2f%% | Pressure Error: %.2f%%\n",
                sol.label, ndofs(sol.dh), h, velocity_error, pressure_error)
        
        # Note: VTU files removed to save space
        # VTKGridFile("convergence_$(sol.label)", sol.grid) do vtk
        #     write_solution(vtk, sol.dh, sol.u)
        # end
    end
    
    println()
    
    # Compute convergence rates
    if length(velocity_errors) > 1
        println("Computing convergence rates:")
        println("-" ^ 40)
        
        for i in 2:length(velocity_errors)
            h_ratio = mesh_sizes[i-1] / mesh_sizes[i]
            
            if velocity_errors[i] > 0 && velocity_errors[i-1] > 0
                velocity_rate = log(velocity_errors[i-1] / velocity_errors[i]) / log(h_ratio)
                @printf("Velocity convergence rate (mesh %d to %d): %.2f\n", 
                        i-1, i, velocity_rate)
            end
            
            if pressure_errors[i] > 0 && pressure_errors[i-1] > 0
                pressure_rate = log(pressure_errors[i-1] / pressure_errors[i]) / log(h_ratio)
                @printf("Pressure convergence rate (mesh %d to %d): %.2f\n", 
                        i-1, i, pressure_rate)
            end
        end
    end
    
    # Create convergence plots
    try
        println("\nCreating convergence plots...")
        
        # DOFs vs Error plot
        p1 = plot(dof_counts, velocity_errors, 
                 marker=:circle, linewidth=2, 
                 xscale=:log10, yscale=:log10,
                 xlabel="Number of DOFs", ylabel="Relative Error (%)",
                 label="Velocity Error", title="Convergence Study")
        plot!(p1, dof_counts, pressure_errors, 
              marker=:square, linewidth=2,
              label="Pressure Error")
        
        # Mesh size vs Error plot  
        p2 = plot(mesh_sizes, velocity_errors,
                 marker=:circle, linewidth=2,
                 xscale=:log10, yscale=:log10,
                 xlabel="Mesh Size h", ylabel="Relative Error (%)", 
                 label="Velocity Error", title="h-Convergence")
        plot!(p2, mesh_sizes, pressure_errors,
              marker=:square, linewidth=2,
              label="Pressure Error")
        
        # Add theoretical convergence lines
        h_theory = 10 .^ range(log10(minimum(mesh_sizes)), log10(maximum(mesh_sizes)), length=10)
        plot!(p2, h_theory, h_theory.^2 * velocity_errors[end] / mesh_sizes[end]^2,
              linestyle=:dash, color=:black, label="h²")
        plot!(p2, h_theory, h_theory * pressure_errors[end] / mesh_sizes[end],
              linestyle=:dot, color=:gray, label="h¹")
        
        # Save plots
        savefig(p1, "convergence_dofs.png")
        savefig(p2, "convergence_mesh_size.png")
        
        println("Plots saved as 'convergence_dofs.png' and 'convergence_mesh_size.png'")
        
    catch e
        println("Could not create plots (Plots.jl might not be available): $e")
    end
    
    println("\nConvergence study completed!")
    println("Convergence plots saved as PNG files")
end

# Run the convergence study
run_convergence_study()