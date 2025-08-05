using Ferrite, FerriteGmsh, Gmsh, Tensors, LinearAlgebra, SparseArrays, WriteVTK

function setup_grid(h = 0.05)
    # Initialize Gmsh
    Gmsh.initialize()
    gmsh.option.set_number("General.Terminal", 1)
    gmsh.option.set_number("General.Verbosity", 2)

    # Geometry parameters
    tank_length = 1.0
    tank_radius = 0.1
    nozzle_length = 0.15
    nozzle_radius = 0.03
    wedge_angle_deg = 45.0
    mesh_size_fine = 0.008
    mesh_size_coarse = 0.025
    
    println("\n📐 Creating 3D Wedge Tank:")
    println("   Tank Length: $(tank_length)")
    println("   Tank Radius: $(tank_radius)")
    println("   Wedge Angle: $(wedge_angle_deg)°")
    
    wedge_angle_rad = deg2rad(wedge_angle_deg)
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("3D_wedge_tank")
    
    # ========================================================================
    # STEP 1: Create 2D Profile with Structured Mesh
    # ========================================================================
    println("\n🔧 Step 1: Creating structured 2D tank profile...")
    
    # Geometry parameters (using your naming convention)
    L = tank_length
    R = tank_radius
    Lin = -nozzle_length
    Rin = nozzle_radius
    lc1 = mesh_size_fine
    lc2 = mesh_size_coarse
    
    # Define 10 points for structured mesh
    p1  = gmsh.model.geo.addPoint(0,   -R,   0, lc1, 1)
    p2  = gmsh.model.geo.addPoint(L,   -R,   0, lc2, 2)
    p3  = gmsh.model.geo.addPoint(L,   -Rin, 0, lc2, 3)
    p4  = gmsh.model.geo.addPoint(L,    Rin, 0, lc2, 4)
    p5  = gmsh.model.geo.addPoint(L,    R,   0, lc2, 5)
    p6  = gmsh.model.geo.addPoint(0,    R,   0, lc1, 6)
    p7  = gmsh.model.geo.addPoint(0,    Rin, 0, lc1, 7)
    p8  = gmsh.model.geo.addPoint(Lin,  Rin, 0, lc1, 8)
    p9  = gmsh.model.geo.addPoint(Lin, -Rin, 0, lc1, 9)
    p10 = gmsh.model.geo.addPoint(0,   -Rin, 0, lc1, 10)

    p11 = gmsh.model.geo.addPoint(0,   0,   0, lc1, 11)
    p12 = gmsh.model.geo.addPoint(L,   0,   0, lc2, 12)
    p13 = gmsh.model.geo.addPoint(Lin, 0,   0, lc1, 13)
    
    # Define 13 lines
    l1 = gmsh.model.geo.addLine(1, 2, 1)
    l2 = gmsh.model.geo.addLine(2, 3, 2)
    l3 = gmsh.model.geo.addLine(3, 12, 3)
    l4 = gmsh.model.geo.addLine(4, 5, 4)
    l5 = gmsh.model.geo.addLine(5, 6, 5)
    l6 = gmsh.model.geo.addLine(6, 7, 6)
    l7 = gmsh.model.geo.addLine(7, 8, 7)
    l8 = gmsh.model.geo.addLine(8, 13, 8)
    l9 = gmsh.model.geo.addLine(9, 10, 9)
    l10 = gmsh.model.geo.addLine(10, 1, 10)
    l11 = gmsh.model.geo.addLine(10, 11, 11)
    l12 = gmsh.model.geo.addLine(3, 10, 12)
    l13 = gmsh.model.geo.addLine(7, 4, 13)

    l14 = gmsh.model.geo.addLine(13, 9, 14)  
    l15 = gmsh.model.geo.addLine(13, 11, 15)
    l16 = gmsh.model.geo.addLine(11, 7, 16)
    l17 = gmsh.model.geo.addLine(11, 12, 17)
    l18 = gmsh.model.geo.addLine(12, 4, 18)
    
    # Define 4 structured surfaces
    loop1 = gmsh.model.geo.addCurveLoop([1, 2, 12, 10], 1)
    loop2 = gmsh.model.geo.addCurveLoop([-12, 3, -17, -11], 2)
    loop3 = gmsh.model.geo.addCurveLoop([13, 4, 5, 6], 3)
    loop4 = gmsh.model.geo.addCurveLoop([9, 11, -15, 14], 4)

    loop5 = gmsh.model.geo.addCurveLoop([17, 18, -13, -16], 5)
    loop6 = gmsh.model.geo.addCurveLoop([15, 16, 7, 8], 6)
    
    surf1 = gmsh.model.geo.addPlaneSurface([1], 1)
    surf2 = gmsh.model.geo.addPlaneSurface([2], 2)
    surf3 = gmsh.model.geo.addPlaneSurface([3], 3)
    surf4 = gmsh.model.geo.addPlaneSurface([4], 4)

    surf5 = gmsh.model.geo.addPlaneSurface([5], 5)
    surf6 = gmsh.model.geo.addPlaneSurface([6], 6)
    
    # Set transfinite curves for structured mesh
    gmsh.model.geo.mesh.setTransfiniteCurve(1, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(12, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, 5)
    gmsh.model.geo.mesh.setTransfiniteCurve(10, 5)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(12, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(13, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, 5)
    gmsh.model.geo.mesh.setTransfiniteCurve(11, 5)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(13, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(5, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(4, 5)
    gmsh.model.geo.mesh.setTransfiniteCurve(6, 5)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(9, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(7, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(11, 5)
    gmsh.model.geo.mesh.setTransfiniteCurve(8, 5)

    gmsh.model.geo.mesh.setTransfiniteCurve(17, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(15, 10)
    gmsh.model.geo.mesh.setTransfiniteCurve(14, 5)
    gmsh.model.geo.mesh.setTransfiniteCurve(16, 5)
    gmsh.model.geo.mesh.setTransfiniteCurve(18, 5)
    
    # Set transfinite surfaces
    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setTransfiniteSurface(2)
    gmsh.model.geo.mesh.setTransfiniteSurface(3)
    gmsh.model.geo.mesh.setTransfiniteSurface(4)

    gmsh.model.geo.mesh.setTransfiniteSurface(5)
    gmsh.model.geo.mesh.setTransfiniteSurface(6)
    
    gmsh.model.geo.synchronize()
    
    # ========================================================================
    # STEP 2: Generate structured 2D mesh
    # ========================================================================
    println("🔧 Step 2: Generating structured 2D mesh...")
    
    # Set recombine for all surfaces to get quadrilaterals
    gmsh.model.mesh.setRecombine(2, 1)
    gmsh.model.mesh.setRecombine(2, 2)
    gmsh.model.mesh.setRecombine(2, 3)
    gmsh.model.mesh.setRecombine(2, 4)
    gmsh.model.mesh.setRecombine(2, 5)
    gmsh.model.mesh.setRecombine(2, 6)

    gmsh.fltk.run()

    
    # ========================================================================
    # STEP 3: Revolve to create 3D wedge with tetrahedra
    # ========================================================================
    println("🔧 Step 3: Revolving structured surfaces to create 3D wedge...")

    num_layers = 8
    println("   Revolving all surfaces with $(num_layers) layers...")

    # ⚠️ Disable recombination options
    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    gmsh.option.setNumber("Mesh.Recombine3DAll", 0)

    # Optional: Use unstructured algorithms for triangles and tets
    gmsh.option.setNumber("Mesh.Algorithm", 6)    # Delaunay triangles
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay tetrahedra

    # Do NOT recombine
    revolve_result = gmsh.model.geo.revolve(
        [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6)],
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        wedge_angle_rad,
        [num_layers],
        [],
        true  # <<<< This disables recombination and avoids hexahedra/prisms
    )

    gmsh.model.geo.synchronize()

   
    inlet_tag = 100

    # Define all revolved surface tags
    revolved_surfaces = [tag for (dim, tag) in revolve_result if dim == 2]

    # Wall surface tags = all revolved surfaces except inlet and internal ones
    other_tags = Set([52, 30, 56, 79, 35, 101, 57])  # Internal or periodic surfaces to exclude
    wall_tags = [tag for tag in revolved_surfaces if tag != inlet_tag && !(tag in other_tags)]

    # Add physical groups
    gmsh.model.addPhysicalGroup(2, [inlet_tag], -1, "inlet")
    gmsh.model.addPhysicalGroup(2, wall_tags, -1, "wall")

    gmsh.model.addPhysicalGroup(2, [1], -1, "Γ1")
    gmsh.model.addPhysicalGroup(2, [35], -1, "Γ35")
    gmsh.model.addPhysicalGroup(2, [2], -1, "Γ2")
    gmsh.model.addPhysicalGroup(2, [57], -1, "Γ57")
    gmsh.model.addPhysicalGroup(2, [3], -1, "Γ3")
    gmsh.model.addPhysicalGroup(2, [101], -1, "Γ101")
    gmsh.model.addPhysicalGroup(2, [4], -1, "Γ4")
    gmsh.model.addPhysicalGroup(2, [79], -1, "Γ79")

    # ========================================================================
    # STEP 4: Mesh generation
    # ========================================================================
    println("🔧 Step 4: Generating tetrahedral 3D mesh...")

    # Set general mesh options for tetrahedra
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)

    

    # Run meshing
    gmsh.model.mesh.generate(3)
    println("   ✅ 3D tetrahedral mesh generated successfully")

    # Print element types to confirm
    element_types = gmsh.model.mesh.getElementTypes()
    println("   📊 Element types found: $(element_types)")
    for elem_type in element_types
        elem_name = gmsh.model.mesh.getElementProperties(elem_type)[1]
        num_elements = length(gmsh.model.mesh.getElementsByType(elem_type)[1])
        println("   - Type $(elem_type) ($(elem_name)): $(num_elements) elements")
    end


    # Show 3D mesh in Gmsh
    println("   🖥️  Showing 3D mesh...")
    gmsh.fltk.run()

    # Save in MSH format for Ferrite
    gmsh.write("Simranjeet/paraview/tank_3D_wedge_hexahedra.msh")


    # Convert to Ferrite grid
    grid = FerriteGmsh.togrid()

    # Save VTK via Ferrite for proper 3D visualization
    VTKGridFile("Simranjeet/paraview/tank_3D_wedge_hexahedra", grid) do vtk
        # Add any cell/node data here if needed
    end

    gmsh.finalize()
    println("\n✅ Done!")
    return grid
end


function setup_fevalues(ipu, ipp, ipg)
    qr = QuadratureRule{RefTetrahedron}(2)
    cvu = CellValues(qr, ipu, ipg)
    cvp = CellValues(qr, ipp, ipg)
    qr_facet = FacetQuadratureRule{RefTetrahedron}(2)
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

    # Combine all boundary facets
    set = union(
        getfacetset(dh.grid, "inlet"),
        getfacetset(dh.grid, "wall")
    )

    # DoF and buffer setup
    range_p = dof_range(dh, :p)
    element_dofs = zeros(Int, ndofs_per_cell(dh))
    element_dofs_p = view(element_dofs, range_p)
    element_coords = zeros(Vec{3}, 4)  # Triangular face → 3 nodes

    Ce = zeros(1, length(range_p))  # Local constraint vector

    for (cell_id, face_id) in set
        Ce .= 0
        getcoordinates!(element_coords, dh.grid, cell_id)
        reinit!(fvp, element_coords, face_id)
        celldofs!(element_dofs, dh, cell_id)

        for qp in 1:getnquadpoints(fvp)
            dΓ = getdetJdV(fvp, qp)
            for i in 1:getnbasefunctions(fvp)
                Ce[1, i] += shape_value(fvp, qp, i) * dΓ
            end
        end

        assemble!(assembler, [1], element_dofs_p, Ce)
    end

    C, _ = finish_assemble(assembler)

    # Normalize constraint to set one pressure DoF as reference
    _, J, V = findnz(C)
    _, constrained_idx = findmax(abs2, V)
    constrained_dof = J[constrained_idx]
    V ./= V[constrained_idx]

    return AffineConstraint(
        constrained_dof,
        [J[i] => -V[i] for i in eachindex(J) if J[i] != constrained_dof],
        0.0,
    )
end


function setup_constraints(dh, fvp)
    ch = ConstraintHandler(dh)

    # Inlet: Dirichlet BC for velocity (along x-direction)
    Γin = getfacetset(dh.grid, "inlet")
    inlet_bc = Dirichlet(:u, Γin, (x, t) -> Vec{3, Float64}((1.0, 0.0, 0.0)))
    add!(ch, inlet_bc)

    # Walls: No-slip (zero velocity)
    Γwall = getfacetset(dh.grid, "wall") 
    wall_bc = Dirichlet(:u, Γwall, (x, t) -> Vec{3, Float64}((0.0, 0.0, 0.0)))
    add!(ch, wall_bc)

    # Periodic BCs
    θ = π / 4  # 45 degrees
    R = Tensor{2,3,Float64}([
        1.0  0.0        0.0;
        0.0  cos(θ)   -sin(θ);
        0.0  sin(θ)    cos(θ)
    ])

    add!(ch, PeriodicDirichlet(
        :u,
        collect_periodic_facets(dh.grid, "Γ35", "Γ1", x -> R ⋅ x),
        R,
        [1, 2, 3]
    ))
    add!(ch, PeriodicDirichlet(
        :u,
        collect_periodic_facets(dh.grid, "Γ57", "Γ2", x -> R ⋅ x),
        R,
        [1, 2, 3]
    ))
    add!(ch, PeriodicDirichlet(
        :u,
        collect_periodic_facets(dh.grid, "Γ79", "Γ3", x -> R ⋅ x),
        R,
        [1, 2, 3]
    ))
    add!(ch, PeriodicDirichlet(
        :u,
        collect_periodic_facets(dh.grid, "Γ101", "Γ4", x -> R ⋅ x),
        R,
        [1, 2, 3]
    ))

    # Mean pressure constraint
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
    range_p = dof_range(dh, :p)
    ndofs_u = length(range_u)
    ndofs_p = length(range_p)

    # For 3D: velocity shape functions are Vec{3}, gradients are 3x3 tensors
    ϕᵤ = Vector{Vec{3, Float64}}(undef, ndofs_u)
    ∇ϕᵤ = Vector{Tensor{2, 3, Float64, 9}}(undef, ndofs_u)
    divϕᵤ = Vector{Float64}(undef, ndofs_u)
    ϕₚ = Vector{Float64}(undef, ndofs_p)

    for cell in CellIterator(dh)
        reinit!(cvu, cell)
        reinit!(cvp, cell)
        ke .= 0
        fe .= 0

        for qp in 1:getnquadpoints(cvu)
            dΩ = getdetJdV(cvu, qp)

            # Collect basis function values
            for i in 1:ndofs_u
                ϕᵤ[i] = shape_value(cvu, qp, i)
                ∇ϕᵤ[i] = shape_gradient(cvu, qp, i)
                divϕᵤ[i] = shape_divergence(cvu, qp, i)
            end
            for i in 1:ndofs_p
                ϕₚ[i] = shape_value(cvp, qp, i)
            end

            # Assemble local stiffness matrix
            # velocity-velocity (∇u : ∇v)
            for (i, I) in pairs(range_u), (j, J) in pairs(range_u)
                ke[I, J] += (∇ϕᵤ[i] ⊡ ∇ϕᵤ[j]) * dΩ
            end

            # velocity-pressure (∇⋅u * p)
            for (i, I) in pairs(range_u), (j, J) in pairs(range_p)
                ke[I, J] += (-divϕᵤ[i] * ϕₚ[j]) * dΩ
            end

            # pressure-velocity (q * ∇⋅v)
            for (i, I) in pairs(range_p), (j, J) in pairs(range_u)
                ke[I, J] += (-divϕᵤ[j] * ϕₚ[i]) * dΩ
            end

            # Right-hand side (source term): simple localized Gaussian bump
            x = spatial_coordinate(cvu, qp, getcoordinates(cell))
            b = exp(-100 * norm(x - Vec{3}((0.75, 0.0, 0.0)))^2)  # e.g. in x-direction
            bv = Vec{3, Float64}((b, 0.0, 0.0))

            for (i, I) in pairs(range_u)
                fe[I] += (ϕᵤ[i] ⋅ bv) * dΩ
            end
        end

        assemble!(assembler, celldofs(cell), ke, fe)
    end

    return K, f
end


function main()
    # Grid
    println("Setting up grid...")
    h = 0.05 # approximate element size
    grid = setup_grid()
    # Interpolations
    ipu = Lagrange{RefTetrahedron, 2}()^3 # quadratic
    ipp = Lagrange{RefTetrahedron, 1}()   # linear
    # Dofs
    dh = setup_dofs(grid, ipu, ipp)
    # FE values
    ipg = Lagrange{RefTetrahedron, 1}() # linear geometric interpolation
    cvu, cvp, fvp = setup_fevalues(ipu, ipp, ipg)
    # Boundary conditions
    ch = setup_constraints(dh, fvp)
    # Global tangent matrix and rhs
    coupling = [true true; true false] # no coupling between pressure test/trial functions
    K = allocate_matrix(dh, ch; coupling = coupling)
    f = zeros(ndofs(dh))
    # Assemble system
    assemble_system!(K, f, dh, cvu, cvp)
    # Apply boundary conditions and solve
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)
    # Export the solution
    VTKGridFile("Simranjeet/paraview/3D_stokes_flow_tank_hexahedra", grid) do vtk
        write_solution(vtk, dh, u)
        Ferrite.write_constraints(vtk, ch)
    end


    return
end

main()