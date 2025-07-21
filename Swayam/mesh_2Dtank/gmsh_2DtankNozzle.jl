# ============================================================================
# 2D Tank-Nozzle Generation using GMSH with Boundary Conditions
# ============================================================================
# Creates a 2D tank with nozzle geometry with proper boundary condition labels for Ferrite.jl
# Based on the 3D tank design but simplified to 2D for easier development and testing

println("Starting 2D Tank-Nozzle Generation with Boundary Conditions...")

# Import GMSH
try
    using Gmsh: gmsh
    println("✅ GMSH imported successfully!")
catch
    using gmsh
    println("✅ GMSH imported successfully!")
end

function create_2d_tank_nozzle_with_bcs()
    println("Starting 2D Tank-Nozzle Generation with Boundary Conditions...")
    
    # Geometry parameters (same as 3D version)
    tank_length = 1.0
    tank_radius = 0.1
    nozzle_length = 0.15
    nozzle_radius = 0.03
    mesh_size_fine = 0.008
    mesh_size_coarse = 0.025
    
    println("\n📐 Creating 2D Tank-Nozzle geometry:")
    println("   Tank Length: $(tank_length)")
    println("   Tank Radius: $(tank_radius)")
    println("   Nozzle Length: $(nozzle_length)")
    println("   Nozzle Radius: $(nozzle_radius)")
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("2D_tank_nozzle_with_BCs")
    
    # ========================================================================
    # STEP 1: Create 2D Tank-Nozzle Profile with Structured Mesh
    # ========================================================================
    println("\n🔧 Step 1: Creating structured 2D tank-nozzle profile...")
    
    # Geometry parameters (using naming convention from 3D version)
    L = tank_length
    R = tank_radius
    Lin = -nozzle_length
    Rin = nozzle_radius
    lc1 = mesh_size_fine
    lc2 = mesh_size_coarse
    
    # Define 10 points for structured mesh (same as 3D version but in 2D)
    p1  = gmsh.model.geo.addPoint(0,   -R,   0, lc1, 1)   # Bottom-left tank
    p2  = gmsh.model.geo.addPoint(L,   -R,   0, lc2, 2)   # Bottom-right tank
    p3  = gmsh.model.geo.addPoint(L,   -Rin, 0, lc2, 3)   # Right-bottom nozzle junction
    p4  = gmsh.model.geo.addPoint(L,    Rin, 0, lc2, 4)   # Right-top nozzle junction
    p5  = gmsh.model.geo.addPoint(L,    R,   0, lc2, 5)   # Top-right tank
    p6  = gmsh.model.geo.addPoint(0,    R,   0, lc1, 6)   # Top-left tank
    p7  = gmsh.model.geo.addPoint(0,    Rin, 0, lc1, 7)   # Left-top nozzle junction
    p8  = gmsh.model.geo.addPoint(Lin,  Rin, 0, lc1, 8)   # Top nozzle end (inlet)
    p9  = gmsh.model.geo.addPoint(Lin, -Rin, 0, lc1, 9)   # Bottom nozzle end (inlet)
    p10 = gmsh.model.geo.addPoint(0,   -Rin, 0, lc1, 10)  # Left-bottom nozzle junction
    
    # Define 13 lines (same connectivity as 3D version)
    l1 = gmsh.model.geo.addLine(1, 2, 1)   # Bottom tank wall
    l2 = gmsh.model.geo.addLine(2, 3, 2)   # Right tank wall (bottom part)
    l3 = gmsh.model.geo.addLine(3, 4, 3)   # Right tank wall (nozzle section)
    l4 = gmsh.model.geo.addLine(4, 5, 4)   # Right tank wall (top part)
    l5 = gmsh.model.geo.addLine(5, 6, 5)   # Top tank wall
    l6 = gmsh.model.geo.addLine(6, 7, 6)   # Left tank wall (top part)
    l7 = gmsh.model.geo.addLine(7, 8, 7)   # Top nozzle wall
    l8 = gmsh.model.geo.addLine(8, 9, 8)   # Inlet opening (boundary line)
    l9 = gmsh.model.geo.addLine(9, 10, 9)  # Bottom nozzle wall
    l10 = gmsh.model.geo.addLine(10, 1, 10) # Left tank wall (bottom part)
    l11 = gmsh.model.geo.addLine(10, 7, 11) # Internal line (vertical)
    l12 = gmsh.model.geo.addLine(3, 10, 12) # Internal line (horizontal bottom)
    l13 = gmsh.model.geo.addLine(7, 4, 13)  # Internal line (horizontal top)
    
    # Define 4 structured surfaces (quadrilateral regions)
    loop1 = gmsh.model.geo.addCurveLoop([1, 2, 12, 10], 1)      # Bottom-right region
    loop2 = gmsh.model.geo.addCurveLoop([-12, 3, -13, -11], 2)  # Center region
    loop3 = gmsh.model.geo.addCurveLoop([13, 4, 5, 6], 3)       # Top-right region
    loop4 = gmsh.model.geo.addCurveLoop([9, 11, 7, 8], 4)       # Nozzle region
    
    surf1 = gmsh.model.geo.addPlaneSurface([1], 1)
    surf2 = gmsh.model.geo.addPlaneSurface([2], 2)
    surf3 = gmsh.model.geo.addPlaneSurface([3], 3)
    surf4 = gmsh.model.geo.addPlaneSurface([4], 4)
    
    # Set transfinite curves for structured mesh
    # We need to ensure opposite sides of each surface have matching node counts
    
    # Surface 1: loop1 = [1, 2, 12, 10] (Bottom-right region)
    # Opposite pairs: (1,12) and (2,10)
    gmsh.model.geo.mesh.setTransfiniteCurve(1, 15)   # Bottom tank
    gmsh.model.geo.mesh.setTransfiniteCurve(12, 15)  # Internal horizontal bottom (match with 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, 8)    # Right tank bottom
    gmsh.model.geo.mesh.setTransfiniteCurve(10, 8)   # Left tank bottom (match with 2)
    
    # Surface 2: loop2 = [-12, 3, -13, -11] (Center region)
    # Opposite pairs: (12,13) and (3,11)
    # 12 already set to 15, so 13 must be 15
    # 3 and 11 must match
    gmsh.model.geo.mesh.setTransfiniteCurve(13, 15)  # Internal horizontal top (match with 12)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, 8)    # Right tank nozzle
    gmsh.model.geo.mesh.setTransfiniteCurve(11, 8)   # Internal vertical (match with 3)
    
    # Surface 3: loop3 = [13, 4, 5, 6] (Top-right region)
    # Opposite pairs: (13,5) and (4,6)
    # 13 already set to 15, so 5 must be 15
    # 4 and 6 must match
    gmsh.model.geo.mesh.setTransfiniteCurve(5, 15)   # Top tank (match with 13)
    gmsh.model.geo.mesh.setTransfiniteCurve(4, 8)    # Right tank top
    gmsh.model.geo.mesh.setTransfiniteCurve(6, 8)    # Left tank top (match with 4)
    
    # Surface 4: loop4 = [9, 11, 7, 8] (Nozzle region)
    # Opposite pairs: (9,7) and (11,8)
    # 11 already set to 8, so 8 must be 8
    # 9 and 7 must match
    gmsh.model.geo.mesh.setTransfiniteCurve(9, 10)   # Bottom nozzle
    gmsh.model.geo.mesh.setTransfiniteCurve(7, 10)   # Top nozzle (match with 9)
    gmsh.model.geo.mesh.setTransfiniteCurve(8, 8)    # Inlet opening (match with 11)
    
    # Set transfinite surfaces for structured quadrilateral mesh
    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setTransfiniteSurface(2)
    gmsh.model.geo.mesh.setTransfiniteSurface(3)
    gmsh.model.geo.mesh.setTransfiniteSurface(4)
    
    # Set recombine for all surfaces to get quadrilaterals instead of triangles
    gmsh.model.geo.mesh.setRecombine(2, 1)
    gmsh.model.geo.mesh.setRecombine(2, 2)
    gmsh.model.geo.mesh.setRecombine(2, 3)
    gmsh.model.geo.mesh.setRecombine(2, 4)
    
    gmsh.model.geo.synchronize()

    # ========================================================================
    # STEP 2: Define Physical Groups for Boundary Conditions
    # ========================================================================
    println("🔧 Step 2: Defining physical groups for boundary conditions...")
    
    # Classify boundary lines:
    # - Wall lines: all solid boundaries (tank walls, nozzle walls)
    # - Inlet line: the nozzle opening where fluid enters
    
    wall_lines = [1, 2, 3, 4, 5, 6, 7, 9, 10]  # All walls except inlet
    inlet_lines = [8]  # Inlet opening at nozzle end
    
    # Create physical groups for boundary conditions
    gmsh.model.addPhysicalGroup(1, wall_lines, -1, "walls")
    gmsh.model.addPhysicalGroup(1, inlet_lines, -1, "inlet")
    
    # Create physical group for the fluid domain (all surfaces)
    gmsh.model.addPhysicalGroup(2, [1, 2, 3, 4], -1, "fluid_domain")
    
    println("   ✅ Created physical groups:")
    println("      - 'walls': $(length(wall_lines)) boundary lines (no-slip/wall BC)")
    println("      - 'inlet': $(length(inlet_lines)) boundary line (inflow BC)")
    println("      - 'fluid_domain': 4 surfaces (fluid domain)")
    
    # ========================================================================
    # STEP 3: Generate 2D Structured Mesh
    # ========================================================================
    println("🔧 Step 3: Generating structured 2D mesh...")
    
    # Set mesh options for structured quadrilateral mesh
    gmsh.option.setNumber("Mesh.RecombineAll", 1)     # Enable recombination globally
    gmsh.option.setNumber("Mesh.Algorithm", 8)         # Frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.ElementOrder", 1)      # Linear elements
    
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    
    # Check element types
    element_types = gmsh.model.mesh.getElementTypes()
    println("   📊 Element types found: $(element_types)")
    
    has_quads = false
    for elem_type in element_types
        elem_name = gmsh.model.mesh.getElementProperties(elem_type)[1]
        num_elements = length(gmsh.model.mesh.getElementsByType(elem_type)[1])
        println("   - Type $(elem_type) ($(elem_name)): $(num_elements) elements")
        
        # Check if we have quadrilateral elements (type 3 = 4-node quadrilateral)
        if elem_type == 3
            has_quads = true
            println("   ✅ Found quadrilateral elements!")
        end
    end
    
    if !has_quads
        println("   ⚠️  Warning: No quadrilateral elements found. You may have triangular elements.")
    end
    
    # ========================================================================
    # STEP 4: Save and Visualize
    # ========================================================================
    println("🔧 Step 4: Saving mesh and displaying results...")
    
    # Save the mesh
    gmsh.write("tank_2D_with_BCs.msh")
    println("   ✅ 2D tank-nozzle mesh saved: tank_2D_with_BCs.msh")
    
    # Show mesh info
    mesh_info = gmsh.model.mesh.getNodes()
    num_nodes = length(mesh_info[1])
    println("   📊 Total nodes: $(num_nodes)")
    
    # Print summary of boundary conditions for Ferrite.jl
    println("\n🏷️  Boundary Condition Labels for Ferrite.jl:")
    println("   - 'walls': Tank and nozzle wall boundaries (no-slip/wall BC)")
    println("   - 'inlet': Nozzle opening boundary (inflow BC)")
    println("   - 'fluid_domain': Interior surfaces (fluid domain)")
    
    # Print physical group information
    println("\n📋 Physical Group Details:")
    physical_groups = gmsh.model.getPhysicalGroups()
    for (dim, tag) in physical_groups
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("   - Dimension $(dim), Tag $(tag), Name '$(name)': $(length(entities)) entities")
    end
    
    # Print mesh verification
    println("\n🔍 Mesh Verification:")
    verify_2d_boundary_conditions(wall_lines, inlet_lines, Lin)
    
    # Show the mesh
    println("   🖥️  Showing 2D tank-nozzle mesh with boundary conditions...")
    gmsh.fltk.run()
    
    gmsh.finalize()
    println("\n✅ Done! 2D tank-nozzle mesh ready for Ferrite.jl")
end

"""
    verify_2d_boundary_conditions(wall_lines, inlet_lines, Lin)

Verify that the 2D boundary conditions are correctly set up.
"""
function verify_2d_boundary_conditions(wall_lines, inlet_lines, Lin)
    println("   Checking boundary condition setup...")
    
    # Check wall lines
    println("   Wall lines ($(length(wall_lines))): $(wall_lines)")
    for line_tag in wall_lines
        bbox = gmsh.model.getBoundingBox(1, line_tag)
        x_min, y_min = bbox[1], bbox[2]
        x_max, y_max = bbox[4], bbox[5]
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        println("      Line $(line_tag): center at ($(round(x_center,digits=4)), $(round(y_center,digits=4)))")
    end
    
    # Check inlet lines
    println("   Inlet lines ($(length(inlet_lines))): $(inlet_lines)")
    for line_tag in inlet_lines
        bbox = gmsh.model.getBoundingBox(1, line_tag)
        x_min, y_min = bbox[1], bbox[2]
        x_max, y_max = bbox[4], bbox[5]
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        expected_x = Lin
        
        println("      Line $(line_tag): center at ($(round(x_center,digits=4)), $(round(y_center,digits=4)))")
        
        if abs(x_center - expected_x) < 0.01
            println("      ✅ Inlet at correct x position: $(round(x_center,digits=4)) ≈ $(expected_x)")
        else
            println("      ⚠️  Inlet at wrong x position: $(round(x_center,digits=4)) ≠ $(expected_x)")
        end
    end
    
    total_boundary_lines = length(wall_lines) + length(inlet_lines)
    println("   ✅ Total boundary lines classified: $(total_boundary_lines)")
    
    return true
end

# Run the function
create_2d_tank_nozzle_with_bcs()
