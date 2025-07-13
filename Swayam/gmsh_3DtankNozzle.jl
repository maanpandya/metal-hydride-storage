# ============================================================================
# 3D Wedge Tank Generation using GMSH Revolution
# ============================================================================
# Creates a 3D wedge tank by revolving a 2D rectangular profile around the x-axis

println("Starting 3D Wedge Tank Generation...")

# Import GMSH
try
    using Gmsh: gmsh
    println("✅ GMSH imported successfully!")
catch
    using gmsh
    println("✅ GMSH imported successfully!")
end

function create_3d_wedge_tank()
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
    
    # Define 13 lines
    l1 = gmsh.model.geo.addLine(1, 2, 1)
    l2 = gmsh.model.geo.addLine(2, 3, 2)
    l3 = gmsh.model.geo.addLine(3, 4, 3)
    l4 = gmsh.model.geo.addLine(4, 5, 4)
    l5 = gmsh.model.geo.addLine(5, 6, 5)
    l6 = gmsh.model.geo.addLine(6, 7, 6)
    l7 = gmsh.model.geo.addLine(7, 8, 7)
    l8 = gmsh.model.geo.addLine(8, 9, 8)
    l9 = gmsh.model.geo.addLine(9, 10, 9)
    l10 = gmsh.model.geo.addLine(10, 1, 10)
    l11 = gmsh.model.geo.addLine(10, 7, 11)
    l12 = gmsh.model.geo.addLine(3, 10, 12)
    l13 = gmsh.model.geo.addLine(7, 4, 13)
    
    # Define 4 structured surfaces
    loop1 = gmsh.model.geo.addCurveLoop([1, 2, 12, 10], 1)
    loop2 = gmsh.model.geo.addCurveLoop([-12, 3, -13, -11], 2)
    loop3 = gmsh.model.geo.addCurveLoop([13, 4, 5, 6], 3)
    loop4 = gmsh.model.geo.addCurveLoop([9, 11, 7, 8], 4)
    
    surf1 = gmsh.model.geo.addPlaneSurface([1], 1)
    surf2 = gmsh.model.geo.addPlaneSurface([2], 2)
    surf3 = gmsh.model.geo.addPlaneSurface([3], 3)
    surf4 = gmsh.model.geo.addPlaneSurface([4], 4)
    
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
    
    # Set transfinite surfaces
    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setTransfiniteSurface(2)
    gmsh.model.geo.mesh.setTransfiniteSurface(3)
    gmsh.model.geo.mesh.setTransfiniteSurface(4)
    
    gmsh.model.geo.synchronize()
    
    # Assign physical groups
    gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4, 5, 6, 7, 9, 10], -1, "wall")
    gmsh.model.addPhysicalGroup(1, [8], -1, "inlet")
    gmsh.model.addPhysicalGroup(2, [1, 2, 3, 4], -1, "omega")
    
    # ========================================================================
    # STEP 2: Generate structured 2D mesh
    # ========================================================================
    println("🔧 Step 2: Generating structured 2D mesh...")
    
    # Set recombine for all surfaces to get quadrilaterals
    gmsh.model.mesh.setRecombine(2, 1)
    gmsh.model.mesh.setRecombine(2, 2)
    gmsh.model.mesh.setRecombine(2, 3)
    gmsh.model.mesh.setRecombine(2, 4)
    
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    
    # Save and show 2D mesh
    gmsh.write("tank_2D_structured.msh")
    println("   ✅ 2D structured mesh saved: tank_2D_structured.msh")
    println("   🖥️  Showing 2D structured mesh - close window to continue...")
    gmsh.fltk.run()  # Show 2D mesh first
    
    # ========================================================================
    # STEP 3: Revolve to create 3D wedge with structured mesh
    # ========================================================================
    println("🔧 Step 3: Revolving structured surfaces to create 3D wedge...")
    
    # Calculate number of layers for revolution
    num_layers = 8  # Number of elements in circumferential direction
    
    println("   Revolving all surfaces with $(num_layers) layers...")
    
    # Method 1: Try revolving with numElements parameter for structured layers
    # This should create structured layers in the circumferential direction
    
    revolve_result = gmsh.model.geo.revolve(
        [(2, 1), (2, 2), (2, 3), (2, 4)],  # All 4 surfaces to revolve
        0.0, 0.0, 0.0,                      # Point on rotation axis
        1.0, 0.0, 0.0,                      # Rotation axis (X-axis)
        wedge_angle_rad,                    # Rotation angle
        [num_layers]                        # Number of structured layers
    )
    
    gmsh.model.geo.synchronize()
    
    # Set transfinite volumes for the created 3D regions
    # This should force structured hexahedral meshing
    volumes = gmsh.model.getEntities(3)
    println("   Found $(length(volumes)) volumes after revolution")
    
    for (dim, tag) in volumes
        try
            # Set the volume to be transfinite (structured)
            gmsh.model.geo.mesh.setTransfiniteVolume(tag)
            println("   Set transfinite volume $(tag)")
        catch e
            println("   Warning: Could not set transfinite for volume $(tag): $(e)")
        end
    end
    
    # Set recombination for volumes to ensure hexahedral elements
    for (dim, tag) in volumes
        try
            gmsh.model.mesh.setRecombine(3, tag)
            println("   Set recombination for volume $(tag)")
        catch e
            println("   Warning: Could not set recombination for volume $(tag): $(e)")
        end
    end
    
    println("   ✅ Revolution completed with structured layers")
    
    # ========================================================================
    # STEP 4: Generate 3D mesh with hexahedral elements
    # ========================================================================
    println("🔧 Step 4: Generating 3D hexahedral mesh...")
    
    # Check what entities we have after extrusion
    volumes = gmsh.model.getEntities(3)
    surfaces = gmsh.model.getEntities(2)
    println("   Found $(length(volumes)) volumes and $(length(surfaces)) surfaces after extrusion")
    
    # The extrusion with recombine=true should already create hexahedral elements
    # But let's ensure the settings are correct
    
    # Set 3D mesh algorithm for structured elements
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)     # Delaunay (works well with structured)
    gmsh.option.setNumber("Mesh.Recombine3DAll", 1)  # Enable 3D recombination
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # Subdivision algorithm
    
    # The mesh should already be created by the extrusion, but let's generate to be sure
    try
        gmsh.model.mesh.generate(3)
        println("   ✅ 3D hexahedral mesh generated successfully")
        
        # Check element types
        element_types = gmsh.model.mesh.getElementTypes()
        println("   📊 Element types found: $(element_types)")
        
        for elem_type in element_types
            elem_name = gmsh.model.mesh.getElementProperties(elem_type)[1]
            num_elements = length(gmsh.model.mesh.getElementsByType(elem_type)[1])
            println("   - Type $(elem_type) ($(elem_name)): $(num_elements) elements")
        end
        
    catch e
        println("   ❌ Error generating 3D mesh: $(e)")
    end
    
    # ========================================================================
    # STEP 5: Save and visualize
    # ========================================================================
    gmsh.write("tank_3D_wedge.msh")
    println("   ✅ 3D wedge mesh saved: tank_3D_wedge.msh")
    
    # Show mesh info
    mesh_info = gmsh.model.mesh.getNodes()
    num_nodes = length(mesh_info[1])
    println("   📊 Total nodes: $(num_nodes)")
    
    # Show 3D mesh
    println("   🖥️  Showing 3D mesh...")
    gmsh.fltk.run()
    
    gmsh.finalize()
    println("\n✅ Done!")
end

# Run the example
create_3d_wedge_tank()