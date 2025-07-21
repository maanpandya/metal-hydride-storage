# ============================================================================
# 3D Half-Tank Generation using GMSH Revolution with Boundary Conditions
# ============================================================================
# Creates a 3D half-tank (45° wedge) by revolving a 2D rectangular profile around the x-axis
# with proper boundary condition labels for Ferrite.jl

# Include verification functions
include("mesh_verification.jl")

println("Starting 3D Wedge Tank Generation with Boundary Conditions...")

# Import GMSH
try
    using Gmsh: gmsh
    println("✅ GMSH imported successfully!")
catch
    using gmsh
    println("✅ GMSH imported successfully!")
end

function create_3d_wedge_tank_with_bcs()
    println("Starting 3D Half-Tank Generation with Boundary Conditions...")
    
    # Geometry parameters
    tank_length = 1.0
    tank_radius = 0.1
    nozzle_length = 0.15
    nozzle_radius = 0.03
    wedge_angle_deg = 45.0  # Half tank - 45 degree wedge
    mesh_size_fine = 0.008
    mesh_size_coarse = 0.025
    
    println("\n📐 Creating 3D Half-Tank (45° wedge):")
    println("   Tank Length: $(tank_length)")
    println("   Tank Radius: $(tank_radius)")
    println("   Wedge Angle: $(wedge_angle_deg)° (half tank)")
    
    wedge_angle_rad = deg2rad(wedge_angle_deg)
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("3D_half_tank_with_BCs")
    
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
    l1 = gmsh.model.geo.addLine(1, 2, 1)   # Bottom wall
    l2 = gmsh.model.geo.addLine(2, 3, 2)   # Right wall (bottom part)
    l3 = gmsh.model.geo.addLine(3, 4, 3)   # Right wall (nozzle)
    l4 = gmsh.model.geo.addLine(4, 5, 4)   # Right wall (top part)
    l5 = gmsh.model.geo.addLine(5, 6, 5)   # Top wall
    l6 = gmsh.model.geo.addLine(6, 7, 6)   # Symmetry line (top)
    l7 = gmsh.model.geo.addLine(7, 8, 7)   # Nozzle wall (top)
    l8 = gmsh.model.geo.addLine(8, 9, 8)   # Inlet opening
    l9 = gmsh.model.geo.addLine(9, 10, 9)  # Nozzle wall (bottom)
    l10 = gmsh.model.geo.addLine(10, 1, 10) # Symmetry line (bottom)
    l11 = gmsh.model.geo.addLine(10, 7, 11) # Internal line
    l12 = gmsh.model.geo.addLine(3, 10, 12) # Internal line
    l13 = gmsh.model.geo.addLine(7, 4, 13)  # Internal line
    
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
    
    # Create physical groups for 2D boundaries (will be used to identify 3D boundaries)
    wall_lines = [1, 2, 3, 4, 5, 6, 7, 9, 10]  # All walls except inlet
    inlet_lines = [8]  # Inlet opening
    
    gmsh.model.addPhysicalGroup(1, wall_lines, -1, "wall_lines_2D")
    gmsh.model.addPhysicalGroup(1, inlet_lines, -1, "inlet_lines_2D")
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
    
    # ========================================================================
    # STEP 3: Revolve to create 3D half-tank with structured mesh
    # ========================================================================
    println("🔧 Step 3: Revolving structured surfaces to create 3D half-tank...")
    
    # Calculate number of layers for revolution (45° wedge)
    num_layers = 8  # Number of elements in circumferential direction for 45°
    
    println("   Revolving all surfaces through $(wedge_angle_deg)° with $(num_layers) layers...")
    println("   This creates a half-tank with two flat cut surfaces for air exposure")
    
    # Set GMSH options for hexahedral meshing BEFORE revolution
    gmsh.option.setNumber("Mesh.RecombineAll", 1)     # Enable recombination globally
    gmsh.option.setNumber("Mesh.Algorithm", 8)         # Frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)       # Delaunay for 3D
    gmsh.option.setNumber("Mesh.Recombine3DAll", 1)    # Enable 3D recombination
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # All-hex subdivision
    
    # Revolve with structured layers
    revolve_result = gmsh.model.geo.revolve(
        [(2, 1), (2, 2), (2, 3), (2, 4)],  # All 4 surfaces to revolve
        0.0, 0.0, 0.0,                      # Point on rotation axis
        1.0, 0.0, 0.0,                      # Rotation axis (X-axis)
        wedge_angle_rad,                    # Rotation angle
        [num_layers],                       # Number of structured layers
        [],                                 # Heights (empty for equal spacing)
        true                               # Recombine to create hexahedra
    )
    
    gmsh.model.geo.synchronize()
    
    # ========================================================================
    # STEP 4: Identify and label boundary surfaces for Ferrite
    # ========================================================================
    println("🔧 Step 4: Identifying and labeling boundary surfaces...")
    
    # Get all surfaces after revolution
    all_surfaces = gmsh.model.getEntities(2)
    println("   Found $(length(all_surfaces)) surfaces after revolution")
    
    # We need to identify different types of surfaces:
    # 1. Wall surfaces (curved tank walls, nozzle walls)
    # 2. Inlet surface (the small opening at the end of the nozzle)
    # 3. Air section surfaces (the flat cut surfaces created by the wedge)
    
    wall_surfaces = []
    inlet_surfaces = []
    air_section_surfaces = []
    
    # Iterate through all surfaces and classify them using detailed analysis
    for (dim, surf_tag) in all_surfaces
        classification = classify_surface_detailed(surf_tag, Lin)
        
        if classification == :inlet
            push!(inlet_surfaces, surf_tag)
        elseif classification == :air_section
            push!(air_section_surfaces, surf_tag)
        else  # classification == :wall
            push!(wall_surfaces, surf_tag)
        end
    end
    
    println("\n📊 Surface classification summary:")
    println("   Wall surfaces: $(length(wall_surfaces))")
    println("   Inlet surfaces: $(length(inlet_surfaces))")  
    println("   Air section surfaces: $(length(air_section_surfaces))")
    
    # Create physical groups for boundary conditions
    if !isempty(wall_surfaces)
        gmsh.model.addPhysicalGroup(2, wall_surfaces, -1, "walls")
        println("   ✅ Created 'walls' physical group with $(length(wall_surfaces)) surfaces")
    end
    
    if !isempty(inlet_surfaces)
        gmsh.model.addPhysicalGroup(2, inlet_surfaces, -1, "inlet")
        println("   ✅ Created 'inlet' physical group with $(length(inlet_surfaces)) surfaces")
    end
    
    if !isempty(air_section_surfaces)
        gmsh.model.addPhysicalGroup(2, air_section_surfaces, -1, "air_section")
        println("   ✅ Created 'air_section' physical group with $(length(air_section_surfaces)) surfaces")
    end
    
    # Create physical group for the volume (fluid domain)
    volumes = gmsh.model.getEntities(3)
    if !isempty(volumes)
        volume_tags = [tag for (dim, tag) in volumes]
        gmsh.model.addPhysicalGroup(3, volume_tags, -1, "fluid_domain")
        println("   ✅ Created 'fluid_domain' physical group with $(length(volume_tags)) volumes")
    end
    
    # ========================================================================
    # VERIFICATION: Check surface classifications
    # ========================================================================
    verify_surface_classification(all_surfaces, wall_surfaces, inlet_surfaces, air_section_surfaces, Lin)
    
    # Set transfinite volumes for structured hexahedral meshing
    for (dim, tag) in volumes
        try
            gmsh.model.geo.mesh.setTransfiniteVolume(tag)
            println("   Set transfinite volume $(tag)")
        catch e
            println("   Warning: Could not set transfinite for volume $(tag): $(e)")
        end
    end
    
    # Ensure recombination for all surfaces and volumes
    for (dim, tag) in all_surfaces
        try
            gmsh.model.mesh.setRecombine(2, tag)
        catch e
            # Some surfaces might not support recombination
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
    
    # ========================================================================
    # STEP 5: Generate 3D mesh with hexahedral elements
    # ========================================================================
    println("🔧 Step 5: Generating 3D hexahedral mesh...")
    
    # Additional settings to ensure hexahedral elements
    gmsh.option.setNumber("Mesh.RecombineAll", 1)           # Global recombination
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1.0)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)           # Linear elements
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)      # Disable high-order
    
    # Generate the mesh step by step
    try
        # First generate 1D mesh (edges)
        gmsh.model.mesh.generate(1)
        println("   ✅ 1D mesh generated")
        
        # Then generate 2D mesh (surfaces)
        gmsh.model.mesh.generate(2)
        println("   ✅ 2D mesh generated")
        
        # Finally generate 3D mesh (volumes)
        gmsh.model.mesh.generate(3)
        println("   ✅ 3D hexahedral mesh generated successfully")
        
        # Check element types to verify we have hexahedra
        element_types = gmsh.model.mesh.getElementTypes()
        println("   📊 Element types found: $(element_types)")
        
        has_hexahedra = false
        for elem_type in element_types
            elem_name = gmsh.model.mesh.getElementProperties(elem_type)[1]
            num_elements = length(gmsh.model.mesh.getElementsByType(elem_type)[1])
            println("   - Type $(elem_type) ($(elem_name)): $(num_elements) elements")
            
            # Check if we have hexahedral elements (type 5 = 8-node hexahedron)
            if elem_type == 5
                has_hexahedra = true
                println("   ✅ Found hexahedral elements!")
            end
        end
        
        if !has_hexahedra
            println("   ⚠️  Warning: No hexahedral elements found. You may have tetrahedral elements.")
            println("   💡 This might be due to complex geometry. Consider simplifying the mesh.")
        end
        
    catch e
        println("   ❌ Error generating 3D mesh: $(e)")
        println("   🔄 Trying fallback mesh generation...")
        
        # Fallback: Try without some strict settings
        gmsh.option.setNumber("Mesh.Algorithm3D", 4)  # Try Frontal algorithm
        try
            gmsh.model.mesh.generate(3)
            println("   ✅ 3D mesh generated with fallback settings")
        catch e2
            println("   ❌ Fallback also failed: $(e2)")
        end
    end
    
    # ========================================================================
    # STEP 6: Save and visualize
    # ========================================================================
    
    # Save debug mesh with color-coded surfaces before generating 3D mesh
    gmsh.write("tank_3D_half_surfaces_debug.msh")
    println("   💾 Debug mesh saved: tank_3D_half_surfaces_debug.msh")
    println("   🎨 In GMSH GUI: Tools > Options > Mesh > Surface labels = 1 to see surface numbers")
    
    gmsh.write("tank_3D_half_with_BCs.msh")
    println("   ✅ 3D half-tank mesh with boundary conditions saved: tank_3D_half_with_BCs.msh")
    
    # Show mesh info
    mesh_info = gmsh.model.mesh.getNodes()
    num_nodes = length(mesh_info[1])
    println("   📊 Total nodes: $(num_nodes)")
    
    # Print summary of boundary conditions for Ferrite
    println("\n🏷️  Boundary Condition Labels for Ferrite.jl:")
    println("   - 'walls': Tank and nozzle wall surfaces (no-slip/wall BC)")
    println("   - 'inlet': Small nozzle opening (inflow BC)")
    println("   - 'air_section': Exposed cut surfaces (outflow/symmetry BC)")
    println("   - 'fluid_domain': Interior volume (fluid domain)")
    
    # Print physical group details and verification instructions
    print_physical_groups()
    print_verification_instructions()
    
    # Show 3D mesh
    println("   🖥️  Showing 3D half-tank mesh with boundary conditions...")
    gmsh.fltk.run()
    
    gmsh.finalize()
    println("\n✅ Done! Half-tank mesh ready for Ferrite.jl")
end

create_3d_wedge_tank_with_bcs()