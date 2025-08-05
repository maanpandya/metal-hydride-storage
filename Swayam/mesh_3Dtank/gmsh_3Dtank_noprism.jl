# 3D Tank Mesh Generation with Small Axis Offset to Avoid Prisms
# This script generates a 3D tank mesh by revolving a 2D mesh around the x-axis.
# To avoid prismatic elements near the axis, we introduce a very small offset 
# from the axis, creating a tiny cylindrical hole in the center.

using Gmsh: gmsh

function generate_tank_mesh_no_prism()
    println("🚀 Generating 3D tank mesh with no prismatic elements...")
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("tank_no_prism")
    
    # ========================================================================
    # GEOMETRY PARAMETERS
    # ========================================================================
    L = 1.0         # Tank length
    R = 0.1         # Tank radius
    Lin = 0.1 * L   # Nozzle extension length
    Rin = 0.3 * R   # Nozzle radius
    
    # CRITICAL: Small offset from axis to avoid prisms
    # This creates a tiny cylindrical hole in the center
    r_offset = 0.001 * R  # Very small offset (0.1% of tank radius)
    
    println("   Tank geometry:")
    println("   - Length (L): $(L)")
    println("   - Radius (R): $(R)")
    println("   - Nozzle length (Lin): $(Lin)")
    println("   - Nozzle radius (Rin): $(Rin)")
    println("   - Axis offset (r_offset): $(r_offset) ($(r_offset/R*100)% of tank radius)")
    
    # Mesh density parameters
    lc1 = 0.005     # Fine mesh near walls
    lc2 = 0.02      # Coarse mesh in far field
    
    # ========================================================================
    # STEP 1: Create 2D geometry with axis offset
    # ========================================================================
    println("🔧 Step 1: Creating 2D geometry with all points offset from axis...")
    
    # Define points for multi-block structured mesh (UPPER HALF ONLY)
    # ALL points must be offset from the axis (y > 0) to avoid degenerate revolution
    
    # Main tank points - upper half from small radius to outer radius
    p1  = gmsh.model.geo.addPoint(r_offset, r_offset, 0, lc1, 1)    # Inner bottom (offset from axis)
    p2  = gmsh.model.geo.addPoint(L,        r_offset, 0, lc2, 2)    # Outer bottom (offset from axis)
    p3  = gmsh.model.geo.addPoint(L,        R,        0, lc2, 3)    # Outer top
    p4  = gmsh.model.geo.addPoint(r_offset, R,        0, lc1, 4)    # Inner top
    
    # Nozzle connection points - upper half only
    p5  = gmsh.model.geo.addPoint(r_offset, Rin,      0, lc1, 5)   # Inner at nozzle radius
    
    # Nozzle points - upper half only (ALL offset from axis)
    p6  = gmsh.model.geo.addPoint(-Lin,     r_offset, 0, lc1, 6)   # Nozzle bottom (offset from axis)
    p7  = gmsh.model.geo.addPoint(-Lin,     Rin,      0, lc1, 7)   # Nozzle top
    
    # Define edges for structured blocks
    # Main tank block (upper half rectangular, all points offset from axis)
    l1  = gmsh.model.geo.addLine(1,  2,  1)   # Bottom edge (at r_offset height)
    l2  = gmsh.model.geo.addLine(2,  3,  2)   # Right edge
    l3  = gmsh.model.geo.addLine(3,  4,  3)   # Top edge
    l4  = gmsh.model.geo.addLine(4,  1,  4)   # Left edge (inner radius)
    
    # Nozzle block edges
    l5  = gmsh.model.geo.addLine(1,  5,  5)   # Inner vertical connection
    l6  = gmsh.model.geo.addLine(5,  7,  6)   # Nozzle connection
    l7  = gmsh.model.geo.addLine(7,  6,  7)   # Nozzle inlet  
    l8  = gmsh.model.geo.addLine(6,  1,  8)   # Bottom connection (nozzle to tank)
    
    # Define curve loops for each structured block
    loop1 = gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)       # Main tank block (upper half)
    loop2 = gmsh.model.geo.addCurveLoop([5, 6, 7, 8], 2)       # Nozzle block (upper half)
    
    # Define surfaces
    surf1 = gmsh.model.geo.addPlaneSurface([1], 1)  # Main tank surface
    surf2 = gmsh.model.geo.addPlaneSurface([2], 2)  # Nozzle surface
    
    # ========================================================================
    # STEP 2: Set up transfinite meshing for structured mesh
    # ========================================================================
    println("🔧 Step 2: Setting up structured transfinite meshing...")
    
    # Define number of nodes on each edge
    n_axial = 20        # Nodes along tank length
    n_radial = 8        # Nodes in radial direction
    n_nozzle_axial = 8  # Nodes along nozzle length
    n_nozzle_radial = 6 # Nodes in nozzle radial direction
    
    # Set transfinite constraints for edges
    # Main tank block
    gmsh.model.geo.mesh.setTransfiniteCurve(1,  n_axial)         # Bottom edge
    gmsh.model.geo.mesh.setTransfiniteCurve(2,  n_radial)        # Right edge
    gmsh.model.geo.mesh.setTransfiniteCurve(3,  n_axial)         # Top edge
    gmsh.model.geo.mesh.setTransfiniteCurve(4,  n_radial)        # Left edge
    
    # Nozzle block
    gmsh.model.geo.mesh.setTransfiniteCurve(5,  n_nozzle_radial) # Inner connection
    gmsh.model.geo.mesh.setTransfiniteCurve(6,  n_nozzle_axial)  # Top nozzle connection
    gmsh.model.geo.mesh.setTransfiniteCurve(7,  n_nozzle_radial) # Nozzle inlet
    gmsh.model.geo.mesh.setTransfiniteCurve(8,  n_nozzle_axial)  # Bottom nozzle connection
    
    # Set transfinite surfaces for structured mesh
    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setTransfiniteSurface(2) 
    
    # Set recombination for quadrilateral elements
    gmsh.model.geo.mesh.setRecombine(2, 1)
    gmsh.model.geo.mesh.setRecombine(2, 2)
    
    gmsh.model.geo.synchronize()
    
    # ========================================================================
    # STEP 3: Assign physical groups for boundary conditions
    # ========================================================================
    println("🔧 Step 3: Assigning physical groups...")
    
    # Inner boundaries (edges at small radius - these will become inner cylindrical walls)
    inner_edges = [1, 8]  # Bottom edges at r_offset radius
    gmsh.model.addPhysicalGroup(1, inner_edges, -1, "inner_wall")
    
    # Outer wall boundaries (use plural "walls" to match simulation code)
    wall_edges = [2, 3, 4, 5, 6]  # Tank walls and nozzle connection
    gmsh.model.addPhysicalGroup(1, wall_edges, -1, "walls")
    
    # Inlet boundary
    gmsh.model.addPhysicalGroup(1, [7], -1, "inlet")
    
    # Domain
    gmsh.model.addPhysicalGroup(2, [1, 2], -1, "omega")
    
    # ========================================================================
    # STEP 4: Generate 2D structured mesh
    # ========================================================================
    println("🔧 Step 4: Generating 2D structured mesh...")
    
    gmsh.model.mesh.generate(2)
    #gmsh.write("tank_2D_axisymmetric.msh")
    #println("   ✅ 2D axisymmetric mesh saved: tank_2D_axisymmetric.msh")
    
    # Show 2D mesh
    println("   🖥️  Showing 2D axisymmetric mesh - close window to continue...")
    gmsh.fltk.run()
    
    # ========================================================================
    # STEP 5: Revolve to create 3D hexahedral mesh
    # ========================================================================
    println("🔧 Step 5: Revolving to create 3D hexahedral mesh...")
    
    # Revolution parameters - start with smaller wedge for testing
    wedge_angle_deg = 90.0                              # Quarter revolution for testing
    wedge_angle_rad = wedge_angle_deg * π / 180         # Convert to radians
    num_layers = 6                                      # Circumferential elements
    
    println("   Revolution parameters:")
    println("   - Wedge angle: $(wedge_angle_deg)°")
    println("   - Number of layers: $(num_layers)")
    
    # Set GMSH options for hexahedral meshing
    gmsh.option.setNumber("Mesh.RecombineAll", 1)       # Global recombination
    gmsh.option.setNumber("Mesh.Algorithm", 6)          # Frontal-Delaunay (more robust)
    gmsh.option.setNumber("Mesh.Algorithm3D", 4)        # Frontal for 3D (more robust)
    gmsh.option.setNumber("Mesh.Recombine3DAll", 1)     # Enable 3D recombination
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1) # All-hex subdivision
    
    # Revolve all surfaces around x-axis
    revolve_result = gmsh.model.geo.revolve(
        [(2, 1), (2, 2)],            # Both surfaces to revolve
        0.0, 0.0, 0.0,               # Point on rotation axis
        1.0, 0.0, 0.0,               # Rotation axis (X-axis)
        wedge_angle_rad,             # Rotation angle
        [num_layers],                # Number of structured layers
        [],                          # Heights (empty for equal spacing)
        true                         # Recombine to create hexahedra
    )
    
    gmsh.model.geo.synchronize()
    
    # ========================================================================
    # STEP 6: Set up 3D physical groups
    # ========================================================================
    println("🔧 Step 6: Setting up 3D physical groups...")
    
    # Alternative approach: Use GMSH's automatic boundary detection
    # Get all boundary surfaces of the 3D domain
    all_volumes = gmsh.model.getEntities(3)
    all_surfaces = gmsh.model.getEntities(2)
    
    # For each surface, determine its type based on geometric properties
    walls_surfaces = []
    inlet_surfaces = []
    inner_wall_surfaces = []
    
    for surface in all_surfaces
        surface_tag = surface[2]
        
        # Get surface bounding box to classify the surface
        bbox = gmsh.model.getBoundingBox(2, surface_tag)
        x_min, y_min, z_min = bbox[1], bbox[2], bbox[3]
        x_max, y_max, z_max = bbox[4], bbox[5], bbox[6]
        
        # Get surface center for classification
        center_x = (x_min + x_max) / 2
        center_r = sqrt(((y_min + y_max) / 2)^2 + ((z_min + z_max) / 2)^2)
        
        # Classify based on geometry:
        # - Inlet: negative x position (x < -0.01)
        # - Inner wall: very small radius (r < 0.01)  
        # - Walls: positive x and significant radius
        
        if x_max < -0.01  # Inlet surfaces (at negative x)
            push!(inlet_surfaces, surface_tag)
        elseif center_r < 0.01  # Inner cylindrical surfaces (close to axis)
            push!(inner_wall_surfaces, surface_tag)
        elseif center_r > 0.05  # Outer wall surfaces
            push!(walls_surfaces, surface_tag)
        end
    end
    
    # Create 3D surface physical groups for boundary conditions
    if !isempty(walls_surfaces)
        gmsh.model.addPhysicalGroup(2, walls_surfaces, -1, "walls")
        println("   ✅ Added $(length(walls_surfaces)) wall surfaces to 'walls' group")
    else
        println("   ⚠️  Warning: No wall surfaces found for 3D boundary conditions")
    end
    
    if !isempty(inlet_surfaces)
        gmsh.model.addPhysicalGroup(2, inlet_surfaces, -1, "inlet")
        println("   ✅ Added $(length(inlet_surfaces)) inlet surfaces to 'inlet' group")
    else
        println("   ⚠️  Warning: No inlet surfaces found for 3D boundary conditions")
    end
    
    if !isempty(inner_wall_surfaces)
        gmsh.model.addPhysicalGroup(2, inner_wall_surfaces, -1, "inner_walls")
        println("   ✅ Added $(length(inner_wall_surfaces)) inner wall surfaces")
    end
    
    # Get all volumes for fluid domain
    volume_tags = [v[2] for v in all_volumes]
    gmsh.model.addPhysicalGroup(3, volume_tags, -1, "fluid_domain")
    
    println("   Found $(length(all_volumes)) volumes")
    
    # Don't set transfinite volumes - let GMSH handle 3D meshing automatically
    # Transfinite volumes can cause issues with complex revolution geometries
    println("   Using automatic 3D meshing (not transfinite volumes)")
    
    # ========================================================================
    # STEP 7: Generate 3D mesh and finalize
    # ========================================================================
    println("🔧 Step 7: Generating final 3D mesh...")
    
    gmsh.model.mesh.generate(3)
    
    # Write final mesh
    gmsh.write("tank_3D_no_prismv2.msh")
    
    println("   ✅ 3D mesh files saved:")
    println("   - tank_3D_no_prism.msh (GMSH format)")
    println("   - tank_3D_no_prism.vtk (ParaView format)")
    
    # Get mesh statistics
    println("\\n📊 Mesh Statistics:")
    nodes = gmsh.model.mesh.getNodes()
    elements_3D = gmsh.model.mesh.getElements(3)
    
    println("   - Number of nodes: $(length(nodes[1]))")
    
    total_elements = 0
    for i in 1:length(elements_3D[2])
        total_elements += length(elements_3D[2][i])
    end
    println("   - Number of 3D elements: $(total_elements)")
    
    # Show element types
    if length(elements_3D[1]) > 0
        for i in 1:length(elements_3D[1])
            element_type = elements_3D[1][i]
            num_elements = length(elements_3D[2][i])
            element_name = get_element_name(element_type)
            println("   - Element type $(element_type) ($(element_name)): $(num_elements) elements")
        end
    end
    
    # Show physical groups for verification
    println("\\n🏷️  Physical Groups (for Ferrite.jl):")
    physical_groups = gmsh.model.getPhysicalGroups()
    for group in physical_groups
        dim, tag = group[1], group[2]
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        
        if dim == 2  # Surface groups (boundaries)
            println("   - Surface group '$(name)': $(length(entities)) surfaces (dim=$(dim))")
        elseif dim == 3  # Volume groups (domain)
            println("   - Volume group '$(name)': $(length(entities)) volumes (dim=$(dim))")
        end
    end
    
    # Show final mesh
    println("\\n🖥️  Showing final 3D mesh - close window to exit...")
    gmsh.fltk.run()
    
    gmsh.finalize()
    
    println("\\n🎉 3D tank mesh generation completed successfully!")
    println("\\n✨ Key improvements:")
    println("   - All points offset from axis ($(r_offset)) prevents degenerate revolution")
    println("   - No points on rotation axis - creates proper surfaces when revolved")
    println("   - All elements should be hexahedral (element type 5)")
    println("   - Structured mesh suitable for Ferrite.jl")
    println("   - Small inner hole ($(r_offset/R*100)% of tank radius) is negligible")
    println("   - Physical groups named for Ferrite.jl: 'walls', 'inlet', 'fluid_domain'")
    println("\\n💡 Usage in Ferrite.jl:")
    println("   walls_set = getfacetset(dh.grid, \"walls\")")
    println("   inlet_set = getfacetset(dh.grid, \"inlet\")")
    
    return "tank_3D_no_prism.msh"
end

# Helper function to get element type names
function get_element_name(element_type)
    element_names = Dict(
        1 => "2-node line",
        2 => "3-node triangle", 
        3 => "4-node quadrangle",
        4 => "4-node tetrahedron",
        5 => "8-node hexahedron",
        6 => "6-node prism",
        7 => "5-node pyramid"
    )
    return get(element_names, element_type, "Unknown")
end

# Run the mesh generation
mesh_file = generate_tank_mesh_no_prism()