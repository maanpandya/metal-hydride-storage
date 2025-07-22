# 2D Tank-Nozzle Generation using GMSH with Boundary Conditions

try
    using Gmsh: gmsh
catch
    using gmsh
end

function create_2d_tank_nozzle_with_bcs()

    OUTPUT_MESH_FILE_NAME = "tank_2D_quad_default.msh" # add the .msh extension

    # Geometry parameters
    tank_length = 1.0
    tank_radius = 0.1
    nozzle_length = 0.15
    nozzle_radius = 0.03
    
    # Mesh subdivision (refinement) parameters
    n_nozzle_vertical = 30 # same as center tank veritical subdivision      
    n_nozzle_horizontal = 30 # unique to nozzle (only used for nozzle)
    n_tank_vertical = 15 # only used in tank top and bottom part (not center)
    n_tank_horizontal = 100 # used for tank center, bottom and top but not nozzle
    
    # Initialize GMSH
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("2D_tank_nozzle_with_BCs")
    
    # Create geometry
    L = tank_length
    R = tank_radius
    Lin = -nozzle_length
    Rin = nozzle_radius
    
    # Define points (mesh size not needed for transfinite mesh)
    p1  = gmsh.model.geo.addPoint(0,   -R,   0, 0, 1)
    p2  = gmsh.model.geo.addPoint(L,   -R,   0, 0, 2)
    p3  = gmsh.model.geo.addPoint(L,   -Rin, 0, 0, 3)
    p4  = gmsh.model.geo.addPoint(L,    Rin, 0, 0, 4)
    p5  = gmsh.model.geo.addPoint(L,    R,   0, 0, 5)
    p6  = gmsh.model.geo.addPoint(0,    R,   0, 0, 6)
    p7  = gmsh.model.geo.addPoint(0,    Rin, 0, 0, 7)
    p8  = gmsh.model.geo.addPoint(Lin,  Rin, 0, 0, 8)
    p9  = gmsh.model.geo.addPoint(Lin, -Rin, 0, 0, 9)
    p10 = gmsh.model.geo.addPoint(0,   -Rin, 0, 0, 10)
    
    # Define lines
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
    
    # Define surfaces
    loop1 = gmsh.model.geo.addCurveLoop([1, 2, 12, 10], 1)
    loop2 = gmsh.model.geo.addCurveLoop([-12, 3, -13, -11], 2)
    loop3 = gmsh.model.geo.addCurveLoop([13, 4, 5, 6], 3)
    loop4 = gmsh.model.geo.addCurveLoop([9, 11, 7, 8], 4)
    
    surf1 = gmsh.model.geo.addPlaneSurface([1], 1)
    surf2 = gmsh.model.geo.addPlaneSurface([2], 2)
    surf3 = gmsh.model.geo.addPlaneSurface([3], 3)
    surf4 = gmsh.model.geo.addPlaneSurface([4], 4)
    
    # Set transfinite curves for structured mesh
    # Surface 1 (bottom-right): horizontal tank lines + vertical tank lines
    gmsh.model.geo.mesh.setTransfiniteCurve(1, n_tank_horizontal)   # Bottom tank
    gmsh.model.geo.mesh.setTransfiniteCurve(12, n_tank_horizontal)  # Internal horizontal bottom
    gmsh.model.geo.mesh.setTransfiniteCurve(2, n_tank_vertical)     # Right tank bottom
    gmsh.model.geo.mesh.setTransfiniteCurve(10, n_tank_vertical)    # Left tank bottom
    
    # Surface 2 (center): horizontal tank lines + vertical nozzle lines
    gmsh.model.geo.mesh.setTransfiniteCurve(13, n_tank_horizontal)     # Internal horizontal top
    gmsh.model.geo.mesh.setTransfiniteCurve(3, n_nozzle_vertical)      # Right tank nozzle section
    gmsh.model.geo.mesh.setTransfiniteCurve(11, n_nozzle_vertical)     # Internal vertical
    
    # Surface 3 (top-right): horizontal tank lines + vertical tank lines
    gmsh.model.geo.mesh.setTransfiniteCurve(5, n_tank_horizontal)   # Top tank
    gmsh.model.geo.mesh.setTransfiniteCurve(4, n_tank_vertical)     # Right tank top
    gmsh.model.geo.mesh.setTransfiniteCurve(6, n_tank_vertical)     # Left tank top
    
    # Surface 4 (nozzle): horizontal nozzle lines + vertical nozzle lines
    gmsh.model.geo.mesh.setTransfiniteCurve(9, n_nozzle_horizontal)    # Bottom nozzle
    gmsh.model.geo.mesh.setTransfiniteCurve(7, n_nozzle_horizontal)    # Top nozzle
    gmsh.model.geo.mesh.setTransfiniteCurve(8, n_nozzle_vertical)      # Inlet opening
    
    # Set transfinite surfaces and recombine
    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setTransfiniteSurface(2)
    gmsh.model.geo.mesh.setTransfiniteSurface(3)
    gmsh.model.geo.mesh.setTransfiniteSurface(4)
    
    gmsh.model.geo.mesh.setRecombine(2, 1)
    gmsh.model.geo.mesh.setRecombine(2, 2)
    gmsh.model.geo.mesh.setRecombine(2, 3)
    gmsh.model.geo.mesh.setRecombine(2, 4)
    
    gmsh.model.geo.synchronize()

    # Define physical groups for boundary conditions
    wall_lines = [1, 2, 3, 4, 5, 6, 7, 9, 10]
    inlet_lines = [8]
    
    gmsh.model.addPhysicalGroup(1, wall_lines, -1, "walls")
    gmsh.model.addPhysicalGroup(1, inlet_lines, -1, "inlet")
    gmsh.model.addPhysicalGroup(2, [1, 2, 3, 4], -1, "fluid_domain")
    
    # Generate mesh
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    
    gmsh.model.mesh.generate(2)
    
    # Save mesh
    gmsh.write(OUTPUT_MESH_FILE_NAME)
    
    # Show mesh
    gmsh.fltk.run()
    
    gmsh.finalize()
end

# Run the function
create_2d_tank_nozzle_with_bcs()
