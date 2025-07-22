# ============================================================================
# Main Script for 3D Tank Mesh Generation
# ============================================================================
# Main entry point for generating 3D tank meshes with boundary conditions

println("🚀 Starting 3D Tank Mesh Generation Pipeline...")

# Include required files
include("gmsh_3Dtanknozzle_with_BCs.jl")
include("mesh_verification.jl")

# Import GMSH
try
    using Gmsh: gmsh
    println("✅ GMSH imported successfully!")
catch
    using gmsh
    println("✅ GMSH imported successfully!")
end

function main()
    println("\n" * "="^60)
    println("         3D TANK MESH GENERATION PIPELINE")
    println("="^60)
    
    create_3d_wedge_tank_with_bcs()
end

# Configuration options
const CONFIG = Dict(
    # Geometry parameters
    :tank_length => 1.0,
    :tank_radius => 0.1,
    :nozzle_length => 0.15,
    :nozzle_radius => 0.03,
    :wedge_angle_deg => 45.0,
    
    # Mesh parameters
    :mesh_size_fine => 0.008,
    :mesh_size_coarse => 0.025,
    :num_layers => 8,
    
    # Output options
    :save_debug_mesh => true,
    :show_gmsh_gui => true,
    :verbose => true
)

"""
    run_with_config(config::Dict)

Run mesh generation with custom configuration.
"""
function run_with_config(config::Dict)
    println("🔧 Running with custom configuration:")
    for (key, value) in config
        println("   $(key): $(value)")
    end
    
    # You can modify the mesh generation function to accept these parameters
    # For now, it uses the hardcoded values in the function
    main()
end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
