# ============================================================================
# Mesh Verification Functions for 3D Tank Generation
# ============================================================================
# Contains functions to verify and test mesh boundary condition labeling

"""
    verify_surface_classification(all_surfaces, wall_surfaces, inlet_surfaces, air_section_surfaces, Lin)

Verify that all surfaces are correctly classified for boundary conditions.
"""
function verify_surface_classification(all_surfaces, wall_surfaces, inlet_surfaces, air_section_surfaces, Lin)
    println("\n🔍 VERIFICATION: Checking surface classifications...")
    
    # Expected surface counts for a 45-degree wedge:
    # - 2 air section surfaces (the two flat cut planes)
    # - 1 inlet surface (the small circular opening)
    # - Remaining surfaces should be walls
    
    println("   Expected for 45° wedge:")
    println("   - Air sections: 2 (flat cut surfaces)")
    println("   - Inlet: 1 (small opening)")
    println("   - Walls: $(length(all_surfaces) - 3) (all other surfaces)")
    
    println("\n   Actual classification:")
    println("   - Air sections: $(length(air_section_surfaces))")
    println("   - Inlet: $(length(inlet_surfaces))")
    println("   - Walls: $(length(wall_surfaces))")
    
    # Detailed verification
    total_classified = length(wall_surfaces) + length(inlet_surfaces) + length(air_section_surfaces)
    total_surfaces = length(all_surfaces)
    
    if total_classified == total_surfaces
        println("   ✅ All $(total_surfaces) surfaces classified!")
    else
        println("   ❌ Missing surfaces! Classified $(total_classified) out of $(total_surfaces)")
    end
    
    # Check if we have the expected inlet
    if length(inlet_surfaces) == 1
        println("   ✅ Found exactly 1 inlet surface (expected)")
    elseif length(inlet_surfaces) == 0
        println("   ⚠️  No inlet surfaces found - check inlet detection logic")
    else
        println("   ⚠️  Found $(length(inlet_surfaces)) inlet surfaces - expected 1")
    end
    
    # Check if we have the expected air sections
    if length(air_section_surfaces) == 2
        println("   ✅ Found exactly 2 air section surfaces (expected for 45° wedge)")
    elseif length(air_section_surfaces) == 0
        println("   ⚠️  No air section surfaces found - check cut plane detection")
    else
        println("   ⚠️  Found $(length(air_section_surfaces)) air section surfaces - expected 2")
    end
    
    # Validate inlet position
    if !isempty(inlet_surfaces)
        for inlet_surf in inlet_surfaces
            bbox = gmsh.model.getBoundingBox(2, inlet_surf)
            x_center = (bbox[1] + bbox[4]) / 2
            expected_x = Lin
            if abs(x_center - expected_x) < 0.01
                println("   ✅ Inlet surface $(inlet_surf) at correct x position: $(round(x_center,digits=4)) ≈ $(expected_x)")
            else
                println("   ⚠️  Inlet surface $(inlet_surf) at wrong x position: $(round(x_center,digits=4)) ≠ $(expected_x)")
            end
        end
    end
    
    # Validate air section planes
    if !isempty(air_section_surfaces)
        y_plane_found = false
        z_plane_found = false
        
        for air_surf in air_section_surfaces
            bbox = gmsh.model.getBoundingBox(2, air_surf)
            y_span = abs(bbox[5] - bbox[2])
            z_span = abs(bbox[6] - bbox[3])
            y_center = (bbox[2] + bbox[5]) / 2
            z_center = (bbox[3] + bbox[6]) / 2
            
            if y_span < 0.01 && abs(y_center) < 0.01
                y_plane_found = true
                println("   ✅ Air section surface $(air_surf) is y=0 plane")
            elseif z_span < 0.01 && abs(z_center) < 0.01
                z_plane_found = true
                println("   ✅ Air section surface $(air_surf) is z=0 plane")
            else
                println("   ⚠️  Air section surface $(air_surf) doesn't appear to be a cut plane")
            end
        end
        
        if y_plane_found && z_plane_found
            println("   ✅ Both y=0 and z=0 cut planes found correctly")
        else
            println("   ⚠️  Missing cut planes: y=0 found: $(y_plane_found), z=0 found: $(z_plane_found)")
        end
    end
    
    return total_classified == total_surfaces
end

"""
    classify_surface_detailed(surf_tag, Lin; flat_tol=1e-6)

Classify a single surface with detailed analysis and logging.
Returns classification type as Symbol: :inlet, :air_section, or :wall
"""
function classify_surface_detailed(surf_tag, Lin; flat_tol=1e-6)
    try
        # Get surface bounding box to help classify it
        bbox = gmsh.model.getBoundingBox(2, surf_tag)
        x_min, y_min, z_min = bbox[1], bbox[2], bbox[3]
        x_max, y_max, z_max = bbox[4], bbox[5], bbox[6]
        
        # Calculate center and dimensions
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        z_center = (z_min + z_max) / 2
        
        # Calculate surface dimensions
        x_span = abs(x_max - x_min)
        y_span = abs(y_max - y_min)
        z_span = abs(z_max - z_min)
        
        println("   Surface $(surf_tag): bbox=[$(round(x_min,digits=4)),$(round(x_max,digits=4))] x [$(round(y_min,digits=4)),$(round(y_max,digits=4))] x [$(round(z_min,digits=4)),$(round(z_max,digits=4))]")
        println("      Center: ($(round(x_center,digits=4)), $(round(y_center,digits=4)), $(round(z_center,digits=4)))")
        println("      Spans: x=$(round(x_span,digits=4)), y=$(round(y_span,digits=4)), z=$(round(z_span,digits=4))")
        
        # Classification logic with detailed checks
        if abs(x_center - Lin) < 1e-3 && x_span < flat_tol
            # Surface is at the inlet location and flat in x-direction
            println("      ✅ INLET: x ≈ $(Lin) and flat in x-direction")
            return :inlet
        
        elseif y_span < flat_tol && abs(y_center) < 1e-3
            # Flat surface in y=0 plane (wedge cut surface)
            println("      ✅ AIR_SECTION: flat in y=0 plane")
            return :air_section
        
        elseif z_span < flat_tol && abs(z_center) < 1e-3
            # Flat surface in z=0 plane (wedge cut surface)
            println("      ✅ AIR_SECTION: flat in z=0 plane")
            return :air_section
        
        # Additional check for inlet surfaces that might be slightly off
        elseif abs(x_center - Lin) < 0.01 && x_span < 0.01
            # Surface is near the inlet location
            println("      ✅ INLET: near x ≈ $(Lin) (relaxed tolerance)")
            return :inlet
        
        # Check for air section surfaces with relaxed tolerance
        elseif y_span < 0.01 && abs(y_center) < 0.01
            println("      ✅ AIR_SECTION: near y=0 plane (relaxed tolerance)")
            return :air_section
        
        elseif z_span < 0.01 && abs(z_center) < 0.01
            println("      ✅ AIR_SECTION: near z=0 plane (relaxed tolerance)")
            return :air_section
        
        # All other surfaces are walls
        else
            println("      ✅ WALL: curved/other surface")
            return :wall
        end
        
    catch e
        # If we can't get bounding box, assume it's a wall
        println("      ❌ WALL: classification failed - $(e)")
        return :wall
    end
end

"""
    print_physical_groups()

Print details of all physical groups in the current GMSH model.
"""
function print_physical_groups()
    println("\n📋 Physical Group Details:")
    physical_groups = gmsh.model.getPhysicalGroups()
    for (dim, tag) in physical_groups
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        entity_type = dim == 1 ? "curves" : (dim == 2 ? "surfaces" : "volumes")
        println("   $(name) ($(dim)D): $(length(entities)) $(entity_type) - tags: $(entities)")
    end
end

"""
    print_verification_instructions()

Print manual verification instructions for the user.
"""
function print_verification_instructions()
    println("\n🔧 Manual Verification Tips:")
    println("   1. Open tank_3D_half_surfaces_debug.msh in GMSH")
    println("   2. Go to Tools > Options > Mesh > General")
    println("   3. Set 'Surface labels' = 1 to see surface numbers")
    println("   4. Check that:")
    println("      - Inlet surfaces are at x ≈ -0.15 (nozzle end)")
    println("      - Air section surfaces are the flat cut planes")
    println("      - Wall surfaces are all curved tank/nozzle surfaces")
    println("   5. Use View > Physical groups to highlight each BC type")
end
