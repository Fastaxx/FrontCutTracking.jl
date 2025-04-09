
"""
    redistribute_markers!(interface::Interface; 
                       target_spacing=:auto, 
                       min_angle_deg=15.0)

Redistribute markers along the interface to maintain uniform spacing and smooth representation.

# Arguments
- `interface`: The Interface object to modify
- `target_spacing`: Desired distance between markers or :auto to calculate from average
- `min_angle_deg`: Minimum angle (in degrees) for adjacent segments

# Returns
The modified interface
"""
function redistribute_markers!(interface::Interface; 
                            target_spacing=:auto,
                            min_angle_deg=15.0)
    
    # Get current positions
    positions = [marker.position for marker in interface.markers]
    n_markers = length(positions)
    
    if n_markers <= 3
        return interface  # Not enough markers to redistribute
    end
    
    # Calculate current segment lengths
    segment_lengths = Float64[]
    if interface.closed
        for i in 1:n_markers
            next_i = i % n_markers + 1
            push!(segment_lengths, norm(positions[next_i] - positions[i]))
        end
    else
        for i in 1:n_markers-1
            push!(segment_lengths, norm(positions[i+1] - positions[i]))
        end
    end
    
    # Determine target spacing
    if target_spacing == :auto
        mean_spacing = mean(segment_lengths)
        target = mean_spacing
    else
        target = float(target_spacing)
    end
    
    # Convert minimum angle to radians
    min_angle = min_angle_deg * π / 180
    
    # Create a new set of positions with proper spacing
    new_positions = Vector{Vector{Float64}}()
    
    # For closed interfaces, parameterize by cumulative distance
    if interface.closed
        # Calculate cumulative distance along the interface
        cumulative_distance = [0.0]
        for i in 1:n_markers
            next_i = i % n_markers + 1
            push!(cumulative_distance, cumulative_distance[end] + norm(positions[next_i] - positions[i]))
        end
        
        # Total length of the interface
        total_length = cumulative_distance[end]
        
        # Calculate number of points in resampled interface
        n_new_points = max(4, round(Int, total_length / target))
        
        # Generate new points evenly spaced by arc length
        for i in 0:n_new_points-1
            t = i * total_length / n_new_points
            
            # Find the segment containing position t
            seg_idx = findlast(d -> d <= t, cumulative_distance)
            next_idx = seg_idx % n_markers + 1
            
            # Interpolate position
            seg_start = positions[seg_idx]
            seg_end = positions[next_idx]
            seg_length = norm(seg_end - seg_start)
            
            if seg_length > 1e-10
                # Linear interpolation parameter
                alpha = (t - cumulative_distance[seg_idx]) / seg_length
                new_pos = seg_start + alpha * (seg_end - seg_start)
                push!(new_positions, new_pos)
            else
                push!(new_positions, seg_start)
            end
        end
    else
        # For open curves, similar approach but keep endpoints fixed
        cumulative_distance = [0.0]
        for i in 1:n_markers-1
            push!(cumulative_distance, cumulative_distance[end] + norm(positions[i+1] - positions[i]))
        end
        
        # Total length
        total_length = cumulative_distance[end]
        
        # Always keep the endpoints
        push!(new_positions, positions[1])
        
        # Calculate number of internal points
        n_internal_points = max(2, round(Int, total_length / target) - 1)
        
        # Generate new internal points
        for i in 1:n_internal_points
            t = i * total_length / (n_internal_points + 1)
            
            # Find the segment containing position t
            seg_idx = findlast(d -> d <= t, cumulative_distance)
            
            # Interpolate position
            seg_start = positions[seg_idx]
            seg_end = positions[seg_idx+1]
            seg_length = norm(seg_end - seg_start)
            
            if seg_length > 1e-10
                alpha = (t - cumulative_distance[seg_idx]) / seg_length
                new_pos = seg_start + alpha * (seg_end - seg_start)
                push!(new_positions, new_pos)
            end
        end
        
        # Add the last endpoint
        push!(new_positions, positions[end])
    end
    
    # Create new markers and connectivity
    new_markers = [Marker(pos, i) for (i, pos) in enumerate(new_positions)]
    
    if interface.closed
        new_connectivity = [(i, i % length(new_markers) + 1) for i in 1:length(new_markers)]
    else
        new_connectivity = [(i, i + 1) for i in 1:length(new_markers)-1]
    end
    
    # Update the interface
    interface.markers = new_markers
    interface.connectivity = new_connectivity
    
    return interface
end

# Add exports for the new functions
export plot_interface_evolution, evolve_by_curvature_flow, simulate_interface_evolution, redistribute_markers!

"""
    refine_by_curvature!(interface::Interface; 
                      curvature_threshold=:auto,
                      relative_threshold=1.5, 
                      max_refinement_level=3,
                      min_segment_length=nothing,
                      smooth_refinement=true)

Adaptively refine an interface by adding markers in regions of high curvature.
"""
function refine_by_curvature!(interface::Interface; 
                           curvature_threshold=:auto,
                           relative_threshold=1.5, 
                           max_refinement_level=3,
                           min_segment_length=nothing,
                           smooth_refinement=true)
    
    # Compute curvature for all markers
    curvatures = compute_curvature(interface)
    abs_curvatures = abs.(curvatures)
    
    # Determine the threshold for refinement
    if curvature_threshold == :auto
        # Use relative thresholding based on mean curvature
        mean_abs_curvature = mean(abs_curvatures)
        threshold = mean_abs_curvature * relative_threshold
    else
        threshold = float(curvature_threshold)
    end
    
    # Create a dictionary mapping marker IDs to their curvature
    curvature_map = Dict(marker.id => abs_curvatures[i] for (i, marker) in enumerate(interface.markers))
    
    # Create a mapping from marker IDs to positions
    position_map = Dict(marker.id => marker.position for marker in interface.markers)
    
    # Track which segments have been refined and their level of refinement
    # Use a non-normalized version of segments to track refinement
    refinement_levels = Dict{Tuple{Int, Int}, Int}()
    
    # Initialize refinement levels for all segments (preserve original directions)
    for conn in interface.connectivity
        refinement_levels[conn] = 0
    end
    
    # First, verify all connections actually exist
    connectivity_check = Set([(conn[1], conn[2]) for conn in interface.connectivity])
    connectivity_check_rev = Set([(conn[2], conn[1]) for conn in interface.connectivity])
    
    # Identify high-curvature markers
    high_curv_markers = findall(c -> c >= threshold, abs_curvatures)
    
    # Store segments to refine as actual connectivity tuples to preserve directionality
    segments_to_refine = Tuple{Int, Int}[]
    
    # For each high curvature marker, add its segments to the refinement list
    for i in high_curv_markers
        marker_id = interface.markers[i].id
        neighbor_ids = get_neighbors(interface, marker_id)
        
        for neighbor_id in neighbor_ids
            # Check which direction the connection exists
            if (marker_id, neighbor_id) in connectivity_check
                push!(segments_to_refine, (marker_id, neighbor_id))
            elseif (neighbor_id, marker_id) in connectivity_check
                push!(segments_to_refine, (neighbor_id, marker_id))
            end
        end
    end
    
    # Sort segments by curvature to refine highest curvature regions first
    sort!(segments_to_refine, by=segment -> max(
        get(curvature_map, segment[1], 0.0),
        get(curvature_map, segment[2], 0.0)
    ), rev=true)
    
    # Process each segment for refinement
    i = 1
    while i <= length(segments_to_refine)
        id1, id2 = segments_to_refine[i]
        
        # Skip if already at maximum refinement level
        if get(refinement_levels, (id1, id2), 0) >= max_refinement_level
            i += 1
            continue
        end
        
        # Get positions of the two endpoints
        pos1 = position_map[id1]
        pos2 = position_map[id2]
        
        # Calculate segment length
        segment_length = norm(pos2 - pos1)
        
        # Check minimum segment length if specified
        if min_segment_length !== nothing && segment_length < 2 * min_segment_length
            i += 1
            continue
        end
        
        # Calculate midpoint
        midpoint = 0.5 * (pos1 + pos2)
        
        # Double check that connection exists (for debugging)
        if !((id1, id2) in connectivity_check || (id2, id1) in connectivity_check_rev)
            println("Warning: Connection between $id1 and $id2 not found in connectivity")
            println("Connections: ", interface.connectivity)
            i += 1
            continue
        end
        
        # Add a new marker at the midpoint using the original connection direction
        try
            new_id = add_marker!(interface, midpoint, (id1, id2))
            
            # Update our maps with the new marker
            position_map[new_id] = midpoint
            
            # Estimate curvature at the new point by averaging neighbors
            curvature_map[new_id] = 0.5 * (
                get(curvature_map, id1, 0.0) + 
                get(curvature_map, id2, 0.0)
            )
            
            # Update refinement levels for the new segments
            new_level = get(refinement_levels, (id1, id2), 0) + 1
            
            # Add the new connections to refinement levels (with original directionality)
            # Connections are now (id1,new_id) and (new_id,id2)
            refinement_levels[(id1, new_id)] = new_level
            refinement_levels[(new_id, id2)] = new_level
            
            # Update our connectivity check sets
            push!(connectivity_check, (id1, new_id))
            push!(connectivity_check, (new_id, id2))
            push!(connectivity_check_rev, (new_id, id1))
            push!(connectivity_check_rev, (id2, new_id))
            
            # If the new marker has high curvature and we're not at max level,
            # add new segments for further refinement
            if get(curvature_map, new_id, 0.0) >= threshold && new_level < max_refinement_level
                push!(segments_to_refine, (id1, new_id))
                push!(segments_to_refine, (new_id, id2))
            end
        catch e
            println("Error refining segment ($id1, $id2): $e")
        end
        
        i += 1
    end
    
    # Apply smoothing if requested
    if smooth_refinement
        # Simple Laplacian smoothing for interior points
        new_positions = Dict{Int, Vector{Float64}}()
        
        for marker in interface.markers
            # Only smooth newly added points
            if marker.id > length(curvatures)
                neighbors = get_neighbors(interface, marker.id)
                if length(neighbors) >= 2
                    # Average neighbor positions
                    avg_pos = zeros(2)
                    for neighbor_id in neighbors
                        avg_pos += position_map[neighbor_id]
                    end
                    avg_pos /= length(neighbors)
                    
                    # Blend original and average position
                    new_pos = 0.7 * marker.position + 0.3 * avg_pos
                    new_positions[marker.id] = new_pos
                end
            end
        end
        
        # Update marker positions
        for (id, pos) in new_positions
            idx = findfirst(m -> m.id == id, interface.markers)
            if idx !== nothing
                interface.markers[idx] = Marker(pos, id)
                position_map[id] = pos
            end
        end
    end
    
    return interface
end
