

"""
    compute_normals(interface::Interface)

Compute the normal vectors at each marker of the interface.
Returns a vector of normal vectors corresponding to each marker.
"""
function compute_normals(interface::Interface)
    n_markers = length(interface.markers)
    normals = Vector{Vector{Float64}}(undef, n_markers)
    
    # Create a dictionary for fast lookup of marker positions by ID
    positions = Dict(marker.id => marker.position for marker in interface.markers)
    
    # For each marker, compute its normal vector
    for (idx, marker) in enumerate(interface.markers)
        # Get the neighbors of the current marker
        neighbor_ids = get_neighbors(interface, marker.id)
        
        if length(neighbor_ids) == 2
            # Regular case: marker has two neighbors
            prev_pos = positions[neighbor_ids[1]]
            next_pos = positions[neighbor_ids[2]]
            
            # Compute tangent as the average direction
            tangent = next_pos - prev_pos
            
            # Rotate 90 degrees to get the normal: [tx, ty] -> [-ty, tx]
            normal = [-tangent[2], tangent[1]]
            
            # For closed curves, ensure the normal points outward
            if interface.closed
                # For a closed curve, the normal should point outward
                # We can check this by comparing with the vector from center to marker
                center = [0.0, 0.0]
                for m in interface.markers
                    center += m.position
                end
                center ./= n_markers
                
                outward = marker.position - center
                
                # If the dot product is negative, flip the normal
                if dot(normal, outward) < 0
                    normal = -normal
                end
            end
        elseif length(neighbor_ids) == 1
            # Endpoint of an open curve
            neighbor_pos = positions[neighbor_ids[1]]
            
            # Find if this is the first or last marker in the line
            is_first_marker = marker.id == interface.markers[1].id
            
            # Compute tangent based on position in the line
            # This ensures consistent normal direction at both endpoints
            if is_first_marker
                # First endpoint - tangent points into the line
                tangent = normalize(neighbor_pos - marker.position)
            else
                # Last endpoint - tangent points out of the line
                tangent = normalize(marker.position - neighbor_pos)
            end
            
            # Rotate 90 degrees to get normal (always to the same side)
            normal = [-tangent[2], tangent[1]]
        else
            # Edge case: isolated marker or more than 2 neighbors
            # Just use a default normal pointing upward
            normal = [0.0, 1.0]
        end
        
        # Normalize to unit length
        if norm(normal) > 0
            normal = normalize(normal)
        end
        
        normals[idx] = normal
    end
    
    return normals
end

"""
    compute_curvature(interface::Interface)

Compute the curvature at each marker of the interface.
Returns a vector of curvature values corresponding to each marker.
Positive curvature indicates bending toward the normal direction.
"""
function compute_curvature(interface::Interface)
    n_markers = length(interface.markers)
    curvatures = Vector{Float64}(undef, n_markers)
    
    # Create a dictionary for fast lookup of marker positions by ID
    positions = Dict(marker.id => marker.position for marker in interface.markers)
    
    # Get normals (needed for sign of curvature)
    normals = compute_normals(interface)
    
    # For each marker, compute the curvature
    for (idx, marker) in enumerate(interface.markers)
        # Get the neighbors of the current marker
        neighbor_ids = get_neighbors(interface, marker.id)
        
        if length(neighbor_ids) == 2
            # Regular case: marker has two neighbors
            prev_pos = positions[neighbor_ids[1]]
            next_pos = positions[neighbor_ids[2]]
            curr_pos = marker.position
            
            # Calculate vectors to neighboring points
            v1 = prev_pos - curr_pos
            v2 = next_pos - curr_pos
            
            # Calculate the Menger curvature using the cross product formula
            # κ = 2 * |v1 × v2| / (|v1| * |v2| * |v1 - v2|)
            
            # For 2D vectors, the cross product magnitude is |v1|*|v2|*sin(θ)
            # For 2D, v1 × v2 = v1[1]*v2[2] - v1[2]*v2[1]
            cross_product = v1[1]*v2[2] - v1[2]*v2[1]
            
            norm_v1 = norm(v1)
            norm_v2 = norm(v2)
            norm_diff = norm(v2 - v1)
            
            if norm_v1 > 1e-10 && norm_v2 > 1e-10 && norm_diff > 1e-10
                unsigned_curvature = 2 * abs(cross_product) / (norm_v1 * norm_v2 * norm_diff)
                
                # Determine sign of curvature using the normal vector
                # Compute bisector of the angle between v1 and v2
                bisector = normalize(normalize(v1) + normalize(v2))
                
                # If bisector points in same general direction as normal, curvature is positive
                curvatures[idx] = dot(bisector, normals[idx]) > 0 ? 
                                 -unsigned_curvature : unsigned_curvature
            else
                curvatures[idx] = 0.0
            end
        elseif length(neighbor_ids) == 1 && !interface.closed
            # Endpoint of an open curve
            
            # For endpoints, we use the curvature of a circle through 
            # this point and two nearby points
            if idx == 1 && n_markers >= 3
                # First point: use first three points
                p1 = marker.position
                p2 = interface.markers[2].position
                p3 = interface.markers[3].position
                curvatures[idx] = circle_curvature(p1, p2, p3)
            elseif idx == n_markers && n_markers >= 3
                # Last point: use last three points
                p1 = interface.markers[n_markers-2].position
                p2 = interface.markers[n_markers-1].position
                p3 = marker.position
                curvatures[idx] = circle_curvature(p1, p2, p3)
            else
                # If not enough points, set curvature to zero
                curvatures[idx] = 0.0
            end
        else
            # Isolated marker or other special case
            curvatures[idx] = 0.0
        end
    end
    
    return curvatures
end

"""
Helper function to compute curvature from three points using circumscribed circle
"""
function circle_curvature(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})
    # Compute side lengths of triangle
    a = norm(p2 - p3)
    b = norm(p1 - p3)
    c = norm(p1 - p2)
    
    # Avoid division by zero
    if a*b*c < 1e-10
        return 0.0
    end
    
    # Compute semiperimeter
    s = (a + b + c) / 2
    
    # Area of triangle using Heron's formula
    area = sqrt(max(0.0, s * (s - a) * (s - b) * (s - c)))
    
    if area < 1e-10
        return 0.0
    end
    
    # Radius of circumscribed circle
    R = (a * b * c) / (4 * area)
    
    # Curvature is 1/R (sign is determined later)
    return 1/R
end