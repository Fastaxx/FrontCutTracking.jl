module FrontCutTracking

using LinearAlgebra
using Graphs
using MetaGraphs
using CairoMakie
using Colors
using Statistics

export Marker, Interface
export create_circle, create_rectangle, create_line
export add_marker!, remove_marker!, redistribute_markers!
export compute_normals, compute_curvature, interface_length
export to_graph, from_graph, get_neighbors
export plot_interface, plot_graph, plot_interface_evolution
export curvature_statistics, mean_curvature, curvature_extrema, curvature_histogram

"""
    Marker(position, id)

A marker on the interface with position (a 2D vector [x,y]) and unique identifier.
"""
struct Marker
    position::Vector{Float64}
    id::Int
    
    # Constructor that enforces 2D coordinates
    function Marker(position::Vector{Float64}, id::Int)
        if length(position) != 2
            throw(ArgumentError("Position must be a 2D vector"))
        end
        new(position, id)
    end
end

"""
    Interface

A 2D interface represented by markers and their connectivity.
"""
mutable struct Interface
    markers::Vector{Marker}
    connectivity::Vector{Tuple{Int,Int}}  # Pairs of marker indices that are connected
    closed::Bool  # Whether the interface is closed (e.g., a circle) or open (e.g., a line)
end

"""
    to_graph(interface::Interface)

Convert an Interface to a MetaGraph where vertices represent markers and edges represent connections.
"""
function to_graph(interface::Interface)
    # Create a graph with vertices equal to number of markers
    g = MetaGraph()
    
    # Add vertices with marker data
    for marker in interface.markers
        add_vertex!(g)
        set_prop!(g, marker.id, :position, marker.position)
        set_prop!(g, marker.id, :id, marker.id)
    end
    
    # Add edges based on connectivity
    for (i, j) in interface.connectivity
        add_edge!(g, i, j)
    end
    
    return g
end

"""
    from_graph(g::AbstractMetaGraph, closed::Bool=false)

Convert a MetaGraph back to an Interface structure.
"""
function from_graph(g::AbstractMetaGraph, closed::Bool=false)
    # Extract markers from graph vertices
    markers = Marker[]
    
    for v in vertices(g)
        if has_prop(g, v, :position) && has_prop(g, v, :id)
            push!(markers, Marker(get_prop(g, v, :position), get_prop(g, v, :id)))
        else
            error("Graph vertex $v is missing required properties")
        end
    end
    
    # Extract connectivity from edges
    connectivity = Tuple{Int,Int}[]
    
    for e in edges(g)
        push!(connectivity, (e.src, e.dst))
    end
    
    return Interface(markers, connectivity, closed)
end

"""
    get_neighbors(interface::Interface, id::Int)

Get the IDs of markers connected to the marker with the given ID.
"""
function get_neighbors(interface::Interface, id::Int)
    neighbors = Int[]
    
    for (i, j) in interface.connectivity
        if i == id
            push!(neighbors, j)
        elseif j == id
            push!(neighbors, i)
        end
    end
    
    return neighbors
end

"""
    create_circle(center::Vector{Float64}, radius::Float64, num_markers::Int)

Create a circular interface with given center, radius and number of markers.
"""
function create_circle(center::Vector{Float64}, radius::Float64, num_markers::Int)
    markers = Vector{Marker}(undef, num_markers)
    
    for i in 1:num_markers
        θ = 2π * (i - 1) / num_markers
        x = center[1] + radius * cos(θ)
        y = center[2] + radius * sin(θ)
        markers[i] = Marker([x, y], i)
    end
    
    # Create connectivity (each marker is connected to the next one, and the last is connected to the first)
    connectivity = [(i, i % num_markers + 1) for i in 1:num_markers]
    
    return Interface(markers, connectivity, true)
end

"""
    create_rectangle(x_min::Float64, y_min::Float64, x_max::Float64, y_max::Float64, markers_per_side::Int)

Create a rectangular interface with specified boundaries and number of markers per side.
"""
function create_rectangle(x_min::Float64, y_min::Float64, x_max::Float64, y_max::Float64, markers_per_side::Int)
    # Total number of markers (4 sides)
    total_markers = 4 * markers_per_side
    markers = Vector{Marker}(undef, total_markers)
    
    # Place markers at corners and along sides
    if markers_per_side == 1
        # Special case: just place markers at the corners
        markers[1] = Marker([x_min, y_min], 1) # bottom-left
        markers[2] = Marker([x_max, y_min], 2) # bottom-right
        markers[3] = Marker([x_max, y_max], 3) # top-right
        markers[4] = Marker([x_min, y_max], 4) # top-left
    else
        # Calculate the step sizes for each side
        dx = (x_max - x_min) / markers_per_side
        dy = (y_max - y_min) / markers_per_side
        
        # Create markers for the bottom side (left to right)
        for i in 1:markers_per_side
            x = x_min + (i - 1) * dx
            markers[i] = Marker([x, y_min], i)
        end
        
        # Create markers for the right side (bottom to top)
        for i in 1:markers_per_side
            y = y_min + (i - 1) * dy
            markers[markers_per_side + i] = Marker([x_max, y], markers_per_side + i)
        end
        
        # Create markers for the top side (right to left)
        for i in 1:markers_per_side
            x = x_max - (i - 1) * dx
            markers[2 * markers_per_side + i] = Marker([x, y_max], 2 * markers_per_side + i)
        end
        
        # Create markers for the left side (top to bottom)
        for i in 1:markers_per_side
            y = y_max - (i - 1) * dy
            markers[3 * markers_per_side + i] = Marker([x_min, y], 3 * markers_per_side + i)
        end
    end
    
    # Create connectivity
    connectivity = [(i, i % total_markers + 1) for i in 1:total_markers]
    
    return Interface(markers, connectivity, true)
end

"""
    create_line(start::Vector{Float64}, finish::Vector{Float64}, num_markers::Int)

Create a line interface from start to finish with specified number of markers.
"""
function create_line(start::Vector{Float64}, finish::Vector{Float64}, num_markers::Int)
    markers = Vector{Marker}(undef, num_markers)
    
    for i in 1:num_markers
        t = (i - 1) / (num_markers - 1)
        x = start[1] + t * (finish[1] - start[1])
        y = start[2] + t * (finish[2] - start[2])
        markers[i] = Marker([x, y], i)
    end
    
    # Create connectivity (each marker is connected to the next one)
    connectivity = [(i, i + 1) for i in 1:(num_markers-1)]
    
    return Interface(markers, connectivity, false)
end

"""
    add_marker!(interface::Interface, position::Vector{Float64}, between_ids::Tuple{Int,Int})

Add a new marker at the specified position between two connected markers.
"""
function add_marker!(interface::Interface, position::Vector{Float64}, between_ids::Tuple{Int,Int})
    id1, id2 = between_ids
    
    # Check if the two markers are connected (directly)
    is_connected = false
    connection_found = nothing
    
    for conn in interface.connectivity
        if (conn[1] == id1 && conn[2] == id2) || (conn[1] == id2 && conn[2] == id1)
            is_connected = true
            connection_found = conn
            break
        end
    end
    
    if !is_connected
        # Debug information
        error("Markers $id1 and $id2 are not connected")
    end
    
    # Create new marker with next available ID
    new_id = maximum([marker.id for marker in interface.markers]) + 1
    new_marker = Marker(position, new_id)
    
    # Add the marker to the list
    push!(interface.markers, new_marker)
    
    # Remove the direct connection that was found
    interface.connectivity = filter(conn -> conn != connection_found, interface.connectivity)
    
    # Add connections in the proper direction (preserve the original orientation)
    if connection_found[1] == id1 && connection_found[2] == id2
        push!(interface.connectivity, (id1, new_id))
        push!(interface.connectivity, (new_id, id2))
    else
        push!(interface.connectivity, (id2, new_id))
        push!(interface.connectivity, (new_id, id1))
    end
    
    return new_id
end

"""
    remove_marker!(interface::Interface, id::Int)

Remove the marker with the specified ID and update connectivity.
"""
function remove_marker!(interface::Interface, id::Int)
    # Get neighbors directly before removing
    neighbors_to_reconnect = get_neighbors(interface, id)
    
    # Find the actual marker to remove
    marker_index = findfirst(m -> m.id == id, interface.markers)
    if isnothing(marker_index)
        error("Marker with ID $id not found")
    end
    
    # Remove the marker
    deleteat!(interface.markers, marker_index)
    
    # Remove all connections to/from this marker
    filter!(conn -> conn[1] != id && conn[2] != id, interface.connectivity)
    
    # If it had exactly two neighbors, reconnect them
    if length(neighbors_to_reconnect) == 2
        push!(interface.connectivity, (neighbors_to_reconnect[1], neighbors_to_reconnect[2]))
    end
    
    return nothing
end

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

"""
    curvature_statistics(interface::Interface; absolute=false)

Calculate statistical measures of interface curvature.

# Arguments
- `interface`: The Interface object
- `absolute`: Whether to use absolute curvature values

# Returns
A named tuple with the following fields:
- `minimum`: Minimum curvature
- `maximum`: Maximum curvature
- `mean`: Mean curvature
- `median`: Median curvature
- `std`: Standard deviation of curvature
- `total`: Total curvature (integral along interface)
- `rms`: Root mean square curvature
"""
function curvature_statistics(interface::Interface; absolute=false)
    curvatures = compute_curvature(interface)
    
    if absolute
        curvatures = abs.(curvatures)
    end
    
    # Calculate segment lengths for weighted statistics
    positions = [marker.position for marker in interface.markers]
    n = length(positions)
    
    # For closed interfaces, connect the last point to the first
    if interface.closed
        segment_lengths = [norm(positions[i % n + 1] - positions[i]) for i in 1:n]
    else
        # For open interfaces, use half-segment lengths for endpoints
        segment_lengths = zeros(n)
        for i in 2:(n-1)
            # For interior points, use average of adjacent segment lengths
            segment_lengths[i] = (norm(positions[i] - positions[i-1]) + 
                                norm(positions[i+1] - positions[i])) / 2
        end
        # For endpoints, use their single adjacent segment length
        if n > 1
            segment_lengths[1] = norm(positions[2] - positions[1]) / 2
            segment_lengths[n] = norm(positions[n] - positions[n-1]) / 2
        end
    end
    
    # Make sure segment_lengths is not zero to avoid division issues
    if sum(segment_lengths) < 1e-10
        segment_lengths = ones(n)
    end
    
    # Normalize weights
    weights = segment_lengths ./ sum(segment_lengths)
    
    # Calculate statistics
    min_curv = minimum(curvatures)
    max_curv = maximum(curvatures)
    mean_curv = sum(curvatures .* weights)
    
    # Sort for median calculation
    sorted_indices = sortperm(curvatures)
    sorted_weights = weights[sorted_indices]
    
    # Calculate weighted median
    cumulative_weight = 0.0
    median_idx = 1
    for (i, w) in enumerate(sorted_weights)
        cumulative_weight += w
        if cumulative_weight >= 0.5
            median_idx = i
            break
        end
    end
    median_curv = curvatures[sorted_indices[median_idx]]
    
    # Calculate weighted standard deviation
    variance = sum(weights .* (curvatures .- mean_curv).^2)
    std_curv = sqrt(variance)
    
    # Calculate total curvature (integral along the interface)
    total_curv = sum(curvatures .* segment_lengths)
    
    # Calculate RMS curvature
    rms_curv = sqrt(sum(weights .* curvatures.^2))
    
    return (
        minimum = min_curv,
        maximum = max_curv,
        mean = mean_curv,
        median = median_curv,
        std = std_curv,
        total = total_curv,
        rms = rms_curv
    )
end

"""
    mean_curvature(interface::Interface; absolute=false)

Compute the mean curvature of an interface.

# Arguments
- `interface`: The Interface object
- `absolute`: Whether to use absolute curvature values

# Returns
The mean curvature value
"""
function mean_curvature(interface::Interface; absolute=false)
    stats = curvature_statistics(interface, absolute=absolute)
    return stats.mean
end

"""
    curvature_extrema(interface::Interface; n=1, min_separation=0)

Find the locations of maximum and minimum curvature along the interface.

# Arguments
- `interface`: The Interface object
- `n`: Number of extrema to find (of each type)
- `min_separation`: Minimum number of markers that must separate extrema points

# Returns
A named tuple with:
- `maxima`: Vector of (index, value, position) tuples for maximum curvature points
- `minima`: Vector of (index, value, position) tuples for minimum curvature points
"""
function curvature_extrema(interface::Interface; n=1, min_separation=0)
    curvatures = compute_curvature(interface)
    positions = [marker.position for marker in interface.markers]
    
    # Find local maxima and minima
    n_markers = length(curvatures)
    is_maximum = falses(n_markers)
    is_minimum = falses(n_markers)
    
    # Handle closed interfaces specially
    if interface.closed
        for i in 1:n_markers
            prev_idx = i == 1 ? n_markers : i - 1
            next_idx = i == n_markers ? 1 : i + 1
            
            if curvatures[i] > curvatures[prev_idx] && curvatures[i] > curvatures[next_idx]
                is_maximum[i] = true
            elseif curvatures[i] < curvatures[prev_idx] && curvatures[i] < curvatures[next_idx]
                is_minimum[i] = true
            end
        end
    else
        # For open interfaces, skip endpoints
        for i in 2:(n_markers-1)
            if curvatures[i] > curvatures[i-1] && curvatures[i] > curvatures[i+1]
                is_maximum[i] = true
            elseif curvatures[i] < curvatures[i-1] && curvatures[i] < curvatures[i+1]
                is_minimum[i] = true
            end
        end
    end
    
    # Get indices of maxima and minima
    max_indices = findall(is_maximum)
    min_indices = findall(is_minimum)
    
    # Sort by curvature value
    sort!(max_indices, by=i -> -curvatures[i])  # Descending order for maxima
    sort!(min_indices, by=i -> curvatures[i])   # Ascending order for minima
    
    # Apply minimum separation if needed
    if min_separation > 0
        max_indices = filter_by_separation(max_indices, min_separation, n_markers, interface.closed)
        min_indices = filter_by_separation(min_indices, min_separation, n_markers, interface.closed)
    end
    
    # Take the top n
    max_indices = max_indices[1:min(n, length(max_indices))]
    min_indices = min_indices[1:min(n, length(min_indices))]
    
    # Prepare return values
    maxima = [(i, curvatures[i], positions[i]) for i in max_indices]
    minima = [(i, curvatures[i], positions[i]) for i in min_indices]
    
    return (maxima=maxima, minima=minima)
end

"""
Helper function to filter points by minimum separation
"""
function filter_by_separation(indices, min_separation, n_markers, is_closed)
    if isempty(indices) || min_separation <= 0
        return indices
    end
    
    filtered = [indices[1]]
    
    for idx in indices[2:end]
        # Check if this index is far enough from all already selected
        is_far_enough = true
        
        for selected in filtered
            # Calculate separation, accounting for closed interfaces
            if is_closed
                separation = min(
                    mod(idx - selected, n_markers),
                    mod(selected - idx, n_markers)
                )
            else
                separation = abs(idx - selected)
            end
            
            if separation < min_separation
                is_far_enough = false
                break
            end
        end
        
        if is_far_enough
            push!(filtered, idx)
        end
    end
    
    return filtered
end

"""
    curvature_histogram(interface::Interface; bins=10, range=:auto, absolute=false)

Compute a histogram of curvature values along the interface.

# Arguments
- `interface`: The Interface object
- `bins`: Number of histogram bins
- `range`: Range for histogram bins, :auto or tuple (min, max)
- `absolute`: Whether to use absolute curvature values

# Returns
A named tuple with:
- `edges`: Bin edges
- `counts`: Bin counts
- `weights`: Length-weighted bin counts (normalized)
"""
function curvature_histogram(interface::Interface; bins=10, range=:auto, absolute=false)
    curvatures = compute_curvature(interface)
    
    if absolute
        curvatures = abs.(curvatures)
    end
    
    # Determine range if auto
    if range === :auto
        min_val, max_val = extrema(curvatures)
        range_padding = 0.05 * (max_val - min_val)
        range = (min_val - range_padding, max_val + range_padding)
    end
    
    # Create bin edges
    bin_width = (range[2] - range[1]) / bins
    edges = [range[1] + i * bin_width for i in 0:bins]
    
    # Calculate segment lengths for weighting
    positions = [marker.position for marker in interface.markers]
    n = length(positions)
    
    if interface.closed
        segment_lengths = [norm(positions[i % n + 1] - positions[i]) for i in 1:n]
    else
        segment_lengths = zeros(n)
        for i in 2:n
            segment_lengths[i-1] += norm(positions[i] - positions[i-1]) / 2
            segment_lengths[i] += norm(positions[i] - positions[i-1]) / 2
        end
    end
    
    # Initialize counts
    counts = zeros(Int, bins)
    weighted_counts = zeros(Float64, bins)
    
    # Fill histogram
    for (i, k) in enumerate(curvatures)
        # Find bin index
        bin_idx = min(bins, max(1, Int(floor((k - range[1]) / bin_width) + 1)))
        counts[bin_idx] += 1
        weighted_counts[bin_idx] += segment_lengths[i]
    end
    
    # Normalize weighted counts
    if sum(segment_lengths) > 0
        weighted_counts ./= sum(segment_lengths)
    end
    
    return (edges=edges, counts=counts, weights=weighted_counts)
end

"""
    interface_length(interface::Interface)

Compute the total length of the interface.
"""
function interface_length(interface::Interface)
    total_length = 0.0
    
    # Find positions of markers by ID
    positions = Dict(marker.id => marker.position for marker in interface.markers)
    
    # Sum the lengths of all segments
    for (id1, id2) in interface.connectivity
        segment_length = norm(positions[id1] - positions[id2])
        total_length += segment_length
    end
    
    return total_length
end

"""
    plot_interface(interface::Interface; 
                  show_ids=false, 
                  marker_size=10, 
                  connection_width=1.5, 
                  marker_color=:blue, 
                  connection_color=:black, 
                  show_normals=false,
                  normal_scale=0.1,
                  normal_color=:red,
                  figure_size=(800, 600),
                  axis_aspect=DataAspect(),
                  show_curvature=false,
                  colormap=:viridis,
                  show_colorbar=true,
                  curvature_scale=:auto)

Plot an interface showing markers and connections.

# Arguments
- `interface`: The Interface object to plot
- `show_ids`: Whether to show marker IDs
- `marker_size`: Size of markers
- `connection_width`: Width of connection lines
- `marker_color`: Color of markers (can be a single color or vector of colors)
- `connection_color`: Color of connection lines
- `show_normals`: Whether to show normal vectors
- `normal_scale`: Scale factor for normal vectors
- `normal_color`: Color of normal vectors
- `figure_size`: Size of the figure in pixels
- `axis_aspect`: Aspect ratio for the axis
- `show_curvature`: Whether to color markers by curvature
- `colormap`: Colormap to use for curvature visualization
- `show_colorbar`: Whether to show a colorbar
- `curvature_scale`: Scale for curvature values (:auto or tuple of (min, max))
"""
function plot_interface(interface::Interface; 
                       show_ids=false, 
                       marker_size=10, 
                       connection_width=1.5, 
                       marker_color=:blue, 
                       connection_color=:black, 
                       show_normals=false,
                       normal_scale=0.1,
                       normal_color=:red,
                       figure_size=(800, 600),
                       axis_aspect=DataAspect(),
                       show_curvature=false,
                       colormap=:viridis,
                       show_colorbar=true,
                       curvature_scale=:auto)
    
    # Create figure
    fig = Figure(size=figure_size)
    ax = Axis(fig[1, 1], aspect=axis_aspect)
    
    # Extract marker positions for plotting
    marker_positions = [marker.position for marker in interface.markers]
    x = [pos[1] for pos in marker_positions]
    y = [pos[2] for pos in marker_positions]
    
    # Handle curvature coloring
    if show_curvature
        curvatures = compute_curvature(interface)
        
        # Determine curvature scale
        if curvature_scale == :auto
            # Use symmetric limits for better visualization
            max_abs_curv = maximum(abs.(curvatures))
            curv_limits = (-max_abs_curv, max_abs_curv)
        else
            curv_limits = curvature_scale
        end
        
        # Create colormap
        cmap = cgrad(colormap)
        
        # Normalize curvature values to [0, 1] for colormap
        norm_curvatures = (clamp.(curvatures, curv_limits[1], curv_limits[2]) .- curv_limits[1]) / 
                        (curv_limits[2] - curv_limits[1])
        
        # Map normalized curvatures to colors
        marker_color = [cmap[nc] for nc in norm_curvatures]
        
        # Add colorbar if requested
        if show_colorbar
            Colorbar(fig[1, 2], colormap=cmap, limits=curv_limits, 
                    label="Curvature", height=Relative(0.8))
        end
    end
    
    # Plot markers
    scatter!(ax, x, y, color=marker_color, markersize=marker_size)
    
    # Plot connections
    positions_by_id = Dict(marker.id => marker.position for marker in interface.markers)
    
    for (id1, id2) in interface.connectivity
        pos1 = positions_by_id[id1]
        pos2 = positions_by_id[id2]
        lines!(ax, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=connection_color, linewidth=connection_width)
    end
    
    # Show marker IDs if requested
    if show_ids
        for marker in interface.markers
            text!(ax, marker.position[1], marker.position[2], 
                  text=string(marker.id), align=(:center, :center), 
                  offset=(0, marker_size/2 + 5), 
                  color=:black, fontsize=12)
        end
    end
    
    # Show normal vectors if requested
    if show_normals
        normals = compute_normals(interface)
        for (i, marker) in enumerate(interface.markers)
            pos = marker.position
            normal = normals[i]
            normal_end = pos + normal_scale * normal
            arrows!(ax, [pos[1]], [pos[2]], 
                    [normal_end[1] - pos[1]], [normal_end[2] - pos[2]], 
                    color=normal_color, arrowsize=10)
        end
    end
    
    # Set axis limits automatically to fit all points with some margin
    min_x, max_x = minimum(x), maximum(x)
    min_y, max_y = minimum(y), maximum(y)
    range_x = max_x - min_x
    range_y = max_y - min_y
    margin = 0.1
    
    xlims!(ax, min_x - margin * range_x, max_x + margin * range_x)
    ylims!(ax, min_y - margin * range_y, max_y + margin * range_y)
    
    # Set labels and title
    title_text = show_curvature ? "Interface Curvature Visualization" : "Interface Visualization"
    ax.title = title_text
    ax.xlabel = "x"
    ax.ylabel = "y"
    
    return fig
end

"""
    plot_graph(interface::Interface; 
              layout=:spring, 
              node_size=30, 
              show_ids=true,
              node_color=:lightblue, 
              edge_color=:gray,
              figure_size=(800, 600))

Plot the graph representation of an interface.

# Arguments
- `interface`: The Interface object to plot
- `layout`: Graph layout algorithm (:spring, :spectral, or :circular)
- `node_size`: Size of nodes
- `show_ids`: Whether to show node IDs
- `node_color`: Color of nodes
- `edge_color`: Color of edges
- `figure_size`: Size of the figure in pixels
"""
function plot_graph(interface::Interface; 
                   layout=:spring, 
                   node_size=30, 
                   show_ids=true,
                   node_color=:lightblue, 
                   edge_color=:gray,
                   figure_size=(800, 600))
    
    # Create a graph
    g = to_graph(interface)
    
    # Create figure
    fig = Figure(size=figure_size)
    ax = Axis(fig[1, 1])
    
    # Determine layout
    if layout == :spring
        positions = spring_layout(g)
    elseif layout == :spectral
        positions = spectral_layout(g)
    elseif layout == :circular
        positions = circular_layout(g)
    elseif layout == :original
        # Use original positions from interface
        positions = Dict(v => get_prop(g, v, :position) for v in vertices(g))
    else
        error("Unsupported layout: $layout")
    end
    
    # Get x, y coordinates
    xs = [positions[v][1] for v in vertices(g)]
    ys = [positions[v][2] for v in vertices(g)]
    
    # Draw edges
    for e in edges(g)
        lines!(ax, [positions[e.src][1], positions[e.dst][1]], 
              [positions[e.src][2], positions[e.dst][2]], 
              color=edge_color, linewidth=2)
    end
    
    # Draw vertices
    scatter!(ax, xs, ys, color=node_color, markersize=node_size)
    
    # Add IDs if requested
    if show_ids
        for v in vertices(g)
            text!(ax, positions[v][1], positions[v][2], 
                 text=string(v), align=(:center, :center), 
                 color=:black, fontsize=12)
        end
    end
    
    # Remove axis decorations for cleaner look
    hidedecorations!(ax)
    hidespines!(ax)
    
    ax.title = "Graph Representation"
    
    return fig
end

"""
    spring_layout(g)

Simple spring layout algorithm for graphs.
"""
function spring_layout(g)
    positions = Dict()
    n = nv(g)
    
    # Initialize with random positions
    for v in vertices(g)
        positions[v] = [rand(), rand()]
    end
    
    # Run force-directed algorithm (simplified version)
    k = 0.2  # optimal distance
    iterations = 50
    
    for _ in 1:iterations
        # Calculate forces
        forces = Dict(v => [0.0, 0.0] for v in vertices(g))
        
        # Repulsive forces between all vertices
        for u in vertices(g)
            for v in vertices(g)
                if u != v
                    d = positions[v] - positions[u]
                    distance = norm(d)
                    if distance > 0
                        f = (k^2 / distance) * normalize(d)
                        forces[v] += f
                        forces[u] -= f
                    end
                end
            end
        end
        
        # Attractive forces between adjacent vertices
        for e in edges(g)
            d = positions[e.dst] - positions[e.src]
            distance = norm(d)
            if distance > 0
                f = (distance^2 / k) * normalize(d)
                forces[e.src] += f
                forces[e.dst] -= f
            end
        end
        
        # Update positions
        for v in vertices(g)
            positions[v] += 0.1 * forces[v]
        end
    end
    
    return positions
end

"""
    spectral_layout(g)

Simple spectral layout algorithm.
"""
function spectral_layout(g)
    # Use eigenvectors of the Laplacian matrix for layout
    n = nv(g)
    L = laplacian_matrix(g)
    # If available, get the second and third smallest eigenvectors
    if n > 2
        F = eigen(Matrix(L))
        # Use the 2nd and 3rd smallest eigenvectors for x, y coordinates
        x = F.vectors[:, 2]
        y = F.vectors[:, 3]
    else
        # For very small graphs, just use circular layout
        return circular_layout(g)
    end
    
    positions = Dict()
    for (i, v) in enumerate(vertices(g))
        positions[v] = [x[i], y[i]]
    end
    
    return positions
end

"""
    circular_layout(g)

Arrange vertices in a circle.
"""
function circular_layout(g)
    positions = Dict()
    n = nv(g)
    
    # Arrange vertices in a circle
    for (i, v) in enumerate(vertices(g))
        angle = 2π * (i - 1) / n
        positions[v] = [cos(angle), sin(angle)]
    end
    
    return positions
end

"""
    plot_interface_evolution(interfaces::Vector{Interface};
                           show_curvature=true,
                           colormap=:viridis,
                           frame_duration=0.1,
                           show_time=true,
                           show_normals=false,
                           normal_scale=0.1,
                           fixed_limits=false,
                           filename="interface_evolution.mp4",
                           fps=15,
                           marker_size=10,
                           connection_width=1.5)

Create an animation showing the evolution of an interface over time.

# Arguments
- `interfaces`: Vector of Interface objects representing the evolution over time
- `show_curvature`: Whether to color markers by curvature
- `colormap`: Colormap to use for curvature visualization
- `frame_duration`: Duration of each frame in seconds
- `show_time`: Whether to show time step counter
- `show_normals`: Whether to display normal vectors
- `normal_scale`: Scale factor for normal vectors
- `fixed_limits`: Whether to keep axis limits fixed across frames
- `filename`: Output file name for saving the animation
- `fps`: Frames per second for the saved animation
- `marker_size`: Size of interface markers
- `connection_width`: Width of connection lines

# Returns
A Figure object containing the animation
"""
function plot_interface_evolution(interfaces::Vector{Interface};
                                 show_curvature=true,
                                 colormap=:viridis,
                                 frame_duration=0.1,
                                 show_time=true,
                                 show_normals=false,
                                 normal_scale=0.1,
                                 fixed_limits=false,
                                 filename="interface_evolution.mp4",
                                 fps=15,
                                 marker_size=10,
                                 connection_width=1.5)
    
    # Create figure with proper size
    fig = Figure(size=(800, 600))
    
    # Create a layout that uses the full figure width
    if show_curvature
        # Use entire figure for the GridLayout
        ax = Axis(fig[1, 1], aspect=DataAspect())
        
        # Find global min/max curvature across all frames for consistent coloring
        all_curvatures = Float64[]
        for interface in interfaces
            append!(all_curvatures, compute_curvature(interface))
        end
        
        # Use symmetric limits for better visualization
        max_abs_curv = maximum(abs.(all_curvatures))
        curv_limits = (-max_abs_curv, max_abs_curv)
        
        # Create colormap
        cmap = cgrad(colormap)
        
        # Add colorbar with appropriate sizing - to position [1, 2]
        cbar = Colorbar(fig[1, 2], colormap=cmap, limits=curv_limits, 
                      label="Curvature", height=Relative(0.6))
        
        # Set proper column width ratios for the whole figure
        colsize!(fig.layout, 1, Relative(0.9))  # Plot gets 90% of width
        colsize!(fig.layout, 2, Relative(0.1))  # Colorbar gets 10% of width
    else
        ax = Axis(fig[1, 1], aspect=DataAspect())
    end
    
    # Set up time display if requested
    time_text = show_time ? Observable("Time: 0") : nothing
    if show_time
        Label(fig[0, 1:2], time_text[])  # Span across columns
    end
    
    # Find global axis limits if fixed_limits is true
    if fixed_limits
        x_min, x_max = Inf, -Inf
        y_min, y_max = Inf, -Inf
        
        for interface in interfaces
            positions = [marker.position for marker in interface.markers]
            x_vals = [pos[1] for pos in positions]
            y_vals = [pos[2] for pos in positions]
            
            x_min = min(x_min, minimum(x_vals))
            x_max = max(x_max, maximum(x_vals))
            y_min = min(y_min, minimum(y_vals))
            y_max = max(y_max, maximum(y_vals))
        end
        
        # Add some padding
        padding_x = 0.1 * (x_max - x_min)
        padding_y = 0.1 * (y_max - y_min)
        
        # Set fixed limits
        xlims!(ax, x_min - padding_x, x_max + padding_x)
        ylims!(ax, y_min - padding_y, y_max + padding_y)
    end
    
    # Create a record of the animation
    record(fig, filename, 1:length(interfaces); framerate=fps) do frame_idx
        # Clear previous frame
        empty!(ax)
        
        # Get the current interface
        interface = interfaces[frame_idx]
        
        # Update time text if showing time
        if show_time
            time_text[] = "Time: $(frame_idx-1)"
        end
        
        # Plot the current interface
        marker_positions = [marker.position for marker in interface.markers]
        x = [pos[1] for pos in marker_positions]
        y = [pos[2] for pos in marker_positions]
        
        # Handle curvature coloring
        marker_colors = :blue
        if show_curvature
            curvatures = compute_curvature(interface)
            
            # Normalize curvature values to [0, 1] for colormap
            norm_curvatures = (clamp.(curvatures, curv_limits[1], curv_limits[2]) .- curv_limits[1]) / 
                            (curv_limits[2] - curv_limits[1])
            
            # Map normalized curvatures to colors
            marker_colors = [cmap[nc] for nc in norm_curvatures]
        end
        
        # Plot markers
        scatter!(ax, x, y, color=marker_colors, markersize=marker_size)
        
        # Plot connections
        positions_by_id = Dict(marker.id => marker.position for marker in interface.markers)
        
        for (id1, id2) in interface.connectivity
            pos1 = positions_by_id[id1]
            pos2 = positions_by_id[id2]
            lines!(ax, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black, linewidth=connection_width)
        end
        
        # Show normal vectors if requested
        if show_normals
            normals = compute_normals(interface)
            for (i, marker) in enumerate(interface.markers)
                pos = marker.position
                normal = normals[i]
                normal_end = pos + normal_scale * normal
                arrows!(ax, [pos[1]], [pos[2]], 
                        [normal_end[1] - pos[1]], [normal_end[2] - pos[2]], 
                        color=:red, arrowsize=10)
            end
        end
        
        # Set axis limits if not fixed
        if !fixed_limits
            autolimits!(ax)
        end
    end
    
    return fig
end
"""
    evolve_by_curvature_flow(interface::Interface, dt::Float64; 
                          motion_factor=1.0, 
                          noise=0.0,
                          constrain_motion=false)

Evolve an interface according to curvature flow, where each point moves in the 
direction of its normal proportional to the curvature.

# Arguments
- `interface`: The Interface object
- `dt`: Time step size
- `motion_factor`: Factor to scale the motion (can be negative)
- `noise`: Amount of random noise to add to the motion
- `constrain_motion`: Whether to constrain motion to prevent self-intersections

# Returns
A new Interface object representing the evolved state
"""
function evolve_by_curvature_flow(interface::Interface, dt::Float64;
                               motion_factor=1.0,
                               noise=0.0,
                               constrain_motion=false)
    
    # Copy the interface
    new_markers = Vector{Marker}(undef, length(interface.markers))
    
    # Compute normals and curvatures
    normals = compute_normals(interface)
    curvatures = compute_curvature(interface)
    
    # Move each marker along its normal proportional to curvature
    for (i, marker) in enumerate(interface.markers)
        # Calculate displacement vector
        displacement = normals[i] * curvatures[i] * motion_factor * dt
        
        # Add random noise if requested
        if noise > 0.0
            noise_vector = [2*rand()-1, 2*rand()-1]
            noise_vector = normalize(noise_vector) * noise * dt
            displacement += noise_vector
        end
        
        # Apply displacement to get new position
        new_position = marker.position + displacement
        
        # Create new marker with updated position
        new_markers[i] = Marker(new_position, marker.id)
    end
    
    # If constraining motion, check for self-intersections and adjust if needed
    if constrain_motion
        # Implement collision avoidance here
        # For now, just a placeholder
    end
    
    # Create new interface with updated markers
    return Interface(new_markers, interface.connectivity, interface.closed)
end

"""
    simulate_interface_evolution(initial_interface::Interface, 
                              n_steps::Int, 
                              dt::Float64;
                              evolution_function=evolve_by_curvature_flow,
                              kwargs...)

Simulate the evolution of an interface over time using a specified evolution function.

# Arguments
- `initial_interface`: The starting interface
- `n_steps`: Number of time steps
- `dt`: Time step size
- `evolution_function`: Function used to evolve interface (defaults to curvature flow)
- `kwargs...`: Additional arguments passed to the evolution function

# Returns
A vector of Interface objects representing the evolution over time
"""
function simulate_interface_evolution(initial_interface::Interface,
                                   n_steps::Int,
                                   dt::Float64;
                                   evolution_function=evolve_by_curvature_flow,
                                   kwargs...)
    
    # Initialize result array with the initial interface
    interfaces = Vector{Interface}(undef, n_steps + 1)
    interfaces[1] = initial_interface
    
    # Evolve interface for n_steps
    for i in 1:n_steps
        interfaces[i+1] = evolution_function(interfaces[i], dt; kwargs...)
    end
    
    return interfaces
end

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

# Add the new function to exports
export refine_by_curvature!

end # module