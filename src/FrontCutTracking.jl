module FrontCutTracking

using LinearAlgebra
using Graphs
using MetaGraphs

export Marker, Interface
export create_circle, create_rectangle, create_line
export add_marker!, remove_marker!, redistribute_markers!
export compute_normals, compute_curvature, interface_length
export to_graph, from_graph, get_neighbors

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
    
    # Check if the two markers are connected
    is_connected = (id1, id2) in interface.connectivity || (id2, id1) in interface.connectivity
    if !is_connected
        error("Markers $id1 and $id2 are not connected")
    end
    
    # Create new marker with next available ID
    new_id = maximum([marker.id for marker in interface.markers]) + 1
    new_marker = Marker(position, new_id)
    
    # Add the marker to the list
    push!(interface.markers, new_marker)
    
    # Remove the direct connection between id1 and id2
    interface.connectivity = filter(conn -> conn != (id1, id2) && conn != (id2, id1), interface.connectivity)
    
    # Add connections to the new marker
    push!(interface.connectivity, (id1, new_id))
    push!(interface.connectivity, (new_id, id2))
    
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
    # Implementation will go here
end

"""
    compute_curvature(interface::Interface)

Compute the curvature at each marker of the interface.
"""
function compute_curvature(interface::Interface)
    # Implementation will go here
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

end # module