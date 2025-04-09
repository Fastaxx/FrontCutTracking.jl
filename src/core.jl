
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