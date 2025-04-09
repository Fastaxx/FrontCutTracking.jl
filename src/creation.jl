

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