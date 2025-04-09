
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
