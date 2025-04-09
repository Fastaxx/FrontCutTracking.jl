
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