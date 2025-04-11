using FrontCutTracking
using CairoMakie

# ---- DECIMATION EXAMPLE ----
# Create a shape with dense markers - a circle with 200 markers
dense_circle = create_circle([0.0, 0.0], 3.0, 200)

# Plot the original dense interface
dense_fig = plot_interface(dense_circle, show_ids=false, 
                          marker_size=4, connection_width=1)
save("dense_circle.png", dense_fig)

# Decimate to reduce markers in low-curvature regions
decimated_circle = deepcopy(dense_circle)
decimate_by_curvature!(decimated_circle, 
                     relative_threshold=0.3, 
                     min_markers=30,
                     smooth_after=true)

# Plot the decimated interface
decimated_fig = plot_interface(decimated_circle, show_ids=false,
                              marker_size=4, connection_width=1)
save("decimated_circle.png", decimated_fig)

# Compare the number of markers
println("Dense circle: $(length(dense_circle.markers)) markers")
println("Decimated circle: $(length(decimated_circle.markers)) markers")

# Create comparison figure
compare_fig = Figure(size=(1000, 400))
ax1 = Axis(compare_fig[1, 1], aspect=DataAspect(), title="Original Interface ($(length(dense_circle.markers)) markers)")
ax2 = Axis(compare_fig[1, 2], aspect=DataAspect(), title="Decimated Interface ($(length(decimated_circle.markers)) markers)")

# Plot original
orig_positions = [marker.position for marker in dense_circle.markers]
orig_x = [pos[1] for pos in orig_positions]
orig_y = [pos[2] for pos in orig_positions]
scatter!(ax1, orig_x, orig_y, markersize=4, color=:blue)
for (id1, id2) in dense_circle.connectivity
    pos1 = dense_circle.markers[findfirst(m -> m.id == id1, dense_circle.markers)].position
    pos2 = dense_circle.markers[findfirst(m -> m.id == id2, dense_circle.markers)].position
    lines!(ax1, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black, linewidth=1)
end

# Plot decimated
dec_positions = [marker.position for marker in decimated_circle.markers]
dec_x = [pos[1] for pos in dec_positions]
dec_y = [pos[2] for pos in dec_positions]
scatter!(ax2, dec_x, dec_y, markersize=6, color=:red)
for (id1, id2) in decimated_circle.connectivity
    pos1 = decimated_circle.markers[findfirst(m -> m.id == id1, decimated_circle.markers)].position
    pos2 = decimated_circle.markers[findfirst(m -> m.id == id2, decimated_circle.markers)].position
    lines!(ax2, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black, linewidth=1.5)
end

save("decimation_comparison.png", compare_fig)

# Example with non-uniform curvature - star shape
t = range(0, 2π, length=200)
positions = [[5 * (1 + 0.3 * sin(5θ)) * cos(θ), 5 * (1 + 0.3 * sin(5θ)) * sin(θ)] for θ in t]
markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
dense_star = Interface(markers, connectivity, true)

# Plot original with curvature coloring
original_star_fig = plot_interface(dense_star, show_curvature=true, 
                                 marker_size=4, colormap=:viridis)
save("dense_star.png", original_star_fig)

# Apply adaptive decimation
decimated_star = deepcopy(dense_star)
decimate_by_curvature!(decimated_star, 
                     relative_threshold=0.2,
                     preserve_feature_points=true,
                     smooth_after=true)

# Plot decimated with curvature coloring
decimated_star_fig = plot_interface(decimated_star, show_curvature=true, 
                                  marker_size=6, colormap=:viridis)
save("decimated_star.png", decimated_star_fig)

# Create comparison figure showing how decimation preserves features
compare_star_fig = Figure(size=(1000, 400))
ax1 = Axis(compare_star_fig[1, 1], aspect=DataAspect(), 
          title="Original Star ($(length(dense_star.markers)) markers)")
ax2 = Axis(compare_star_fig[1, 2], aspect=DataAspect(), 
          title="Decimated Star ($(length(decimated_star.markers)) markers)")

# Plot original
"""
Helper function to plot an interface on an axis
"""
function plot_on_axis!(ax, interface; marker_color=:blue, marker_size=5, line_color=:black, line_width=1)
    positions = [marker.position for marker in interface.markers]
    x = [pos[1] for pos in positions]
    y = [pos[2] for pos in positions]
    
    # Plot the markers
    scatter!(ax, x, y, color=marker_color, markersize=marker_size)
    
    # Plot the connections
    for (id1, id2) in interface.connectivity
        pos1 = interface.markers[findfirst(m -> m.id == id1, interface.markers)].position
        pos2 = interface.markers[findfirst(m -> m.id == id2, interface.markers)].position
        lines!(ax, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=line_color, linewidth=line_width)
    end
    
    # Set axis limits based on the data
    padding = 0.1 * (maximum(x) - minimum(x))
    xlims!(ax, minimum(x) - padding, maximum(x) + padding)
    ylims!(ax, minimum(y) - padding, maximum(y) + padding)
end

plot_on_axis!(ax1, dense_star)
plot_on_axis!(ax2, decimated_star)

save("star_decimation_comparison.png", compare_star_fig)