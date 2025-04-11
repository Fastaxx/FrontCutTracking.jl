using FrontCutTracking
using CairoMakie

# ---- LENGTH-SCALE ANALYSIS EXAMPLES ----

# Create a shape with multi-scale features - a star with small perturbations
t = range(0, 2π, length=150)
# Combine large-scale star shape with small-scale noise
positions = [
    [5 * (1 + 0.3 * sin(5θ) + 0.03 * sin(30θ)) * cos(θ), 
     5 * (1 + 0.3 * sin(5θ) + 0.03 * sin(30θ)) * sin(θ)] 
    for θ in t
]
markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
multi_scale_shape = Interface(markers, connectivity, true)

# Visualize the original shape with curvature coloring
multi_scale_fig = plot_interface(multi_scale_shape, show_curvature=true, 
                               colormap=:viridis, marker_size=6)
save("multi_scale_shape.png", multi_scale_fig)

# Compute length scale spectrum
spectrum = length_scale_spectrum(multi_scale_shape, max_scales=20, 
                               feature_measure=:curvature_extrema)

# Plot the length scale spectrum
scale_fig = Figure(size=(800, 500))
ax = Axis(scale_fig[1, 1], 
         title="Length Scale Spectrum",
         xlabel="Scale (Length)", 
         ylabel="Feature Significance", 
         xscale=log10)

# Plot significance curve
lines!(ax, spectrum.scales, spectrum.significance, linewidth=3, color=:blue)
scatter!(ax, spectrum.scales, spectrum.significance, markersize=10, color=:blue)

# Extract dominant scales
dominant_scales = dominant_length_scales(multi_scale_shape, max_scales=15, threshold=0.15)
println("Dominant length scales: $dominant_scales")

# Mark dominant scales on the plot
for scale in dominant_scales
    # Find closest point in the spectrum
    idx = argmin(abs.(spectrum.scales .- scale))
    scatter!(ax, [spectrum.scales[idx]], [spectrum.significance[idx]], 
            markersize=15, color=:red, marker=:star5)
end

save("length_scale_spectrum.png", scale_fig)

# Decompose the shape at multiple scales
decomposition = length_scale_decomposition(multi_scale_shape, n_scales=3)

# Create a visualization showing the different scales
decomp_fig = Figure(size=(1200, 400))
titles = ["Fine Scale", "Medium Scale", "Coarse Scale"]

# Plot each scale component
for (i, interface) in enumerate(decomposition)
    ax = Axis(decomp_fig[1, i], aspect=DataAspect(), title=titles[i])
    
    # Plot the interface
    positions = [marker.position for marker in interface.markers]
    x = [pos[1] for pos in positions]
    y = [pos[2] for pos in positions]
    
    # Plot the markers and connections
    scatter!(ax, x, y, markersize=6, color=:blue)
    
    for (id1, id2) in interface.connectivity
        pos1 = interface.markers[findfirst(m -> m.id == id1, interface.markers)].position
        pos2 = interface.markers[findfirst(m -> m.id == id2, interface.markers)].position
        lines!(ax, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black)
    end
    
    # Set consistent limits across all plots
    xlims!(ax, -7, 7)
    ylims!(ax, -7, 7)
end

save("scale_decomposition.png", decomp_fig)

# Example of filtering at a specific scale
filtered_shape = filter_by_scale(multi_scale_shape, dominant_scales[1])

# Compare original vs filtered (specific scale)
compare_fig = Figure(size=(900, 400))
ax1 = Axis(compare_fig[1, 1], aspect=DataAspect(), title="Original Shape")
ax2 = Axis(compare_fig[1, 2], aspect=DataAspect(), 
          title="Filtered Shape (Scale: $(round(dominant_scales[1], digits=2)))")

# Helper function to plot an interface on an axis
function plot_interface_on_axis!(ax, interface; show_curvature=false)
    positions = [marker.position for marker in interface.markers]
    x = [pos[1] for pos in positions]
    y = [pos[2] for pos in positions]
    
    # Plot markers
    if show_curvature
        curvatures = compute_curvature(interface)
        scatter!(ax, x, y, markersize=6, color=curvatures)
    else
        scatter!(ax, x, y, markersize=6, color=:blue)
    end
    
    # Plot connections
    for (id1, id2) in interface.connectivity
        pos1 = interface.markers[findfirst(m -> m.id == id1, interface.markers)].position
        pos2 = interface.markers[findfirst(m -> m.id == id2, interface.markers)].position
        lines!(ax, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black)
    end
    
    # Set consistent limits
    xlims!(ax, -7, 7)
    ylims!(ax, -7, 7)
end

# Plot original and filtered shapes
plot_interface_on_axis!(ax1, multi_scale_shape, show_curvature=true)
plot_interface_on_axis!(ax2, filtered_shape, show_curvature=true)

save("scale_filtering_comparison.png", compare_fig)