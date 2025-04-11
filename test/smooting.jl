using FrontCutTracking
using CairoMakie

# ---- SMOOTHING EXAMPLES ----

# Create a noisy interface - a circle with random perturbations
radius = 1.0
num_points = 50
noise_level = 0.1

t = range(0, 2π, length=num_points+1)[1:num_points]
noisy_positions = [[radius * cos(θ) + noise_level * (2rand() - 1), 
                   radius * sin(θ) + noise_level * (2rand() - 1)] for θ in t]
noisy_markers = [Marker(pos, i) for (i, pos) in enumerate(noisy_positions)]
noisy_connectivity = [(i, i % length(noisy_markers) + 1) for i in 1:length(noisy_markers)]
noisy_circle = Interface(noisy_markers, noisy_connectivity, true)

# Plot the original noisy interface
noisy_fig = plot_interface(noisy_circle, show_curvature=true, colormap=:viridis)
save("noisy_circle.png", noisy_fig)

# Apply basic Laplacian smoothing
smoothed_circle1 = deepcopy(noisy_circle)
smooth_interface!(smoothed_circle1, iterations=10, lambda=0.5)
laplacian_fig = plot_interface(smoothed_circle1, show_curvature=true, colormap=:viridis)
save("laplacian_smoothed_circle.png", laplacian_fig)

# Apply curvature-weighted smoothing
smoothed_circle2 = deepcopy(noisy_circle)
smooth_interface_curvature_weighted!(smoothed_circle2, iterations=10, lambda=0.3, curvature_weight=0.7)
curvature_smoothed_fig = plot_interface(smoothed_circle2, show_curvature=true, colormap=:viridis)
save("curvature_weighted_smoothed_circle.png", curvature_smoothed_fig)

# Apply Taubin smoothing (prevents shrinkage)
smoothed_circle3 = deepcopy(noisy_circle)
taubin_smooth_interface!(smoothed_circle3, iterations=10)
taubin_fig = plot_interface(smoothed_circle3, show_curvature=true, colormap=:viridis)
save("taubin_smoothed_circle.png", taubin_fig)

# Compare all methods in a single figure
comparison_fig = Figure(size=(1600, 800))
ax1 = Axis(comparison_fig[1, 1], aspect=DataAspect(), title="Original Noisy Interface")
ax2 = Axis(comparison_fig[1, 2], aspect=DataAspect(), title="Laplacian Smoothing")
ax3 = Axis(comparison_fig[2, 1], aspect=DataAspect(), title="Curvature-Weighted Smoothing")
ax4 = Axis(comparison_fig[2, 2], aspect=DataAspect(), title="Taubin Smoothing")

# Helper function to plot an interface on an axis
function plot_on_axis!(ax, interface)
    positions = [marker.position for marker in interface.markers]
    x = [pos[1] for pos in positions]
    y = [pos[2] for pos in positions]
    
    # Plot the markers
    scatter!(ax, x, y, color=:blue, markersize=5)
    
    # Plot the connections
    for (id1, id2) in interface.connectivity
        pos1 = interface.markers[findfirst(m -> m.id == id1, interface.markers)].position
        pos2 = interface.markers[findfirst(m -> m.id == id2, interface.markers)].position
        lines!(ax, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black)
    end
    
    # Plot the perfect circle as reference
    θ = range(0, 2π, length=100)
    perfect_x = radius * cos.(θ)
    perfect_y = radius * sin.(θ)
    lines!(ax, perfect_x, perfect_y, color=:red, linestyle=:dash, linewidth=1)
    
    # Set consistent limits
    limits = (-1.2, 1.2, -1.2, 1.2)
    xlims!(ax, limits[1], limits[2])
    ylims!(ax, limits[3], limits[4])
end

# Plot all interfaces
plot_on_axis!(ax1, noisy_circle)
plot_on_axis!(ax2, smoothed_circle1)
plot_on_axis!(ax3, smoothed_circle2)
plot_on_axis!(ax4, smoothed_circle3)

save("smoothing_comparison.png", comparison_fig)

# Example with an open curve (sine wave with noise)
x_vals = range(0, 2π, length=30)
noisy_sine_positions = [[x, sin(x) + noise_level * (2rand() - 1)] for x in x_vals]
noisy_sine_markers = [Marker(pos, i) for (i, pos) in enumerate(noisy_sine_positions)]
noisy_sine_connectivity = [(i, i+1) for i in 1:length(noisy_sine_markers)-1]
noisy_sine = Interface(noisy_sine_markers, noisy_sine_connectivity, false)

# Apply different smoothing techniques
smooth_sine1 = deepcopy(noisy_sine)
smooth_sine2 = deepcopy(noisy_sine)
smooth_sine3 = deepcopy(noisy_sine)

# Standard Laplacian with endpoint preservation
smooth_interface!(smooth_sine1, iterations=5, lambda=0.5, preserve_endpoints=true)

# Curvature-weighted with endpoint preservation
smooth_interface_curvature_weighted!(smooth_sine2, iterations=5, lambda=0.3, 
                                   curvature_weight=0.6, preserve_endpoints=true)

# Taubin smoothing with endpoint preservation
taubin_smooth_interface!(smooth_sine3, iterations=10, preserve_endpoints=true)

# Create comparison figure
sine_fig = Figure(size=(1600, 400))
ax1 = Axis(sine_fig[1, 1], aspect=DataAspect(), title="Noisy Sine Wave")
ax2 = Axis(sine_fig[1, 2], aspect=DataAspect(), title="Laplacian Smoothing")
ax3 = Axis(sine_fig[1, 3], aspect=DataAspect(), title="Curvature-Weighted")
ax4 = Axis(sine_fig[1, 4], aspect=DataAspect(), title="Taubin Smoothing")

# Helper to plot sine curve
function plot_sine!(ax, interface)
    positions = [marker.position for marker in interface.markers]
    x = [pos[1] for pos in positions]
    y = [pos[2] for pos in positions]
    
    # Plot the markers
    scatter!(ax, x, y, color=:blue, markersize=5)
    
    # Plot the connections
    for (id1, id2) in interface.connectivity
        pos1 = interface.markers[findfirst(m -> m.id == id1, interface.markers)].position
        pos2 = interface.markers[findfirst(m -> m.id == id2, interface.markers)].position
        lines!(ax, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black)
    end
    
    # Plot perfect sine as reference
    ref_x = range(0, 2π, length=100)
    ref_y = sin.(ref_x)
    lines!(ax, ref_x, ref_y, color=:red, linestyle=:dash, linewidth=1)
    
    # Set consistent limits
    xlims!(ax, 0, 2π)
    ylims!(ax, -1.5, 1.5)
end

plot_sine!(ax1, noisy_sine)
plot_sine!(ax2, smooth_sine1)
plot_sine!(ax3, smooth_sine2)
plot_sine!(ax4, smooth_sine3)

save("sine_smoothing_comparison.png", sine_fig)