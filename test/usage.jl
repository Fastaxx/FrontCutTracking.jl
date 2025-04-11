using FrontCutTracking
using CairoMakie

# Create a circular interface
center = [0.0, 0.0]
radius = 1.0
num_markers = 20
interface = create_circle(center, radius, num_markers)

# Basic interface visualization
fig1 = plot_interface(interface)
save("interface.png", fig1)

# Interface with IDs and normals
fig2 = plot_interface(interface, show_ids=true, show_normals=true, 
                     marker_color=:orange, normal_color=:green)
save("interface_with_ids_and_normals.png", fig2)

# Graph representation
fig3 = plot_graph(interface, layout=:spring)
save("interface_graph_spring.png", fig3)

fig4 = plot_graph(interface, layout=:circular, 
                  node_color=:skyblue, edge_color=:navy)
save("interface_graph_circular.png", fig4)

# Evolution visualization (example with a circle expanding over time)
interfaces = [create_circle([0.0, 0.0], r, 20) for r in 1.0:0.5:3.0]
fig5 = plot_interface_evolution(interfaces)
save("interface_evolution.png", fig5)

# Example rectangular interface
x_min, y_min = 0.0, 0.0
x_max, y_max = 2.0, 1.0
markers_per_side = 5
        
rectangle = create_rectangle(x_min, y_min, x_max, y_max, markers_per_side)

# Rectangle with IDs and normals
fig6 = plot_interface(rectangle, show_ids=true, show_normals=true, 
                     marker_color=:orange, normal_color=:green)
save("rectangle_interface_with_ids_and_normals.png", fig6)

# Example Line Interface
start = [0.0, 0.0]
finish = [1.0, 1.0]
num_markers = 5
        
line = create_line(start, finish, num_markers)

# Line with IDs and normals
fig7 = plot_interface(line, show_ids=true, show_normals=true, 
                     marker_color=:orange, normal_color=:green)
save("line_interface_with_ids_and_normals.png", fig7)

# ---- CURVATURE VISUALIZATION EXAMPLES ----

# Create a more complex shape for interesting curvature
t = range(0, 2π, length=50)
markers = [Marker([16*sin(θ)^3, 13*cos(θ) - 5*cos(2θ) - 2*cos(3θ) - cos(4θ)], i) 
        for (i, θ) in enumerate(t)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
heart = Interface(markers, connectivity, true)

# Basic curvature visualization for circle
fig8 = plot_interface(interface, show_curvature=true)
save("circle_curvature.png", fig8)

# Curvature visualization with normals
fig9 = plot_interface(interface, show_curvature=true, show_normals=true, 
                    colormap=:coolwarm, normal_color=:black)
save("circle_curvature_with_normals.png", fig9)

# Curvature visualization for rectangle
fig10 = plot_interface(rectangle, show_curvature=true, colormap=:plasma)
save("rectangle_curvature.png", fig10)

# Curvature visualization for heart shape
fig12 = plot_interface(heart, show_curvature=true, colormap=:viridis, 
                     show_normals=true, normal_scale=0.05, normal_color=:white)
save("heart_curvature.png", fig12)

# Create a curved line instead of straight line for better visualization
t = range(0, π, length=20)
line_positions = [[t[i], sin(t[i])] for i in 1:length(t)]
line_markers = [Marker(pos, i) for (i, pos) in enumerate(line_positions)]
line_connectivity = [(i, i+1) for i in 1:(length(line_positions)-1)]
curved_line = Interface(line_markers, line_connectivity, false)

# Curvature visualization for curved line
fig11 = plot_interface(curved_line, show_curvature=true, colormap=:turbo, 
                     show_normals=true, normal_scale=0.05, normal_color=:black)
save("line_curvature.png", fig11)

# Add examples for curvature statistics

# Calculate curvature statistics for a circle
circle = create_circle([0.0, 0.0], 2.0, 50)
stats = curvature_statistics(circle)
println("Circle curvature statistics:")
println("  Mean: $(stats.mean), should be close to 0.5 (1/radius)")
println("  Min: $(stats.minimum)")
println("  Max: $(stats.maximum)")
println("  Total: $(stats.total), should be close to 2π")

# Find curvature extrema for a more complex shape
t = range(0, 2π, length=100)
# Create a shape with varying curvature (cardioid)
positions = [[2(1+cos(θ))*cos(θ), 2(1+cos(θ))*sin(θ)] for θ in t]
markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
cardioid = Interface(markers, connectivity, true)

# Compute and visualize curvature
fig_cardioid = plot_interface(cardioid, show_curvature=true, colormap=:viridis)
save("cardioid_curvature.png", fig_cardioid)

# Find extrema of curvature
extrema = curvature_extrema(cardioid, n=3, min_separation=5)
println("Cardioid curvature extrema:")
println("  Maxima: $(extrema.maxima)")
println("  Minima: $(extrema.minima)")

# Create curvature histogram
hist = curvature_histogram(cardioid, bins=20)
# Visualize histogram
fig_hist = Figure(size=(600, 400))
ax = Axis(fig_hist[1, 1], xlabel="Curvature", ylabel="Frequency")
barplot!(ax, (hist.edges[1:end-1] .+ hist.edges[2:end]) ./ 2, hist.counts)
save("cardioid_curvature_histogram.png", fig_hist)

# Compare mean curvature across different shapes
println("\nMean curvature comparison:")
println("  Circle (r=2): $(mean_curvature(circle))")
println("  Cardioid: $(mean_curvature(cardioid))")
println("  Rectangle: $(mean_curvature(create_rectangle(0.0, 0.0, 2.0, 1.0, 5)))")

# ---- INTERFACE EVOLUTION EXAMPLES ----
# Create an initial interface - let's use a simple circle
initial_interface = create_circle([0.0, 0.0], 5.0, 50)

# Simulate evolution - circle under inward curvature flow
# This should collapse the circle at a constant rate
interfaces = simulate_interface_evolution(
    initial_interface,
    50,          # 50 time steps 
    0.1,         # dt = 0.1
    motion_factor=-0.5  # negative to move inward
)

# Create the animation
fig = plot_interface_evolution(
    interfaces,
    show_curvature=true,
    show_normals=true,
    normal_scale=0.2,
    fixed_limits=true,
    filename="circle_collapse.mp4",
    fps=10
)

# Let's also create a more complex example - star shape evolving
t = range(0, 2π, length=100)
positions = [[5 * (1 + 0.3 * sin(5θ)) * cos(θ), 5 * (1 + 0.3 * sin(5θ)) * sin(θ)] for θ in t]
markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
star = Interface(markers, connectivity, true)

# Simulate evolution with curvature flow and redistribution
star_interfaces = Vector{Interface}()
push!(star_interfaces, star)
current = star

for i in 1:500
    global current
    # Evolve with curvature flow
    current = evolve_by_curvature_flow(current, 0.01, motion_factor=-1.0)
    
    # Redistribute markers every 10 steps
    if i % 20 == 0
        redistribute_markers!(current)
    end
    
    push!(star_interfaces, current)
end

# Create animation
fig2 = plot_interface_evolution(
    star_interfaces,
    show_curvature=true,
    show_normals=false,
    fixed_limits=true,
    filename="star_evolution.mp4",
    fps=15
)

# ---- CURVATURE REFINEMENT EXAMPLE ----
# Create a shape with varying curvature - a cardioid
t = range(0, 2π, length=40)
cardioid_points = [[2 * (1 - cos(θ)) * cos(θ), 2 * (1 - cos(θ)) * sin(θ)] for θ in t]
markers = [Marker(pos, i) for (i, pos) in enumerate(cardioid_points)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
cardioid = Interface(markers, connectivity, true)

# Plot the original interface with curvature coloring
fig1 = plot_interface(cardioid, show_curvature=true, show_ids=false,
                    marker_size=6, colormap=:viridis)
save("cardioid_original.png", fig1)

fig1b = plot_graph(cardioid, layout=:circular, node_color=:skyblue, edge_color=:navy)
save("cardioid_graph_spring.png", fig1b)

# Adaptive refinement based on curvature
refined_cardioid = deepcopy(cardioid)
refine_by_curvature!(refined_cardioid, relative_threshold=1.2, max_refinement_level=2,smooth_refinement=true)

# Plot the refined interface
fig2 = plot_interface(refined_cardioid, show_curvature=true, show_ids=false,
                    marker_size=6, colormap=:viridis)
save("cardioid_refined.png", fig2)

# Compare the number of markers
println("Original interface: $(length(cardioid.markers)) markers")
println("Refined interface: $(length(refined_cardioid.markers)) markers")

# Create a figure showing both original and refined for comparison
fig3 = Figure(size=(1000, 400))

ax1 = Axis(fig3[1, 1], aspect=DataAspect(), title="Original Interface")
ax2 = Axis(fig3[1, 2], aspect=DataAspect(), title="Refined Interface")

# Plot original
orig_positions = [marker.position for marker in cardioid.markers]
orig_x = [pos[1] for pos in orig_positions]
orig_y = [pos[2] for pos in orig_positions]
scatter!(ax1, orig_x, orig_y, markersize=6, color=:blue)
for (id1, id2) in cardioid.connectivity
    pos1 = cardioid.markers[findfirst(m -> m.id == id1, cardioid.markers)].position
    pos2 = cardioid.markers[findfirst(m -> m.id == id2, cardioid.markers)].position
    lines!(ax1, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black)
end

# Plot refined
ref_positions = [marker.position for marker in refined_cardioid.markers]
ref_x = [pos[1] for pos in ref_positions]
ref_y = [pos[2] for pos in ref_positions]
scatter!(ax2, ref_x, ref_y, markersize=6, color=:red)
for (id1, id2) in refined_cardioid.connectivity
    pos1 = refined_cardioid.markers[findfirst(m -> m.id == id1, refined_cardioid.markers)].position
    pos2 = refined_cardioid.markers[findfirst(m -> m.id == id2, refined_cardioid.markers)].position
    lines!(ax2, [pos1[1], pos2[1]], [pos1[2], pos2[2]], color=:black)
end

save("comparison.png", fig3)

