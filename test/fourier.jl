using FrontCutTracking
using CairoMakie

# Create a shape with interesting features for analysis
t = range(0, 2π, length=200)
# Simple flower shape with 5 petals and small noise
positions = [
    [3 * (1 + 0.3 * cos(5θ) + 0.02 * sin(15θ)) * cos(θ), 
     3 * (1 + 0.3 * cos(5θ) + 0.02 * sin(15θ)) * sin(θ)] 
    for θ in t
]
markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
flower = Interface(markers, connectivity, true)

# Visualize the original shape
flower_fig = plot_interface(flower, show_curvature=true, 
                           colormap=:viridis, marker_size=4)
save("flower_shape.png", flower_fig)

# Compute standard Fourier descriptors
descriptors = compute_fourier_descriptors(flower)

# Plot the Fourier spectrum
power_spectrum = fourier_power_spectrum(flower)

spectrum_fig = Figure(size=(900, 400))
ax1 = Axis(spectrum_fig[1, 1], 
         title="Fourier Power Spectrum",
         xlabel="Frequency", 
         ylabel="Power (log scale)",
         yscale=log10)

# Plot power spectrum (use only first 50 frequencies for clarity)
max_freq = 50
hist!(ax1, power_spectrum.frequencies[1:max_freq], 
     power_spectrum.power[1:max_freq])

# Plot cumulative power
ax2 = Axis(spectrum_fig[1, 2], 
         title="Cumulative Power",
         xlabel="Frequency", 
         ylabel="Cumulative Power (%)")

lines!(ax2, power_spectrum.frequencies[1:max_freq], 
      power_spectrum.cumulative_power[1:max_freq] * 100,
      linewidth=2, color=:blue)

# Mark 95% power threshold
threshold_idx = findfirst(x -> x >= 0.95, power_spectrum.cumulative_power)
if !isnothing(threshold_idx) && threshold_idx <= max_freq
    lines!(ax2, [0, threshold_idx], [95, 95], linestyle=:dash, color=:red)
    lines!(ax2, [threshold_idx, threshold_idx], [0, 95], linestyle=:dash, color=:red)
    text!(ax2, "95% at freq $(threshold_idx-1)", 
         position=(threshold_idx + 2, 95), 
         fontsize=12)
end

save("fourier_spectrum.png", spectrum_fig)

# Identify dominant frequencies
dom_freqs = dominant_frequencies(flower, threshold=0.99)
println("Dominant frequencies: $dom_freqs")

# Create reconstructions with different frequency cutoffs
cutoffs = [5, 10, 20, 50]
reconstructions = [fourier_filter(flower, cutoff, filter_type=:lowpass) for cutoff in cutoffs]

# Visualize reconstructions
recon_fig = Figure(size=(1200, 800))

# Original shape
ax_orig = Axis(recon_fig[1, 1], aspect=DataAspect(), title="Original Shape")
pos = [marker.position for marker in flower.markers]
x = [p[1] for p in pos]
y = [p[2] for p in pos]
lines!(ax_orig, x, y, color=:black)
scatter!(ax_orig, x, y, markersize=2, color=:blue)

# Reconstructions with different frequency cutoffs
for (i, (recon, cutoff)) in enumerate(zip(reconstructions, cutoffs))
    row, col = divrem(i, 2) .+ (1, 2)
    ax = Axis(recon_fig[row, col], aspect=DataAspect(), 
             title="Reconstruction (f ≤ $cutoff)")
    
    pos = [marker.position for marker in recon.markers]
    x = [p[1] for p in pos]
    y = [p[2] for p in pos]
    lines!(ax, x, y, color=:black)
    scatter!(ax, x, y, markersize=2, color=:red)
end

#save("fourier_reconstructions.png", recon_fig)

# Elliptic Fourier descriptors example
efd = elliptic_fourier_descriptors(flower, num_harmonics=30)

# Reconstruct with different numbers of harmonics
harmonic_counts = [1, 3, 5, 10, 20]
efd_reconstructions = [
    reconstruct_from_elliptic_fourier(efd, 200, h) 
    for h in harmonic_counts
]

# Visualize elliptic Fourier reconstructions
efd_fig = Figure(size=(1200, 800))

# Original shape
ax_orig = Axis(efd_fig[1, 1], aspect=DataAspect(), title="Original Shape")
pos = [marker.position for marker in flower.markers]
x = [p[1] for p in pos]
y = [p[2] for p in pos]
lines!(ax_orig, x, y, color=:black)
scatter!(ax_orig, x, y, markersize=2, color=:blue)

# Reconstructions with different harmonic counts
for (i, (recon, h)) in enumerate(zip(efd_reconstructions, harmonic_counts))
    row, col = divrem(i, 3) .+ (1, 2)
    ax = Axis(efd_fig[row, col], aspect=DataAspect(), 
             title="EFD Reconstruction ($h harmonics)")
    
    pos = [marker.position for marker in recon.markers]
    x = [p[1] for p in pos]
    y = [p[2] for p in pos]
    lines!(ax, x, y, color=:black)
    scatter!(ax, x, y, markersize=2, color=:red)
end

save("elliptic_fourier_reconstructions.png", efd_fig)

# Shape comparison using Fourier descriptors
# Create three shapes: a circle, a square and a 5-pointed star
circle = create_circle([0.0, 0.0], 3.0, 100)

square_points = []
for side in 1:4
    θ1 = (side - 1) * π/2
    θ2 = side * π/2
    side_points = [
        [3 * cos(θ), 3 * sin(θ)] for θ in range(θ1, θ2, length=25)
    ]
    append!(square_points, side_points)
end
square_markers = [Marker(pos, i) for (i, pos) in enumerate(square_points)]
square_connectivity = [(i, i % length(square_markers) + 1) for i in 1:length(square_markers)]
square = Interface(square_markers, square_connectivity, true)

star_t = range(0, 2π, length=100)
star_points = [
    [3 * (1 + 0.5 * cos(5θ - π/2)) * cos(θ), 
     3 * (1 + 0.5 * cos(5θ - π/2)) * sin(θ)] 
    for θ in star_t
]
star_markers = [Marker(pos, i) for (i, pos) in enumerate(star_points)]
star_connectivity = [(i, i % length(star_markers) + 1) for i in 1:length(star_markers)]
star = Interface(star_markers, star_connectivity, true)

# Compare shapes
circle_square_similarity = fourier_shape_similarity(circle, square)
circle_star_similarity = fourier_shape_similarity(circle, star)
square_star_similarity = fourier_shape_similarity(square, star)

println("Shape similarity comparison (lower = more similar):")
println("  Circle to Square: $circle_square_similarity")
println("  Circle to Star: $circle_star_similarity")
println("  Square to Star: $square_star_similarity")

# Visualize the shapes for comparison
compare_fig = Figure(size=(1200, 400))
ax1 = Axis(compare_fig[1, 1], aspect=DataAspect(), title="Circle")
ax2 = Axis(compare_fig[1, 2], aspect=DataAspect(), title="Square")
ax3 = Axis(compare_fig[1, 3], aspect=DataAspect(), title="Star")

for (ax, shape) in zip([ax1, ax2, ax3], [circle, square, star])
    pos = [marker.position for marker in shape.markers]
    x = [p[1] for p in pos]
    y = [p[2] for p in pos]
    lines!(ax, x, y, color=:black, linewidth=2)
    xlims!(ax, -4, 4)
    ylims!(ax, -4, 4)
end

save("shape_comparison.png", compare_fig)

# High-frequency filtering example (remove noise)
# Create a noisy circle
noisy_t = range(0, 2π, length=200)
noisy_circle_points = [
    [3 * (1 + 0.1 * sin(20θ) * rand()) * cos(θ), 
     3 * (1 + 0.1 * sin(20θ) * rand()) * sin(θ)] 
    for θ in noisy_t
]
noisy_markers = [Marker(pos, i) for (i, pos) in enumerate(noisy_circle_points)]
noisy_connectivity = [(i, i % length(noisy_markers) + 1) for i in 1:length(noisy_markers)]
noisy_circle = Interface(noisy_markers, noisy_connectivity, true)

# Apply low-pass filter to remove high-frequency noise
smoothed_circle = fourier_filter(noisy_circle, 5, filter_type=:lowpass)

# Visualize noise removal
noise_fig = Figure(size=(900, 400))
ax1 = Axis(noise_fig[1, 1], aspect=DataAspect(), title="Noisy Circle")
ax2 = Axis(noise_fig[1, 2], aspect=DataAspect(), title="Smoothed Circle")

for (ax, shape) in zip([ax1, ax2], [noisy_circle, smoothed_circle])
    pos = [marker.position for marker in shape.markers]
    x = [p[1] for p in pos]
    y = [p[2] for p in pos]
    lines!(ax, x, y, color=:black, linewidth=1.5)
    scatter!(ax, x, y, markersize=2, color=:blue)
    xlims!(ax, -4, 4)
    ylims!(ax, -4, 4)
end

save("fourier_smoothing.png", noise_fig)