using FrontCutTracking
using CairoMakie
using LinearAlgebra
using Statistics
using SpecialFunctions  # For computing elliptic integrals

"""
Run convergence tests on scalar geometric properties like length and curvature.
Tests how these properties converge to their exact values as resolution increases.
"""
function run_convergence_tests()
    # Define range of mesh resolutions to test
    resolutions = [10, 20, 40, 80, 160, 320, 640]
    
    # -------------- CIRCLE CONVERGENCE TESTS --------------
    
    # Circle properties - known analytical values
    radius = 2.0
    center = [0.0, 0.0]
    exact_length = 2π * radius
    exact_curvature = 1.0 / radius
    exact_total_curvature = 2π  # Total curvature of any simple closed curve = 2π
    
    # Initialize results arrays
    length_errors = Float64[]
    mean_curvature_errors = Float64[]
    min_curvature_errors = Float64[]
    max_curvature_errors = Float64[]
    total_curvature_errors = Float64[]
    
    # Run tests for different resolutions
    for num_markers in resolutions
        circle = create_circle(center, radius, num_markers)
        
        # Calculate interface length
        measured_length = interface_length(circle)
        push!(length_errors, abs(measured_length - exact_length) / exact_length)
        
        # Calculate curvature statistics
        stats = curvature_statistics(circle)
        push!(mean_curvature_errors, abs(stats.mean - exact_curvature) / exact_curvature)
        push!(min_curvature_errors, abs(stats.minimum - exact_curvature) / exact_curvature)
        push!(max_curvature_errors, abs(stats.maximum - exact_curvature) / exact_curvature)
        push!(total_curvature_errors, abs(stats.total - exact_total_curvature) / exact_total_curvature)
        
        println("Circle with $num_markers markers:")
        println("  Length error: $(length_errors[end] * 100)%")
        println("  Mean curvature error: $(mean_curvature_errors[end] * 100)%")
        println("  Min/Max curvature error: $(min_curvature_errors[end] * 100)% / $(max_curvature_errors[end] * 100)%")
        println("  Total curvature error: $(total_curvature_errors[end] * 100)%")
    end
    
    # -------------- ELLIPSE CONVERGENCE TESTS --------------
    
    # Ellipse properties
    a, b = 3.0, 1.0  # Semi-major and semi-minor axes
    
    # Compute exact ellipse perimeter using elliptic integral
    # P = 4a * E(e) where e = sqrt(1 - b²/a²)
    e = sqrt(1 - (b/a)^2)
    exact_ellipse_length = 4 * a * elliptic_e(e)
    
    ellipse_length_errors = Float64[]
    ellipse_curvature_l1_errors = Float64[]
    ellipse_curvature_l2_errors = Float64[]
    ellipse_curvature_linf_errors = Float64[]
    
    for num_markers in resolutions
        # Create ellipse points
        t = range(0, 2π, length=num_markers+1)[1:num_markers]
        positions = [[a * cos(θ), b * sin(θ)] for θ in t]
        markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
        connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
        ellipse = Interface(markers, connectivity, true)
        
        # Calculate lengths
        measured_length = interface_length(ellipse)
        push!(ellipse_length_errors, abs(measured_length - exact_ellipse_length) / exact_ellipse_length)
        
        # Calculate numerical curvature
        curvatures = compute_curvature(ellipse)
        
        # Calculate analytical curvature for comparison
        analytical_curvatures = Float64[]
        for θ in t
            κ = (a * b) / (a^2 * sin(θ)^2 + b^2 * cos(θ)^2)^(3/2)
            push!(analytical_curvatures, κ)
        end
        
        # Compute different error norms
        differences = abs.(curvatures .- analytical_curvatures)
        rel_differences = differences ./ analytical_curvatures
        
        # L1 norm (mean absolute error)
        push!(ellipse_curvature_l1_errors, mean(differences))
        
        # L2 norm (root mean squared error)
        push!(ellipse_curvature_l2_errors, sqrt(mean(differences.^2)))
        
        # L∞ norm (maximum error)
        push!(ellipse_curvature_linf_errors, maximum(differences))
        
        println("Ellipse with $num_markers markers:")
        println("  Length error: $(ellipse_length_errors[end] * 100)%")
        println("  Curvature L1 error: $(ellipse_curvature_l1_errors[end])")
        println("  Curvature L2 error: $(ellipse_curvature_l2_errors[end])")
        println("  Curvature L∞ error: $(ellipse_curvature_linf_errors[end])")
    end
    
    # -------------- VISUALIZE CONVERGENCE RESULTS --------------
    
    # Create convergence plot for circle
    fig = Figure(size=(1200, 800))
    ax1 = Axis(fig[1, 1], title="Circle Property Convergence",
             xlabel="Number of Markers", ylabel="Relative Error (L∞ norm)",
             xscale=log10, yscale=log10,
             limits=(nothing, (1e-16, 1e0)))  # Set y-axis limits to show smaller errors
    
    lines!(ax1, resolutions, length_errors, label="Length", linewidth=2, color=:blue)
    scatter!(ax1, resolutions, length_errors, markersize=8, color=:blue)
    lines!(ax1, resolutions, mean_curvature_errors, label="Curvature", linewidth=2, color=:red) 
    #lines!(ax1, resolutions, min_curvature_errors, label="Min Curvature", linewidth=2, color=:green, linestyle=:dash)
    #lines!(ax1, resolutions, max_curvature_errors, label="Max Curvature", linewidth=2, color=:orange, linestyle=:dash)
    
    # Add reference lines showing O(h²) and O(h) convergence
    ref_h = resolutions[1] ./ resolutions
    ref_h2 = (resolutions[1] ./ resolutions).^2
    
    # Scale reference lines to match the data
    h_scale = length_errors[1] / ref_h[1] * 0.5
    h2_scale = length_errors[1] / ref_h2[1] * 0.1
    
    lines!(ax1, resolutions, h_scale .* ref_h, label="O(h)", linestyle=:dash, color=:black)
    lines!(ax1, resolutions, h2_scale .* ref_h2, label="O(h²)", linestyle=:dot, color=:black)
    
    axislegend(ax1, position=:lb)
    
    # Create convergence plot for ellipse
    ax2 = Axis(fig[1, 2], title="Ellipse Length Convergence",
             xlabel="Number of Markers", ylabel="Relative Error (L∞ norm)",
             xscale=log10, yscale=log10)
    
    lines!(ax2, resolutions, ellipse_length_errors, label="Length", linewidth=2, color=:blue)
    scatter!(ax2, resolutions, ellipse_length_errors, markersize=8, color=:blue)
    
    # Add reference lines
    h_scale_e = ellipse_length_errors[1] / ref_h[1] * 0.5
    h2_scale_e = ellipse_length_errors[1] / ref_h2[1] * 0.1
    
    lines!(ax2, resolutions, h_scale_e .* ref_h, label="O(h)", linestyle=:dash, color=:black)
    lines!(ax2, resolutions, h2_scale_e .* ref_h2, label="O(h²)", linestyle=:dot, color=:black)
    
    axislegend(ax2, position=:lb)
    
    # Create convergence plot for ellipse curvature with different norms
    ax3 = Axis(fig[2, 1:2], title="Ellipse Curvature Error Convergence",
             xlabel="Number of Markers", ylabel="Curvature Error (different norms)",
             xscale=log10, yscale=log10)
    
    lines!(ax3, resolutions, ellipse_curvature_l1_errors, label="L¹ Error", linewidth=2, color=:blue)
    scatter!(ax3, resolutions, ellipse_curvature_l1_errors, markersize=8, color=:blue)
    lines!(ax3, resolutions, ellipse_curvature_l2_errors, label="L² Error", linewidth=2, color=:red)
    scatter!(ax3, resolutions, ellipse_curvature_l2_errors, markersize=8, color=:red)
    lines!(ax3, resolutions, ellipse_curvature_linf_errors, label="L∞ Error", linewidth=2, color=:green)
    scatter!(ax3, resolutions, ellipse_curvature_linf_errors, markersize=8, color=:green)
    
    # Add reference lines
    h_scale_c = ellipse_curvature_l2_errors[1] / ref_h[1] * 0.5
    h2_scale_c = ellipse_curvature_l2_errors[1] / ref_h2[1] * 0.1
    
    lines!(ax3, resolutions, h_scale_c .* ref_h, label="O(h)", linestyle=:dash, color=:black)
    lines!(ax3, resolutions, h2_scale_c .* ref_h2, label="O(h²)", linestyle=:dot, color=:black)
    
    axislegend(ax3, position=:lb)
    
    save("convergence_results.png", fig)
    return fig
end

# Helper function to compute elliptic integral of the second kind
# This gives the exact perimeter of an ellipse
function elliptic_e(k)
    # For an ellipse with semi-major axis a and semi-minor axis b
    # e = sqrt(1 - b²/a²) is the eccentricity
    # The complete elliptic integral of the second kind E(e) gives the perimeter as 4a*E(e)
    
    # Implementation using SpecialFunctions.jl or approximation
    if isdefined(Main, :SpecialFunctions) && isdefined(Main.SpecialFunctions, :ellipe)
        return SpecialFunctions.ellipe(k^2)
    else
        # Approximate using series expansion (first few terms)
        m = k^2
        return π/2 * (1 - m/4 - 3*m^2/64 - 5*m^3/256)
    end
end

# Run the tests
run_convergence_tests()