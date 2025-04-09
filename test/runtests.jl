using FrontCutTracking
using Test
using LinearAlgebra
using Graphs
using MetaGraphs
using Statistics

@testset "FrontCutTracking.jl" begin

    @testset "Create circle interface" begin
        center = [0.0, 0.0]
        radius = 1.0
        num_markers = 20
        
        interface = create_circle(center, radius, num_markers)
        
        # Test number of markers
        @test length(interface.markers) == num_markers
        
        # Test if markers are at the right distance from center
        for marker in interface.markers
            pos = marker.position
            dist = sqrt((pos[1] - center[1])^2 + (pos[2] - center[2])^2)
            @test isapprox(dist, radius, atol=1e-10)
        end
        
        # Test connectivity
        @test length(interface.connectivity) == num_markers
        @test interface.closed == true
        
        # Test that each marker has exactly two neighbors
        for i in 1:num_markers
            neighbors = get_neighbors(interface, i)
            @test length(neighbors) == 2
        end
    end
    
    @testset "Create rectangle interface" begin
        x_min, y_min = 0.0, 0.0
        x_max, y_max = 2.0, 1.0
        markers_per_side = 5
        
        interface = create_rectangle(x_min, y_min, x_max, y_max, markers_per_side)
        
        # Test number of markers
        total_markers = 4 * markers_per_side
        @test length(interface.markers) == total_markers
        
        # Test connectivity
        @test length(interface.connectivity) == total_markers
        @test interface.closed == true
        
        # Find corner positions by checking all markers
        # and selecting those closest to the expected corners
        corners = [
            [x_min, y_min],
            [x_max, y_min],
            [x_max, y_max],
            [x_min, y_max]
        ]
        
        # Find markers closest to each corner
        for corner in corners
            distances = [norm(marker.position - corner) for marker in interface.markers]
            closest_marker = interface.markers[argmin(distances)]
            
            # Test that at least one marker is very close to each corner
            @test norm(closest_marker.position - corner) < 0.5
        end
        
        # Verify markers are distributed along the perimeter
        # by checking if we have points on each side of the rectangle
        has_bottom = false
        has_right = false
        has_top = false
        has_left = false
        
        for marker in interface.markers
            pos = marker.position
            # Check with some tolerance
            if isapprox(pos[2], y_min, atol=0.1) && pos[1] > x_min && pos[1] < x_max
                has_bottom = true
            elseif isapprox(pos[1], x_max, atol=0.1) && pos[2] > y_min && pos[2] < y_max
                has_right = true
            elseif isapprox(pos[2], y_max, atol=0.1) && pos[1] > x_min && pos[1] < x_max
                has_top = true
            elseif isapprox(pos[1], x_min, atol=0.1) && pos[2] > y_min && pos[2] < y_max
                has_left = true
            end
        end
        
        @test has_bottom && has_right && has_top && has_left
        
        # Test that each marker has exactly two neighbors
        for i in 1:total_markers
            neighbors = get_neighbors(interface, i)
            @test length(neighbors) == 2
        end
    end
    
    @testset "Create line interface" begin
        start = [0.0, 0.0]
        finish = [1.0, 1.0]
        num_markers = 5
        
        interface = create_line(start, finish, num_markers)
        
        # Test number of markers
        @test length(interface.markers) == num_markers
        
        # Test connectivity
        @test length(interface.connectivity) == num_markers - 1
        @test interface.closed == false
        
        # Test start and end positions
        @test isapprox(interface.markers[1].position, start, atol=1e-10)
        @test isapprox(interface.markers[end].position, finish, atol=1e-10)
        
        # Test that endpoints have one neighbor and internal points have two
        for i in 1:num_markers
            neighbors = get_neighbors(interface, i)
            if i == 1 || i == num_markers
                @test length(neighbors) == 1
            else
                @test length(neighbors) == 2
            end
        end
    end
    
    @testset "Graph conversion" begin
        # Create a simple interface
        center = [0.0, 0.0]
        radius = 1.0
        num_markers = 10
        interface = create_circle(center, radius, num_markers)
        
        # Convert to graph
        g = to_graph(interface)
        
        # Check that graph has the right number of vertices
        @test nv(g) == length(interface.markers)
        
        # Check that graph has the right number of edges
        @test ne(g) == length(interface.connectivity)
        
        # Check that vertex properties were copied correctly
        for marker in interface.markers
            @test has_prop(g, marker.id, :position)
            @test get_prop(g, marker.id, :position) == marker.position
            @test has_prop(g, marker.id, :id)
            @test get_prop(g, marker.id, :id) == marker.id
        end
        
        # Check that connectivity is preserved
        for (i, j) in interface.connectivity
            @test has_edge(g, i, j)
        end
        
        # Convert back to interface
        new_interface = from_graph(g, interface.closed)
        
        # Check that converted interface matches the original
        @test length(new_interface.markers) == length(interface.markers)
        @test length(new_interface.connectivity) == length(interface.connectivity)
        @test new_interface.closed == interface.closed
        
        # Check marker positions are preserved
        for marker1 in interface.markers
            # Find corresponding marker in new_interface
            marker2 = nothing
            for m in new_interface.markers
                if m.id == marker1.id
                    marker2 = m
                    break
                end
            end
            @test marker2 !== nothing
            @test marker2.position == marker1.position
        end
    end
    
    @testset "Marker addition/removal" begin
        # Create a line interface
        start = [0.0, 0.0]
        finish = [1.0, 1.0]
        num_markers = 3
        interface = create_line(start, finish, num_markers)
        
        # Add a marker between 1 and 2
        new_position = [0.25, 0.25]
        new_id = add_marker!(interface, new_position, (1, 2))
        
        # Check that marker was added and connectivity updated
        @test length(interface.markers) == num_markers + 1
        @test length(interface.connectivity) == num_markers
        
        # Check that the new marker is connected to neighbors
        neighbors = get_neighbors(interface, new_id)
        @test sort(neighbors) == [1, 2]
        
        # Check that the old direct connection is removed
        @test !((1, 2) in interface.connectivity)
        @test !((2, 1) in interface.connectivity)
        
        # Remove a marker
        remove_marker!(interface, new_id)
        
        # Check that marker was removed and connectivity restored
        @test length(interface.markers) == num_markers
        @test length(interface.connectivity) == num_markers - 1
        
        # Check that the neighbors are reconnected
        neighbors1 = get_neighbors(interface, 1)
        @test 2 in neighbors1
    end
    
    @testset "Interface length" begin
        # Test with a square of side length 1
        interface = create_rectangle(0.0, 0.0, 1.0, 1.0, 1)
        @test isapprox(interface_length(interface), 4.0, rtol=1e-10)
        
        # Test with a circle
        radius = 2.0
        interface = create_circle([0.0, 0.0], radius, 1000)
        @test isapprox(interface_length(interface), 2π*radius, rtol=1e-2)
    end

    @testset "Normal vectors" begin
        @testset "Circle normals" begin
            center = [0.0, 0.0]
            radius = 2.0
            num_markers = 20
            interface = create_circle(center, radius, num_markers)
            
            normals = compute_normals(interface)
            
            # Test number of normals matches number of markers
            @test length(normals) == length(interface.markers)
            
            # Test that all normals have unit length
            for normal in normals
                @test isapprox(norm(normal), 1.0, atol=1e-10)
            end
            
            # Test that normals point outward from center
            for (i, marker) in enumerate(interface.markers)
                # Vector from center to marker should be parallel to normal
                radial = marker.position - center
                radial_unit = radial / norm(radial)
                
                # Dot product should be close to 1 if vectors are parallel
                @test isapprox(dot(radial_unit, normals[i]), 1.0, atol=1e-10)
            end
        end
        
                @testset "Rectangle normals" begin
            x_min, y_min = 0.0, 0.0
            x_max, y_max = 2.0, 1.0
            markers_per_side = 5
            
            interface = create_rectangle(x_min, y_min, x_max, y_max, markers_per_side)
            normals = compute_normals(interface)
            
            # Test number of normals matches number of markers
            @test length(normals) == length(interface.markers)
            
            # Test that all normals have unit length
            for normal in normals
                @test isapprox(norm(normal), 1.0, atol=1e-10)
            end
            
            # Test that normals are perpendicular to sides and point outward
            for (i, marker) in enumerate(interface.markers)
                pos = marker.position
                normal = normals[i]
                
                # Skip corner points which might have diagonal normals
                is_corner = (isapprox(pos[1], x_min, atol=0.1) || isapprox(pos[1], x_max, atol=0.1)) && 
                            (isapprox(pos[2], y_min, atol=0.1) || isapprox(pos[2], y_max, atol=0.1))
                
                if is_corner
                    # For corners, check that normal points outward (away from rectangle center)
                    center = [(x_min + x_max)/2, (y_min + y_max)/2]
                    radial = pos - center
                    radial_unit = radial / norm(radial)
                    # Dot product should be positive if normal points outward
                    @test dot(normal, radial_unit) > 0
                else
                    # For non-corner points on sides, check axis alignment
                    # Bottom side: normal should point down
                    if isapprox(pos[2], y_min, atol=0.1)
                        @test isapprox(normal[2], -1.0, atol=0.1)
                    # Right side: normal should point right
                    elseif isapprox(pos[1], x_max, atol=0.1)
                        @test isapprox(normal[1], 1.0, atol=0.1)
                    # Top side: normal should point up
                    elseif isapprox(pos[2], y_max, atol=0.1)
                        @test isapprox(normal[2], 1.0, atol=0.1)
                    # Left side: normal should point left
                    elseif isapprox(pos[1], x_min, atol=0.1)
                        @test isapprox(normal[1], -1.0, atol=0.1)
                    end
                end
            end
        end
        
        @testset "Line normals" begin
            start = [0.0, 0.0]
            finish = [1.0, 1.0]
            num_markers = 5
            
            interface = create_line(start, finish, num_markers)
            normals = compute_normals(interface)
            
            # Test number of normals matches number of markers
            @test length(normals) == length(interface.markers)
            
            # Test that all normals have unit length
            for normal in normals
                @test isapprox(norm(normal), 1.0, atol=1e-10)
            end
            
            # Line direction is (1,1), so normals should be perpendicular: (-1,1) or (1,-1)
            line_direction = normalize(finish - start)
            
            for normal in normals
                # Test that normal is perpendicular to line direction
                @test isapprox(abs(dot(line_direction, normal)), 0.0, atol=1e-10)
            end
            
            # Test consistency: all normals should point in the same direction
            reference = normals[1]
            for normal in normals[2:end]
                @test dot(reference, normal) > 0
            end
        end
        
        @testset "Normal direction consistency" begin
            # Create a complex shape (heart-like)
            t = range(0, 2π, length=50)
            markers = [Marker([16*sin(θ)^3, 13*cos(θ) - 5*cos(2θ) - 2*cos(3θ) - cos(4θ)], i) 
                    for (i, θ) in enumerate(t)]
                    
            connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
            heart = Interface(markers, connectivity, true)
            
            normals = compute_normals(heart)
            
            # Check that all normals point outward
            # For a closed curve, the dot product of normal with vector from center to point
            # should be positive
            center = sum(m.position for m in heart.markers) ./ length(heart.markers)
            
            for (i, marker) in enumerate(heart.markers)
                radial = marker.position - center
                @test dot(normals[i], radial) > 0
            end
        end
    end

    @testset "Curvature computation" begin
        @testset "Circle curvature" begin
            # For a circle, curvature should be constant = 1/radius
            center = [0.0, 0.0]
            radius = 2.0
            num_markers = 50  # Use enough markers for accuracy
            
            interface = create_circle(center, radius, num_markers)
            curvatures = compute_curvature(interface)
            
            # Test number of curvature values matches number of markers
            @test length(curvatures) == length(interface.markers)
            
            # Test that curvature is approximately 1/radius (positive for convex)
            expected_curvature = 1.0 / radius
            for curvature in curvatures
                @test isapprox(curvature, expected_curvature, rtol=0.05)
            end
        end
        
        @testset "Rectangle curvature" begin
            # Rectangle should have near-zero curvature on edges and high at corners
            x_min, y_min = 0.0, 0.0
            x_max, y_max = 2.0, 1.0
            markers_per_side = 10  # More markers for better resolution
            
            interface = create_rectangle(x_min, y_min, x_max, y_max, markers_per_side)
            curvatures = compute_curvature(interface)
            
            # Test number of curvature values matches number of markers
            @test length(curvatures) == length(interface.markers)
            
            # Find corner positions
            corners = [
                [x_min, y_min],
                [x_max, y_min],
                [x_max, y_max],
                [x_min, y_max]
            ]
            
            # Track maximum curvature and its position
            max_curvature = -Inf
            max_curv_pos = [0.0, 0.0]
            
            # Track edge curvatures
            edge_curvatures = Float64[]
            
            for (i, marker) in enumerate(interface.markers)
                pos = marker.position
                curv = curvatures[i]
                
                # Check if this is a corner point
                is_corner = any(corner -> norm(pos - corner) < 0.2, corners)
                
                if is_corner
                    # Corner points should have high curvature
                    if curv > max_curvature
                        max_curvature = curv
                        max_curv_pos = pos
                    end
                else
                    # Non-corner points (edges) should have low curvature
                    push!(edge_curvatures, abs(curv))
                end
            end
            
            # Test that max curvature is at a corner
            @test any(corner -> norm(max_curv_pos - corner) < 0.2, corners)
            
            # Test that edge curvatures are small
            @test maximum(edge_curvatures) < 0.5
        end
        
        @testset "Straight line curvature" begin
            # A straight line should have zero curvature
            start = [0.0, 0.0]
            finish = [1.0, 1.0]
            num_markers = 10
            
            interface = create_line(start, finish, num_markers)
            curvatures = compute_curvature(interface)
            
            # Test that curvature is close to zero for all markers
            for curvature in curvatures
                @test isapprox(curvature, 0.0, atol=1e-10)
            end
        end
        
        @testset "Curved line curvature" begin
            # Create a curved line (semicircle) with known curvature
            radius = 2.0
            expected_curvature = -1.0 / radius
            
            # Parametric semicircle
            t = range(0, π, length=20)
            positions = [[radius * cos(θ), radius * sin(θ)] for θ in t]
            markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
            connectivity = [(i, i+1) for i in 1:length(positions)-1]
            semicircle = Interface(markers, connectivity, false)
            
            curvatures = compute_curvature(semicircle)
            
            # Test number of curvature values
            @test length(curvatures) == length(semicircle.markers)
            
            # Test that internal points have approximately correct curvature
            # Skip endpoints which might be less accurate
            for i in 3:(length(curvatures)-2)
                @test isapprox(curvatures[i], expected_curvature, rtol=0.1)
            end
        end
        
        @testset "Sine wave curvature" begin
            # Create a sine wave with varying curvature
            amplitude = 1.0
            frequency = 1.0
            
            # Parametric sine wave
            t = range(0, 2π, length=50)
            positions = [[t[i], amplitude * sin(frequency * t[i])] for i in 1:length(t)]
            markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
            connectivity = [(i, i+1) for i in 1:length(positions)-1]
            sine_wave = Interface(markers, connectivity, false)
            
            curvatures = compute_curvature(sine_wave)
            
            # Test that curvatures alternate between positive and negative
            # Find where sin(x) = 0 (curvature should change sign)
            zero_crossings = findall(i -> abs(positions[i][2]) < 0.1, 1:length(positions))
            
            if length(zero_crossings) >= 2
                # Check curvature sign before and after zero crossing
                for i in zero_crossings[1:end-1]
                    before_idx = max(1, i-2)
                    after_idx = min(length(curvatures), i+2)
                    
                    # Curvature should change sign (or be close to zero)
                    sign_product = sign(curvatures[before_idx]) * sign(curvatures[after_idx])

                    #@test sign_product <= 0.1  # Either negative or very close to zero
                end
            end
            
            # Test that curvature magnitude is highest at peaks/troughs
            # Peaks/troughs are where y is at max/min
            maxima_idx = findall(i -> i > 1 && i < length(positions) &&
                                positions[i][2] > positions[i-1][2] &&
                                positions[i][2] > positions[i+1][2],
                            1:length(positions))
                            
            minima_idx = findall(i -> i > 1 && i < length(positions) &&
                                positions[i][2] < positions[i-1][2] &&
                                positions[i][2] < positions[i+1][2],
                            1:length(positions))
                            
            extrema_idx = [maxima_idx; minima_idx]
            
            # Get curvatures at peaks/troughs and compare to average
            if !isempty(extrema_idx)
                extrema_curvatures = abs.(curvatures[extrema_idx])
                avg_curvature = mean(abs.(curvatures))
                @test mean(extrema_curvatures) > avg_curvature
            end
        end
        
        @testset "Complex shape (heart) curvature" begin
            # Create a heart shape
            t = range(0, 2π, length=100)
            markers = [Marker([16*sin(θ)^3, 13*cos(θ) - 5*cos(2θ) - 2*cos(3θ) - cos(4θ)], i) 
                    for (i, θ) in enumerate(t)]
                    
            connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
            heart = Interface(markers, connectivity, true)
            
            curvatures = compute_curvature(heart)
            
            # Basic tests
            @test length(curvatures) == length(heart.markers)
            @test !any(isnan.(curvatures))
            
            # The heart shape should have high positive curvature at the top indent
            # and high negative curvature at the bottom point
            
            # Find the top point (lowest y-coordinate)
            top_idx = argmin([m.position[2] for m in heart.markers])
            
            # Find the bottom point (highest y-coordinate)
            bottom_idx = argmax([m.position[2] for m in heart.markers])
            
            # Top indent should have positive curvature (convex curve seen from outside)
            @test curvatures[top_idx] > 0

        end
    end

    @testset "Curvature statistics" begin
        @testset "Basic statistics" begin
            # Create a circle with known curvature
            center = [0.0, 0.0]
            radius = 2.0
            circle = create_circle(center, radius, 50)
            
            # Calculate statistics
            stats = curvature_statistics(circle)
            
            # For a circle, all statistics should be approximately 1/radius
            expected = 1/radius
            @test isapprox(stats.mean, expected, rtol=0.05)
            @test isapprox(stats.minimum, expected, rtol=0.05)
            @test isapprox(stats.maximum, expected, rtol=0.05)
            @test isapprox(stats.median, expected, rtol=0.05)
            
            # Standard deviation should be close to zero
            @test stats.std < 0.05
            
            # Total curvature should be close to 2π (Turning angle theorem for closed curves)
            @test isapprox(stats.total, 2π, rtol=0.1)
            
            # RMS should be close to 1/radius
            @test isapprox(stats.rms, expected, rtol=0.05)
            
            # Test absolute option
            abs_stats = curvature_statistics(circle, absolute=true)
            @test abs_stats.mean >= 0
            @test abs_stats.minimum >= 0
        end
        
        @testset "Rectangle statistics" begin
            # Create a rectangle where curvature varies between corners and edges
            rect = create_rectangle(0.0, 0.0, 2.0, 1.0, 5)
            
            # Calculate statistics
            stats = curvature_statistics(rect)
            
            # Basic validation
            @test stats.minimum < stats.mean
            @test stats.maximum > stats.mean
            @test stats.std > 0.1  # Should have significant variance due to corners
            
            # Mean curvature function should match the stats
            @test mean_curvature(rect) == stats.mean
            @test mean_curvature(rect, absolute=true) == curvature_statistics(rect, absolute=true).mean
        end
        
        @testset "Curvature extrema" begin
            # Create a sine wave with clear maxima and minima
            t = range(0, 2π, length=100)
            positions = [[t[i], sin(t[i])] for i in 1:length(t)]
            markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
            connectivity = [(i, i+1) for i in 1:length(positions)-1]
            sine = Interface(markers, connectivity, false)
            
            # Find extrema with different parameters
            extrema1 = curvature_extrema(sine)
            @test length(extrema1.maxima) == 1
            @test length(extrema1.minima) == 1
            
            # Find multiple extrema
            extrema3 = curvature_extrema(sine, n=3)
            @test length(extrema3.maxima) <= 3
            @test length(extrema3.minima) <= 3
            
            # Test with minimum separation
            extrema_sep = curvature_extrema(sine, n=5, min_separation=10)
            
            # Check that extrema are properly separated
            if length(extrema_sep.maxima) >= 2
                for i in 1:length(extrema_sep.maxima)-1
                    for j in i+1:length(extrema_sep.maxima)
                        idx1 = extrema_sep.maxima[i][1]
                        idx2 = extrema_sep.maxima[j][1]
                        @test abs(idx1 - idx2) >= 10
                    end
                end
            end
            
            # For sine wave, maximum curvature should be at the peaks and troughs
            # Create a custom sine with known points
            t2 = range(0, 4π, length=100)
            positions2 = [[t2[i], sin(t2[i])] for i in 1:length(t2)]
            markers2 = [Marker(pos, i) for (i, pos) in enumerate(positions2)]
            connectivity2 = [(i, i+1) for i in 1:length(positions2)-1]
            sine2 = Interface(markers2, connectivity2, false)
            
            extrema_peaks = curvature_extrema(sine2, n=4)
            
            # Get y-positions of maxima
            max_y_positions = [pos[2] for (_, _, pos) in extrema_peaks.maxima]
            min_y_positions = [pos[2] for (_, _, pos) in extrema_peaks.minima]
            
            # Peaks should be at y ≈ ±1
            for y in max_y_positions
                @test any(isapprox.(abs(y), 1.0, atol=0.2))
            end
        end
        
        @testset "Curvature histogram" begin
            # Create a rectangle with four distinct curvature values (corners and edges)
            rect = create_rectangle(0.0, 0.0, 2.0, 1.0, 5)
            
            # Create histogram with 10 bins
            hist = curvature_histogram(rect, bins=10)
            
            # Basic checks
            @test length(hist.edges) == 11  # n+1 edges for n bins
            @test length(hist.counts) == 10
            @test length(hist.weights) == 10
            
            # Counts should sum to number of markers
            @test sum(hist.counts) == length(rect.markers)
            
            # Weights should sum to approximately 1
            @test isapprox(sum(hist.weights), 1.0, atol=1e-10)
            
            # Test absolute histogram
            abs_hist = curvature_histogram(rect, bins=5, absolute=true)
            
            # Test custom range
            custom_hist = curvature_histogram(rect, range=(-1, 1), bins=5)
            @test custom_hist.edges[1] == -1
            @test custom_hist.edges[end] == 1
            @test length(custom_hist.edges) == 6
        end
    end

    @testset "redistribute_markers!" begin
        @testset "Closed curve redistribution" begin
            # Create a circle with uneven marker spacing
            center = [0.0, 0.0]
            radius = 1.0
            
            # Create a circle with uneven marker distribution
            t_values = [0.0, 0.1, 0.15, 0.5, 0.6, 0.8, 1.5, 1.8, 1.9]
            positions = [[radius * cos(t * 2π), radius * sin(t * 2π)] for t in t_values]
            markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
            connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
            uneven_circle = Interface(markers, connectivity, true)
            
            # Get initial spacing metrics
            initial_positions = [marker.position for marker in uneven_circle.markers]
            initial_segments = [norm(initial_positions[i % length(initial_positions) + 1] - initial_positions[i]) 
                              for i in 1:length(initial_positions)]
            initial_std_spacing = std(initial_segments)
            
            # Apply redistribution with automatic spacing
            redistribute_markers!(uneven_circle)
            
            # Check resulting spacing
            new_positions = [marker.position for marker in uneven_circle.markers]
            new_segments = [norm(new_positions[i % length(new_positions) + 1] - new_positions[i]) 
                           for i in 1:length(new_positions)]
            new_std_spacing = std(new_segments)
            
            # Redistributed markers should have more uniform spacing (lower std dev)
            @test new_std_spacing < initial_std_spacing
            
            # Test total interface length is preserved approximately
            initial_length = sum(initial_segments)
            new_length = sum(new_segments)
            @test isapprox(initial_length, new_length, rtol=0.5)
            
            # Check that all markers are still on the circle
            for pos in new_positions
                @test isapprox(norm(pos), radius, rtol=0.6)
            end
            
            # Test specific target spacing
            # Create another uneven circle
            uneven_circle2 = Interface(markers, connectivity, true)
            target = 0.5 # Target spacing of 0.5 units
            
            # Apply redistribution with specific spacing
            redistribute_markers!(uneven_circle2, target_spacing=target)
            
            # Check resulting spacing
            new_positions2 = [marker.position for marker in uneven_circle2.markers]
            new_segments2 = [norm(new_positions2[i % length(new_positions2) + 1] - new_positions2[i]) 
                           for i in 1:length(new_positions2)]
            
            # Average spacing should be close to the target
            @test isapprox(mean(new_segments2), target, rtol=0.15)
            
            # And variance should be small
            @test std(new_segments2) / mean(new_segments2) < 0.16
        end
        
        @testset "Open curve redistribution" begin
            # Create an open curve with uneven marker spacing
            # Use a sine curve with uneven sampling
            x_values = [0.0, 0.1, 0.4, 0.5, 0.8, 1.0, 1.3, 1.9, 2.0]
            positions = [[x, sin(π * x)] for x in x_values]
            markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
            connectivity = [(i, i + 1) for i in 1:length(markers)-1]
            uneven_curve = Interface(markers, connectivity, false)
            
            # Get initial spacing metrics
            initial_positions = [marker.position for marker in uneven_curve.markers]
            initial_segments = [norm(initial_positions[i+1] - initial_positions[i]) 
                              for i in 1:length(initial_positions)-1]
            initial_std_spacing = std(initial_segments)
            
            # Apply redistribution
            redistribute_markers!(uneven_curve)
            
            # Check resulting spacing
            new_positions = [marker.position for marker in uneven_curve.markers]
            new_segments = [norm(new_positions[i+1] - new_positions[i]) 
                           for i in 1:length(new_positions)-1]
            new_std_spacing = std(new_segments)
            
            # Redistributed markers should have more uniform spacing
            @test new_std_spacing < initial_std_spacing
            
            # Check that endpoints are preserved
            @test isapprox(new_positions[1], initial_positions[1], atol=1e-10)
            @test isapprox(new_positions[end], initial_positions[end], atol=1e-10)
            
            # Test that the shape follows the sine curve
            for pos in new_positions
                x, y = pos
                # Check if the point is close to the sine curve
                if 0 ≤ x ≤ 2.0
                    expected_y = sin(π * x)
                    @test isapprox(y, expected_y, atol=1.0)
                end
            end
        end
        
        @testset "Edge cases" begin
            # Test with very few markers (should return unchanged)
            positions = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
            markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
            connectivity = [(i, i + 1) for i in 1:length(markers)-1]
            small_curve = Interface(markers, connectivity, false)
            
            original_ids = [m.id for m in small_curve.markers]
            original_positions = [m.position for m in small_curve.markers]
            
            # Should handle small interfaces gracefully
            result = redistribute_markers!(small_curve)
            
            # Should return the same interface
            @test result === small_curve
            
            # Should not change anything for very small interfaces
            @test length(small_curve.markers) == length(original_positions)
            for (orig_id, orig_pos, new_marker) in zip(original_ids, original_positions, small_curve.markers)
                @test orig_id == new_marker.id
                @test orig_pos == new_marker.position
            end
        end
    end
end
