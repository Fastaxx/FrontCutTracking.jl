using FrontCutTracking
using Test
using LinearAlgebra
using Graphs
using MetaGraphs

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
end