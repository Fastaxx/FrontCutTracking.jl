
"""
    evolve_by_curvature_flow(interface::Interface, dt::Float64; 
                          motion_factor=1.0, 
                          noise=0.0,
                          constrain_motion=false)

Evolve an interface according to curvature flow, where each point moves in the 
direction of its normal proportional to the curvature.

# Arguments
- `interface`: The Interface object
- `dt`: Time step size
- `motion_factor`: Factor to scale the motion (can be negative)
- `noise`: Amount of random noise to add to the motion
- `constrain_motion`: Whether to constrain motion to prevent self-intersections

# Returns
A new Interface object representing the evolved state
"""
function evolve_by_curvature_flow(interface::Interface, dt::Float64;
                               motion_factor=1.0,
                               noise=0.0,
                               constrain_motion=false)
    
    # Copy the interface
    new_markers = Vector{Marker}(undef, length(interface.markers))
    
    # Compute normals and curvatures
    normals = compute_normals(interface)
    curvatures = compute_curvature(interface)
    
    # Move each marker along its normal proportional to curvature
    for (i, marker) in enumerate(interface.markers)
        # Calculate displacement vector
        displacement = normals[i] * curvatures[i] * motion_factor * dt
        
        # Add random noise if requested
        if noise > 0.0
            noise_vector = [2*rand()-1, 2*rand()-1]
            noise_vector = normalize(noise_vector) * noise * dt
            displacement += noise_vector
        end
        
        # Apply displacement to get new position
        new_position = marker.position + displacement
        
        # Create new marker with updated position
        new_markers[i] = Marker(new_position, marker.id)
    end
    
    # If constraining motion, check for self-intersections and adjust if needed
    if constrain_motion
        # Implement collision avoidance here
        # For now, just a placeholder
    end
    
    # Create new interface with updated markers
    return Interface(new_markers, interface.connectivity, interface.closed)
end

"""
    simulate_interface_evolution(initial_interface::Interface, 
                              n_steps::Int, 
                              dt::Float64;
                              evolution_function=evolve_by_curvature_flow,
                              kwargs...)

Simulate the evolution of an interface over time using a specified evolution function.

# Arguments
- `initial_interface`: The starting interface
- `n_steps`: Number of time steps
- `dt`: Time step size
- `evolution_function`: Function used to evolve interface (defaults to curvature flow)
- `kwargs...`: Additional arguments passed to the evolution function

# Returns
A vector of Interface objects representing the evolution over time
"""
function simulate_interface_evolution(initial_interface::Interface,
                                   n_steps::Int,
                                   dt::Float64;
                                   evolution_function=evolve_by_curvature_flow,
                                   kwargs...)
    
    # Initialize result array with the initial interface
    interfaces = Vector{Interface}(undef, n_steps + 1)
    interfaces[1] = initial_interface
    
    # Evolve interface for n_steps
    for i in 1:n_steps
        interfaces[i+1] = evolution_function(interfaces[i], dt; kwargs...)
    end
    
    return interfaces
end
