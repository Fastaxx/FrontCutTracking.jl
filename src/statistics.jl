"""
    curvature_histogram(interface::Interface; bins=10, range=:auto, absolute=false)

Compute a histogram of curvature values along the interface.

# Arguments
- `interface`: The Interface object
- `bins`: Number of histogram bins
- `range`: Range for histogram bins, :auto or tuple (min, max)
- `absolute`: Whether to use absolute curvature values

# Returns
A named tuple with:
- `edges`: Bin edges
- `counts`: Bin counts
- `weights`: Length-weighted bin counts (normalized)
"""
function curvature_histogram(interface::Interface; bins=10, range=:auto, absolute=false)
    curvatures = compute_curvature(interface)
    
    if absolute
        curvatures = abs.(curvatures)
    end
    
    # Determine range if auto
    if range === :auto
        min_val, max_val = extrema(curvatures)
        range_padding = 0.05 * (max_val - min_val)
        range = (min_val - range_padding, max_val + range_padding)
    end
    
    # Create bin edges
    bin_width = (range[2] - range[1]) / bins
    edges = [range[1] + i * bin_width for i in 0:bins]
    
    # Calculate segment lengths for weighting
    positions = [marker.position for marker in interface.markers]
    n = length(positions)
    
    if interface.closed
        segment_lengths = [norm(positions[i % n + 1] - positions[i]) for i in 1:n]
    else
        segment_lengths = zeros(n)
        for i in 2:n
            segment_lengths[i-1] += norm(positions[i] - positions[i-1]) / 2
            segment_lengths[i] += norm(positions[i] - positions[i-1]) / 2
        end
    end
    
    # Initialize counts
    counts = zeros(Int, bins)
    weighted_counts = zeros(Float64, bins)
    
    # Fill histogram
    for (i, k) in enumerate(curvatures)
        # Find bin index
        bin_idx = min(bins, max(1, Int(floor((k - range[1]) / bin_width) + 1)))
        counts[bin_idx] += 1
        weighted_counts[bin_idx] += segment_lengths[i]
    end
    
    # Normalize weighted counts
    if sum(segment_lengths) > 0
        weighted_counts ./= sum(segment_lengths)
    end
    
    return (edges=edges, counts=counts, weights=weighted_counts)
end


"""
    curvature_statistics(interface::Interface; absolute=false)

Calculate statistical measures of interface curvature.

# Arguments
- `interface`: The Interface object
- `absolute`: Whether to use absolute curvature values

# Returns
A named tuple with the following fields:
- `minimum`: Minimum curvature
- `maximum`: Maximum curvature
- `mean`: Mean curvature
- `median`: Median curvature
- `std`: Standard deviation of curvature
- `total`: Total curvature (integral along interface)
- `rms`: Root mean square curvature
"""
function curvature_statistics(interface::Interface; absolute=false)
    curvatures = compute_curvature(interface)
    
    if absolute
        curvatures = abs.(curvatures)
    end
    
    # Calculate segment lengths for weighted statistics
    positions = [marker.position for marker in interface.markers]
    n = length(positions)
    
    # For closed interfaces, connect the last point to the first
    if interface.closed
        segment_lengths = [norm(positions[i % n + 1] - positions[i]) for i in 1:n]
    else
        # For open interfaces, use half-segment lengths for endpoints
        segment_lengths = zeros(n)
        for i in 2:(n-1)
            # For interior points, use average of adjacent segment lengths
            segment_lengths[i] = (norm(positions[i] - positions[i-1]) + 
                                norm(positions[i+1] - positions[i])) / 2
        end
        # For endpoints, use their single adjacent segment length
        if n > 1
            segment_lengths[1] = norm(positions[2] - positions[1]) / 2
            segment_lengths[n] = norm(positions[n] - positions[n-1]) / 2
        end
    end
    
    # Make sure segment_lengths is not zero to avoid division issues
    if sum(segment_lengths) < 1e-10
        segment_lengths = ones(n)
    end
    
    # Normalize weights
    weights = segment_lengths ./ sum(segment_lengths)
    
    # Calculate statistics
    min_curv = minimum(curvatures)
    max_curv = maximum(curvatures)
    mean_curv = sum(curvatures .* weights)
    
    # Sort for median calculation
    sorted_indices = sortperm(curvatures)
    sorted_weights = weights[sorted_indices]
    
    # Calculate weighted median
    cumulative_weight = 0.0
    median_idx = 1
    for (i, w) in enumerate(sorted_weights)
        cumulative_weight += w
        if cumulative_weight >= 0.5
            median_idx = i
            break
        end
    end
    median_curv = curvatures[sorted_indices[median_idx]]
    
    # Calculate weighted standard deviation
    variance = sum(weights .* (curvatures .- mean_curv).^2)
    std_curv = sqrt(variance)
    
    # Calculate total curvature (integral along the interface)
    total_curv = sum(curvatures .* segment_lengths)
    
    # Calculate RMS curvature
    rms_curv = sqrt(sum(weights .* curvatures.^2))
    
    return (
        minimum = min_curv,
        maximum = max_curv,
        mean = mean_curv,
        median = median_curv,
        std = std_curv,
        total = total_curv,
        rms = rms_curv
    )
end

"""
    mean_curvature(interface::Interface; absolute=false)

Compute the mean curvature of an interface.

# Arguments
- `interface`: The Interface object
- `absolute`: Whether to use absolute curvature values

# Returns
The mean curvature value
"""
function mean_curvature(interface::Interface; absolute=false)
    stats = curvature_statistics(interface, absolute=absolute)
    return stats.mean
end

"""
    curvature_extrema(interface::Interface; n=1, min_separation=0)

Find the locations of maximum and minimum curvature along the interface.

# Arguments
- `interface`: The Interface object
- `n`: Number of extrema to find (of each type)
- `min_separation`: Minimum number of markers that must separate extrema points

# Returns
A named tuple with:
- `maxima`: Vector of (index, value, position) tuples for maximum curvature points
- `minima`: Vector of (index, value, position) tuples for minimum curvature points
"""
function curvature_extrema(interface::Interface; n=1, min_separation=0)
    curvatures = compute_curvature(interface)
    positions = [marker.position for marker in interface.markers]
    
    # Find local maxima and minima
    n_markers = length(curvatures)
    is_maximum = falses(n_markers)
    is_minimum = falses(n_markers)
    
    # Handle closed interfaces specially
    if interface.closed
        for i in 1:n_markers
            prev_idx = i == 1 ? n_markers : i - 1
            next_idx = i == n_markers ? 1 : i + 1
            
            if curvatures[i] > curvatures[prev_idx] && curvatures[i] > curvatures[next_idx]
                is_maximum[i] = true
            elseif curvatures[i] < curvatures[prev_idx] && curvatures[i] < curvatures[next_idx]
                is_minimum[i] = true
            end
        end
    else
        # For open interfaces, skip endpoints
        for i in 2:(n_markers-1)
            if curvatures[i] > curvatures[i-1] && curvatures[i] > curvatures[i+1]
                is_maximum[i] = true
            elseif curvatures[i] < curvatures[i-1] && curvatures[i] < curvatures[i+1]
                is_minimum[i] = true
            end
        end
    end
    
    # Get indices of maxima and minima
    max_indices = findall(is_maximum)
    min_indices = findall(is_minimum)
    
    # Sort by curvature value
    sort!(max_indices, by=i -> -curvatures[i])  # Descending order for maxima
    sort!(min_indices, by=i -> curvatures[i])   # Ascending order for minima
    
    # Apply minimum separation if needed
    if min_separation > 0
        max_indices = filter_by_separation(max_indices, min_separation, n_markers, interface.closed)
        min_indices = filter_by_separation(min_indices, min_separation, n_markers, interface.closed)
    end
    
    # Take the top n
    max_indices = max_indices[1:min(n, length(max_indices))]
    min_indices = min_indices[1:min(n, length(min_indices))]
    
    # Prepare return values
    maxima = [(i, curvatures[i], positions[i]) for i in max_indices]
    minima = [(i, curvatures[i], positions[i]) for i in min_indices]
    
    return (maxima=maxima, minima=minima)
end

"""
    length_scale_spectrum(interface::Interface; 
                        max_scales=10, 
                        scale_factor=1.2,
                        feature_measure=:curvature_extrema)

Compute a spectrum of feature significance across different length scales.

# Arguments
- `interface`: The Interface object to analyze
- `max_scales`: Maximum number of scales to analyze
- `scale_factor`: Factor by which smoothing increases between scales
- `feature_measure`: Method to quantify features (:curvature_extrema, :curvature_variance, or :fourier)

# Returns
A NamedTuple with:
- `scales`: Vector of length scales (smoothing levels)  
- `significance`: Vector of feature significance values at each scale
- `features`: Vector of feature counts or measures at each scale
"""
function length_scale_spectrum(interface::Interface; 
                             max_scales=10, 
                             scale_factor=1.2,
                             feature_measure=:curvature_extrema)
    
    # Initialize storage for results
    scales = Float64[]
    significance = Float64[]
    features = []
    
    # Create a copy of the interface to progressively smooth
    current = deepcopy(interface)
    
    # Get perimeter of interface for scale normalization
    perimeter = interface_length(current)
    base_scale = perimeter / length(current.markers)
    
    # Calculate initial feature measure before any smoothing
    initial_features = compute_feature_measure(current, feature_measure)
    push!(scales, base_scale)
    push!(significance, 1.0)  # Initial significance is 1.0 (100%)
    push!(features, initial_features)
    
    # Progressive smoothing and feature analysis
    for i in 1:max_scales
        # Calculate smoothing parameters for this scale
        smoothing_scale = base_scale * (scale_factor^i)
        iterations = max(1, round(Int, sqrt(i)))
        lambda = min(0.5, 0.2 * i / max_scales)
        
        # Apply smoothing at this scale
        smooth_interface!(current, iterations=iterations, lambda=lambda)
        
        # Measure features at this scale
        current_features = compute_feature_measure(current, feature_measure)
        
        # Calculate relative significance compared to original
        sig = feature_significance(current_features, initial_features, feature_measure)
        
        # Store results
        push!(scales, smoothing_scale)
        push!(significance, sig)
        push!(features, current_features)
    end
    
    return (scales=scales, significance=significance, features=features)
end

"""
    dominant_length_scales(interface::Interface; 
                         max_scales=15, 
                         threshold=0.2,
                         feature_measure=:curvature_extrema)

Identify the dominant length scales in an interface.

# Arguments
- `interface`: The Interface object to analyze
- `max_scales`: Maximum number of scales to analyze
- `threshold`: Threshold for significance changes to identify dominant scales
- `feature_measure`: Method to quantify features

# Returns
A Vector of dominant length scales
"""
function dominant_length_scales(interface::Interface; 
                              max_scales=15, 
                              threshold=0.2,
                              feature_measure=:curvature_extrema)
    
    # Compute the length scale spectrum
    spectrum = length_scale_spectrum(interface, 
                                    max_scales=max_scales, 
                                    feature_measure=feature_measure)
    
    # Identify scales where significant changes occur
    dominant = Float64[]
    
    # Look for significant drops in the significance curve
    for i in 2:length(spectrum.scales)
        delta = spectrum.significance[i-1] - spectrum.significance[i]
        if delta > threshold
            # This scale represents a significant change
            push!(dominant, spectrum.scales[i])
        end
    end
    
    return dominant
end

"""
    filter_by_scale(interface::Interface, scale::Float64)

Filter an interface to retain only features above a certain scale.

# Arguments
- `interface`: The Interface object to filter
- `scale`: The length scale threshold

# Returns
A new Interface with small-scale features removed
"""
function filter_by_scale(interface::Interface, scale::Float64)
    # Calculate number of smoothing iterations based on scale
    # This is a heuristic conversion from scale to smoothing parameters
    perimeter = interface_length(interface)
    base_scale = perimeter / length(interface.markers)
    
    # Convert scale to relative value
    relative_scale = scale / base_scale
    
    # Calculate smoothing parameters
    iterations = max(1, round(Int, log(relative_scale)))
    lambda = min(0.5, 0.3 * log(relative_scale))
    
    # Create a copy and apply appropriate smoothing
    filtered = deepcopy(interface)
    smooth_interface!(filtered, iterations=iterations, lambda=lambda)
    
    return filtered
end

"""
    length_scale_decomposition(interface::Interface; 
                             n_scales=3,
                             feature_measure=:curvature_extrema)

Decompose an interface into multiple components at different scales.

# Arguments
- `interface`: The Interface object to decompose
- `n_scales`: Number of scale levels to extract
- `feature_measure`: Method to choose scales

# Returns
A Vector of Interfaces representing different scale components
"""
function length_scale_decomposition(interface::Interface; 
                                  n_scales=3,
                                  feature_measure=:curvature_extrema)
    
    # First identify dominant scales
    scales = dominant_length_scales(interface, 
                                  max_scales=max(10, 2*n_scales),
                                  feature_measure=feature_measure)
    
    # If we don't have enough dominant scales, create evenly spaced scales
    if length(scales) < n_scales
        perimeter = interface_length(interface)
        base_scale = perimeter / length(interface.markers)
        max_scale = base_scale * 10  # Arbitrary maximum scale
        
        # Create logarithmically spaced scales
        scales = [base_scale * (max_scale/base_scale)^(i/(n_scales-1)) for i in 0:(n_scales-1)]
    else
        # Take the most significant n_scales
        scales = scales[1:min(n_scales, length(scales))]
    end
    
    # Create filtered interfaces at each scale
    decomposition = [filter_by_scale(interface, scale) for scale in scales]
    
    return decomposition
end

# Helper function to compute feature measures based on chosen method
function compute_feature_measure(interface::Interface, feature_measure::Symbol)
    if feature_measure == :curvature_extrema
        # Count significant extrema points
        extrema = curvature_extrema(interface, n=length(interface.markers) ÷ 4)
        return (maxima=length(extrema.maxima), minima=length(extrema.minima))
    
    elseif feature_measure == :curvature_variance
        # Measure variance in curvature
        curvatures = compute_curvature(interface)
        return var(curvatures)
    
    elseif feature_measure == :fourier
        # Compute simple Fourier decomposition (for closed curves)
        if interface.closed
            positions = [marker.position for marker in interface.markers]
            x = [pos[1] for pos in positions]
            y = [pos[2] for pos in positions]
            
            # Compute FFT (simplified for this use case)
            fx = abs.(fft(x))
            fy = abs.(fft(y))
            
            # Power spectrum (simplified)
            power = fx.^2 + fy.^2
            return power[2:div(end,2)]  # Return non-redundant frequencies
        else
            # For open curves, fall back to curvature variance
            curvatures = compute_curvature(interface)
            return var(curvatures)
        end
    else
        error("Unknown feature measure: $feature_measure")
    end
end

# Helper function to calculate feature significance
function feature_significance(current_features, initial_features, feature_measure::Symbol)
    if feature_measure == :curvature_extrema
        # Ratio of extrema counts
        total_current = current_features.maxima + current_features.minima
        total_initial = initial_features.maxima + initial_features.minima
        
        # Avoid division by zero
        if total_initial == 0
            return 0.0
        else
            return total_current / total_initial
        end
    
    elseif feature_measure == :curvature_variance
        # Ratio of variances
        if initial_features == 0
            return 0.0
        else
            return current_features / initial_features
        end
    
    elseif feature_measure == :fourier
        # For Fourier, use ratio of high-frequency power
        if length(current_features) >= 3 && length(initial_features) >= 3
            # Use upper half of frequencies
            high_idx = div(length(current_features), 2):length(current_features)
            
            current_power = sum(current_features[high_idx])
            initial_power = sum(initial_features[high_idx])
            
            if initial_power == 0
                return 0.0
            else
                return current_power / initial_power
            end
        else
            return 0.0
        end
    else
        return 0.0
    end
end
