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