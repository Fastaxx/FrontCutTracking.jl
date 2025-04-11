"""
    compute_fourier_descriptors(interface::Interface; 
                              num_descriptors=nothing, 
                              normalize=true)

Compute Fourier descriptors for an interface shape.

# Arguments
- `interface`: The Interface object to analyze
- `num_descriptors`: Number of descriptors to return (defaults to half of markers)
- `normalize`: Whether to normalize descriptors by the first coefficient for scale invariance

# Returns
A NamedTuple with:
- `coefficients`: Complex Fourier coefficients
- `magnitudes`: Magnitude of each coefficient
- `phases`: Phase angles of each coefficient
- `normalized`: Whether the descriptors are normalized
"""
function compute_fourier_descriptors(interface::Interface; 
                                   num_descriptors=nothing, 
                                   normalize=true)
    # Get marker positions
    positions = [marker.position for marker in interface.markers]
    
    # Convert to complex representation for easier Fourier analysis
    complex_points = [Complex(pos[1], pos[2]) for pos in positions]
    
    # Compute FFT
    coeffs = fft(complex_points)
    
    # Set default number of descriptors if not specified
    if isnothing(num_descriptors)
        num_descriptors = div(length(coeffs), 2)
    end
    
    # Limit to desired number of descriptors
    coeffs = coeffs[1:min(num_descriptors, length(coeffs))]
    
    # Normalize if requested (for scale/rotation/translation invariance)
    if normalize && length(coeffs) > 1
        # First coefficient represents the centroid
        centroid = coeffs[1] / length(complex_points)
        
        # Center the shape (translation invariance)
        coeffs[1] = 0
        
        # Normalize by the magnitude of the second coefficient (scale invariance)
        if length(coeffs) > 1 && abs(coeffs[2]) > 0
            normalization_factor = abs(coeffs[2])
            coeffs = coeffs ./ normalization_factor
        end
    end
    
    # Calculate magnitudes and phases
    magnitudes = abs.(coeffs)
    phases = angle.(coeffs)
    
    return (
        coefficients = coeffs,
        magnitudes = magnitudes,
        phases = phases,
        normalized = normalize
    )
end

"""
    reconstruct_from_fourier(descriptors::NamedTuple, 
                          num_points=100,
                          freq_range=nothing)

Reconstruct an interface from its Fourier descriptors.

# Arguments
- `descriptors`: Fourier descriptors from compute_fourier_descriptors
- `num_points`: Number of points in the reconstructed interface
- `freq_range`: Range of frequencies to include in reconstruction (e.g., 1:10)

# Returns
A new Interface object reconstructed from the descriptors
"""
function reconstruct_from_fourier(descriptors::NamedTuple, 
                               num_points=100,
                               freq_range=nothing)
    
    coeffs = descriptors.coefficients
    
    # Use all coefficients if frequency range not specified
    if isnothing(freq_range)
        freq_range = 1:length(coeffs)
    end
    
    # Limit to available coefficients
    freq_range = freq_range[freq_range .<= length(coeffs)]
    
    # Create full spectrum for ifft (with zeros for unused frequencies)
    full_spectrum = zeros(Complex{Float64}, num_points)
    
    # Fill in the specified frequency components
    for (i, freq) in enumerate(freq_range)
        if freq <= length(coeffs)
            idx = freq
            full_spectrum[idx] = coeffs[idx]
            
            # Add conjugate for proper reconstruction (if not DC component)
            if freq > 1 && idx < num_points - idx + 2
                full_spectrum[num_points - idx + 2] = conj(coeffs[idx])
            end
        end
    end
    
    # Inverse FFT to get reconstructed points
    reconstructed = ifft(full_spectrum)
    
    # Convert back to [x, y] positions
    positions = [[real(z), imag(z)] for z in reconstructed]
    
    # Create markers and connectivity
    markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
    connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
    
    # Create new interface
    return Interface(markers, connectivity, true)
end

"""
    fourier_power_spectrum(interface::Interface; normalize=true)

Compute the power spectrum of a shape using Fourier analysis.

# Arguments
- `interface`: The Interface object to analyze
- `normalize`: Whether to normalize the spectrum

# Returns
A NamedTuple with:
- `frequencies`: Vector of frequencies
- `power`: Power at each frequency
- `cumulative_power`: Cumulative power distribution
"""
function fourier_power_spectrum(interface::Interface; normalize=true)
    # Compute Fourier descriptors
    descriptors = compute_fourier_descriptors(interface, normalize=normalize)
    
    # Calculate power (squared magnitude)
    power = abs2.(descriptors.coefficients)
    
    # Create frequency vector
    frequencies = 0:(length(power)-1)
    
    # Calculate cumulative power
    cum_power = cumsum(power) ./ sum(power)
    
    return (
        frequencies = frequencies,
        power = power,
        cumulative_power = cum_power
    )
end

"""
    dominant_frequencies(interface::Interface; 
                       threshold=0.95, 
                       max_frequencies=10)

Identify dominant frequencies in a shape.

# Arguments
- `interface`: The Interface object to analyze
- `threshold`: Cumulative power threshold (0-1)
- `max_frequencies`: Maximum number of frequencies to return

# Returns
A Vector of dominant frequency indices
"""
function dominant_frequencies(interface::Interface; 
                            threshold=0.95, 
                            max_frequencies=10)
    
    # Compute power spectrum
    spectrum = fourier_power_spectrum(interface)
    
    # Sort frequencies by power (descending)
    sorted_indices = sortperm(spectrum.power, rev=true)
    
    # Get cumulative power for sorted frequencies
    sorted_power = spectrum.power[sorted_indices]
    cum_power = cumsum(sorted_power) / sum(sorted_power)
    
    # Find frequencies that capture threshold of power
    threshold_idx = findfirst(x -> x >= threshold, cum_power)
    
    # Limit by threshold and max_frequencies
    if isnothing(threshold_idx)
        threshold_idx = length(cum_power)
    end
    
    num_freqs = min(threshold_idx, max_frequencies)
    
    # Return the dominant frequency indices
    return sorted_indices[1:num_freqs]
end

"""
    fourier_filter(interface::Interface, 
                 cutoff_frequency::Int; 
                 filter_type=:lowpass,
                 num_points=nothing)

Filter an interface shape using Fourier coefficients.

# Arguments
- `interface`: The Interface object to filter
- `cutoff_frequency`: Frequency cutoff point
- `filter_type`: Type of filter: :lowpass, :highpass, or :bandpass
- `num_points`: Number of points in the filtered interface (defaults to original count)

# Returns
A new Interface object with filtered frequency components
"""
function fourier_filter(interface::Interface, 
                      cutoff_frequency::Int; 
                      filter_type=:lowpass,
                      band_end=nothing,
                      num_points=nothing)
    
    if isnothing(num_points)
        num_points = length(interface.markers)
    end
    
    # Compute descriptors
    descriptors = compute_fourier_descriptors(interface, num_descriptors=div(num_points, 2) + 1)
    
    # Define frequency range based on filter type
    if filter_type == :lowpass
        freq_range = 1:cutoff_frequency
    elseif filter_type == :highpass
        freq_range = cutoff_frequency:length(descriptors.coefficients)
    elseif filter_type == :bandpass
        if isnothing(band_end)
            error("band_end must be specified for bandpass filter")
        end
        freq_range = cutoff_frequency:band_end
    else
        error("Unknown filter type: $filter_type")
    end
    
    # Reconstruct from filtered frequencies
    return reconstruct_from_fourier(descriptors, num_points, freq_range)
end

"""
    fourier_shape_similarity(interface1::Interface, interface2::Interface; 
                           num_descriptors=20, 
                           method=:euclidean)

Calculate similarity between two shapes using Fourier descriptors.

# Arguments
- `interface1`: First Interface object
- `interface2`: Second Interface object
- `num_descriptors`: Number of descriptors to use in comparison
- `method`: Distance method (:euclidean, :manhattan, or :cosine)

# Returns
A similarity score (lower values indicate more similar shapes)
"""
function fourier_shape_similarity(interface1::Interface, interface2::Interface; 
                                num_descriptors=20, 
                                method=:euclidean)
    
    # Compute normalized descriptors for both shapes
    desc1 = compute_fourier_descriptors(interface1, num_descriptors=num_descriptors, normalize=true)
    desc2 = compute_fourier_descriptors(interface2, num_descriptors=num_descriptors, normalize=true)
    
    # Use magnitudes for comparison (phase-invariant)
    mag1 = desc1.magnitudes
    mag2 = desc2.magnitudes
    
    # Ensure same length
    min_len = min(length(mag1), length(mag2))
    mag1 = mag1[1:min_len]
    mag2 = mag2[1:min_len]
    
    # Calculate distance based on chosen method
    if method == :euclidean
        return norm(mag1 - mag2)
    elseif method == :manhattan
        return sum(abs.(mag1 - mag2))
    elseif method == :cosine
        return 1.0 - dot(mag1, mag2) / (norm(mag1) * norm(mag2))
    else
        error("Unknown distance method: $method")
    end
end

"""
    elliptic_fourier_descriptors(interface::Interface; 
                               num_harmonics=10,
                               normalize=true)

Compute Elliptic Fourier Descriptors for a closed curve.
These provide a more intuitive description of shape contours.

# Arguments
- `interface`: The Interface object to analyze
- `num_harmonics`: Number of harmonics to compute
- `normalize`: Whether to normalize for size/rotation invariance

# Returns
A NamedTuple with coefficient matrices A, B, C, and D
"""
function elliptic_fourier_descriptors(interface::Interface; 
                                    num_harmonics=10,
                                    normalize=true)
    # Get marker positions
    positions = [marker.position for marker in interface.markers]
    
    # Calculate path segments (differences between consecutive points)
    Δxs = diff([last(positions)[1]; [pos[1] for pos in positions]])
    Δys = diff([last(positions)[2]; [pos[2] for pos in positions]])
    
    # Calculate path lengths (Euclidean distance for each segment)
    Δts = sqrt.(Δxs.^2 + Δys.^2)
    
    # Cumulative path length
    t = [0; cumsum(Δts)]
    T = t[end]  # Total path length
    
    # Initialize coefficient matrices
    A = zeros(num_harmonics)
    B = zeros(num_harmonics)
    C = zeros(num_harmonics)
    D = zeros(num_harmonics)
    
    # Compute coefficients for each harmonic
    for n in 1:num_harmonics
        # Calculate coefficients based on elliptic Fourier analysis
        A[n] = (T / (2 * n^2 * π^2)) * sum(Δxs .* (cos.(2 * n * π * t[2:end] / T) - cos.(2 * n * π * t[1:end-1] / T)) ./ Δts)
        B[n] = (T / (2 * n^2 * π^2)) * sum(Δxs .* (sin.(2 * n * π * t[2:end] / T) - sin.(2 * n * π * t[1:end-1] / T)) ./ Δts)
        C[n] = (T / (2 * n^2 * π^2)) * sum(Δys .* (cos.(2 * n * π * t[2:end] / T) - cos.(2 * n * π * t[1:end-1] / T)) ./ Δts)
        D[n] = (T / (2 * n^2 * π^2)) * sum(Δys .* (sin.(2 * n * π * t[2:end] / T) - sin.(2 * n * π * t[1:end-1] / T)) ./ Δts)
    end
    
    # Add DC components (average x and y)
    A₀ = sum(Δxs .* (t[2:end] - t[1:end-1])) / T
    C₀ = sum(Δys .* (t[2:end] - t[1:end-1])) / T
    
    # Normalize if requested
    if normalize && num_harmonics > 0
        # Rotation invariance based on first harmonic
        θ = 0.5 * atan(2 * (A[1] * B[1] + C[1] * D[1]), A[1]^2 + C[1]^2 - B[1]^2 - D[1]^2)
        
        # Scale factor for size invariance
        scale = sqrt(A[1]^2 + B[1]^2 + C[1]^2 + D[1]^2)
        
        # Rotation matrix
        cosθ, sinθ = cos(θ), sin(θ)
        
        # Apply normalization
        for n in 1:num_harmonics
            An, Bn = A[n], B[n]
            Cn, Dn = C[n], D[n]
            
            A[n] = (An * cosθ + Bn * sinθ) / scale
            B[n] = (Bn * cosθ - An * sinθ) / scale
            C[n] = (Cn * cosθ + Dn * sinθ) / scale
            D[n] = (Dn * cosθ - Cn * sinθ) / scale
        end
        
        # Reset the transformed offset (for translation invariance)
        A₀, C₀ = 0.0, 0.0
    end
    
    return (
        A₀ = A₀,
        C₀ = C₀,
        A = A,
        B = B,
        C = C,
        D = D
    )
end

"""
    reconstruct_from_elliptic_fourier(efd::NamedTuple, num_points=100, num_harmonics=nothing)

Reconstruct a shape from elliptic Fourier descriptors.

# Arguments
- `efd`: Elliptic Fourier Descriptors from elliptic_fourier_descriptors
- `num_points`: Number of points in reconstruction
- `num_harmonics`: Number of harmonics to use (defaults to all)

# Returns
A new Interface object reconstructed from the descriptors
"""
function reconstruct_from_elliptic_fourier(efd::NamedTuple, num_points=100, num_harmonics=nothing)
    # Use all harmonics if not specified
    if isnothing(num_harmonics)
        num_harmonics = length(efd.A)
    end
    
    # Limit to available harmonics
    num_harmonics = min(num_harmonics, length(efd.A))
    
    # Parametric angle values
    t = range(0, 2π, length=num_points+1)[1:num_points]
    
    # Initialize coordinates
    x = zeros(num_points)
    y = zeros(num_points)
    
    # Add DC components
    x .+= efd.A₀
    y .+= efd.C₀
    
    # Add each harmonic contribution
    for n in 1:num_harmonics
        x .+= efd.A[n] .* cos.(n .* t) .+ efd.B[n] .* sin.(n .* t)
        y .+= efd.C[n] .* cos.(n .* t) .+ efd.D[n] .* sin.(n .* t)
    end
    
    # Create markers and connectivity
    positions = [[x[i], y[i]] for i in 1:num_points]
    markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
    connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
    
    # Create reconstructed interface
    return Interface(markers, connectivity, true)
end