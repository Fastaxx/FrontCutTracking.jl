module FrontCutTracking

using LinearAlgebra
using Graphs
using MetaGraphs
using CairoMakie
using Colors
using Statistics
using FFTW

# Include all necessary modules
include("core.jl")
include("creation.jl")
include("geometry.jl")
include("graph.jl")
include("vizualization.jl")
include("evolution.jl")
include("modification.jl")
include("statistics.jl")
include("analysis.jl")

# Export all necessary functions and types
export Marker, Interface
export create_circle, create_rectangle, create_line
export add_marker!, remove_marker!, redistribute_markers!
export compute_normals, compute_curvature, interface_length
export to_graph, from_graph, get_neighbors
export plot_interface, plot_graph, plot_interface_evolution
export curvature_statistics, mean_curvature, curvature_extrema, curvature_histogram
export refine_by_curvature!, decimate_by_curvature!
export smooth_interface!, smooth_interface_curvature_weighted!, taubin_smooth_interface!
export length_scale_spectrum, dominant_length_scales, filter_by_scale, length_scale_decomposition
export compute_fourier_descriptors, reconstruct_from_fourier
export fourier_power_spectrum, dominant_frequencies
export fourier_filter, fourier_shape_similarity
export elliptic_fourier_descriptors, reconstruct_from_elliptic_fourier

end # module