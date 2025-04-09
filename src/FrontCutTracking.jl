module FrontCutTracking

using LinearAlgebra
using Graphs
using MetaGraphs
using CairoMakie
using Colors
using Statistics

# Include all necessary modules
include("core.jl")
include("creation.jl")
include("geometry.jl")
include("graph.jl")
include("vizualization.jl")
include("evolution.jl")
include("modification.jl")
include("statistics.jl")

# Export all necessary functions and types
export Marker, Interface
export create_circle, create_rectangle, create_line
export add_marker!, remove_marker!, redistribute_markers!
export compute_normals, compute_curvature, interface_length
export to_graph, from_graph, get_neighbors
export plot_interface, plot_graph, plot_interface_evolution
export curvature_statistics, mean_curvature, curvature_extrema, curvature_histogram
export refine_by_curvature!

end # module