# FrontCutTracking.jl

[![Build Status](https://github.com/Fastaxx/FrontCutTracking.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Fastaxx/FrontCutTracking.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for creating, analyzing, and evolving geometric interfaces with a focus on curvature-based operations.

## Features

- **Interface Creation**: Easily create circles, rectangles, lines, and custom shapes
- **Curvature Analysis**: Compute and visualize curvature along interfaces
- **Normal Vectors**: Calculate normal vectors at each marker point
- **Interface Evolution**: Simulate curvature-driven flow and other evolution mechanisms
- **Adaptive Refinement**: Automatically add markers in regions of high curvature
- **Marker Redistribution**: Maintain uniform spacing along evolving interfaces
- **Graph Conversion**: Transform interfaces to graph representations for topological analysis
- **Rich Visualizations**: Plot interfaces with customizable styles and animations

## Installation

```julia
using Pkg
Pkg.add("FrontCutTracking")
```

## Quick Start

```julia
using FrontCutTracking
using CairoMakie

# Create a circular interface
circle = create_circle([0.0, 0.0], 1.0, 40)

# Visualize the interface with curvature coloring and normals
fig = plot_interface(circle, show_curvature=true, show_normals=true, 
                   colormap=:viridis, normal_scale=0.1)

# Simulate evolution under curvature flow
interfaces = simulate_interface_evolution(
    circle,
    50,          # 50 time steps 
    0.1,         # dt = 0.1
    evolution_function=evolve_by_curvature_flow,
    motion_factor=-0.5  # negative for inward flow
)

# Create an animation
plot_interface_evolution(
    interfaces,
    show_curvature=true,
    show_normals=true,
    fixed_limits=true,
    filename="circle_collapse.mp4"
)
```

## Showcase

### Interface Types
![Interface Types](/img/interface_with_ids_and_normals.png)

![Graph Representation](/img/interface_graph_circular.png)

### Curvature Visualization
![Curvature Visualization](/img/line_curvature.png)

### Interface Evolution (Curvature Flow)
![Evolution Animation](/img/star_evolution.gif)

### Adaptive Refinement
![Adaptive Refinement](/img/comparison.png)

## Key Functionality

### Creating Interfaces

```julia
# Create standard shapes
circle = create_circle([0.0, 0.0], 1.0, 40)
rectangle = create_rectangle(0.0, 0.0, 2.0, 1.0, 5)
line = create_line([0.0, 0.0], [1.0, 1.0], 10)

# Create custom shapes
t = range(0, 2π, length=50)
positions = [[16*sin(θ)^3, 13*cos(θ) - 5*cos(2θ) - 2*cos(3θ) - cos(4θ)] for θ in t]
markers = [Marker(pos, i) for (i, pos) in enumerate(positions)]
connectivity = [(i, i % length(markers) + 1) for i in 1:length(markers)]
heart = Interface(markers, connectivity, true)
```

### Analyzing Curvature

```julia
# Calculate curvature along the interface
curvatures = compute_curvature(heart)

# Get curvature statistics
stats = curvature_statistics(heart)
println("Mean curvature: $(stats.mean)")
println("Maximum curvature: $(stats.maximum)")

# Find curvature extrema
extrema = curvature_extrema(heart, n=3)
```

### Interface Evolution and Refinement

```julia
# Evolve an interface using curvature flow
evolved = evolve_by_curvature_flow(heart, 0.1, motion_factor=0.5)

# Redistribute markers for more uniform spacing
redistribute_markers!(evolved, target_spacing=0.2)

# Adaptively refine in high-curvature regions
refine_by_curvature!(evolved, relative_threshold=1.2)
```

## Documentation

For detailed documentation and examples, visit our documentation site (coming soon).

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.