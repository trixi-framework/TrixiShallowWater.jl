# Tutorials

The tutorials on these pages will guide you through solving shallow water equations using TrixiShallowWater.jl. 
Before diving into the TrixiShallowWater.jl specific tutorials, we recommend familiarizing yourself with 
the core concepts of Trixi.jl by completing the recommended Trixi tutorials below.

## Prerequisites: Trixi.jl Tutorials

To get the most out of the TrixiShallowWater tutorials, we recommend to get familiar with the [Trixi.jl
documentation](https://trixi-framework.github.io/Trixi.jl/stable/) and first complete the following Trixi.jl tutorials:

- [First steps in Trixi.jl](https://trixi-framework.github.io/Trixi.jl/stable/tutorials/first_steps/getting_started/)
- [Introduction to DG Methods](https://trixi-framework.github.io/Trixi.jl/stable/tutorials/scalar_linear_advection_1d/)
- [DGSEM with flux differencing](https://trixi-framework.github.io/Trixi.jl/stable/tutorials/DGSEM_FluxDiff/)
- [Shock capturing with flux differencing and stage limiter](https://trixi-framework.github.io/Trixi.jl/stable/tutorials/shock_capturing/)
- [Non-periodic boundaries](https://trixi-framework.github.io/Trixi.jl/stable/tutorials/non_periodic_boundaries/)
- [Explicit time stepping](https://trixi-framework.github.io/Trixi.jl/stable/tutorials/time_stepping/)

## TrixiShallowWater.jl Tutorials

Once you're comfortable with Trixi.jl, you can start exploring the TrixiShallowWater tutorials. 
These tutorials focus on solving shallow water equations, including wetting and drying scenarios, 
using the features from TrixiShallowWater.jl. Each tutorial is designed to be self-contained and includes:
- A detailed explanation of the problem setup.
- Step-by-step instructions for configuring and running the simulation.
- Visualization and analysis of the results.

### Available Tutorials

- [Dam Break over Triangular Bottom Topography](@ref Dam-break-over-triangular-bottom-topography):
   Simulate a dam break scenario with triangular bottom topography, including wetting and drying effects.
- [Okushiri Tsunami](@ref Okushiri-Tsunami): 
   Learn how to generate a mesh file from real bathymetry data to replicate the Okushiri tsunami experiment
   and visualize the results with Paraview.

## Getting Help

If you encounter any issues or have questions while working through the tutorials, feel free to:
- Check the [Trixi.jl documentation](https://trixi-framework.github.io/Trixi.jl/stable/).
- Open an issue on the [TrixiShallowWater.jl GitHub repository](https://github.com/trixi-framework/TrixiShallowWater.jl/issues).