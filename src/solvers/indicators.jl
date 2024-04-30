# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

"""
    IndicatorHennemannGassnerShallowWater(equations::AbstractEquations, basis;
                                          alpha_max=0.5,
                                          alpha_min=0.001,
                                          alpha_smooth=true,
                                          variable)

Modified version of the [`Trixi.IndicatorHennemannGassner`](@extref)
indicator used for shock-capturing for shallow water equations. After
the element-wise values for the blending factors are computed an additional check
is made to see if the element is partially wet. In this case, partially wet elements
are set to use the pure finite volume scheme that is guaranteed to be well-balanced
for this wet/dry transition state of the flow regime.

See also [`Trixi.VolumeIntegralShockCapturingHG`](@extref).

## References

- Hennemann, Gassner (2020)
  "A provably entropy stable subcell shock capturing approach for high order split form DG"
  [arXiv: 2008.12044](https://arxiv.org/abs/2008.12044)
"""
struct IndicatorHennemannGassnerShallowWater{RealT <: Real, Variable, Cache} <:
       Trixi.AbstractIndicator
    alpha_max::RealT
    alpha_min::RealT
    alpha_smooth::Bool
    variable::Variable
    cache::Cache
end

# this method is used when the indicator is constructed as for shock-capturing volume integrals
# of the shallow water equations
# It modifies the shock-capturing indicator to use full FV method in dry elements
# or partially dry elements containing a wet/dry transition.
function IndicatorHennemannGassnerShallowWater(equations::Union{Trixi.AbstractShallowWaterEquations,
                                                                AbstractShallowWaterMultiLayerEquations},
                                               basis;
                                               alpha_max = 0.5,
                                               alpha_min = 0.001,
                                               alpha_smooth = true,
                                               variable)
    alpha_max, alpha_min = promote(alpha_max, alpha_min)
    cache = Trixi.create_cache(IndicatorHennemannGassner, equations, basis)
    IndicatorHennemannGassnerShallowWater{typeof(alpha_max), typeof(variable),
                                          typeof(cache)}(alpha_max, alpha_min,
                                                         alpha_smooth, variable, cache)
end

function Base.show(io::IO, indicator::IndicatorHennemannGassnerShallowWater)
    @nospecialize indicator # reduce precompilation time

    print(io, "IndicatorHennemannGassnerShallowWater(")
    print(io, indicator.variable)
    print(io, ", alpha_max=", indicator.alpha_max)
    print(io, ", alpha_min=", indicator.alpha_min)
    print(io, ", alpha_smooth=", indicator.alpha_smooth)
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain",
                   indicator::IndicatorHennemannGassnerShallowWater)
    @nospecialize indicator # reduce precompilation time

    if get(io, :compact, false)
        show(io, indicator)
    else
        setup = [
            "indicator variable" => indicator.variable,
            "max. α" => indicator.alpha_max,
            "min. α" => indicator.alpha_min,
            "smooth α" => (indicator.alpha_smooth ? "yes" : "no"),
        ]
        Trixi.summary_box(io, "IndicatorHennemannGassnerShallowWater", setup)
    end
end

include("indicators_1d.jl")
include("indicators_2d.jl")
end # @muladd
