# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# Infiltration models
abstract type InfiltrationModel{RealT} end

"""
    HortonModel(f_0, f_c, k)

The Horton infiltration model.
Input parameters are the initial infiltration capacity `f_0`, the final infiltration capacity `f_c`,
and the rate of capacity decrease `k`.
This model assumes an exponential relation to compute the the maximum 
infiltration rate. 

For more details, see the papers:
- Robert E. Horton (1933)
  "The role of infiltration in the hydrologic cycle"
  [DOI:10.1029/TR014i001p00446](https://doi.org/10.1029/TR014i001p00446)
- J. Fernández-Pato, D. Caviedes-Voullième, P. García-Navarro (2016)
  "Rainfall/runoff simulation with 2D full shallow water equations: Sensitivity analysis and
  calibration of infiltration parameters"
  [DOI: 10.1016/j.jhydrol.2016.03.021](http://dx.doi.org/10.1016/j.jhydrol.2016.03.021)

!!! warning "Experimental code"
    This is an experimental feature and may change in future releases
"""
struct HortonModel{RealT} <: InfiltrationModel{RealT}
    f_0::RealT  # initial infiltration capacity (m / s)
    f_c::RealT  # final infiltration capacity (m / s)
    k::RealT    # rate of capacity decrease (1 / s)
end

# Compute the soil infiltration rate for the Horton model
function infiltration_rate(x, t, infiltration_model::HortonModel)
    (; f_0, f_c, k) = infiltration_model
    return f_c + (f_0 - f_c) * exp(-k * t)
end

"""
    GreenAmptModel(K_s, psi, dtheta)

The Green-Ampt infiltration model.
Input parameters are the saturated hydraulic conductivity `K_s`, the average wetting front suction 
head `psi`, and the difference between soil porosity and volumetric water content `dtheta`.
This model requires a nonlinear solve for the  which requires a nonlinear solve for the 
cumulative infiltration F(t). 

For more details, see the papers:
- W. Heber Green and G. A. Ampt (1911)
  "Studies on Soil Physics"
  [DOI:10.1017/S0021859600001441](https://doi.org/10.1017/S0021859600001441)
- J. Fernández-Pato, D. Caviedes-Voullième, P. García-Navarro (2016)
  "Rainfall/runoff simulation with 2D full shallow water equations: Sensitivity analysis and
  calibration of infiltration parameters"
  [DOI: 10.1016/j.jhydrol.2016.03.021](http://dx.doi.org/10.1016/j.jhydrol.2016.03.021)

!!! warning "Experimental code"
    This is an experimental feature and may change in future releases
"""
struct GreenAmptModel{RealT} <: InfiltrationModel{RealT}
    K_s::RealT    # saturated hydraulic conductivity (m / s)
    psi::RealT    # average wetting front suction head (m)
    dtheta::RealT # difference between soil porosity and volumetric water content
end

# Compute the soil infiltration rate for the Green-Ampt model
function infiltration_rate(x, t, infiltration_model::GreenAmptModel)
    (; K_s, psi, dtheta) = infiltration_model

    # Determine cumulative infiltration F(t) iteratively using Steffensen's method.
    fx = ZeroProblem(F -> F - psi * dtheta * log(1 + F / (psi * dtheta)) - K_s * t, 0.0)
    F = solve(fx, Order2())

    return K_s + K_s * psi * dtheta / F
end

"""
    SourceTermsRain(precipitation_rate, infiltration_model, equations)

Source term for precipitation and soil infiltration.

Required input arguments are the `precipitation_rate` in m/s and an `infiltration_model` that accounts
for the soil infiltration. The source term is then computed as the difference between precipitation
and infiltration.

The precipitation rate can be provided either as a constant value or a function of space and time.

The `infiltration_model` can be set to the [`HortonModel`](@ref) and the [`GreenAmptModel`](@ref).

!!! warning "Experimental code"
    This is an experimental feature and may change in future releases
"""
struct SourceTermsRain{F <: Function, InfiltrationT <: InfiltrationModel}
    precipitation_rate::F
    infiltration_model::InfiltrationT
end

function SourceTermsRain(precipitation_rate::RealT,
                         infiltration_model::InfiltrationModel{RealT},
                         equations) where {RealT <: Real}
    # Convert function output to the correct type
    precipitation_rate = convert(RealT, precipitation_rate)

    return SourceTermsRain((x, t) -> precipitation_rate, infiltration_model)
end

function SourceTermsRain(precipitation_rate::Function,
                         infiltration_model::InfiltrationModel{RealT},
                         equations) where {RealT <: Real}
    # Check if the function output is of the correct type
    if !(typeof(precipitation_rate(one(RealT), one(RealT))) == RealT)
        throw(ArgumentError("Source term functions must return a value of type $(RealT)"))
    end

    return SourceTermsRain((x, t) -> precipitation_rate(x, t), infiltration_model)
end
end # @muladd
