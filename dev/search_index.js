var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TrixiShallowWater","category":"page"},{"location":"#TrixiShallowWater","page":"Home","title":"TrixiShallowWater","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TrixiShallowWater.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TrixiShallowWater]","category":"page"},{"location":"#TrixiShallowWater.flux_hll_chen_noelle","page":"Home","title":"TrixiShallowWater.flux_hll_chen_noelle","text":"flux_hll_chen_noelle = FluxHLL(min_max_speed_chen_noelle)\n\nAn instance of FluxHLL specific to the shallow water equations that uses the wave speed estimates from min_max_speed_chen_noelle. This HLL flux is guaranteed to have zero numerical mass flux out of a \"dry\" element, maintain positivity of the water height, and satisfy an entropy inequality.\n\nFor complete details see Section 2.4 of the following reference\n\nGuoxian Chen and Sebastian Noelle (2017) A new hydrostatic reconstruction scheme based on subcell reconstructions DOI: 10.1137/15M1053074\n\n\n\n\n\n","category":"constant"},{"location":"#TrixiShallowWater.IndicatorHennemannGassnerShallowWater","page":"Home","title":"TrixiShallowWater.IndicatorHennemannGassnerShallowWater","text":"IndicatorHennemannGassnerShallowWater(equations::AbstractEquations, basis;\n                                      alpha_max=0.5,\n                                      alpha_min=0.001,\n                                      alpha_smooth=true,\n                                      variable)\n\nModified version of the IndicatorHennemannGassner indicator used for shock-capturing for shallow water equations. After the element-wise values for the blending factors are computed an additional check is made to see if the element is partially wet. In this case, partially wet elements are set to use the pure finite volume scheme that is guaranteed to be well-balanced for this wet/dry transition state of the flow regime.\n\nSee also VolumeIntegralShockCapturingHG.\n\nReferences\n\nHennemann, Gassner (2020) \"A provably entropy stable subcell shock capturing approach for high order split form DG\" arXiv: 2008.12044\n\n\n\n\n\n","category":"type"},{"location":"#TrixiShallowWater.PositivityPreservingLimiterShallowWater","page":"Home","title":"TrixiShallowWater.PositivityPreservingLimiterShallowWater","text":"PositivityPreservingLimiterShallowWater(; variables)\n\nThe limiter is specifically designed for the shallow water equations. It is applied to all scalar variables in their given order using the defined threshold_limiter from the ShallowWaterEquationsWetDry1D struct or the ShallowWaterEquationsWetDry2D struct to determine the minimal acceptable values. The order of the variables is important and might have a strong influence on the robustness.\n\nAs opposed to the standard version of the PositivityPreservingLimiterZhangShu, nodes with a water height below the threshold_limiter are treated in a special way. To avoid numerical problems caused by velocities close to zero, the velocity is cut off, such that the node can be identified as \"dry\". The special feature of the ShallowWaterEquationsWetDry used here is that the bottom topography is stored as an additional quantity in the solution vector u. However, the value of the bottom topography should not be changed. That is why, it is not limited.\n\nAfter the limiting process is applied to all degrees of freedom, for safety reasons, the threshold_limiter is applied again on all the DG nodes in order to avoid water height below. In the case where the cell mean value is below the threshold before applying the limiter, there could still be dry nodes afterwards due to the logic of the limiter.\n\nThis fully-discrete positivity-preserving limiter is based on the work of\n\nZhang, Shu (2011) Maximum-principle-satisfying and positivity-preserving high-order schemes for conservation laws: survey and new developments doi: 10.1098/rspa.2011.0153\n\n\n\n\n\n","category":"type"},{"location":"#TrixiShallowWater.ShallowWaterEquationsWetDry1D","page":"Home","title":"TrixiShallowWater.ShallowWaterEquationsWetDry1D","text":"ShallowWaterEquationsWetDry1D(; gravity, H0 = 0, threshold_limiter = nothing threshold_wet = nothing)\n\nShallow water equations (SWE) in one space dimension. The equations are given by\n\nbeginaligned\n  fracpartial hpartial t + fracpartialpartial x(h v) = 0 \n    fracpartialpartial t(h v) + fracpartialpartial xleft(h v^2 + fracg2h^2right)\n    + g h fracpartial bpartial x = 0\nendaligned\n\nThe unknown quantities of the SWE are the water height h and the velocity v. The gravitational constant is denoted by g and the (possibly) variable bottom topography function b(x). Conservative variable water height h is measured from the bottom topography b, therefore one also defines the total water height as H = h + b.\n\nThe additional quantity H_0 is also available to store a reference value for the total water height that is useful to set initial conditions or test the \"lake-at-rest\" well-balancedness.\n\nAlso, there are two thresholds which prevent numerical problems as well as instabilities. Both of them do not have to be passed, as default values are defined within the struct. The first one, threshold_limiter, is used in PositivityPreservingLimiterShallowWater on the water height, as a (small) shift on the initial condition and cutoff before the next time step. The second one, threshold_wet, is applied on the water height to define when the flow is \"wet\" before calculating the numerical flux.\n\nThe bottom topography function b(x) is set inside the initial condition routine for a particular problem setup. To test the conservative form of the SWE one can set the bottom topography variable b to zero.\n\nIn addition to the unknowns, Trixi.jl currently stores the bottom topography values at the approximation points despite being fixed in time. This is done for convenience of computing the bottom topography gradients on the fly during the approximation as well as computing auxiliary quantities like the total water height H or the entropy variables. This affects the implementation and use of these equations in various ways:\n\nThe flux values corresponding to the bottom topography must be zero.\nThe bottom topography values must be included when defining initial conditions, boundary conditions or source terms.\nAnalysisCallback analyzes this variable.\nTrixi.jl's visualization tools will visualize the bottom topography by default.\n\nReferences for the SWE are many but a good introduction is available in Chapter 13 of the book:\n\nRandall J. LeVeque (2002) Finite Volume Methods for Hyperbolic Problems DOI: 10.1017/CBO9780511791253\n\n\n\n\n\n","category":"type"},{"location":"#TrixiShallowWater.ShallowWaterEquationsWetDry2D","page":"Home","title":"TrixiShallowWater.ShallowWaterEquationsWetDry2D","text":"ShallowWaterEquationsWetDry2D(; gravity, H0 = 0, threshold_limiter = nothing, threshold_wet = nothing)\n\nShallow water equations (SWE) in two space dimensions. The equations are given by\n\nbeginaligned\n  fracpartial hpartial t + fracpartialpartial x(h v_1)\n    + fracpartialpartial y(h v_2) = 0 \n    fracpartialpartial t(h v_1) + fracpartialpartial xleft(h v_1^2 + fracg2h^2right)\n    + fracpartialpartial y(h v_1 v_2) + g h fracpartial bpartial x = 0 \n    fracpartialpartial t(h v_2) + fracpartialpartial x(h v_1 v_2)\n    + fracpartialpartial yleft(h v_2^2 + fracg2h^2right) + g h fracpartial bpartial y = 0\nendaligned\n\nThe unknown quantities of the SWE are the water height h and the velocities mathbfv = (v_1 v_2)^T. The gravitational constant is denoted by g and the (possibly) variable bottom topography function b(xy). Conservative variable water height h is measured from the bottom topography b, therefore one also defines the total water height as H = h + b.\n\nThe additional quantity H_0 is also available to store a reference value for the total water height that is useful to set initial conditions or test the \"lake-at-rest\" well-balancedness.\n\nAlso, there are two thresholds which prevent numerical problems as well as instabilities. Both of them do not have to be passed, as default values are defined within the struct. The first one, threshold_limiter, is used in PositivityPreservingLimiterShallowWater on the water height, as a (small) shift on the initial condition and cutoff before the next time step. The second one, threshold_wet, is applied on the water height to define when the flow is \"wet\" before calculating the numerical flux.\n\nThe bottom topography function b(xy) is set inside the initial condition routine for a particular problem setup. To test the conservative form of the SWE one can set the bottom topography variable b to zero.\n\nIn addition to the unknowns, TrixiShallowWater.jl currently stores the bottom topography values at the approximation points despite being fixed in time. This is done for convenience of computing the bottom topography gradients on the fly during the approximation as well as computing auxiliary quantities like the total water height H or the entropy variables. This affects the implementation and use of these equations in various ways:\n\nThe flux values corresponding to the bottom topography must be zero.\nThe bottom topography values must be included when defining initial conditions, boundary conditions or source terms.\nAnalysisCallback analyzes this variable.\nTrixi.jl's visualization tools will visualize the bottom topography by default.\n\nReferences for the SWE are many but a good introduction is available in Chapter 13 of the book:\n\nRandall J. LeVeque (2002) Finite Volume Methods for Hyperbolic Problems DOI: 10.1017/CBO9780511791253\n\n\n\n\n\n","category":"type"},{"location":"#Trixi.boundary_condition_slip_wall-Tuple{Any, AbstractVector, Any, Any, Any, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.boundary_condition_slip_wall","text":"boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function,\n                             equations::ShallowWaterEquationsWetDry2D)\n\nCreate a boundary state by reflecting the normal velocity component and keep the tangential velocity component unchanged. The boundary water height is taken from the internal value. For details see Section 9.2.5 of the book:\n\nEleuterio F. Toro (2001) Shock-Capturing Methods for Free-Surface Shallow Flows 1st edition ISBN 0471987662\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.boundary_condition_slip_wall-Tuple{Any, Any, Any, Any, Any, Any, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.boundary_condition_slip_wall","text":"boundary_condition_slip_wall(u_inner, orientation_or_normal, x, t, surface_flux_function,\n                              equations::ShallowWaterEquationsWetDry1D)\n\nCreate a boundary state by reflecting the normal velocity component and keep the tangential velocity component unchanged. The boundary water height is taken from the internal value.\n\nFor details see Section 9.2.5 of the book:\n\nEleuterio F. Toro (2001) Shock-Capturing Methods for Free-Surface Shallow Flows 1st edition ISBN 0471987662\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.boundary_condition_slip_wall-Tuple{Any, Any, Any, Any, Any, Any, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.boundary_condition_slip_wall","text":"boundary_condition_slip_wall(u_inner, orientation, direction, x, t,\n                             surface_flux_function, equations::ShallowWaterEquationsWetDry2D)\n\nShould be used together with TreeMesh.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_fjordholm_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.flux_fjordholm_etal","text":"flux_fjordholm_etal(u_ll, u_rr, orientation,\n                    equations::ShallowWaterEquationsWetDry1D)\n\nTotal energy conservative (mathematical entropy for shallow water equations). When the bottom topography is nonzero this should only be used as a surface flux otherwise the scheme will not be well-balanced. For well-balancedness in the volume flux use flux_wintermeyer_etal.\n\nDetails are available in Eq. (4.1) in the paper:\n\nUlrik S. Fjordholm, Siddhartha Mishr and Eitan Tadmor (2011) Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography DOI: 10.1016/j.jcp.2011.03.042\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_fjordholm_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.flux_fjordholm_etal","text":"flux_fjordholm_etal(u_ll, u_rr, orientation_or_normal_direction,\n                    equations::ShallowWaterEquationsWetDry2D)\n\nTotal energy conservative (mathematical entropy for shallow water equations). When the bottom topography is nonzero this should only be used as a surface flux otherwise the scheme will not be well-balanced. For well-balancedness in the volume flux use flux_wintermeyer_etal.\n\nDetails are available in Eq. (4.1) in the paper:\n\nUlrik S. Fjordholm, Siddhartha Mishr and Eitan Tadmor (2011) Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography DOI: 10.1016/j.jcp.2011.03.042\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_audusse_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.flux_nonconservative_audusse_etal","text":"flux_nonconservative_audusse_etal(u_ll, u_rr, orientation::Integer,\n                                  equations::ShallowWaterEquationsWetDry1D)\n\nNon-symmetric two-point surface flux that discretizes the nonconservative (source) term. The discretization uses the hydrostatic_reconstruction_audusse_etal on the conservative variables.\n\nThis hydrostatic reconstruction ensures that the finite volume numerical fluxes remain well-balanced for discontinuous bottom topographies ShallowWaterEquationsWetDry1D. Should be used together with FluxHydrostaticReconstruction and hydrostatic_reconstruction_audusse_etal in the surface flux to ensure consistency.\n\nFurther details on the hydrostatic reconstruction and its motivation can be found in\n\nEmmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004) A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows DOI: 10.1137/S1064827503431090\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_audusse_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.flux_nonconservative_audusse_etal","text":"flux_nonconservative_audusse_etal(u_ll, u_rr, orientation::Integer,\n                                  equations::ShallowWaterEquationsWetDry2D)\nflux_nonconservative_audusse_etal(u_ll, u_rr,\n                                  normal_direction_ll     ::AbstractVector,\n                                  normal_direction_average::AbstractVector,\n                                  equations::ShallowWaterEquationsWetDry2D)\n\nNon-symmetric two-point surface flux that discretizes the nonconservative (source) term. The discretization uses the hydrostatic_reconstruction_audusse_etal on the conservative variables.\n\nThis hydrostatic reconstruction ensures that the finite volume numerical fluxes remain well-balanced for discontinuous bottom topographies ShallowWaterEquationsWetDry2D. Should be used together with FluxHydrostaticReconstruction and hydrostatic_reconstruction_audusse_etal in the surface flux to ensure consistency.\n\nFurther details for the hydrostatic reconstruction and its motivation can be found in\n\nEmmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004) A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows DOI: 10.1137/S1064827503431090\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_ersing_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.flux_nonconservative_ersing_etal","text":"flux_nonconservative_ersing_etal(u_ll, u_rr, orientation::Integer,\n                                 equations::ShallowWaterEquationsWetDry1D)\n\nNon-symmetric path-conservative two-point volume flux discretizing the nonconservative (source) term that contains the gradient of the bottom topography ShallowWaterEquationsWetDry1D.\n\nThis is a modified version of flux_nonconservative_wintermeyer_etal that gives entropy  conservation and well-balancedness in both the volume and surface when combined with  flux_wintermeyer_etal.\n\nFor further details see:\n\nPatrick Ersing, Andrew R. Winters (2023) An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on  curvilinear meshes DOI: 10.48550/arXiv.2306.12699\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_ersing_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.flux_nonconservative_ersing_etal","text":"flux_nonconservative_ersing_etal(u_ll, u_rr, orientation::Integer,\n                                 equations::ShallowWaterEquationsWetDry2D)\nflux_nonconservative_ersing_etal(u_ll, u_rr,\n                                 normal_direction_ll::AbstractVector,\n                                 normal_direction_average::AbstractVector,\n                                 equations::ShallowWaterEquationsWetDry2D)\n\nwarning: Experimental code\nThis numerical flux is experimental and may change in any future release.\n\nNon-symmetric path-conservative two-point volume flux discretizing the nonconservative (source) term that contains the gradient of the bottom topography ShallowWaterEquationsWetDry2D.\n\nOn curvilinear meshes, this nonconservative flux depends on both the contravariant vector (normal direction) at the current node and the averaged one. This is different from numerical fluxes used to discretize conservative terms.\n\nThis is a modified version of flux_nonconservative_wintermeyer_etal that gives entropy  conservation and well-balancedness in both the volume and surface when combined with  flux_wintermeyer_etal.\n\nFor further details see:\n\nPatrick Ersing, Andrew R. Winters (2023) An entropy stable discontinuous Galerkin method for the two-layer shallow water equations on  curvilinear meshes DOI: 10.48550/arXiv.2306.12699\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_fjordholm_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.flux_nonconservative_fjordholm_etal","text":"flux_nonconservative_fjordholm_etal(u_ll, u_rr, orientation::Integer,\n                                    equations::ShallowWaterEquationsWetDry1D)\n\nNon-symmetric two-point surface flux discretizing the nonconservative (source) term of that contains the gradient of the bottom topography ShallowWaterEquationsWetDry1D.\n\nThis contains additional terms compared to flux_nonconservative_wintermeyer_etal that account for possible discontinuities in the bottom topography function. Thus, this flux should be used in general at interfaces. For flux differencing volume terms, flux_nonconservative_wintermeyer_etal is analytically equivalent but slightly cheaper.\n\nFurther details for the original finite volume formulation are available in\n\nUlrik S. Fjordholm, Siddhartha Mishr and Eitan Tadmor (2011) Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography DOI: 10.1016/j.jcp.2011.03.042\n\nand for curvilinear 2D case in the paper:\n\nNiklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017) An entropy stable nodal discontinuous Galerkin method for the two dimensional shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry DOI: 10.1016/j.jcp.2017.03.036\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_fjordholm_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.flux_nonconservative_fjordholm_etal","text":"flux_nonconservative_fjordholm_etal(u_ll, u_rr, orientation::Integer,\n                                    equations::ShallowWaterEquationsWetDry2D)\nflux_nonconservative_fjordholm_etal(u_ll, u_rr,\n                                    normal_direction_ll     ::AbstractVector,\n                                    normal_direction_average::AbstractVector,\n                                    equations::ShallowWaterEquationsWetDry2D)\n\nNon-symmetric two-point surface flux discretizing the nonconservative (source) term of that contains the gradient of the bottom topography ShallowWaterEquationsWetDry2D.\n\nOn curvilinear meshes, this nonconservative flux depends on both the contravariant vector (normal direction) at the current node and the averaged one. This is different from numerical fluxes used to discretize conservative terms.\n\nThis contains additional terms compared to flux_nonconservative_wintermeyer_etal that account for possible discontinuities in the bottom topography function. Thus, this flux should be used in general at interfaces. For flux differencing volume terms, flux_nonconservative_wintermeyer_etal is analytically equivalent but slightly cheaper.\n\nFurther details for the original finite volume formulation are available in\n\nUlrik S. Fjordholm, Siddhartha Mishr and Eitan Tadmor (2011) Well-balanced and energy stable schemes for the shallow water equations with discontinuous topography DOI: 10.1016/j.jcp.2011.03.042\n\nand for curvilinear 2D case in the paper:\n\nNiklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017) An entropy stable nodal discontinuous Galerkin method for the two dimensional shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry DOI: 10.1016/j.jcp.2017.03.036\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_wintermeyer_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.flux_nonconservative_wintermeyer_etal","text":"flux_nonconservative_wintermeyer_etal(u_ll, u_rr, orientation::Integer,\n                                      equations::ShallowWaterEquationsWetDry1D)\n\nNon-symmetric two-point volume flux discretizing the nonconservative (source) term that contains the gradient of the bottom topography ShallowWaterEquationsWetDry1D.\n\nFurther details are available in the paper:#include(\"numerical_fluxes.jl\")\n\nNiklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017) An entropy stable nodal discontinuous Galerkin method for the two dimensional shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry DOI: 10.1016/j.jcp.2017.03.036\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_nonconservative_wintermeyer_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.flux_nonconservative_wintermeyer_etal","text":"flux_nonconservative_wintermeyer_etal(u_ll, u_rr, orientation::Integer,\n                                      equations::ShallowWaterEquationsWetDry2D)\nflux_nonconservative_wintermeyer_etal(u_ll, u_rr,\n                                      normal_direction_ll     ::AbstractVector,\n                                      normal_direction_average::AbstractVector,\n                                      equations::ShallowWaterEquationsWetDry2D)\n\nNon-symmetric two-point volume flux discretizing the nonconservative (source) term that contains the gradient of the bottom topography ShallowWaterEquationsWetDry2D.\n\nOn curvilinear meshes, this nonconservative flux depends on both the contravariant vector (normal direction) at the current node and the averaged one. This is different from numerical fluxes used to discretize conservative terms.\n\nFurther details are available in the paper:\n\nNiklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017) An entropy stable nodal discontinuous Galerkin method for the two dimensional shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry DOI: 10.1016/j.jcp.2017.03.036\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_wintermeyer_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.flux_wintermeyer_etal","text":"flux_wintermeyer_etal(u_ll, u_rr, orientation,\n                      equations::ShallowWaterEquationsWetDry1D)\n\nTotal energy conservative (mathematical entropy for shallow water equations) split form. When the bottom topography is nonzero this scheme will be well-balanced when used as a volume_flux. The surface_flux should still use, e.g., flux_fjordholm_etal.\n\nFurther details are available in Theorem 1 of the paper:\n\nNiklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017) An entropy stable nodal discontinuous Galerkin method for the two dimensional shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry DOI: 10.1016/j.jcp.2017.03.036\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.flux_wintermeyer_etal-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.flux_wintermeyer_etal","text":"flux_wintermeyer_etal(u_ll, u_rr, orientation_or_normal_direction,\n                      equations::ShallowWaterEquationsWetDry2D)\n\nTotal energy conservative (mathematical entropy for shallow water equations) split form. When the bottom topography is nonzero this scheme will be well-balanced when used as a volume_flux. The surface_flux should still use, e.g., flux_fjordholm_etal.\n\nFurther details are available in Theorem 1 of the paper:\n\nNiklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and David A. Kopriva (2017) An entropy stable nodal discontinuous Galerkin method for the two dimensional shallow water equations on unstructured curvilinear meshes with discontinuous bathymetry DOI: 10.1016/j.jcp.2017.03.036\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.hydrostatic_reconstruction_audusse_etal-Tuple{Any, Any, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.hydrostatic_reconstruction_audusse_etal","text":"hydrostatic_reconstruction_audusse_etal(u_ll, u_rr, orientation::Integer,\n                                        equations::ShallowWaterEquationsWetDry1D)\n\nA particular type of hydrostatic reconstruction on the water height to guarantee well-balancedness for a general bottom topography ShallowWaterEquationsWetDry1D. The reconstructed solution states u_ll_star and u_rr_star variables are then used to evaluate the surface numerical flux at the interface. Use in combination with the generic numerical flux routine FluxHydrostaticReconstruction.\n\nFurther details on this hydrostatic reconstruction and its motivation can be found in\n\nEmmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004) A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows DOI: 10.1137/S1064827503431090\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.hydrostatic_reconstruction_audusse_etal-Tuple{Any, Any, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.hydrostatic_reconstruction_audusse_etal","text":"hydrostatic_reconstruction_audusse_etal(u_ll, u_rr, orientation_or_normal_direction,\n                                        equations::ShallowWaterEquationsWetDry2D)\n\nA particular type of hydrostatic reconstruction on the water height to guarantee well-balancedness for a general bottom topography ShallowWaterEquationsWetDry2D. The reconstructed solution states u_ll_star and u_rr_star variables are used to evaluate the surface numerical flux at the interface. Use in combination with the generic numerical flux routine FluxHydrostaticReconstruction.\n\nFurther details for the hydrostatic reconstruction and its motivation can be found in\n\nEmmanuel Audusse, François Bouchut, Marie-Odile Bristeau, Rupert Klein, and Benoit Perthame (2004) A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows DOI: 10.1137/S1064827503431090\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.initial_condition_convergence_test-Tuple{Any, Any, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.initial_condition_convergence_test","text":"initial_condition_convergence_test(x, t, equations::ShallowWaterEquationsWetDry2D)\n\nA smooth initial condition used for convergence tests in combination with source_terms_convergence_test (and BoundaryConditionDirichlet(initial_condition_convergence_test) in non-periodic domains).\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.initial_condition_weak_blast_wave-Tuple{Any, Any, ShallowWaterEquationsWetDry1D}","page":"Home","title":"Trixi.initial_condition_weak_blast_wave","text":"initial_condition_weak_blast_wave(x, t, equations::ShallowWaterEquationsWetDry1D)\n\nA weak blast wave discontinuity useful for testing, e.g., total energy conservation. Note for the shallow water equations to the total energy acts as a mathematical entropy function.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.initial_condition_weak_blast_wave-Tuple{Any, Any, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.initial_condition_weak_blast_wave","text":"initial_condition_weak_blast_wave(x, t, equations::ShallowWaterEquationsWetDry2D)\n\nA weak blast wave discontinuity useful for testing, e.g., total energy conservation. Note for the shallow water equations to the total energy acts as a mathematical entropy function.\n\n\n\n\n\n","category":"method"},{"location":"#Trixi.source_terms_convergence_test-Tuple{Any, Any, Any, ShallowWaterEquationsWetDry2D}","page":"Home","title":"Trixi.source_terms_convergence_test","text":"source_terms_convergence_test(u, x, t, equations::ShallowWaterEquationsWetDry2D)\n\nSource terms used for convergence tests in combination with initial_condition_convergence_test (and BoundaryConditionDirichlet(initial_condition_convergence_test) in non-periodic domains).\n\nThis manufactured solution source term is specifically designed for the bottom topography function b(x,y) = 2 + 0.5 * sin(sqrt(2)*pi*x) + 0.5 * sin(sqrt(2)*pi*y) as defined in initial_condition_convergence_test.\n\n\n\n\n\n","category":"method"},{"location":"#TrixiShallowWater.flux_nonconservative_chen_noelle-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"TrixiShallowWater.flux_nonconservative_chen_noelle","text":"flux_nonconservative_chen_noelle(u_ll, u_rr,\n                                 orientation::Integer,\n                                 equations::ShallowWaterEquationsWetDry1D)\n\nNon-symmetric two-point surface flux that discretizes the nonconservative (source) term. The discretization uses the hydrostatic_reconstruction_chen_noelle on the conservative variables.\n\nShould be used together with FluxHydrostaticReconstruction and hydrostatic_reconstruction_chen_noelle in the surface flux to ensure consistency.\n\nFurther details on the hydrostatic reconstruction and its motivation can be found in\n\nGuoxian Chen and Sebastian Noelle (2017) A new hydrostatic reconstruction scheme based on subcell reconstructions DOI:10.1137/15M1053074\n\n\n\n\n\n","category":"method"},{"location":"#TrixiShallowWater.flux_nonconservative_chen_noelle-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"TrixiShallowWater.flux_nonconservative_chen_noelle","text":"flux_nonconservative_chen_noelle(u_ll, u_rr,\n                                 orientation::Integer,\n                                 equations::ShallowWaterEquationsWetDry2D)\nflux_nonconservative_chen_noelle(u_ll, u_rr,\n                                 normal_direction_ll      ::AbstractVector,\n                                 normal_direction_average ::AbstractVector,\n                                 equations::ShallowWaterEquationsWetDry2D)\n\nNon-symmetric two-point surface flux that discretizes the nonconservative (source) term. The discretization uses the hydrostatic_reconstruction_chen_noelle on the conservative variables.\n\nShould be used together with FluxHydrostaticReconstruction and hydrostatic_reconstruction_chen_noelle in the surface flux to ensure consistency.\n\nFurther details on the hydrostatic reconstruction and its motivation can be found in\n\nGuoxian Chen and Sebastian Noelle (2017) A new hydrostatic reconstruction scheme based on subcell reconstructions DOI:10.1137/15M1053074\n\n\n\n\n\n","category":"method"},{"location":"#TrixiShallowWater.hydrostatic_reconstruction_chen_noelle-Tuple{Any, Any, ShallowWaterEquationsWetDry1D}","page":"Home","title":"TrixiShallowWater.hydrostatic_reconstruction_chen_noelle","text":"hydrostatic_reconstruction_chen_noelle(u_ll, u_rr, orientation::Integer,\n                                       equations::ShallowWaterEquationsWetDry1D)\n\nA particular type of hydrostatic reconstruction of the water height to guarantee well-balancedness for a general bottom topography of the ShallowWaterEquationsWetDry1D. The reconstructed solution states u_ll_star and u_rr_star variables are used to evaluate the surface numerical flux at the interface. The key idea is a linear reconstruction of the bottom and water height at the interfaces using subcells. Use in combination with the generic numerical flux routine FluxHydrostaticReconstruction.\n\nFurther details on this hydrostatic reconstruction and its motivation can be found in\n\nGuoxian Chen and Sebastian Noelle (2017) A new hydrostatic reconstruction scheme based on subcell reconstructions DOI:10.1137/15M1053074\n\n\n\n\n\n","category":"method"},{"location":"#TrixiShallowWater.hydrostatic_reconstruction_chen_noelle-Tuple{Any, Any, ShallowWaterEquationsWetDry2D}","page":"Home","title":"TrixiShallowWater.hydrostatic_reconstruction_chen_noelle","text":"hydrostatic_reconstruction_chen_noelle(u_ll, u_rr, orientation::Integer,\n                                       equations::ShallowWaterEquationsWetDry2D)\n\nA particular type of hydrostatic reconstruction of the water height to guarantee well-balancedness for a general bottom topography of the ShallowWaterEquationsWetDry2D. The reconstructed solution states u_ll_star and u_rr_star variables are then used to evaluate the surface numerical flux at the interface. The key idea is a linear reconstruction of the bottom and water height at the interfaces using subcells. Use in combination with the generic numerical flux routine FluxHydrostaticReconstruction.\n\nFurther details on this hydrostatic reconstruction and its motivation can be found in\n\nGuoxian Chen and Sebastian Noelle (2017) A new hydrostatic reconstruction scheme based on subcell reconstructions DOI:10.1137/15M1053074\n\n\n\n\n\n","category":"method"},{"location":"#TrixiShallowWater.min_max_speed_chen_noelle-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry1D}","page":"Home","title":"TrixiShallowWater.min_max_speed_chen_noelle","text":"min_max_speed_chen_noelle(u_ll, u_rr, orientation::Integer,\n                          equations::ShallowWaterEquations1D)\n\nThe approximated speeds for the HLL type numerical flux used by Chen and Noelle for their hydrostatic reconstruction. As they state in the paper, these speeds are chosen for the numerical flux to ensure positivity and to satisfy an entropy inequality.\n\nFurther details on this hydrostatic reconstruction and its motivation can be found in\n\nGuoxian Chen and Sebastian Noelle (2017) A new hydrostatic reconstruction scheme based on subcell reconstructions DOI:10.1137/15M1053074\n\n\n\n\n\n","category":"method"},{"location":"#TrixiShallowWater.min_max_speed_chen_noelle-Tuple{Any, Any, Integer, ShallowWaterEquationsWetDry2D}","page":"Home","title":"TrixiShallowWater.min_max_speed_chen_noelle","text":"min_max_speed_chen_noelle(u_ll, u_rr, orientation::Integer,\n                          equations::ShallowWaterEquations2D)\nmin_max_speed_chen_noelle(u_ll, u_rr, normal_direction::AbstractVector,\n                          equations::ShallowWaterEquations2D)\n\nSpecial estimate of the minimal and maximal wave speed of the shallow water equations for the left and right states u_ll, u_rr. These approximate speeds are used for the HLL-type numerical flux flux_hll_chen_noelle. These wave speed estimates together with a particular hydrostatic reconstruction technique guarantee that the numerical flux is positive and satisfies an entropy inequality.\n\nFurther details on this hydrostatic reconstruction and its motivation can be found in the reference below. The definition of the wave speeds are given in Equation (2.20).\n\nGuoxian Chen and Sebastian Noelle (2017) A new hydrostatic reconstruction scheme based on subcell reconstructions DOI:10.1137/15M1053074\n\n\n\n\n\n","category":"method"}]
}
