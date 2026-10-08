function _insulation_evaporation(conv, atmos, area_convection, insulation_temperature, air_temperature,
                                  insulation_wetness, insulation_test;
                                  gas_fractions, smoothing::SmoothingStrategy=HardBound())
    # Always evaluate the evaporation expression and gate by a numeric mask. With
    # an early `|| return 0.0u"W"` the function's return type becomes a Union of
    # {computed Quantity, literal 0.0u"W"}, which Enzyme's reverse-mode flags as
    # EnzymeRuntimeActivityError. `safe_step` collapses to the exact ternary mask
    # under HardBound and to a smooth Heaviside under SmoothBound.
    evap_pars_ins = AnimalEvaporationParameters(;
        skin_wetness = insulation_wetness,
        eye_fraction = 0.0,
        bare_skin_fraction = 1.0,
    )
    mass_ins = TransferCoefficients(conv.mass_transfer_coefficient.combined, conv.mass_transfer_coefficient.combined, conv.mass_transfer_coefficient.forced)
    flow = evaporation(evap_pars_ins, mass_ins, atmos, area_convection, insulation_temperature, air_temperature;
                       gas_fractions).evaporation_heat_flow
    # scale=0.01 for wetness (typical 1e-3–5e-2 fraction);
    # scale=1e-6u"m" for insulation_test — well below any nonzero physical value
    # (≥3e-6 m for the thinnest realistic fur), so the Heaviside is essentially
    # a step except in the structurally-zero "naked" case.
    mask = safe_step(smoothing, insulation_wetness; scale=0.01) *
           safe_step(smoothing, insulation_test;    scale=1.0e-6u"m")
    return mask * flow
end

function _insulation_conductivity(insulation, insulation_pars, side, insulation_temperature, skin_temperature, σ;
                                  smoothing::SmoothingStrategy=HardBound())
    (; longwave_depth_fraction) = insulation_pars
    insulation_temp_mean = insulation_temperature * 0.7 + skin_temperature * 0.3
    air_k = dry_air_properties(insulation_temp_mean).thermal_conductivity
    side_depth  = _side_value(insulation.fibres, side).depth
    side_fibres = _side_value(insulation_pars, side)
    effective_conductivity = insulation_thermal_conductivity(side_fibres, air_k, side_depth; smoothing).effective_conductivity
    absorption_coefficient = _side_value(insulation.absorption_coefficients, side)
    approx_radiant_temperature = skin_temperature * (1 - longwave_depth_fraction) + insulation_temperature * longwave_depth_fraction
    radiative_conductivity = (16 * σ * approx_radiant_temperature^3) / (3 * absorption_coefficient)
    return (; insulation_conductivity = effective_conductivity + radiative_conductivity, effective_conductivity)
end

# Type-stable selector for the dorsal/ventral field of a `BodyRegionValues`
# / `DorsalVentral` / `InsulationParameters` etc. Dispatching on the
# `BodySide` singleton lets the compiler pick the right field at compile
# time and erase the wrong branch entirely.
_side_value(x, ::Dorsal) = x.dorsal
_side_value(x, ::Ventral) = x.ventral

function _radiation_coefficients(area, view_factors, ϵ, σ, radiant_temperature, env_temps)
    T = env_temps
    F = view_factors
    coeff(vf, t) = area * vf * 4 * ϵ * σ * ((radiant_temperature + t) / 2)^3
    RadiationCoeffs(coeff(F.sky, T.sky), coeff(F.bush, T.bush),
                    coeff(F.vegetation, T.vegetation), coeff(F.ground, T.ground))
end

# Longwave exchange with the parts a part sees, over the fractions of its view they
# take (the `neighbours` of `silhouette_factors`), a tuple of `(; fraction, temperature)`.
_add_neighbour_radiation(base, area, ::Tuple{}, ϵ, σ, radiant_temperature) = base
function _add_neighbour_radiation(base, area, neighbours::Tuple, ϵ, σ, radiant_temperature)
    flow = base
    for nb in neighbours
        coeff = area * nb.fraction * 4 * ϵ * σ * ((radiant_temperature + nb.temperature) / 2)^3
        flow += coeff * (radiant_temperature - nb.temperature)
    end
    return flow
end

# The same exchange as the sums of the coefficients and of the coefficients times the
# neighbour temperatures, for the bare-skin fixed point.
_neighbour_coefficients(area, ::Tuple{}, ϵ, σ, surface_temperature) =
    (; csum = zero(area * 4 * ϵ * σ * surface_temperature^3),
       cTsum = zero(area * 4 * ϵ * σ * surface_temperature^3 * surface_temperature))
function _neighbour_coefficients(area, neighbours::Tuple, ϵ, σ, surface_temperature)
    csum = zero(area * 4 * ϵ * σ * surface_temperature^3)
    cTsum = zero(area * 4 * ϵ * σ * surface_temperature^3 * surface_temperature)
    for nb in neighbours
        c = area * nb.fraction * 4 * ϵ * σ * ((surface_temperature + nb.temperature) / 2)^3
        csum += c
        cTsum += c * nb.temperature
    end
    return (; csum, cTsum)
end
