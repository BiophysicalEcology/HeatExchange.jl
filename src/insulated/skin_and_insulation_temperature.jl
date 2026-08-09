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

function _insulation_conductivity(insulation, insulation_pars, side, insulation_temperature, skin_temperature, longwave_depth_fraction, σ;
                                  smoothing::SmoothingStrategy=HardBound())
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
@inline _side_value(x, ::Dorsal)  = x.dorsal
@inline _side_value(x, ::Ventral) = x.ventral

# Update only the matching side's conductivity field — dispatch picks the
# right setproperties call at compile time.
@inline _set_side_conductivity(c, ::Dorsal,  v) = setproperties(c, (; dorsal  = v))
@inline _set_side_conductivity(c, ::Ventral, v) = setproperties(c, (; ventral = v))

function _radiation_coefficients(area, view_factors, ϵ, σ, radiant_temperature, env_temps)
    T = env_temps
    F = view_factors
    coeff(vf, t) = area * vf * 4 * ϵ * σ * ((radiant_temperature + t) / 2)^3
    RadiationCoeffs(coeff(F.sky, T.sky), coeff(F.bush, T.bush),
                    coeff(F.vegetation, T.vegetation), coeff(F.ground, T.ground))
end

# Inter-part (neighbour) surface exchange — the lumped radiative component of the
# blocked-solid-angle term. A part whose hemisphere is partly occluded by a sibling
# part exchanges longwave with that sibling's outer surface over the neighbour view
# fraction, instead of with the sky/ground it can't see through the sibling. Same
# linearised radiative-conductance form as `_radiation_coefficients`, with the
# neighbour's surface temperature as the far side. The fractions come from
# `view_partition`, where sky + ground + Σ neighbours = 1, so the blocked solid angle
# (and hence its energy) is not lost — it becomes this term. `neighbours` is a tuple
# of `(; fraction, temperature)`.
#
# Dispatched so the empty (free-standing / single-part) case is the *identity* on the
# running radiation total: `radiation_heat_flow` then compiles to exactly the pre-
# coupling expression, adding no new operations on the active/differentiated path
# (the nested-Enzyme IPOPT Hessian is intolerant of even an added `+ zero`). Only a
# part that genuinely has neighbours pays the extra term.
@inline _add_neighbour_radiation(base, area, ::Tuple{}, ϵ, σ, radiant_temperature) = base
@inline function _add_neighbour_radiation(base, area, neighbours::Tuple, ϵ, σ, radiant_temperature)
    flow = base
    for nb in neighbours
        coeff = area * nb.fraction * 4 * ϵ * σ * ((radiant_temperature + nb.temperature) / 2)^3
        flow += coeff * (radiant_temperature - nb.temperature)
    end
    return flow
end

# Neighbour exchange as a linear (coefficient, coefficient·temperature) pair, for the
# bare-skin fixed-point that solves skin temperature as `(Σ cᵢTᵢ) / (Σ cᵢ)`. `csum`
# folds into the denominator, `cTsum` into the numerator, exactly like the sky/ground
# radiation coefficients on that path. Empty → exact zeros (the free-standing case).
@inline _neighbour_coefficients(area, ::Tuple{}, ϵ, σ, surface_temperature) =
    (; csum = zero(area * 4 * ϵ * σ * surface_temperature^3),
       cTsum = zero(area * 4 * ϵ * σ * surface_temperature^3 * surface_temperature))
@inline function _neighbour_coefficients(area, neighbours::Tuple, ϵ, σ, surface_temperature)
    csum = zero(area * 4 * ϵ * σ * surface_temperature^3)
    cTsum = zero(area * 4 * ϵ * σ * surface_temperature^3 * surface_temperature)
    for nb in neighbours
        c = area * nb.fraction * 4 * ϵ * σ * ((surface_temperature + nb.temperature) / 2)^3
        csum += c
        cTsum += c * nb.temperature
    end
    return (; csum, cTsum)
end

# Read the (possibly absent) neighbour exchange list from a packed `environment_vars`.
# Callers that don't set it (every single-body path, the NLP path) get `()`.
@inline _neighbours(environment_vars) =
    hasproperty(environment_vars, :neighbours) ? environment_vars.neighbours : ()

"""
    solve_temperatures(; body, insulation_pars, insulation, geometry_vars, environment_vars, traits, temperature_tolerance, skin_temperature, insulation_temperature)

Solve a part's skin and insulation-surface temperatures at an imposed core temperature.

Insulated parts (`insulation.insulation_test > 0`) root-find the two temperatures on exactly
the residuals the multipart NLP uses — `surface_balance` (the surface energy balance; the
metabolic and respiration terms cancel out of it) and `residual_skin_temperature` — via the
shared `solve_part_heat_balance` primitive (`_solve_temperatures_insulated`). So the
rule-based and NLP paths share one surface-physics implementation (this replaced the former
hand-rolled iterative `solve_with_insulation!`). The bare-skin case (`insulation_test ≤ 0`)
has no insulation shell — the insulated formulation's `log(r_insulation / r_skin)` factors
are singular there — so it keeps its own direct-to-skin balance (`solve_without_insulation!`).

# Returns
NamedTuple with `insulation_temperature`, `skin_temperature`, `flows::HeatFlows`,
`insulation_conductivity`, `tolerance`, `success`, `ntry`.
"""
function solve_temperatures(;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation::InsulationProperties,
    geometry_vars::GeometryVariables,
    environment_vars::NamedTuple,
    traits::NamedTuple,
    temperature_tolerance,
    skin_temperature,
    insulation_temperature,
    geometry=_part_geometry(body),
    smoothing::SmoothingStrategy=HardBound(),
)
    if u"m"(insulation.insulation_test) > zero(u"m"(insulation.insulation_test))
        return _solve_temperatures_insulated(; body, insulation_pars, insulation, geometry_vars,
            environment_vars, traits, temperature_tolerance, skin_temperature,
            insulation_temperature, geometry, smoothing)
    else
        return solve_without_insulation!(body, geometry_vars, environment_vars, traits,
            temperature_tolerance, skin_temperature, insulation_temperature; geometry, smoothing)
    end
end

# Insulated surface solve — root-find skin & insulation on the NLP residuals via the shared
# `solve_part_heat_balance` primitive (replaces the former iterative `solve_with_insulation!`).
function _solve_temperatures_insulated(;
    body::AbstractBody,
    insulation_pars::InsulationParameters,
    insulation::InsulationProperties,
    geometry_vars::GeometryVariables,
    environment_vars::NamedTuple,
    traits::NamedTuple,
    temperature_tolerance,
    skin_temperature,
    insulation_temperature,
    geometry=_part_geometry(body),
    smoothing::SmoothingStrategy=HardBound(),
)
    core_temperature = traits.core_temperature
    # Respiration and metabolic heat cancel out of `surface_balance`, so these two never
    # affect the root — a default (stripped) respiration and a unit metabolic probe just
    # keep `solve_part_heat_balance` on its normal path.
    resp_pars = stripparams(RespirationParameters())
    metabolic_probe = 1.0u"W"

    function surface_residuals(skin, insulation_surface)
        b = solve_part_heat_balance(
            core_temperature, skin, insulation_surface, metabolic_probe;
            body, geometry, insulation_pars, insulation, geometry_vars,
            environment_vars, traits, resp_pars,
            k_flesh = traits.flesh_conductivity, pant = resp_pars.pant,
            skin_wetness = traits.skin_wetness, smoothing,
        )
        surface_balance = b.residual_energy_balance - b.residual_internal_conduction
        return (; surface_balance, residual_skin = b.residual_skin_temperature, balance = b)
    end

    sol = _newton_surface(surface_residuals, skin_temperature, insulation_temperature)
    b = sol.balance

    # Effective insulation conductivity at the converged temperatures (a diagnostic the
    # compartment/coupled solves report).
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    (; insulation_conductivity) = _insulation_conductivity(
        insulation, insulation_pars, geometry_vars.side, sol.insulation, sol.skin,
        geometry_vars.longwave_depth_fraction, σ; smoothing)

    flows = HeatFlows(
        b.convection_heat_flow,
        b.conduction_heat_flow,
        b.net_metabolic_heat_internal,
        b.skin_evaporation_heat_flow,
        b.insulation_evaporation_heat_flow,
        b.radiation_heat_flow,
        b.solar_heat_flow,
        b.sky_radiation_flow,
        b.bush_radiation_flow,
        b.vegetation_radiation_flow,
        b.ground_radiation_flow,
    )
    return (;
        insulation_temperature = sol.insulation,
        skin_temperature = sol.skin,
        flows,
        insulation_conductivity,
        tolerance = temperature_tolerance,
        success = sol.success,
        ntry = sol.iters,
    )
end

# Damped 2×2 Newton for the per-part surface solve. Unknowns are the skin and
# insulation-surface temperatures; residuals are `surface_balance` (W) and `residual_skin`
# (K), each nondimensionalised by its own unit so the Jacobian is dimensionless. The
# Jacobian is a finite difference — this path is not differentiated (the NLP uses
# `part_surface_residuals` directly), so a numeric Jacobian keeps the solve self-contained.
# Each step is backtracked until the residual norm strictly decreases and stays finite,
# which keeps the iterate out of the unphysical (negative-temperature → NaN) region the raw
# Newton step can overshoot into; a singular/degenerate Jacobian falls back to a bounded
# fixed-point nudge (skin toward its own balance, insulation toward energy balance).
function _newton_surface(f, skin0, insulation0; tol=1.0e-9, maxiter=200)
    Tunit = oneunit(skin0)
    scale1(x) = x / oneunit(x)
    resnorm(rr) = hypot(scale1(rr.surface_balance), scale1(rr.residual_skin))
    s = skin0 / Tunit
    n = insulation0 / Tunit
    r = f(s * Tunit, n * Tunit)
    nr = resnorm(r)
    success = false
    iters = 0
    for iter in 1:maxiter
        iters = iter
        if nr < tol
            success = true
            break
        end
        R1 = scale1(r.surface_balance)
        R2 = scale1(r.residual_skin)
        h = 1.0e-4
        rs = f((s + h) * Tunit, n * Tunit)
        rn = f(s * Tunit, (n + h) * Tunit)
        J11 = (scale1(rs.surface_balance) - R1) / h
        J21 = (scale1(rs.residual_skin)   - R2) / h
        J12 = (scale1(rn.surface_balance) - R1) / h
        J22 = (scale1(rn.residual_skin)   - R2) / h
        det = J11 * J22 - J12 * J21
        if !isfinite(det) || abs(det) < 1.0e-12
            # Degenerate Jacobian: bounded fixed-point nudge. `residual_skin = skin −
            # skin_from_balance`, so `s -= R2` is the skin fixed point; nudge insulation
            # down-gradient of the energy balance.
            ds = clamp(R2, -1.0, 1.0)
            dn = clamp(R1, -1.0, 1.0)
        else
            ds = ( J22 * R1 - J12 * R2) / det
            dn = (-J21 * R1 + J11 * R2) / det
        end
        # Backtracking line search: shrink the step until the residual norm decreases and
        # both residuals are finite (rejects overshoots into the NaN region).
        α = 1.0
        snew, nnew, rnew, nrnew = s, n, r, nr
        for _ in 1:40
            snew = s - α * ds
            nnew = n - α * dn
            rnew = f(snew * Tunit, nnew * Tunit)
            nrnew = resnorm(rnew)
            (isfinite(nrnew) && nrnew < nr) && break
            α *= 0.5
        end
        s, n, r, nr = snew, nnew, rnew, nrnew
    end
    return (; skin = s * Tunit, insulation = n * Tunit, balance = r.balance, success, iters)
end

# --- Bare-skin surface solve (no insulation shell). Kept separate because the insulated
# formulation's log(r_insulation/r_skin) conductance factors are singular at zero insulation.
function solve_without_insulation!(
    body::AbstractBody, geometry_vars::GeometryVariables, environment_vars::NamedTuple, traits::NamedTuple, temperature_tolerance, skin_temperature, insulation_temperature;
    geometry=_part_geometry(body),
    smoothing::SmoothingStrategy=HardBound(),
)
    (;
        temperature,
        view_factors,
        atmos,
        fluid,
        solar_flow,
        gas_fractions,
        convection_enhancement,
    ) = environment_vars
    T = temperature
    F = view_factors
    air_temperature = T.air
    neighbours = _neighbours(environment_vars)
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; core_temperature, flesh_conductivity, ϵ_body, skin_wetness, bare_skin_fraction, eye_fraction) = traits
    tolerance = temperature_tolerance
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0
    volume = flesh_volume(body)
    total_area = geometry.total_area
    area_evaporation = total_area
    area_convection = total_area #* (1 - conduction_fraction)
    r_skin = skin_radius(body)
    insulation_evaporation_heat_flow = 0.0u"W"
    conduction_flow = 0.0u"W"

    while true
        ntry += 1
        for i in 1:20

            # Evaporative heat loss
            conv = convection(;
                body,
                area=area_convection,
                air_temperature,
                surface_temperature=skin_temperature,
                wind_speed,
                atmospheric_pressure,
                fluid,
                gas_fractions,
                convection_enhancement,
                characteristic_dim=geometry.characteristic_dim,
                smoothing,
            )
            heat_transfer_coefficient = conv.heat_transfer_coefficient.combined
            evap_pars_local = AnimalEvaporationParameters(;
                skin_wetness,
                eye_fraction,
                bare_skin_fraction,
            )
            atmos_local = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
            skin_evaporation_flow = evaporation(
                evap_pars_local,
                conv.mass_transfer_coefficient,
                atmos_local,
                area_evaporation,
                skin_temperature,
                air_temperature;
                gas_fractions,
            ).evaporation_heat_flow

            # Radiation coefficients for radiant exchange
            _rc = _radiation_coefficients(area_convection, F, ϵ_body, σ, skin_temperature, T)
            sky_radiation_coeff        = _rc.sky
            bush_radiation_coeff       = _rc.bush
            vegetation_radiation_coeff = _rc.vegetation
            ground_radiation_coeff     = _rc.ground
            # Neighbour exchange as a linear (Σ cᵢ, Σ cᵢTᵢ) pair; zero when free-standing.
            nb = _neighbour_coefficients(area_convection, neighbours, ϵ_body, σ, skin_temperature)
            skin_temperature1 =
                ((4.0 * flesh_conductivity * volume) / (r_skin^2) * core_temperature) - skin_evaporation_flow +
                heat_transfer_coefficient * area_convection * T.air +
                solar_flow
            skin_temperature2 =
                sky_radiation_coeff * T.sky + bush_radiation_coeff * T.bush + vegetation_radiation_coeff * T.vegetation + ground_radiation_coeff * T.ground +
                nb.cTsum
            skin_temperature3 =
                ((4.0 * flesh_conductivity * volume) / (r_skin^2)) +
                heat_transfer_coefficient * area_convection +
                sky_radiation_coeff +
                bush_radiation_coeff +
                vegetation_radiation_coeff +
                ground_radiation_coeff +
                nb.csum

            skin_temperature_calc = (skin_temperature1 + skin_temperature2) / skin_temperature3

            sky_radiation_flow = sky_radiation_coeff * (skin_temperature_calc - T.sky)
            bush_radiation_flow = bush_radiation_coeff * (skin_temperature_calc - T.bush)
            vegetation_radiation_flow = vegetation_radiation_coeff * (skin_temperature_calc - T.vegetation)
            ground_radiation_flow = ground_radiation_coeff * (skin_temperature_calc - T.ground)
            neighbour_radiation_flow = nb.csum * skin_temperature_calc - nb.cTsum

            longwave_flow = sky_radiation_flow + bush_radiation_flow + vegetation_radiation_flow + ground_radiation_flow + neighbour_radiation_flow
            convection_flow = heat_transfer_coefficient * area_convection * (skin_temperature_calc - T.air)

            # Build flows (net_metabolic updated on success)
            flows = HeatFlows(
                convection_flow,
                conduction_flow,
                0.0u"W",  # net_metabolic placeholder
                skin_evaporation_flow,
                insulation_evaporation_heat_flow,
                longwave_flow,
                solar_flow,
                sky_radiation_flow,
                bush_radiation_flow,
                vegetation_radiation_flow,
                ground_radiation_flow,
            )

            Δskin_temperature = abs(skin_temperature - skin_temperature_calc)

            if Δskin_temperature < tolerance
                #TODO why is this not shape-specific?
                net_metabolic = (4 * flesh_conductivity * volume / r_skin^2) * (core_temperature - skin_temperature_calc)
                flows = setproperties(flows; net_metabolic)
                return (;
                    insulation_temperature,
                    skin_temperature=skin_temperature_calc,
                    flows,
                    insulation_conductivity=nothing,
                    tolerance,
                    success=true,
                    ntry,
                )
            else
                skin_temperature = skin_temperature_calc
                insulation_temperature = skin_temperature_calc
                ntry += 1
                if ntry == 101
                    if tolerance <= 0.001u"K"
                        tolerance = 0.01u"K"
                        ntry = 0
                    else
                        return (;
                            insulation_temperature,
                            skin_temperature=skin_temperature_calc,
                            flows,
                            insulation_conductivity=nothing,
                            tolerance,
                            success=false,
                            ntry,
                        )
                    end
                end
            end
        end
    end
end
