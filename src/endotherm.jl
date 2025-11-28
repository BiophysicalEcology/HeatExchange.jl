"""
    ellipsoid(; posture, mass, density, core_temperature,
               fur_depth, fur_conductivity, oxygen_extraction_efficiency, stress_factor,
               air_temperature, wind_speed, relative_humidity, Q10,
               minimum_metabolic_rate=missing, P_atmos, metabolic_factor=1)

Ellipsoid endotherm model.

Implements the model from Porter & Kearney (2009) *PNAS* 106:19666–19672,  
with rough water loss estimates. Valid for conditions without solar radiation and  
where air, ground, and sky temperatures are equal.

# Arguments
- `posture`: Ratio of long to short axis of a prolate ellipsoid.
- `mass`: Body mass (mass).
- `density`: Body density (mass/volume).
- `core_temperature`: Core temperature.
- `fur_depth`: Fur depth (length).
- `fur_conductivity`: Fur conductivity (energy/time/length/temperature).
- `oxygen_extraction_efficiency`: Oxygen extraction efficiency (fraction).
- `stress_factor`: Fraction of basal metabolism at which evaporative water loss begins.
- `air_temperature`: Air temperature.
- `wind_speed`: Wind speed (length/time).
- `relative_humidity`: Relative humidity (%).
- `P_atmos`: Atmospheric_pressue (pressure).
- `Q10`: Temperature dependence factor for metabolism.
- `minimum_metabolic_rate`: Optional user-specified minimum metabolic rate (energy/time).
- `f_O2`: Fractional oxygen in the atmosphere (-).
- `metabolic_multiplier`: Scaling factor for the mouse–elephant basal rate.

# Returns
A `NamedTuple` with:
- `skin_temperature`, `upper_critical_air_temperature`, `lower_critical_air_temperature`,
- `Q_gen`, `Q_gen_final`, `Q_respiration`, `Q_evap`,
- `O2_consumption_rate`, `respiratory_water_loss_rate`, `total_water_loss_rate`,
- `basal_metabolic_rate_fraction`, `fractional_mass_loss`.

# References
Porter, W. P., & Kearney, M. R. (2009). *Size, shape, and the thermal niche of endotherms*.  
PNAS, 106(46), 19666–19672.
"""
function ellipsoid_endotherm(;
    posture,
    mass,
    density,
    core_temperature,
    fur_depth,
    fur_conductivity,
    oxygen_extraction_efficiency,
    stress_factor,
    air_temperature,
    wind_speed,
    relative_humidity,
    P_atmos,
    Q10,
    minimum_metabolic_rate=missing,
    metabolic_multiplier=1,
    lethal_desiccation,
    f_O2,
)

    # refinition of variables and parameters to match Porter & Kearney 2009 formulae
    T_f = u"K"(air_temperature) # fluid temperature
    T_c = u"K"(core_temperature) # core body temperature
    k_fur = fur_conductivity
    v = wind_speed
    Q_gen_min = minimum_metabolic_rate

    # avoid divide-by-zero
    posture = posture == 1 ? 1.01 : posture

    # estimate basal metabolism if not provided
    # TODO make a set of allometric models of metabolic rate
    if isnothing(minimum_metabolic_rate) || minimum_metabolic_rate === missing
        mouse_elephant = (10^(-1.462 + 0.675 * log10(ustrip(u"g", mass))) * metabolic_multiplier)u"W"
        Q_gen_min = mouse_elephant * Q10^((ustrip(u"°C", core_temperature) - 37) / 10)
    end

    # constants
    a_coef = 0.6
    b_coef = 0.5
    c_p_air = 1005.8u"J/kg/K"

    # body geometry
    V = mass / density # volume
    b = ((3 * V) / (4 * π * posture))^(1 / 3)
    c = b
    a = b * posture

    k_b = (0.5 + (6.14 * ustrip(u"m", b)) + 0.439)u"W/m/K" # thermal conductivity of body
    numerator = a^2 * b^2 * c^2
    denominator = a^2 * b^2 + a^2 * c^2 + b^2 * c^2
    S2 = numerator / denominator # Eq. 5
    R_b = S2 / (2 * k_b * V) # resistance of body

    a_o = b * posture + fur_depth
    b_o = b + fur_depth
    c_o = c + fur_depth

    e_o = sqrt(a_o^2 - c_o^2) / a_o
    # outer area of ellipsoid # TODO use general function
    A_o = 2 * π * b_o^2 + 2 * π * ((a_o * b_o) / e_o) * asin(e_o)
    R_ins = (b_o - b) / (k_fur * A_o) # insulation radius

    # air properties
    (; μ, k_air, ρ_air) = dry_air_properties(T_f, P_atmos)

    V = (4 / 3) * π * a * b * c # recompute volume
    L_c = V^(1 / 3) # characteristic dimension
    e = sqrt(a^2 - c^2) / a # eccentricity
    A = 2 * π * b^2 + 2 * π * ((a * b) / e) * asin(e) # area of skin

    Re = ρ_air * v * L_c / μ # Reynold's number
    Pr = (μ * c_p_air) / k_air # Prandtl number

    q′′′_numerator = 2 * A * k_b * k_air * (2 + a_coef * (Re^b_coef) * Pr^(1 / 3)) * (T_c - T_f)
    q′′′_denominator = 2 * k_b * L_c * V + A * S2 * k_air * (2 + a_coef * Re^b_coef * Pr^(1 / 3))
    q′′′ = q′′′_numerator / q′′′_denominator
    T_s = T_c - (q′′′ * S2) / (2 * k_b) # skin temperature, Eq. 6

    Gr = abs((ρ_air^2) * (1 / T_f) * Unitful.gn * (L_c^3) * (T_s - T_f) / μ^2) # Grashof number
    Nu_free = 2 + 0.6 * (Gr^0.25) * (Pr^(1 / 3))
    Nu_forced = 0.37 * Re^0.6
    Nu_total = (Nu_free^3 + Nu_forced^3)^(1 / 3) #  Eq. 24
    h_cv = Nu_total * k_air / L_c # Eq. 25
    R_cv = 1 / (h_cv * A_o)
    R_rad = u"K/W"(1 / (4 * A_o * 0.95 * Unitful.σ * T_f^3)) # Eq. 39
    R_total = R_b + R_ins + (R_cv * R_rad) / (R_cv + R_rad)

    upper_critical_air_temperature = T_c - (Q_gen_min * stress_factor * R_total)
    lower_critical_air_temperature = T_c - Q_gen_min * R_total

    Q_gen_required = (T_c - T_f) / R_total
    Q_gen_final = max(Q_gen_required, Q_gen_min)

    # TODO make general function to convert metabolic rate to O2
    O2_consumption_rate = Q_gen_final / 20.1u"J/mL"
    ρ_vap_c = wet_air_properties(T_c, 1.0, P_atmos).ρ_vap
    ρ_vap_f = wet_air_properties(T_f, relative_humidity, P_atmos).ρ_vap

    respiratory_water_loss_rate = (O2_consumption_rate / f_O2 / oxygen_extraction_efficiency) *
                                  (ρ_vap_c - ρ_vap_f)
    latent_heat = enthalpy_of_vaporisation(T_f)
    Q_respiration = u"W"(respiratory_water_loss_rate * latent_heat)

    basal_metabolic_rate_fraction = Q_gen_final / Q_gen_min

    Q_evap = max(Q_respiration, (-Q_gen_final) + Q_gen_min)
    total_water_loss_rate = max((u"J/hr"(Q_evap) / latent_heat), 0.0u"kg/hr")
    fractional_mass_loss = total_water_loss_rate / mass

    return (;
        skin_temperature=u"°C"(T_s),
        upper_critical_air_temperature=u"°C"(upper_critical_air_temperature),
        lower_critical_air_temperature=u"°C"(lower_critical_air_temperature),
        Q_gen_required,
        Q_gen_final,
        Q_respiration,
        Q_evap,
        O2_consumption_rate=u"mL/hr"(O2_consumption_rate),
        respiratory_water_loss_rate=u"g/hr"(respiratory_water_loss_rate),
        total_water_loss_rate=u"g/hr"(total_water_loss_rate),
        basal_metabolic_rate_fraction,
        fractional_mass_loss,
    )
end

function update_T_insulation!(T_insulation, T_insulation_calc, ΔT_insulation, solpro)
    if solpro == 1
        # FIRST SOLUTION PROCEDURE: SET T_insulation GUESS TO THE CALCULATED T_insulation
        T_insulation = T_insulation_calc

    else
        if solpro == 2
            # SECOND SOLUTION PROCEDURE: SET T_insulation TO AVERAGE OF PREVIOUS AND CALCULATED
            T_insulation = (T_insulation_calc + T_insulation) / 2

        else
            # FINAL SOLUTION PROCEDURE: INCREMENTALLY ADJUST T_insulation
            if (T_insulation - T_insulation_calc) < 0.0u"K"
                # T_insulation < T_insulation_calc → increase T_insulation
                if ΔT_insulation > 3.5u"K"
                    T_insulation = T_insulation + 0.5u"K"
                end
                if (ΔT_insulation > 1.0u"K") && (ΔT_insulation < 3.5u"K")
                    T_insulation = T_insulation + 0.05u"K"
                end
                if (ΔT_insulation > 0.1u"K") && (ΔT_insulation < 1.0u"K")
                    T_insulation = T_insulation + 0.05u"K"
                end
                if (ΔT_insulation > 0.01u"K") && (ΔT_insulation < 0.1u"K")
                    T_insulation = T_insulation + 0.005u"K"
                end
                if (ΔT_insulation > 0.0u"K") && (ΔT_insulation < 0.01u"K")
                    T_insulation = T_insulation + 0.0001u"K"
                end
                if (ΔT_insulation > 0.0u"K") && (ΔT_insulation < 0.001u"K")
                    T_insulation = T_insulation + 0.00001u"K"
                end

            else
                # T_insulation > T_insulation_calc → decrease T_insulation
                if ΔT_insulation > 3.5u"K"
                    T_insulation = T_insulation - 0.5u"K"
                end
                if (ΔT_insulation > 1.0u"K") && (ΔT_insulation < 3.5u"K")
                    T_insulation = T_insulation - 0.05u"K"
                end
                if (ΔT_insulation > 0.1u"K") && (ΔT_insulation < 1.0u"K")
                    T_insulation = T_insulation - 0.05u"K"
                end
                if (ΔT_insulation > 0.01u"K") && (ΔT_insulation < 0.1u"K")
                    T_insulation = T_insulation - 0.005u"K"
                end
                if (ΔT_insulation > 0.001u"K") && (ΔT_insulation < 0.01u"K")
                    T_insulation = T_insulation - 0.0001u"K"
                end
                if (ΔT_insulation > 0.0u"K") && (ΔT_insulation < 0.001u"K")
                    T_insulation = T_insulation - 0.00001u"K"
                end
            end
        end
    end
    return T_insulation
end

function solve_without_insulation!(T_skin, T_insulation,
    geometry_out, insulation_pars, insulation_out, geom_vars, env_vars, traits,
    tolerance_simusol
)
    (; side, cd, conduction_fraction,
        longwave_depth_fraction, convection_enhancement) = geom_vars
    (; fluid_type, T_air, T_substrate, T_bush, T_vegetation, T_lower, T_sky, T_conduction,
        rh, wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2) = env_vars
    (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness, insulation_wetness, bare_skin_fraction,
        eye_fraction, insulation_conductivity) = traits

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0

    area_evaporation = get_evaporation_area(geometry_out)
    area_total = get_total_area(geometry_out)
    area_convection = area_total * (1 - conduction_fraction)
    volume = geometry_out.geometry.volume
    r_skin = get_r_skin(geometry_out)
    Q_evap_insulation = 0.0u"W"
    Q_conduction = 0.0u"W"

    while true
        ntry += 1
        for i in 1:20

            # Evaporative heat loss
            (; hc, hd, hd_free) = convection(geometry_out, area_convection, T_air,
                T_insulation, wind_speed, P_atmos, fluid_type, fO2, fCO2, fN2)
            Q_evap_skin = evaporation(; T_surface=T_skin, wetness=skin_wetness, area=area_evaporation, hd,
                hd_free, eye_fraction, bare_fraction=bare_skin_fraction, T_air, 
                rh, P_atmos, fO2, fCO2, fN2).Q_evap

            # Q_rad variables for radiant exchange
            Q_rad1 = area_convection * (F_sky * 4.0 * ϵ_body * σ * ((T_skin + T_sky) / 2)^3)
            Q_rad2 = area_convection * (F_bush * 4.0 * ϵ_body * σ * ((T_skin + T_bush) / 2)^3)
            Q_rad3 = area_convection * (F_vegetation * 4.0 * ϵ_body * σ * ((T_skin + T_vegetation) / 2)^3)
            Q_rad4 = area_convection * (F_ground * 4.0 * ϵ_body * σ * ((T_skin + T_lower) / 2)^3)

            T_skin1 = ((4.0 * k_flesh * volume) / (r_skin^2) * T_core) - Q_evap_skin +
                      hc * area_convection * T_air + Q_solar
            T_skin2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower
            T_skin3 = ((4.0 * k_flesh * volume) / (r_skin^2)) + hc * area_convection + Q_rad1 +
                      Q_rad2 + Q_rad3 + Q_rad4

            T_skin_calc = (T_skin1 + T_skin2) / T_skin3

            Q_rad_sky = Q_rad1 * (T_skin_calc - T_sky)
            Q_rad_bush = Q_rad2 * (T_skin_calc - T_bush)
            Q_rad_vegetation = Q_rad3 * (T_skin_calc - T_vegetation)
            Q_rad_ground = Q_rad4 * (T_skin_calc - T_lower)

            Q_radiation = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
            Q_convection = hc * area_convection * (T_skin_calc - T_air)
            Q_env = Q_radiation + Q_convection - Q_solar

            ΔT_skin = abs(T_skin - T_skin_calc)

            if ΔT_skin < tolerance_simusol
                Q_gen_net = (4 * k_flesh * volume / r_skin^2) * (T_core - T_skin_calc)
                success = true
                return (; T_insulation, T_skin=T_skin_calc, Q_convection, Q_conduction,
                    Q_gen_net, Q_evap_skin, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation,
                    Q_rad_ground, Q_evap_insulation, success, ntry)
            else
                T_skin = T_skin_calc
                T_insulation = T_skin_calc
                ntry += 1

                if ntry == 101
                    if tolerance_simusol <= 0.001u"K"
                        tolerance_simusol = 0.01u"K"
                        ntry = 0
                    else
                        success = false
                        Q_gen_net = 0.0u"W"
                        return (; T_insulation, T_skin_mean, Q_convection, Q_conduction, Q_gen_net,
                            Q_evap_skin, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
                            Q_evap_insulation, success, ntry)
                    end
                end
            end
        end
    end
end

function solve_with_insulation!(T_skin, T_insulation,
    geometry_out, insulation_pars, insulation_out, geom_vars, env_vars, traits,
    tolerance_simusol
)
    (; side, cd, conduction_fraction,
        longwave_depth_fraction, convection_enhancement) = geom_vars
    (; fluid_type, T_air, T_substrate, T_bush, T_vegetation, T_lower, T_sky, T_conduction, rh,
        wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2) = env_vars
    (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness, insulation_wetness, bare_skin_fraction,
        eye_fraction, insulation_conductivity) = traits

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    area_evaporation = get_evaporation_area(geometry_out)
    area_total = get_total_area(geometry_out)
    area_convection = area_total * (1 - conduction_fraction)

    ntry = 0
    solct = 0
    solpro = 1
    success = true
    Q_gen_net = 0.0u"W"

    while ntry < 20
        ntry += 1
        for i in 1:20
            # Evaporative heat loss
            # first from the skin
            (; hc, hd, hd_free) = convection(geometry_out, area_convection, T_air,
                T_insulation, wind_speed, P_atmos, fluid_type, fO2, fCO2, fN2)
            Q_evap_skin = evaporation(; T_surface=T_skin, wetness=skin_wetness, area=area_evaporation, hd,
                hd_free, eye_fraction, bare_fraction=bare_skin_fraction, T_air, rh, P_atmos, 
                fO2, fCO2, fN2).Q_evap
            # second from insulation
            if insulation_wetness > 0
                Q_evap_insulation = evaporation(; T_surface=T_insulation, wetness=insulation_wetness,
                 area=area_convection, hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2).Q_evap
            else
                Q_evap_insulation = 0.0u"W"
            end

            # Radiation properties
            (; effective_conductivities, absorption_coefficients) = insulation_out =
                insulation_properties(;
                    insulation=insulation_pars,
                    insulation_temperature=T_insulation * 0.7 + T_skin * 0.3,
                    ventral_fraction=0.0,
                )
            absorption_coefficient = absorption_coefficients[side+1]
            k_eff = effective_conductivities[side+1]

            # Effective fur conductivity
            if !isnothing(insulation_conductivity)
                k_insulation = insulation_conductivity
            else
                T_rad_approx = T_skin * (1 - longwave_depth_fraction) +
                               T_insulation * longwave_depth_fraction
                k_rad = (16 * σ * T_rad_approx^3) / (3 * absorption_coefficient)
                k_insulation = k_eff + k_rad
            end

            (; T_radiant, T_ins_compressed, cd1, cd2, cd3, dv1, dv2, dv3, dv4) =
                radiant_temperature(; body=geometry_out, insulation=insulation_out, insulation_pars,
                    Q_evap=Q_evap_skin, T_core, T_skin, T_conduction, T_insulation, 
                    k_flesh, k_fat, k_insulation, cd, longwave_depth_fraction, conduction_fraction)

            # Radiative heat fluxes
            Q_rad1 = area_convection * F_sky * 4 * ϵ_body * σ * ((T_radiant + T_sky) / 2)^3
            Q_rad2 = area_convection * F_bush * 4 * ϵ_body * σ * ((T_radiant + T_bush) / 2)^3
            Q_rad3 = area_convection * F_vegetation * 4 * ϵ_body * σ * ((T_radiant + T_vegetation) / 2)^3
            Q_rad4 = area_convection * F_ground * 4 * ϵ_body * σ * ((T_radiant + T_lower) / 2)^3

            if conduction_fraction < 1
                # These calculations are for when there is less than 100% conduction.
                # The term Q_evap_insulation is included for heat lost due to evaporation from
                # the insulation surface
                (; T_insulation_calc, T_radiant2) =
                    insulation_radiant_temperature(; body=geometry_out, insulation=insulation_out,
                        insulation_pars, T_core, T_ins_compressed, T_air, T_sky, T_lower, T_vegetation, T_bush,
                        T_conduction, area_convection, hc, cd, k_insulation, Q_solar, Q_evap_insulation,
                        Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, dv1, dv2, dv3, dv4,
                        longwave_depth_fraction, conduction_fraction)

                Q_rad_sky = Q_rad1 * (T_radiant2 - T_sky)
                Q_rad_bush = Q_rad2 * (T_radiant2 - T_bush)
                Q_rad_vegetation = Q_rad3 * (T_radiant2 - T_vegetation)
                Q_rad_ground = Q_rad4 * (T_radiant2 - T_lower)
                Q_radiation = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
                Q_convection = hc * area_convection * (T_insulation_calc - T_air)
                Q_conduction = cd * (T_ins_compressed - T_conduction)
            else
                (; cf1, T_ins_compressed) =
                    compressed_radiant_temperature(; body=geometry_out, insulation=insulation_out,
                        insulation_pars, k_flesh, k_fat, T_core, T_conduction, cd)
                Q_radiation = 0.0u"W"
                Q_convection = 0.0u"W"
                Q_evap_insulation = 0.0u"W"
                Q_solar = 0.0u"W"
                Q_conduction = cd * (T_ins_compressed - T_conduction)
            end
            Q_env = Q_radiation + Q_convection + Q_conduction + Q_evap_insulation - Q_solar
            T_skin_mean, T_skin_calc1 =
                mean_skin_temperature(; body=geometry_out, insulation=insulation_out,
                    insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc,
                    T_ins_compressed, cd1, cd2, cd3, conduction_fraction)

            ΔT_insulation = abs(T_insulation - T_insulation_calc)
            ΔT_skin = abs(T_skin - T_skin_mean)

            # FIRST convergence test (TFADIFF)
            if ΔT_insulation < tolerance_simusol
                # Next check TSKIN convergence
                if ΔT_skin < tolerance_simusol
                    Q_gen_net = net_metabolic_heat(; body = geometry_out, T_core, T_skin, k_flesh, k_fat)
                    success = true
                    return return (; T_insulation, T_skin=T_skin_mean, Q_convection, Q_conduction,
                        Q_gen_net, Q_evap_skin, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation,
                        Q_rad_ground, Q_evap_insulation, success, ntry)
                else
                    # Not converged, restart iteration
                    if ntry < 20
                        T_skin = T_skin_calc1
                        continue
                    else
                        success_ref[] = 0
                        Q_gen_net = compute_qgennet(
                            ipt, tc, tskcalcav,
                            rflesh, rskin, ak1, ak2, vol,
                            ssqg, bs, bg
                        )
                        return (; T_insulation, T_skin_mean, Q_convection, Q_conduction, Q_gen_net,
                            Q_evap_skin, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
                            Q_evap_insulation, success, ntry)
                    end
                end

            else
                # No ΔT_insulation convergence → update T_insulation
                T_insulation = update_T_insulation!(T_insulation, T_insulation_calc, ΔT_insulation, solpro)
            end
            # update T_skin
            T_skin = T_skin_mean
            solct += 1

            # fallback if stuck (≈ FORTRAN SOLCT, SOLPRO logic)
            if solct ≥ 100
                if solpro != 3
                    solct = 0
                    solpro += 1
                else
                    # Didn't converge → relax tolerance or fail
                    if tolerance_simusol <= 0.001u"K"
                        tolerance_simusol = 0.01u"K"
                        solct = 0
                        solpro = 1
                    else
                        success = false
                        Q_gen_net = 0.0u"W"
                        return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net,
                            Q_evap_skin, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
                            Q_evap_insulation, success, ntry)
                    end
                end
            end
        end
    end
    return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net,
        Q_evap_skin, Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation, Q_rad_ground,
        Q_evap_insulation, success, ntry)
end

radiant_temperature(; body::AbstractBody, insulation, insulation_pars, Q_evap, T_core, T_skin, 
 T_conduction, T_insulation, k_flesh, k_fat, k_insulation, cd, longwave_depth_fraction, 
 conduction_fraction) = radiant_temperature(shape(body), body, insulation, insulation_pars, 
 Q_evap, T_core, T_skin, T_conduction, T_insulation, k_flesh, k_fat, k_insulation, cd, 
 longwave_depth_fraction, conduction_fraction)

function radiant_temperature(shape::Union{Cylinder, Plate}, body, insulation, insulation_pars, 
    Q_evap, T_core, T_skin, T_conduction, T_insulation, k_flesh, k_fat, k_insulation,
     cd, longwave_depth_fraction, conduction_fraction)

    volume = body.geometry.volume
    k_compressed = insulation.insulation_conductivity_compressed
    r_skin = get_r_skin(body)
    r_insulation = get_r_insulation(body)
    r_compressed = r_skin + insulation.insulation_depth_compressed
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation.insulation_depth_compressed
    length = body.geometry.length.length
    compression_fraction =
        (conduction_fraction * 2 * π * k_compressed * length) /
        log(r_compressed / r_skin)

    T_ins_compressed = conduction_fraction > 0 ?
        (compression_fraction * T_skin + cd * T_conduction) / (cd + compression_fraction) :
        0.0u"K"

    cd1 = (k_compressed / log(r_compressed / r_skin)) * conduction_fraction +
          (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)

    cd2 = (k_compressed / log(r_compressed / r_skin)) * conduction_fraction
    cd3 = (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)

    dv1 = 1 +
        ((2 * π * fibre_length * r_flesh^2 * cd1) / (4 * k_flesh * volume)) +
        ((2 * π * fibre_length * r_flesh^2 * cd1) / (2 * k_fat * volume)) *
        log(r_skin / r_flesh)

    dv2 = Q_evap * ((r_flesh^2 * cd1) / (4 * k_flesh * volume)) +
          Q_evap * ((r_flesh^2 * cd1) / (2 * k_fat * volume)) *
          log(r_skin / r_flesh)

    dv3 = ((2 * π * fibre_length) / dv1) *
        (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) *
        r_flesh^2 / (2 * volume)

    dv4 = longwave_depth_fraction < 1 ?
        cd2 + (k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction) :
        1.0

    T_radiant = longwave_depth_fraction < 1 ?
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (T_insulation *
            ((k_insulation / log(r_insulation / r_radiation)) *
            (1 - conduction_fraction))) / dv4 : T_insulation
    return (; T_radiant, T_ins_compressed, cd1, cd2, cd3, dv1, dv2, dv3, dv4)
end

function radiant_temperature(shape::Sphere, body, insulation, insulation_pars, Q_evap, T_core, T_skin,
     T_conduction, T_insulation, k_flesh, k_fat, k_insulation, cd, longwave_depth_fraction, 
     conduction_fraction)

    volume = body.geometry.volume
    k_compressed = insulation.insulation_conductivity_compressed
    r_skin = get_r_skin(body)
    r_insulation = get_r_insulation(body)
    r_compressed = r_skin + insulation.insulation_depth_compressed
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation.insulation_depth_compressed

    compression_fraction =
        (conduction_fraction * 4 * π * k_compressed * r_compressed * r_skin) /
        (r_compressed - r_skin)

    T_ins_compressed = conduction_fraction > 0 ?
        (compression_fraction * T_skin + cd * T_conduction) / (cd + compression_fraction) :
        0.0u"K"

    cd1 = ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction +
          ((k_insulation * r_insulation) / (r_insulation - r_skin)) * (1 - conduction_fraction)

    cd2 = ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction
    cd3 = ((k_insulation * r_insulation) / (r_insulation - r_skin)) * (1 - conduction_fraction)

    dv1 = 1 +
        ((4 * π * r_skin * r_flesh^2 * cd1) / (6 * k_flesh * volume)) +
        ((4 * π * r_skin * r_flesh^3 * cd1) / (3 * k_fat * volume)) *
            ((r_skin - r_flesh) / (r_flesh * r_skin))

    dv2 = Q_evap * ((r_flesh^2 * cd1) / (6 * k_flesh * volume)) +
          Q_evap * ((r_flesh^3 * cd1) / (3 * k_fat * volume)) *
          ((r_skin - r_flesh) / (r_flesh * r_skin))

    dv3 =
        ((4 * π * r_skin) / dv1) *
        (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) *
        r_flesh^3 / (3 * volume * r_radiation)

    dv4 = longwave_depth_fraction < 1 ?
        cd2 + ((k_insulation * r_insulation) / (r_insulation - r_radiation)) * (1 - conduction_fraction) :
        1.0

    T_radiant = longwave_depth_fraction < 1 ?
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (T_insulation *
            ((k_insulation * r_insulation) / (r_insulation - r_radiation) *
            (1 - conduction_fraction))) / dv4 : T_insulation
    return (; T_radiant, T_ins_compressed, cd1, cd2, cd3, dv1, dv2, dv3, dv4)
end

function radiant_temperature(::Ellipsoid, body, insulation, insulation_pars, Q_evap,
     T_core, T_skin, T_conduction, T_insulation, k_flesh, k_fat, k_insulation, cd, 
      longwave_depth_fraction, conduction_fraction)

    volume = body.geometry.volume
    k_compressed = insulation.insulation_conductivity_compressed

    insulation_depth = insulation.insulation_depths[1]
    a_semi_major = body.geometry.length.a_semi_major
    b_semi_minor = body.geometry.length.b_semi_minor
    c_semi_minor = body.geometry.length.c_semi_minor
    bl_compressed = b_semi_minor + insulation_pars.insulation_depth_compressed
    fat =  body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg = (a_square * b_square * c_square) / (a_square * b_square + a_square * c_square + b_square * c_square)

    bg = min(b_semi_minor, b_semi_minor_flesh)
    bs = b_semi_minor
    bl = b_semi_minor + insulation_depth
    br = bs + longwave_depth_fraction * insulation_depth

    compression_fraction =
        (conduction_fraction * 3 * k_compressed * volume * bl_compressed * bs) /
        ((sqrt(3 * ssqg))^3 * (bl_compressed - bs))

    T_ins_compressed = conduction_fraction > 0.0 ?
        (compression_fraction * T_skin + cd * T_conduction) / (cd + compression_fraction) :
        0.0u"K"

    cd1 = ((k_compressed * bl_compressed) / (bl_compressed - bs)) * conduction_fraction +
          ((k_insulation * bl) / (bl - bs)) * (1 - conduction_fraction)

    cd2 = ((k_compressed * bl_compressed) / (bl_compressed - bs)) * conduction_fraction
    cd3 = ((k_insulation * bl) / (bl - bs)) * (1 - conduction_fraction)
    dv1 = 1 + (3 * bs * ssqg * cd1) / (2 * k_flesh * (sqrt(3 * ssqg)^3)) +
        (bs * cd1) / k_fat * ((bs - bg) / (bs * bg))

    dv2 = Q_evap * ((ssqg * cd1) / (2 * k_flesh * volume)) +
        Q_evap * ((sqrt(3 * ssqg)^3 * cd1) / (3 * k_fat * volume)) *
        ((bs - bg) / (bs * bg))
    
    dv3 = (bs / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) / br

    dv4 = longwave_depth_fraction < 1 ?
        cd2 + ((k_insulation * bl) / (bl - br)) * (1 - conduction_fraction) :
        1.0

    T_radiant = longwave_depth_fraction < 1 ?
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (T_insulation *
            ((k_insulation * bl) / (bl - br) *
            (1 - conduction_fraction))) / dv4 : T_insulation
    return (; T_radiant, T_ins_compressed, cd1, cd2, cd3, dv1, dv2, dv3, dv4)
end

insulation_radiant_temperature(; body::AbstractBody, insulation, insulation_pars, T_core, 
    T_ins_compressed, T_air, T_sky, T_lower, T_vegetation, T_bush, T_conduction, area_convection, hc, cd,
     k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, dv1,
      dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction) = 
      insulation_radiant_temperature(shape(body), body, insulation, insulation_pars, T_core,
      T_ins_compressed, T_air, T_sky, T_lower, T_vegetation, T_bush, T_conduction, area_convection, hc, cd,
       k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, dv1,
        dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
        
function insulation_radiant_temperature(shape::Union{Cylinder,Plate}, body, insulation, insulation_pars,
     T_core, T_ins_compressed, T_air, T_sky, T_lower, T_vegetation, T_bush, T_conduction, area_convection, hc,
      cd, k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, 
      dv1, dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
    r_skin = get_r_skin(body)
    r_insulation = get_r_insulation(body)
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation.insulation_depth_compressed
    length = body.geometry.length.length
    if longwave_depth_fraction < 1
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower - 
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 = ((2 * π * LEN) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_conduction - 
            Q_evap_insulation + Q_solar
        T_ins4 = (2 * π * length * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * 
            (((k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction)) / dv4) 
            + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) + (T_insulation_calc *
             ((k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower
        T_ins2 = ((2 * π * length) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_conduction -
         Q_evap_insulation + Q_solar
        T_ins4 = (2 * π * length * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_ins_compressed
    end
    return(; T_insulation_calc, T_radiant2)
end

function insulation_radiant_temperature(shape::Sphere, body, insulation, insulation_pars, T_core, T_insulation,
     T_ins_compressed, T_air, T_sky, T_lower, T_vegetation, T_bush, T_conduction, area_convection,
      hc, cd, k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4,
       cd1, cd2, cd3, dv1, dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
    r_skin = get_r_skin(body)
    r_insulation = get_r_insulation(body)
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation.insulation_depth_compressed
    if longwave_depth_fraction < 1
        T_ins1 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower -
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_conduction - 
         Q_evap_insulation + Q_solar
        T_ins4 = (4 * π * r_skin * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) *
         ((((k_insulation * r_insulation) / (r_insulation - r_radiation)) * 
         (1 - conduction_fraction)) / dv4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) +
         (T_insulation_calc * (((k_insulation * r_insulation) / (r_insulation - r_radiation)) * 
         (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower
        T_ins3 = hc * area_convection * T_ins_compressed - cd * T_ins_compressed + cd * T_conduction - 
         Q_evap_insulation + Q_solar
        T_ins4 = (4 * π * r_skin * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end
    return(; T_insulation_calc, T_radiant2)    
end

function insulation_radiant_temperature(shape::Ellipsoid, body, insulation, insulation_pars, 
    T_core, T_ins_compressed, T_air, T_sky, T_lower, T_vegetation, T_bush, T_conduction, area_convection, 
    hc, cd, k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, 
    dv1, dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
    
    volume = body.geometry.volume

    insulation_depth = insulation.insulation_depths[1]
    a_semi_major = body.geometry.length.a_semi_major
    b_semi_minor = body.geometry.length.b_semi_minor
    c_semi_minor = body.geometry.length.c_semi_minor

    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg = (a_square * b_square * c_square) / (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bl = b_semi_minor + insulation_depth
    br = bs + longwave_depth_fraction * insulation_depth

    if longwave_depth_fraction < 1
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower - 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 = ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) * 
         (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + 
         cd * T_conduction - Q_evap_insulation + Q_solar
        T_ins4 = (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) + 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((((k_insulation * bl) / (bl - br)) * 
         (1 - conduction_fraction)) / dv4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) + 
         (T_insulation_calc * (((k_insulation * bl) / (bl - br)) * (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) * 
         (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_lower
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_conduction - 
         Q_evap_insulation + Q_solar
        T_ins4 = (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) + 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end
    return (; T_insulation_calc, T_radiant2)
end

compressed_radiant_temperature(; body::AbstractBody, insulation, insulation_pars, k_flesh, k_fat, T_core, T_conduction, cd) = compressed_radiant_temperature(shape(body), body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_conduction, cd)

function compressed_radiant_temperature(shape::Union{Cylinder,Plate}, body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_conduction, cd)
    k_compressed = insulation.insulation_conductivity_compressed
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    r_compressed = r_skin + insulation.insulation_depth_compressed
    cf1 = (2 * π * k_compressed * length) / (log(r_compressed / r_skin))
    dv5 = 1 + ((cf1 * r_flesh^2) / (4 * k_flesh * volume)) + ((cf1 * r_flesh^2) / (2 * k_fat * volume)) * log(r_skin / RFLESH)
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_conduction
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end

function compressed_radiant_temperature(shape::Sphere, body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_conduction, cd)
    k_compressed = insulation.insulation_conductivity_compressed
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    r_compressed = r_skin + insulation.insulation_depth_compres
    cf1 = (4 * π * k_compressed * r_compressed) / (r_compressed - r_skin)
    dv5 = 1 + ((cf1 * r_flesh^2.) / (6 * k_flesh * volume)) + ((cf1 * r_flesh^3) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_skin - r_flesh))
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_conduction
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end

function compressed_radiant_temperature(shape::Ellipsoid, body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_conduction, cd)
    volume = body.geometry.volume
    k_compressed = insulation.insulation_conductivity_compressed

    insulation_depth = insulation.insulation_depths[1]
    a_semi_major = body.geometry.length.a_semi_major
    b_semi_minor = body.geometry.length.b_semi_minor
    c_semi_minor = body.geometry.length.c_semi_minor

    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg = (a_square * b_square * c_square) / (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bl = b_semi_minor + insulation_depth
    bl_compressed = b_semi_minor + insulation_pars.insulation_depth_compressed
    bg = min(b_semi_minor, b_semi_minor_flesh)
    
    cf1 = (3 * k_compressed * volume * bl_compressed * bs) / ((((3 * ssqg)^0.5)^3) * (bl - bs))
    dv5 = 1 + ((cf1 * ssqg) / (2 * k_flesh * volume)) + ((cf1 * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) * ((bs - bg) / (bs * bg))
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_conduction
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end

mean_skin_temperature(; body::AbstractBody, insulation, insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc, T_ins_compressed, cd1, cd2, cd3, conduction_fraction) = mean_skin_temperature(shape(body), body, insulation, insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc, T_ins_compressed, cd1, cd2, cd3, conduction_fraction)

function mean_skin_temperature(shape::Union{Cylinder,Plate}, body, insulation, insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc, T_ins_compressed, cd1, cd2, cd3, conduction_fraction)
    volume = body.geometry.volume
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    if conduction_fraction < 1
        T_skin_calc1 = T_core - (((Q_env + Q_evap_skin) * r_flesh^2) / (4 * k_flesh * volume)) - (((Q_env + Q_evap_skin) * r_flesh^2) / (2 * k_fat * volume)) * log(r_skin / r_flesh)
        T_skin_calc2 = ((Q_env * r_flesh^2) / (2 * cd1 * volume)) + ((T_ins_compressed * cd2) / cd1) + ((T_insulation_calc * cd3) / cd1)
    else
        T_skin_calc1 = T_core - ((Q_env * r_flesh^2) / (4 * k_flesh * volume)) - ((Q_env * r_flesh^2.) / (2 * k_fat * volume)) * log(r_skin / r_flesh)
        T_skin_calc2 = (((Q_env * r_flesh^2) / (2 * k_compressed * volume)) * log(r_compressed / r_skin)) + T_ins_compressed
    end
    T_skin_mean = (T_skin_calc1 + T_skin_calc2) / 2
    return (; T_skin_mean, T_skin_calc1)
end

function mean_skin_temperature(shape::Sphere, body, insulation, insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc, T_ins_compressed, cd1, cd2, cd3, conduction_fraction)
    volume = body.geometry.volume
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    if conduction_fraction < 1
        T_skin_calc1 = T_core - (((Q_env + Q_evap_skin) * r_flesh^2) / (6 * k_flesh * volume)) - (((Q_env + Q_evap_skin) * r_flesh^3.) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_skin * r_flesh))
        T_skin_calc2 = ((Q_env * r_flesh^3) / (3 * cd1 * volume * r_skin)) + ((T_ins_compressed * cd2) / cd1) + ((T_insulation_calc * cd3) / cd1)
    else
        T_skin_calc1 = T_core - ((Q_env * r_flesh^2) / (6 * k_flesh * volume)) - ((Q_env * r_flesh^3) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_skin * r_flesh))
        T_skin_calc2 = (((Q_env * r_flesh^3) / (3 * k_compressed * volume)) * ((r_compressed - r_skin) / (r_compressed * r_skin))) + T_ins_compressed
    end
    T_skin_mean = (T_skin_calc1 + T_skin_calc2) / 2
    return (; T_skin_mean, T_skin_calc1)
end

function mean_skin_temperature(shape::Ellipsoid, body, insulation, insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc, T_ins_compressed, cd1, cd2, cd3, conduction_fraction)
    volume = body.geometry.volume
    k_compressed = insulation.insulation_conductivity_compressed
    a_semi_major = body.geometry.length.a_semi_major
    b_semi_minor = body.geometry.length.b_semi_minor
    c_semi_minor = body.geometry.length.c_semi_minor
    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg = (a_square * b_square * c_square) / (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bl_compressed = b_semi_minor + insulation_pars.insulation_depth_compressed
    bg = min(b_semi_minor, b_semi_minor_flesh)

    if conduction_fraction < 1
        T_skin_calc1 = T_core - (((Q_env + Q_evap_skin) * ssqg) / (2 * k_flesh * volume)) - (((Q_env + Q_evap_skin) * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) * ((bs - bg) / (bs * bg))
        T_skin_calc2 = ((Q_env * (((3 * ssqg)^0.5)^3)) / (3 * cd1 * volume * bs)) + ((T_ins_compressed * cd2) / cd1) + ((T_insulation_calc * cd3) / cd1)
    else
        T_skin_calc1 = T_core - ((Q_env * ssqg) / (2 * k_flesh * volume)) - ((Q_env * (((3 * ssqg)^0.5)^3)) / (3 * k_fat * volume)) * ((bs - bg) / (bs * bg))
        T_skin_calc2 = (((Q_env * (((3 * ssqg)^0.5)^3)) / (3 * k_compressed * volume)) * ((bl_compressed - bs) / (bl_compressed * bs))) + T_ins_compressed
    end
        T_skin_mean = (T_skin_calc1 + T_skin_calc2) / 2
    return (; T_skin_mean, T_skin_calc1)
end


net_metabolic_heat(; body::AbstractBody, T_core, T_skin, k_flesh, k_fat) = 
    net_metabolic_heat(shape(body), body, T_core, T_skin, k_flesh, k_fat)

net_metabolic_heat = function(shape::Union{Cylinder,Plate}, body, T_core, T_skin, k_flesh, k_fat)
    volume = body.geometry.volume
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)    
    Q_gen_net = (T_core - T_skin) / ((r_flesh ^ 2 / (4 * k_flesh * volume)) + 
        ((r_flesh^2 / (2 * k_fat * volume)) * log(r_skin / r_flesh)))
    return Q_gen_net
end
net_metabolic_heat = function(shape::Sphere, body, T_core, T_skin, k_flesh, k_fat)
    volume = body.geometry.volume
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    Q_gen_net = (T_core - T_skin) / ((r_flesh ^ 2 / (6 * k_flesh * volume)) + 
        ((r_flesh^3 / (4 * k_fat * volume)) * ((r_skin - r_flesh)/(r_flesh * r_skin))))
    return Q_gen_net
end
net_metabolic_heat = function(shape::Ellipsoid, body, T_core, T_skin, k_flesh, k_fat)
    volume = body.geometry.volume
    a_semi_major = body.geometry.length.a_semi_major
    b_semi_minor = body.geometry.length.b_semi_minor
    c_semi_minor = body.geometry.length.c_semi_minor
    fat = body.geometry.length.fat
    a_semi_major_flesh = a_semi_major - fat
    b_semi_minor_flesh = b_semi_minor - fat
    c_semi_minor_flesh = c_semi_minor - fat

    a_square = min(a_semi_major_flesh^2, a_semi_major^2)
    b_square = min(b_semi_minor_flesh^2, b_semi_minor^2)
    c_square = min(c_semi_minor_flesh^2, c_semi_minor^2)

    ssqg = (a_square * b_square * c_square) / (a_square * b_square + a_square * c_square + b_square * c_square)

    bs = b_semi_minor
    bg = min(b_semi_minor, b_semi_minor_flesh)

    Q_gen_net = (T_core - T_skin) / ((ssqg / (2 * k_flesh * volume)) + 
        (((((3 * ssqg)^0.5)^3) / (3 * k_fat * volume)) * ((bs - bg) / (bg * bs))))
    return Q_gen_net
end
