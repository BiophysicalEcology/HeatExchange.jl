"""
    ellipsoid(; posture, mass, density, core_temperature,
               insulation_depth, k_insulation, oxygen_extraction_efficiency, stress_factor,
               air_temperature, wind_speed, relative_humidity, q10,
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
- `insulation_depth`: Insulation depth (length).
- `k_insulation`: Insulation conductivity (energy/time/length/temperature).
- `oxygen_extraction_efficiency`: Oxygen extraction efficiency (fraction).
- `stress_factor`: Fraction of basal metabolism at which evaporative water loss begins.
- `air_temperature`: Air temperature.
- `wind_speed`: Wind speed (length/time).
- `relative_humidity`: Relative humidity (%).
- `P_atmos`: Atmospheric_pressue (pressure).
- `q10`: Temperature dependence factor for metabolism.
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
    mrate_equation::MetabolicRateEquation = Kleiber(),
    posture,
    mass,
    density,
    core_temperature,
    insulation_depth,
    k_insulation,
    oxygen_extraction_efficiency,
    stress_factor,
    air_temperature,
    wind_speed,
    relative_humidity,
    P_atmos,
    q10,
    minimum_metabolic_rate=missing,
    metabolic_multiplier=1,
    lethal_desiccation,
    f_O2,
)

    # refinition of variables and parameters to match Porter & Kearney 2009 formulae
    T_f = air_temperature |> u"K" # fluid temperature
    T_c = core_temperature |> u"K" # core body temperature
    v = wind_speed
    Q_gen_min = minimum_metabolic_rate

    # avoid divide-by-zero
    posture = posture == 1 ? 1.01 : posture

    # estimate basal metabolism if not provided
    if isnothing(minimum_metabolic_rate) || minimum_metabolic_rate === missing
        allometric_estimate = metabolic_rate(mrate_equation, mass).Q_metab * metabolic_multiplier
        Q_gen_min = allometric_estimate * q10^((ustrip(u"°C", core_temperature) - 37) / 10)
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

    a_o = b * posture + insulation_depth
    b_o = b + insulation_depth
    c_o = c + insulation_depth

    e_o = sqrt(a_o^2 - c_o^2) / a_o
    # outer area of ellipsoid # TODO use general function
    A_o = 2 * π * b_o^2 + 2 * π * ((a_o * b_o) / e_o) * asin(e_o)
    R_ins = (b_o - b) / (k_insulation * A_o) # insulation radius

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
    O2_consumption_rate = u"ml/hr"(Q_gen_final / 20.1u"J/mL")
    ρ_vap_c = wet_air_properties(T_c, 1.0, P_atmos).ρ_vap
    ρ_vap_f = wet_air_properties(T_f, relative_humidity, P_atmos).ρ_vap

    respiratory_water_loss_rate = u"g/hr"((O2_consumption_rate / f_O2 / oxygen_extraction_efficiency) *
                                  (ρ_vap_c - ρ_vap_f))
    latent_heat = enthalpy_of_vaporisation(T_f)
    Q_respiration = u"W"(respiratory_water_loss_rate * latent_heat)

    basal_metabolic_rate_fraction = Q_gen_final / Q_gen_min

    Q_evap = -Q_gen_required + Q_gen_min
    total_water_loss_rate = u"g/hr"(max((u"J/hr"(Q_evap) / latent_heat), 0.0u"kg/hr"))
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

Base.@kwdef struct EndoModelPars{TM,TR,RE,TP,TO} <: AbstractModelParameters
    thermoregulation_mode::TM = Param(1)
    thermoregulate::TR =        Param(true)
    respire::RE =               Param(trie)
    torpor::TP =                Param(false)
    tolerance::TO =             Param(0.001u"K")
end

function endotherm(; model_pars, bodyshape, body_pars, integument_pars, insulation_pars, physio_pars, thermoreg_pars, environmental_pars, organism_vars, environmental_vars)
    shade = 0.0 # TODO how to deal with shade?
    # unpack parameters and variables
    (; thermoregulation_mode, thermoregulate, respire, torpor, tolerance) = model_pars
    (; mass, ρ_flesh, ρ_fat, fat_fraction, shape_b, shape_c) = body_pars
    (; α_body_dorsal, α_body_ventral, ϵ_body_dorsal, ϵ_body_ventral, F_sky,
        F_ground, F_vegetation, F_bush, eye_fraction, skin_wetness, bare_skin_fraction,
        conduction_fraction, ventral_fraction) = integument_pars
    (; insulation_conductivity_dorsal, insulation_conductivity_ventral, 
        fibre_diameter_dorsal,  fibre_diameter_ventral,  fibre_length_dorsal,
        fibre_length_ventral, insulation_depth_dorsal, insulation_depth_ventral,
        fibre_density_dorsal, fibre_density_ventral,  insulation_reflectance_dorsal,   
        insulation_reflectance_ventral, insulation_depth_compressed, fibre_conductivity, 
        longwave_depth_fraction, insulation_wetness) = insulation_pars
    (; Q_minimum, q10, k_flesh, k_fat, fO2_extract, rq, Δ_breath, rh_exit) = physio_pars
    (; insulation_step, shape_b_step, T_core_target, T_core_max, T_core_min, T_core_step,
        k_flesh_step, k_flesh_max, pant, pant_step, pant_multiplier, pant_max,
        skin_wetness_step, skin_wetness_max) = thermoreg_pars
    (; α_ground, ϵ_ground, ϵ_sky, elevation, fluid, fN2, fO2, fCO2, convection_enhancement) = environmental_pars
    (; T_core, T_skin, T_insulation, T_lung, ψ_org) = organism_vars
    (; T_air, T_air_reference, T_sky, T_ground, T_substrate, T_bush, T_vegetation, rh, wind_speed, P_atmos,
        zenith_angle, k_substrate, solar_radiation, direct_radiation, diffuse_radiation)= environmental_vars
    ϵ_body = ϵ_body_dorsal # TODO use both dorsal and ventral
    Q_minimum_ref = Q_minimum
    T_core_ref = T_core_target
    insulation_depth_dorsal_ref = insulation_depth_dorsal
    insulation_depth_ventral_ref = insulation_depth_ventral
    T_air_reference = T_air
    F_ground_reference = F_ground
    F_sky_reference = F_sky
    F_vegetation_reference = F_vegetation
    F_bush_reference = F_bush
    Q_gen = 0.0u"W"
    pant_cost = 0.0u"W"

    while Q_gen < Q_minimum * 0.995
        insulation_temperature = T_insulation * 0.7 + T_skin * 0.3
        insulation_out = insulation_properties(; insulation=insulation_pars, insulation_temperature, ventral_fraction)
        fibre_diameters = insulation_out.fibre_diameters # fur diameter array, mean, dorsal, ventral (m)
        fibre_densities = insulation_out.fibre_densities # fur density array, mean, dorsal, ventral (1/m2)
        insulation_depths = insulation_out.insulation_depths # fur depth array, mean, dorsal, ventral (m)

        fat = Fat(fat_fraction, ρ_fat)
        fur = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
        composite_insulation = (fat, fur)

        # geometric properties
        geometry_pars = Body(bodyshape, CompositeInsulation(fur, fat))

        # solar radiation
        normal_radiation = zenith_angle < 90u"°" ? solar_radiation / cos(zenith_angle) : solar_radiation

        α_body_dorsal = 1 - insulation_reflectance_dorsal #solar absorptivity of dorsal fur (fractional, 0-1)
        α_body_ventral = 1 - insulation_reflectance_ventral # solar absorptivity of ventral fur (fractional, 0-1)

        # correct F_sky for vegetaion overhead
        F_vegetation = F_sky_reference * shade
        F_sky = F_sky_reference - F_vegetation
        F_ground = F_ground_reference

        area_silhouette = silhouette_area(geometry_pars, 0u"°")
        area_total = get_total_area(geometry_pars) # total area, m2
        area_conduction = area_total * conduction_fraction # area of skin for convection/evaporation (total skin area - hair area), m2
        #area_evaporation = get_evaporation_area(geometry_pars)
        #area_convection = area_total * (1 - conduction_fraction)

        Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate = solar_out = solar(α_body_dorsal, α_body_ventral, area_silhouette, area_total, area_conduction, F_ground, F_sky, α_ground, shade, normal_radiation, direct_radiation, diffuse_radiation)
        #Q_dorsal = Q_direct + Q_solar_sky
        Q_ventral = Q_solar_substrate
        #Q_diffuse = Q_solar_sky + Q_solar_substrate

        # # convection
        # (; Q_conv, hc, hd, Sh,
        #     Q_free, Q_forc,
        #     hc_free, hc_forc,
        #     Sh_free, Sh_forc,
        #     hd_free, hd_forc) = convection(geometry_pars, area_convection, T_air, T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)

        # set infrared environment
        T_vegetation = T_air_reference # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (deg C)

        simulsol_out = Vector{NamedTuple}(undef, 2)
        for side in 1:2

            # Calculating solar intensity entering fur. This will depend on whether we are calculating the fur
            # temperature for the dorsal side or the ventral side. The dorsal side will have solar inputs from
            # the direct beam hitting the silhouette area as well as diffuse solar scattered from the sky.
            # The ventral side will have diffuse solar scattered off the substrate.
            # Resetting config factors and solar depending on whether the dorsal side (side=1) or ventral side (side=2) is being estimated.
            if Q_solar > 0.0u"W"
                if side == 1
                    F_sky = F_sky_reference * 2.0 # proportion of upward view that is sky
                    F_vegetation = F_vegetation_reference * 2.0 # proportion of upward view that is vegetation (shade)
                    F_ground = 0.0
                    F_bush = 0.0
                    Q_solar = 2.0 * Q_direct + ((Q_solar_sky / F_sky_reference) * F_sky) # direct x 2 because assuming sun in both directions, and unadjusting Q_solar_sky for config factor imposed in SOLAR_ENDO and back to new larger one in both directions
                else
                    F_sky = 0.0
                    F_vegetation = 0.0
                    F_ground = F_ground_reference * 2.0
                    F_bush = F_bush_reference * 2.0
                    Q_solar = (Q_ventral / (1.0 - F_sky_reference - F_vegetation_reference)) * (1.0 - (2.0 * conduction_fraction)) # unadjust by config factor imposed in SOLAR_ENDO to have it coming in both directions, but also cutting down according to fractional area conducting to ground (in both directions)
                end
            else
                Q_solar = 0.0u"W"
                if side == 1
                    F_sky = F_sky_reference * 2.0
                    F_vegetation = F_vegetation_reference * 2.0
                    F_ground = 0.0
                    F_bush = 0.0
                else
                    F_sky = 0.0
                    F_vegetation = 0.0
                    F_ground = F_ground_reference * 2.0
                    F_bush = F_bush_reference
                end
            end

            # set fur depth and conductivity
            # index for effective_conductivities, is the average (1), front/dorsal (2), back/ventral(3) of the body part
            if Q_solar > 0.0u"W" || (insulation_depths[2] != insulation_depths[3])
                if side == 1
                    insulation = Fur(insulation_depths[2], fibre_diameters[2], fibre_densities[2])
                else
                    insulation = Fur(insulation_depths[3], fibre_diameters[3], fibre_densities[3])
                end
            else
                insulation = Fur(insulation_depths[1], fibre_diameters[1], fibre_densities[1])
            end

            if side == 1
                insulation_conductivity = insulation_conductivity_dorsal
            else
                insulation_conductivity = insulation_conductivity_ventral                
            end
            composite_insulation = (fat, fur)
            geometry_pars = Body(bodyshape, CompositeInsulation(insulation, fat))
            # volume, accounting for any subcutaneous fat
            #volume = geometry_pars.geometry.volume - fat_fraction * mass / ρ_fat

            A_total = get_total_area(geometry_pars)
            r_skin = get_r_skin(geometry_pars) # body radius (including fat), m
            #r_flesh = get_r_flesh(geometry_pars) # body radius flesh only (no fat), m
            r_insulation = r_skin + insulation.thickness # body radius including fur, m
            #diameter = 2 * r_insulation # diameter, m
            #r_radiation = r_skin + (longwave_depth_fraction * insulation.thickness) # effective radiation radius, m
            if geometry_pars.shape isa Cylinder || geometry_pars.shape isa Sphere
                r_compressed = r_skin + insulation_depth_compressed
            else
                r_compressed = r_insulation # Note that this value is never used if conduction not being modeled, but need to have a value for the calculations
            end

            #LEN = ALENTH # length, m

            # Calculating the "cd" variable: Qcond = cd(Tskin-Tsub), where cd = Conduction area*ksub/subdepth
            if side == 2 # doing ventral side, add conduction
                A_conduction = A_total * (conduction_fraction * 2)
                cd = (A_conduction * k_substrate) / 0.025u"m" # assume conduction happens from 2.5 cm depth
            else  # doing dorsal side, no conduction. No need to adjust areas used for convection.
                A_conduction = 0.0u"m^2"
                cd = 0.0u"W/K"
            end

            # package up inputs
            geom_vars = (; side, cd, conduction_fraction, longwave_depth_fraction, convection_enhancement)
            env_vars = (; fluid, T_air, T_ground, T_bush, T_vegetation, T_sky, T_substrate,
                rh, wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2)
            traits = (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness,
                insulation_wetness, bare_skin_fraction, eye_fraction, insulation_conductivity)

            simulsol_out[side,] = simulsol(; T_skin, T_insulation, tolerance, geometry_pars, insulation_pars, insulation_out, geom_vars,
                env_vars, traits)
        end

        T_skin_max = max(simulsol_out[1].T_skin, simulsol_out[2].T_skin)

        # zbrent and respfun

        # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(F_sky+F_vegetation))+(ventral value*F_ground)
        gen_d = simulsol_out[1].Q_gen_net
        gen_v = simulsol_out[2].Q_gen_net
        dmult = F_sky_reference + F_vegetation_reference
        vmult = 1 - dmult # assume that reflectivity of veg below equals reflectivity of soil so vmult left as 1 - dmult
        x = gen_d * dmult + gen_v * vmult # weighted estimate of metabolic heat generation
        Q_sum = x

        # reset configuration factors
        F_bush = F_bush_reference # nearby bush
        F_sky = F_sky_reference # sky
        F_ground = F_ground_reference # ground
        F_vegetation = F_vegetation_reference # vegetation

        # lung temperature and temperature of exhaled air
        T_skin = (simulsol_out[1].T_skin + simulsol_out[2].T_skin) * 0.5
        T_insulation = (simulsol_out[1].T_insulation + simulsol_out[2].T_insulation) * 0.5
        T_lung = (T_core + T_skin) * 0.5 # average of skin and core
        T_air_exit = min(T_air + Δ_breath, T_lung) # temperature of exhaled air, deg C

        if respire
            # now guess for metabolic rate that balances the heat budget while allowing metabolic rate
            # to remain at or above Q_basal, via root-finder ZBRENT
            Q_min = Q_minimum
            Q_m1 = Q_minimum * (-2.)
            Q_m2 = Q_minimum * 10.
            if T_skin_max >= T_core
                Q_m2 = Q_minimum * 1.01
            end
            #tolerance = Q_minimum * tolerance_zbrent

            Q_gen = find_zero(x -> respiration_endotherm(x; T_air_reference=T_air, fO2, fN2, fCO2, P_atmos, Q_min, rq,
                    T_lung, mass, fO2_extract, rh, rh_exit, T_air_exit, pant, Q_sum).balance, (Q_m1, Q_m2), Bisection())

            zbrent_out = respiration_endotherm(Q_gen; T_air_reference=T_air, fO2, fN2, fCO2, P_atmos, Q_min, rq,
                T_lung, mass, fO2_extract, rh, rh_exit, T_air_exit, pant, Q_sum)

            Q_gen = zbrent_out.Q_gen # Q_gen_net
        else
            Q_gen = Q_sum
        end
        # store last-step values
        shape_b_last = shape_b
        k_flesh_last = k_flesh
        T_core_last = T_core
        pant_last = pant
        skin_wetness_last = skin_wetness

        if thermoregulate

            if (insulation_depth_dorsal > insulation_depth_dorsal_ref) &&
               (insulation_depth_ventral > insulation_depth_ventral_ref)
                insulation_depth_dorsal = max(insulation_depth_dorsal_ref, insulation_depth_dorsal -
                                                                           insulation_step * fibre_length_dorsal)
                insulation_depth_ventral = max(insulation_depth_ventral_ref, insulation_depth_ventral -
                                                                             insulation_step * fibre_length_ventral)
            else
                insulation_depth_dorsal = insulation_depth_dorsal_ref
                insulation_depth_ventral = insulation_depth_ventral_ref
                if shape_b < shape_b_max
                    shape_b += shape_b_step
                else
                    shape_b = shape_b_max

                    if k_flesh < k_flesh_max
                        k_flesh += k_flesh_step
                    else
                        k_flesh = k_flesh_max

                        if T_core < T_core_max
                            T_core += T_core_step
                            q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                            if (thermoregulation_mode >= 2) && (pant < pant_max)
                                pant += pant_step
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref

                                if thermoregulation_mode == 3
                                    skin_wetness += skin_wetness_step
                                    if skin_wetness > skin_wetness_max
                                        skin_wetness = skin_wetness_max
                                    end
                                end
                            end

                            Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                        else
                            T_core = T_core_max
                            q10mult = q10^((ustrip(u"K", (T_core - T_core_ref))) / 10)

                            if pant < pant_max
                                pant += pant_step
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref
                                Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                                if thermoregulation_mode == 3
                                    skin_wetness += skin_wetness_step
                                    if skin_wetness > skin_wetness_max
                                        skin_wetness = skin_wetness_max
                                    end
                                end
                            else
                                pant = pant_max
                                pant_cost = ((pant - 1) / (pant_max + 1e-6 - 1)) *
                                            (pant_multiplier - 1) * Q_minimum_ref

                                Q_minimum = (Q_minimum_ref + pant_cost) * q10mult

                                skin_wetness += skin_wetness_step
                                if (skin_wetness > skin_wetness_max) || (skin_wetness_step <= 0)
                                    skin_wetness = skin_wetness_max
                                    return
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return(; T_core_last, Q_gen, shape_b_last, k_flesh_last, pant_last, skin_wetness_last)
end

function simulsol(; T_skin, T_insulation, tolerance, geometry_pars, insulation_pars,
    insulation_out, geom_vars, env_vars, traits)

    insulation_test = u"m"(insulation_out.insulation_test)

    success = true

    if insulation_test > 0.0u"m"
        return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net, Q_evap_skin,
            Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation,
            Q_rad_ground, Q_evap_insulation, success, ntry) = solve_with_insulation!(T_skin, T_insulation,
            geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
            tolerance)
    else
        return (; T_insulation, T_skin, Q_convection, Q_conduction, Q_gen_net, Q_evap_skin,
            Q_solar, Q_rad_sky, Q_rad_bush, Q_rad_vegetation,
            Q_rad_ground, Q_evap_insulation, success, ntry) = solve_without_insulation!(T_skin, T_insulation,
            geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
            tolerance)        
    end
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
    geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
    tolerance
)
    (; side, cd, conduction_fraction,
        longwave_depth_fraction, convection_enhancement) = geom_vars
    (; fluid, T_air, T_bush, T_vegetation, T_ground, T_sky, T_substrate,
        rh, wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2) = env_vars
    (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness, insulation_wetness, bare_skin_fraction,
        eye_fraction, insulation_conductivity) = traits

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0

    area_evaporation = get_evaporation_area(geometry_pars)
    area_total = get_total_area(geometry_pars)
    area_convection = area_total * (1 - conduction_fraction)
    volume = geometry_pars.geometry.volume
    r_skin = get_r_skin(geometry_pars)
    Q_evap_insulation = 0.0u"W"
    Q_conduction = 0.0u"W"

    while true
        ntry += 1
        for i in 1:20

            # Evaporative heat loss
            (; hc, hd, hd_free) = convection(geometry_pars, area_convection, T_air,
                T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)
            Q_evap_skin = evaporation(; T_surface=T_skin, wetness=skin_wetness, area=area_evaporation, hd,
                hd_free, eye_fraction, bare_fraction=bare_skin_fraction, T_air, 
                rh, P_atmos, fO2, fCO2, fN2).Q_evap

            # Q_rad variables for radiant exchange
            Q_rad1 = area_convection * (F_sky * 4.0 * ϵ_body * σ * ((T_skin + T_sky) / 2)^3)
            Q_rad2 = area_convection * (F_bush * 4.0 * ϵ_body * σ * ((T_skin + T_bush) / 2)^3)
            Q_rad3 = area_convection * (F_vegetation * 4.0 * ϵ_body * σ * ((T_skin + T_vegetation) / 2)^3)
            Q_rad4 = area_convection * (F_ground * 4.0 * ϵ_body * σ * ((T_skin + T_ground) / 2)^3)

            T_skin1 = ((4.0 * k_flesh * volume) / (r_skin^2) * T_core) - Q_evap_skin +
                      hc * area_convection * T_air + Q_solar
            T_skin2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
            T_skin3 = ((4.0 * k_flesh * volume) / (r_skin^2)) + hc * area_convection + Q_rad1 +
                      Q_rad2 + Q_rad3 + Q_rad4

            T_skin_calc = (T_skin1 + T_skin2) / T_skin3

            Q_rad_sky = Q_rad1 * (T_skin_calc - T_sky)
            Q_rad_bush = Q_rad2 * (T_skin_calc - T_bush)
            Q_rad_vegetation = Q_rad3 * (T_skin_calc - T_vegetation)
            Q_rad_ground = Q_rad4 * (T_skin_calc - T_ground)

            Q_radiation = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
            Q_convection = hc * area_convection * (T_skin_calc - T_air)
            Q_env = Q_radiation + Q_convection - Q_solar

            ΔT_skin = abs(T_skin - T_skin_calc)

            if ΔT_skin < tolerance
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
                    if tolerance <= 0.001u"K"
                        tolerance = 0.01u"K"
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
    geometry_pars, insulation_pars, insulation_out, geom_vars, env_vars, traits,
    tolerance
)
    (; side, cd, conduction_fraction,
        longwave_depth_fraction, convection_enhancement) = geom_vars
    (; fluid, T_air, T_substrate, T_bush, T_vegetation, T_ground, T_sky, rh,
        wind_speed, P_atmos, F_sky, F_ground, F_bush, F_vegetation, Q_solar, fO2, fCO2, fN2) = env_vars
    (; T_core, k_flesh, k_fat, ϵ_body, skin_wetness, insulation_wetness, bare_skin_fraction,
        eye_fraction, insulation_conductivity) = traits

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    area_evaporation = get_evaporation_area(geometry_pars)
    area_total = get_total_area(geometry_pars)
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
            (; hc, hd, hd_free) = convection(geometry_pars, area_convection, T_air,
                T_insulation, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)
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

            # Effective insulation conductivity
            if !isnothing(insulation_conductivity)
                k_insulation = insulation_conductivity
            else
                T_rad_approx = T_skin * (1 - longwave_depth_fraction) +
                               T_insulation * longwave_depth_fraction
                k_rad = (16 * σ * T_rad_approx^3) / (3 * absorption_coefficient)
                k_insulation = k_eff + k_rad
            end

            (; T_radiant, T_ins_compressed, cd1, cd2, cd3, dv1, dv2, dv3, dv4) =
                radiant_temperature(; body=geometry_pars, insulation=insulation_out, insulation_pars,
                    Q_evap=Q_evap_skin, T_core, T_skin, T_substrate, T_insulation, 
                    k_flesh, k_fat, k_insulation, cd, longwave_depth_fraction, conduction_fraction)

            # Radiative heat fluxes
            Q_rad1 = area_convection * F_sky * 4 * ϵ_body * σ * ((T_radiant + T_sky) / 2)^3
            Q_rad2 = area_convection * F_bush * 4 * ϵ_body * σ * ((T_radiant + T_bush) / 2)^3
            Q_rad3 = area_convection * F_vegetation * 4 * ϵ_body * σ * ((T_radiant + T_vegetation) / 2)^3
            Q_rad4 = area_convection * F_ground * 4 * ϵ_body * σ * ((T_radiant + T_ground) / 2)^3

            if conduction_fraction < 1
                # These calculations are for when there is less than 100% conduction.
                # The term Q_evap_insulation is included for heat lost due to evaporation from
                # the insulation surface
                (; T_insulation_calc, T_radiant2) =
                    insulation_radiant_temperature(; body=geometry_pars, insulation=insulation_out,
                        insulation_pars, T_core, T_ins_compressed, T_air, T_sky, T_ground, T_vegetation, T_bush,
                        T_substrate, area_convection, hc, cd, k_insulation, Q_solar, Q_evap_insulation,
                        Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, dv1, dv2, dv3, dv4,
                        longwave_depth_fraction, conduction_fraction)

                Q_rad_sky = Q_rad1 * (T_radiant2 - T_sky)
                Q_rad_bush = Q_rad2 * (T_radiant2 - T_bush)
                Q_rad_vegetation = Q_rad3 * (T_radiant2 - T_vegetation)
                Q_rad_ground = Q_rad4 * (T_radiant2 - T_ground)
                Q_radiation = Q_rad_sky + Q_rad_bush + Q_rad_vegetation + Q_rad_ground
                Q_convection = hc * area_convection * (T_insulation_calc - T_air)
                Q_conduction = cd * (T_ins_compressed - T_substrate)
            else
                (; cf1, T_ins_compressed) =
                    compressed_radiant_temperature(; body=geometry_pars, insulation=insulation_out,
                        insulation_pars, k_flesh, k_fat, T_core, T_substrate, cd)
                Q_radiation = 0.0u"W"
                Q_convection = 0.0u"W"
                Q_evap_insulation = 0.0u"W"
                Q_solar = 0.0u"W"
                Q_conduction = cd * (T_ins_compressed - T_substrate)
            end
            Q_env = Q_radiation + Q_convection + Q_conduction + Q_evap_insulation - Q_solar
            T_skin_mean, T_skin_calc1 =
                mean_skin_temperature(; body=geometry_pars, insulation=insulation_out,
                    insulation_pars, Q_env, Q_evap_skin, k_flesh, k_fat, T_core, T_insulation_calc,
                    T_ins_compressed, cd1, cd2, cd3, conduction_fraction)

            ΔT_insulation = abs(T_insulation - T_insulation_calc)
            ΔT_skin = abs(T_skin - T_skin_mean)

            # FIRST convergence test (TFADIFF)
            if ΔT_insulation < tolerance
                # Next check TSKIN convergence
                if ΔT_skin < tolerance
                    Q_gen_net = net_metabolic_heat(; body = geometry_pars, T_core, T_skin, k_flesh, k_fat)
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
                    if tolerance <= 0.001u"K"
                        tolerance = 0.01u"K"
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
 T_ground, T_insulation, k_flesh, k_fat, k_insulation, cd, longwave_depth_fraction, 
 conduction_fraction) = radiant_temperature(shape(body), body, insulation, insulation_pars, 
 Q_evap, T_core, T_skin, T_ground, T_insulation, k_flesh, k_fat, k_insulation, cd, 
 longwave_depth_fraction, conduction_fraction)

function radiant_temperature(shape::Union{Cylinder, Plate}, body, insulation, insulation_pars, 
    Q_evap, T_core, T_skin, T_ground, T_insulation, k_flesh, k_fat, k_insulation,
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
        (compression_fraction * T_skin + cd * T_substrate) / (cd + compression_fraction) :
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
     T_ground, T_insulation, k_flesh, k_fat, k_insulation, cd, longwave_depth_fraction, 
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
        (compression_fraction * T_skin + cd * T_substrate) / (cd + compression_fraction) :
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
     T_core, T_skin, T_ground, T_insulation, k_flesh, k_fat, k_insulation, cd, 
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
        (compression_fraction * T_skin + cd * T_substrate) / (cd + compression_fraction) :
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
    T_ins_compressed, T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate, area_convection, hc, cd,
     k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, dv1,
      dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction) = 
      insulation_radiant_temperature(shape(body), body, insulation, insulation_pars, T_core,
      T_ins_compressed, T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate, area_convection, hc, cd,
       k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, dv1,
        dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
        
function insulation_radiant_temperature(shape::Union{Cylinder,Plate}, body, insulation, insulation_pars,
     T_core, T_ins_compressed, T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate, area_convection, hc,
      cd, k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4, cd1, cd2, cd3, 
      dv1, dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
    r_skin = get_r_skin(body)
    r_insulation = get_r_insulation(body)
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation.insulation_depth_compressed
    length = body.geometry.length.length
    if longwave_depth_fraction < 1
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground - 
            (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 = ((2 * π * LEN) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate - 
            Q_evap_insulation + Q_solar
        T_ins4 = (2 * π * length * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * 
            (((k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction)) / dv4) 
            + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) + (T_insulation_calc *
             ((k_insulation / log(r_insulation / r_radiation)) * (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
        T_ins2 = ((2 * π * length) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate -
         Q_evap_insulation + Q_solar
        T_ins4 = (2 * π * length * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_ins_compressed
    end
    return(; T_insulation_calc, T_radiant2)
end

function insulation_radiant_temperature(shape::Sphere, body, insulation, insulation_pars, T_core, T_insulation,
     T_ins_compressed, T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate, area_convection,
      hc, cd, k_insulation, Q_solar, Q_evap_insulation, Q_rad1, Q_rad2, Q_rad3, Q_rad4,
       cd1, cd2, cd3, dv1, dv2, dv3, dv4, longwave_depth_fraction, conduction_fraction)
    r_skin = get_r_skin(body)
    r_insulation = get_r_insulation(body)
    r_radiation = r_skin + insulation_pars.longwave_depth_fraction * insulation.insulation_depth_compressed
    if longwave_depth_fraction < 1
        T_ins1 = ((4 * π * r_skin) / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground -
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate - 
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
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
        T_ins3 = hc * area_convection * T_ins_compressed - cd * T_ins_compressed + cd * T_substrate - 
         Q_evap_insulation + Q_solar
        T_ins4 = (4 * π * r_skin * cd3) / dv1 + (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end
    return(; T_insulation_calc, T_radiant2)    
end

function insulation_radiant_temperature(shape::Ellipsoid, body, insulation, insulation_pars, 
    T_core, T_ins_compressed, T_air, T_sky, T_ground, T_vegetation, T_bush, T_substrate, area_convection, 
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
        T_ins1 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground - 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((dv3 / dv4) + ((T_ins_compressed * cd2) / dv4))
        T_ins2 = ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) * 
         (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + 
         cd * T_substrate - Q_evap_insulation + Q_solar
        T_ins4 = (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) + 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) * ((((k_insulation * bl) / (bl - br)) * 
         (1 - conduction_fraction)) / dv4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = u"K"(dv3 / dv4 + ((T_ins_compressed * cd2) / dv4) + 
         (T_insulation_calc * (((k_insulation * bl) / (bl - br)) * (1 - conduction_fraction))) / dv4)
    else
        T_ins1 = ((3 * volume * bs) / ((((3 * ssqg)^0.5)^3) * dv1)) * 
         (T_core * cd1 - dv2 - T_ins_compressed * cd2)
        T_ins2 = Q_rad1 * T_sky + Q_rad2 * T_bush + Q_rad3 * T_vegetation + Q_rad4 * T_ground
        T_ins3 = hc * area_convection * T_air - cd * T_ins_compressed + cd * T_substrate - 
         Q_evap_insulation + Q_solar
        T_ins4 = (3 * volume * bs * cd3) / ((((3 * ssqg)^0.5)^3) * dv1) + 
         (Q_rad1 + Q_rad2 + Q_rad3 + Q_rad4) + hc * area_convection
        T_insulation_calc = u"K"((T_ins1 + T_ins2 + T_ins3) / T_ins4)
        T_radiant2 = T_insulation_calc
    end
    return (; T_insulation_calc, T_radiant2)
end

compressed_radiant_temperature(; body::AbstractBody, insulation, insulation_pars, k_flesh, k_fat, T_core, T_substrate, cd) = compressed_radiant_temperature(shape(body), body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_substrate, cd)

function compressed_radiant_temperature(shape::Union{Cylinder,Plate}, body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_substrate, cd)
    k_compressed = insulation.insulation_conductivity_compressed
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    r_compressed = r_skin + insulation.insulation_depth_compressed
    cf1 = (2 * π * k_compressed * length) / (log(r_compressed / r_skin))
    dv5 = 1 + ((cf1 * r_flesh^2) / (4 * k_flesh * volume)) + ((cf1 * r_flesh^2) / (2 * k_fat * volume)) * log(r_skin / RFLESH)
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_substrate
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end

function compressed_radiant_temperature(shape::Sphere, body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_substrate, cd)
    k_compressed = insulation.insulation_conductivity_compressed
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    r_compressed = r_skin + insulation.insulation_depth_compres
    cf1 = (4 * π * k_compressed * r_compressed) / (r_compressed - r_skin)
    dv5 = 1 + ((cf1 * r_flesh^2.) / (6 * k_flesh * volume)) + ((cf1 * r_flesh^3) / (3 * k_fat * volume)) * ((r_skin - r_flesh) / (r_skin - r_flesh))
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_substrate
    T_ins_compressed_calc2 = cd + cf1 / dv5
    T_ins_compressed = T_ins_compressed_calc1 / T_ins_compressed_calc2
    return (; cf1, T_ins_compressed)
end

function compressed_radiant_temperature(shape::Ellipsoid, body, insulation, insulation_pars, k_flesh, k_fat, T_core, T_substrate, cd)
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
    T_ins_compressed_calc1 = (cf1 / dv5) * T_core + cd * T_substrate
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

function net_metabolic_heat(shape::Union{Cylinder,Plate}, body, T_core, T_skin, k_flesh, k_fat)
    volume = body.geometry.volume
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)    
    Q_gen_net = (T_core - T_skin) / ((r_flesh ^ 2 / (4 * k_flesh * volume)) + 
        ((r_flesh^2 / (2 * k_fat * volume)) * log(r_skin / r_flesh)))
    return Q_gen_net
end
function net_metabolic_heat(shape::Sphere, body, T_core, T_skin, k_flesh, k_fat)
    volume = body.geometry.volume
    r_skin = get_r_skin(body)
    r_flesh = get_r_flesh(body)
    Q_gen_net = (T_core - T_skin) / ((r_flesh ^ 2 / (6 * k_flesh * volume)) + 
        ((r_flesh^3 / (4 * k_fat * volume)) * ((r_skin - r_flesh)/(r_flesh * r_skin))))
    return Q_gen_net
end
function net_metabolic_heat(shape::Ellipsoid, body, T_core, T_skin, k_flesh, k_fat)
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

function respiration_endotherm(x; T_air_reference, fO2, fN2, fCO2, P_atmos, Q_min, rq,
     T_lung, mass, fO2_extract, rh, rh_exit, T_air_exit, pant, Q_sum)

    T_air = T_air_reference
    gen = x
    # ref_fO2 = 0.2095
    # ref_fN2 = 0.7902
    # fO2 = ref_fO2
    # fN2 = ref_fN2
    # fCO2 = ref_fCO2

    # p_O2 = ref_fO2
    # p_N2 = ref_fN2
    # p_CO2 = ref_fCO2

    # # Allow user override if gas concentrations different from normal, e.g. burrow
    # p_O2 = (p_O2 > fO2 || p_O2 < fO2) ? fO2 : ref_fO2
    # p_N2 = (p_N2 > fN2 || p_N2 < fN2) ? fN2 : ref_fN2
    # p_CO2 = (p_CO2 > fCO2 || p_CO2 < fCO2) ? fCO2 : ref_fCO2

    # # Normalise if the sum is not exactly 1.0
    # tot = p_O2 + p_N2 + p_CO2
    # if tot != 1.0
    #     p_O2 = 1.0 - (p_N2 + p_CO2)
    # end

    # adjust O2 to ensure sum to 1
    if fO2 + fCO2 + fN2 != 1
        fO2 = 1 - (fN2 + fCO2)
    end

    resp_gen = max(gen, Q_min)

    # respgen is assumed to be in J/s
    # o2stp will be in L/s
    # litres of O₂/S @ STP: (data from Kleiber, 1961)
    if rq ≥ 1.0
        # carbohydrate metabolism
        # carbohydrates ≈ 4193 cal/g
        # L/s = (J/s) * (cal/J) * (kcal/cal) * (L O2 / kcal)
        V_O2_STP = resp_gen * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/5.057u"kcal")
    end

    if rq ≤ 0.7
        # fat metabolism; fats ≈ 9400 cal/g
        V_O2_STP = resp_gen * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.7u"kcal")
    else
        # protein metabolism (RQ ≈ 0.8); ≈ 4300 cal/g
        V_O2_STP = resp_gen * (1u"cal"/4.185u"J") * (1u"kcal"/1000u"cal") * (1u"L"/4.5u"kcal")
    end

    P_O2 = P_atmos * fO2
    V_O2 = (V_O2_STP * P_O2 / 273.15u"K") * (T_lung / P_O2)
    #n = PV/RT (ideal gas law: number of moles from press,vol,temp)
    J_O2 = uconvert(u"mol/s", P_atmos * V_O2 / (Unitful.R * T_lung)) # mol O2 consumed
    # moles/s of O2, N2, dry air at entrance [air flow = f(O2 consumption)]
    J_O2_in = J_O2 / fO2_extract # actual oxygen flow in (moles/s), accounting for efficiency of extraction
    J_N2_in = J_O2_in * (fN2 / fO2) #  actual nitrogen flow in (moles/s), accounting for efficiency of extraction
    V_air = V_O2 / fO2 # air flow
    V_CO2 = fCO2 * V_air #O2 flow
    J_CO2_in = P_atmos * V_CO2 / (Unitful.R * T_lung)
    J_air_in = (J_O2_in + J_N2_in + J_CO2_in) * pant
    V_air = uconvert(u"m^3/s", (J_air_in * Unitful.R * 273.15u"K" / 101325u"Pa")) # air volume @ stp (m3/s)
    # computing the vapor pressure at saturation for the subsequent calculation of 
    # actual moles of water based on actual relative humidity
    P_vap_sat = vapour_pressure(T_air)
    J_H2O_in = J_air_in * (P_vap_sat * rh) / (P_atmos - P_vap_sat * rh)
    # moles at exit
    J_O2_out = J_O2_in - J_O2 # remove consumed oxygen from the total
    J_N2_out = J_N2_in
    J_CO2_out = rq * J_O2 + J_CO2_in
    # total moles of air at exit will be approximately the same as at entrance, since 
    # the moles of O2 removed = approx. the # moles of co2 added
    J_air_out = (J_O2_out + J_N2_out + J_CO2_out) * pant
    # calculate vapour pressue of air at exit
    wet_air_out = wet_air_properties(T_air_exit, rh_exit, P_atmos; fO2, fCO2, fN2)
    P_vap_exit = wet_air_out.P_vap
    J_H2O_out = J_air_out * (P_vap_exit / (P_atmos - P_vap_exit))
    # enthalpy = U2-U1, internal energy only, i.e. lat. heat of vap. only involved, since assume 
    # P,T,V constant, so not significant flow energy, PV. (H = U + PV)

    # moles/s lost by breathing:
    J_evap = J_H2O_out - J_H2O_in
    # grams/s lost by breathing = moles lost * gram molecular weight of water:
    m_resp = J_evap * 18u"g/mol"
    
    # putting a cap on water loss for small animals in very cold conditions
    # by assuming they will seek more moderate conditions if they exceed
    # this cap. this will improve stability for solution method.
    # based on data from w.r. welch. 1980. evaporative water loss from
    # endotherms in thermally and hygically complex environments: an
    # empirical approach for interspecific comparisons.
    # j. comp. physiol. 139: 135-143. maximum value recorded for prairie
    # dogs was 0.6 g/(kg-h) = 1.667 x 10**-4 g/(kg-s)
    # highest recorded rate was a resting deer mouse at 8 g/kg-h =
    # 2.22e-03*
    # (edwards & haines 1978. j. comp. physiol. 128: 177-184 in welch 1980)
    # for a 0.01 kg animal, the max. rate would be 1.67 x 10^-6 g/s
    m_resp = min(m_resp, (2.22E-03 * ustrip(u"kg", mass) * 15)u"g/s")


    # get latent heat of vapourisation and compute heat exchange due to respiration
    #L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(T_lung) - 273.15))J / kg # from wet_air_properties
    L_v = enthalpy_of_vaporisation(T_lung)
    # heat loss by breathing (J/s)=(J/kg)*(kg/s)
    Q_resp = uconvert(u"W", L_v * m_resp)

    # note that there is no recovery of heat or moisture assumed in the nose
    Q_net_check = x - Q_resp
    balance = Q_net_check - Q_sum
    return (; balance, Q_resp, m_resp, Q_gen = x, V_air, J_air_in, J_air_out, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out, J_CO2_in, J_CO2_out)
end
