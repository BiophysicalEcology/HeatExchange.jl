# Biophysics

"""
    conduction(; conduction_area, L, surface_temperature, substrate_temperature, substrate_conductivity)
    conduction(conduction_area, L, surface_temperature, substrate_temperature, substrate_conductivity)

Calculate conductive heat transfer between organism and substrate.

# Arguments / Keywords
- `conduction_area`: Contact area with substrate
- `L`: Conduction path length (depth into substrate)
- `surface_temperature`: Organism surface temperature
- `substrate_temperature`: Substrate temperature
- `substrate_conductivity`: Thermal conductivity of substrate

# Returns
- Heat flow from organism to substrate (W). Positive = heat loss.
"""
function conduction(;
    conduction_area=0.0u"m^2",
    L=0.025u"m",
    surface_temperature=u"K"(25.0u"°C"),
    substrate_temperature=u"K"(10.0u"°C"),
    substrate_conductivity=0.1u"W/m/K",
)
    return conduction(conduction_area, L, surface_temperature, substrate_temperature, substrate_conductivity)
end
function conduction(conduction_area, L, surface_temperature, substrate_temperature, substrate_conductivity)
    conduction_area * (substrate_conductivity / L) * (surface_temperature - substrate_temperature)
end

"""
    solar(body, absorptivities, view_factors, solar_conditions; conduction_fraction=0.0)
    solar(body, absorptivities, view_factors, solar_conditions, silhouette_area, conduction_area)

Calculate solar energy balance.

# Arguments
- `body::AbstractBody`: Organism body with geometry
- `absorptivities::Absorptivities`: Body and ground absorptivities
- `view_factors::ViewFactors`: View factors to sky and ground
- `solar_conditions::SolarConditions`: Solar radiation conditions
- `conduction_fraction`: Fraction of body in contact with ground (default 0.0)

# Returns
NamedTuple with solar_flow, direct_flow, solar_sky_flow, solar_substrate_flow
"""
function solar(
    body::AbstractBody,
    absorptivities::Absorptivities,
    view_factors::ViewFactors,
    solar_conditions::SolarConditions;
    conduction_fraction=0.0,
)
    total_area = BiophysicalGeometry.total_area(body)
    silhouette_area = BiophysicalGeometry.silhouette_area(body, solar_conditions.zenith_angle)
    conduction_area = total_area * conduction_fraction
    return solar(body, absorptivities, view_factors, solar_conditions, silhouette_area, conduction_area)
end
function solar(
    body::AbstractBody,
    absorptivities::Absorptivities,
    view_factors::ViewFactors,
    solar_conditions::SolarConditions,
    silhouette_area,
    conduction_area,
)
    body_absorptivity = absorptivities.body
    α_d, α_v, α_g = body_absorptivity.dorsal, body_absorptivity.ventral, absorptivities.ground
    F = view_factors
    (; zenith_angle, global_radiation, diffuse_fraction, shade) = solar_conditions

    total_area = BiophysicalGeometry.total_area(body)

    direct_radiation = global_radiation * (1 - diffuse_fraction)
    diffuse_radiation = global_radiation * diffuse_fraction
    beam_radiation =
        zenith_angle < 90u"°" ? direct_radiation / cos(zenith_angle) : direct_radiation
    direct_flow = α_d * silhouette_area * beam_radiation * (1 - shade)
    solar_sky_flow = α_d * F.sky * total_area * diffuse_radiation * (1 - shade)
    solar_substrate_flow =
        α_v *
        F.ground *
        (total_area - conduction_area) *
        (1 - α_g) *
        global_radiation *
        (1 - shade)
    solar_flow = (direct_flow + solar_substrate_flow + solar_sky_flow)

    return (; solar_flow, direct_flow, solar_sky_flow, solar_substrate_flow)
end

"""
    radiation_in(body, view_factors, emissivities, environmental_temps; conduction_fraction=0.0)

Calculate incoming longwave radiation from environment.

# Arguments
- `body::AbstractBody`: Organism body with geometry
- `view_factors::ViewFactors`: View factors to sky and ground
- `emissivities::Emissivities`: Body and environment emissivities
- `environmental_temps::EnvironmentTemperatures`: Environmental temperatures
- `conduction_fraction`: Fraction of body in contact with ground (default 0.0)

# Returns
NamedTuple with longwave_flow_in, longwave_sky_flow, longwave_substrate_flow
"""
function radiation_in(
    body::AbstractBody,
    view_factors::ViewFactors,
    emissivities::Emissivities,
    environmental_temps::EnvironmentTemperatures;
    conduction_fraction=0.0,
)
    F = view_factors
    body_emissivity = emissivities.body
    ϵ_d, ϵ_v, ϵ_g, ϵ_s = body_emissivity.dorsal, body_emissivity.ventral, emissivities.ground, emissivities.sky
    T = environmental_temps

    total_area = BiophysicalGeometry.total_area(body)
    conduction_area = total_area * conduction_fraction

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    longwave_sky_flow = ϵ_d * F.sky * total_area * ϵ_s * σ * T.sky^4
    longwave_substrate_flow = ϵ_v * F.ground * (total_area - conduction_area) * ϵ_g * σ * T.ground^4
    longwave_flow_in = longwave_sky_flow + longwave_substrate_flow
    return (; longwave_flow_in, longwave_sky_flow, longwave_substrate_flow)
end

"""
    radiation_out(body, view_factors, emissivities, conduction_fraction, dorsal_temperature, ventral_temperature)

Calculate outgoing longwave radiation from organism.

# Arguments
- `body::AbstractBody`: Organism body with geometry
- `view_factors::ViewFactors`: View factors to sky and ground
- `emissivities::Emissivities`: Body emissivities
- `conduction_fraction`: Fraction of body in contact with ground
- `dorsal_temperature`: Dorsal surface temperature
- `ventral_temperature`: Ventral surface temperature

# Returns
NamedTuple with longwave_flow_out, longwave_to_sky_flow, longwave_to_substrate_flow
"""
function radiation_out(
    body::AbstractBody,
    view_factors::ViewFactors,
    emissivities::Emissivities,
    conduction_fraction,
    dorsal_temperature,
    ventral_temperature,
)
    F = view_factors
    ϵ_d, ϵ_v = emissivities.body.dorsal, emissivities.body.ventral

    total_area = BiophysicalGeometry.total_area(body)
    conduction_area = total_area * conduction_fraction

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    longwave_to_sky_flow = total_area * F.sky * ϵ_d * σ * dorsal_temperature^4
    longwave_to_substrate_flow = (total_area - conduction_area) * F.ground * ϵ_v * σ * ventral_temperature^4
    longwave_flow_out = longwave_to_sky_flow + longwave_to_substrate_flow
    return (; longwave_flow_out, longwave_to_sky_flow, longwave_to_substrate_flow)
end

"""
    convection(; body, area, air_temperature, surface_temperature, wind_speed, atmospheric_pressure, fluid, gas_fractions, convection_enhancement)

Calculate combined free and forced convective heat transfer for an organism.

# Keywords
- `body`: Organism body with geometry (must have `shape` and `geometry.characteristic_dimension`)
- `area`: Surface area for convection
- `air_temperature`: Air temperature
- `surface_temperature`: Surface temperature
- `wind_speed`: Wind speed
- `atmospheric_pressure`: Atmospheric pressure
- `fluid`: Fluid type (0 = air, 1 = water)
- `gas_fractions::GasFractions`: Gas fractions for air properties (default: `GasFractions()`)
- `convection_enhancement`: Enhancement factor for forced convection (default: 1.0)

# Returns
NamedTuple with:
- `heat_flow`: Total convective heat loss
- `heat::TransferCoefficients`: Heat transfer coefficients (combined, free, forced)
- `mass::TransferCoefficients`: Mass transfer coefficients (combined, free, forced)
"""
function convection(;
    body,
    area,
    air_temperature,
    surface_temperature,
    wind_speed,
    atmospheric_pressure,
    fluid,
    gas_fractions::GasFractions=GasFractions(),
    convection_enhancement=1.0,
)
    thermal_expansion_coefficient = 1 / air_temperature
    characteristic_dimension = body.geometry.characteristic_dimension
    dry_air_out = dry_air_properties(air_temperature, atmospheric_pressure; gas_fractions)
    vapour_diffusivity = dry_air_out.vapour_diffusivity
    # checking to see if the fluid is water, not air
    if fluid == 1
        water_prop_out = water_properties(air_temperature)
        cp_fluid = water_prop_out.specific_heat
        fluid_density = water_prop_out.density
        fluid_conductivity = water_prop_out.thermal_conductivity
        dynamic_viscosity = water_prop_out.dynamic_viscosity
    else
        cp_fluid = 1.0057E+3u"J/K/kg"
        fluid_density = dry_air_out.density
        fluid_conductivity = dry_air_out.thermal_conductivity
        dynamic_viscosity = dry_air_out.dynamic_viscosity
    end

    # free convection
    prandtl_number = cp_fluid * dynamic_viscosity / u"J/s/K/m"(fluid_conductivity)
    if fluid == 1
        #  water; no meaning
        schmidt_number = 1
    else
        schmidt_number = dynamic_viscosity / (fluid_density * vapour_diffusivity)
    end
    temperature_difference = surface_temperature - air_temperature
    if temperature_difference <= 0.0u"K" # stability check - avoiding zero
        temperature_difference = temperature_difference + 0.00001u"K"
    end
    grashof_number = abs(((fluid_density^2) * thermal_expansion_coefficient * Unitful.gn * (characteristic_dimension^3) * temperature_difference) / (dynamic_viscosity^2))
    reynolds_number = fluid_density * wind_speed * characteristic_dimension / dynamic_viscosity
    free_nusselt_number = nusselt_free(body.shape, grashof_number, prandtl_number)
    free_heat_transfer_coefficient = (free_nusselt_number * fluid_conductivity) / characteristic_dimension # heat transfer coefficient, free
    # calculating the Sherwood number from the Colburn analogy
    # Bird, Stewart & Lightfoot, 1960. Transport Phenomena. Wiley.
    free_sherwood_number = free_nusselt_number * (schmidt_number / prandtl_number)^(1 / 3) # Sherwood number, free
    # calculating the mass transfer coefficient from the Sherwood number
    free_mass_transfer_coefficient = free_sherwood_number * vapour_diffusivity / characteristic_dimension # mass transfer coefficient, free
    free_convection_flow = free_heat_transfer_coefficient * area * (surface_temperature - air_temperature) # free convective heat loss at surface
    # forced convection
    forced_nusselt_number = nusselt_forced(body.shape, reynolds_number) * convection_enhancement
    # forced convection for object
    forced_heat_transfer_coefficient = forced_nusselt_number * fluid_conductivity / characteristic_dimension # heat transfer coefficient, forced
    forced_sherwood_number = forced_nusselt_number * (schmidt_number / prandtl_number)^(1 / 3) # Sherwood number, forced
    forced_mass_transfer_coefficient = forced_sherwood_number * vapour_diffusivity / characteristic_dimension # mass transfer coefficient
    forced_convection_flow = forced_mass_transfer_coefficient * area * (surface_temperature - air_temperature) # forced convective heat transfer
    # combined free and forced convection
    # using Bird, Stewart & Lightfoot's mixed convection formula (p. 445, Transport Phenomena, 2002)
    combined_nusselt_number = (free_nusselt_number^3 + forced_nusselt_number^3)^(1 / 3)
    combined_heat_transfer_coefficient = combined_nusselt_number * (fluid_conductivity / characteristic_dimension) # mixed convection heat transfer
    convection_flow = combined_heat_transfer_coefficient * area * (surface_temperature - air_temperature) # total convective heat loss
    combined_sherwood_number = combined_nusselt_number * (schmidt_number / prandtl_number)^(1 / 3) # Sherwood number, combined
    combined_mass_transfer_coefficient = combined_sherwood_number * vapour_diffusivity / characteristic_dimension # mass transfer coefficient, combined

    heat = TransferCoefficients(combined_heat_transfer_coefficient, free_heat_transfer_coefficient, forced_heat_transfer_coefficient)
    mass = TransferCoefficients(combined_mass_transfer_coefficient, free_mass_transfer_coefficient, forced_mass_transfer_coefficient)
    heat_flow = convection_flow

    return (; heat_flow, heat, mass)
end

"""
    nusselt_free(shape, grashof_number, prandtl_number)

Calculate the Nusselt number for free convection based on body shape.

# Arguments
- `shape`: Body shape (Cylinder, Plate, Sphere, Ellipsoid, or animal-specific shapes)
- `grashof_number`: Grashof number (dimensionless)
- `prandtl_number`: Prandtl number (dimensionless)

# Returns
- Free convection Nusselt number (dimensionless)
"""
function nusselt_free(shape::Union{Cylinder,DesertIguana,LeopardFrog}, grashof_number, prandtl_number)
    #  free convection for a cylinder
    #  from p.334 Kreith (1965): Mc Adam's 1954 recommended coordinates
    rayleigh_number = grashof_number * prandtl_number
    if rayleigh_number < 1.0E-05
        nusselt_number = 0.4
    else
        if rayleigh_number < 0.1
            nusselt_number = 0.976 * rayleigh_number^0.0784
        else
            if rayleigh_number <= 100
                nusselt_number = 1.1173 * rayleigh_number^0.1344
            else
                if rayleigh_number < 10000.0
                    nusselt_number = 0.7455 * rayleigh_number^0.2167
                else
                    if rayleigh_number < 1.0E+09
                        nusselt_number = 0.5168 * rayleigh_number^0.2501
                    else
                        if rayleigh_number < 1.0E+12
                            nusselt_number = 0.5168 * rayleigh_number^0.2501
                        end
                    end
                end
            end
        end
    end
    return nusselt_number
end
function nusselt_free(shape::Plate, grashof_number, prandtl_number)
    rayleigh_number = grashof_number * prandtl_number
    nusselt_number = 0.13 * rayleigh_number ^ (1 / 3) # Gates 1980 eq. 9.77
    return nusselt_number
end
function nusselt_free(shape::Union{Sphere,Ellipsoid}, grashof_number, prandtl_number)
    #  sphere free convection
    #  from p.413 Bird et all (1960) Transport Phenomena
    rayleigh_number = (grashof_number ^ (1 / 4)) * (prandtl_number ^ (1 / 3))
    nusselt_number = 2.0 + 0.60 * rayleigh_number
    if rayleigh_number ≥ 200.0
        value = (grashof_number^(1/4)) * (prandtl_number^(1/3))
        println("rayleigh_number = $rayleigh_number, (grashof_number^(1/4)) * (prandtl_number^(1/3)) = $value")
    end
    return nusselt_number
end

"""
    nusselt_forced(shape, reynolds_number)

Calculate the Nusselt number for forced convection based on body shape.

# Arguments
- `shape`: Body shape (Cylinder, Plate, Sphere, Ellipsoid, or animal-specific shapes)
- `reynolds_number`: Reynolds number (dimensionless)

# Returns
- Forced convection Nusselt number (dimensionless)
"""
function nusselt_forced(shape::Cylinder, reynolds_number)
    #  forced convection of a cylinder
    #  adjusting Nusselt-Reynolds correlation for Reynolds number (p. 260 McAdams, 1954)
    if reynolds_number < 4
        nusselt_number = 0.891 * reynolds_number^0.33
    else
        if reynolds_number < 40
            nusselt_number = 0.821 * reynolds_number^0.385
        else
            if reynolds_number < 4000.0
                nusselt_number = 0.615 * reynolds_number^0.466
            else
                if reynolds_number < 40000.0
                    nusselt_number = 0.174 * reynolds_number^0.618
                else
                    if reynolds_number < 400000.0
                        nusselt_number = 0.0239 * reynolds_number^0.805
                    end
                end
            end
        end
    end
    return nusselt_number
end
function nusselt_forced(shape::Plate, reynolds_number)
    # forced convection of a plate
    0.032 * reynolds_number ^ 0.8
end
function nusselt_forced(shape::Union{Ellipsoid,Sphere,DesertIguana,LeopardFrog}, reynolds_number)
    #  forced convection of a sphere
    0.35 * reynolds_number ^ 0.6 # from McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
end

"""
    evaporation(evap_pars, mass, atmos, area, surface_temperature, air_temperature; kw...)

Compute surface evaporation based on mass transfer coefficient, wetness fractions,
and vapor density gradient between surface and air.

# Arguments
- `evap_pars::EvaporationParameters`: Evaporation parameters (wetness, eye_fraction, bare_skin_fraction)
- `mass::TransferCoefficients`: Mass transfer coefficients (combined, free, forced)
- `atmos::AtmosphericConditions`: Atmospheric conditions (relative_humidity, atmospheric_pressure)
- `area`: Total surface area for evaporation
- `surface_temperature`: Surface temperature
- `air_temperature`: Air temperature

# Keywords
- `water_potential`: Body water potential (J/kg), default 0.0
- `gas_fractions::GasFractions`: Gas fractions for air properties

# Returns
NamedTuple with evaporation_heat_flow, m_cut, m_eyes
"""
function evaporation(
    evap_pars::EvaporationParameters,
    mass::TransferCoefficients,
    atmos::AtmosphericConditions,
    area,
    surface_temperature,
    air_temperature;
    water_potential=0.0u"J/kg",
    gas_fractions::GasFractions=GasFractions(),
)
    (; skin_wetness, eye_fraction, bare_skin_fraction) = evap_pars
    (; relative_humidity, atmospheric_pressure) = atmos

    # effective areas for evaporation, partitioned into eye, insulated and bare
    effective_area_eye = area * eye_fraction
    effective_area_insulated = (area - effective_area_eye) * skin_wetness * (1 - bare_skin_fraction)
    effective_area_bare = (area - effective_area_eye) * skin_wetness * bare_skin_fraction

    # get vapour density at surface based on water potential of body
    molar_mass_water = (u"kg"(1u"molH₂O")) / 1u"mol"
    rh_surf = exp(water_potential / (Unitful.R / molar_mass_water * surface_temperature))
    wet_air_out = wet_air_properties(surface_temperature, rh_surf, atmospheric_pressure; gas_fractions)
    surface_vapour_density = wet_air_out.vapour_density

    # get air vapour density
    wet_air_out = wet_air_properties(air_temperature, relative_humidity, atmospheric_pressure; gas_fractions)
    air_vapour_density = wet_air_out.vapour_density

    # mass of water lost
    m_eyes = mass.combined * effective_area_eye * (surface_vapour_density - air_vapour_density) # forced + free
    m_cut_bare = mass.combined * effective_area_bare * (surface_vapour_density - air_vapour_density) # forced + free
    m_cut_insulated = mass.free * effective_area_insulated * (surface_vapour_density - air_vapour_density) # free
    m_cut = m_cut_bare + m_cut_insulated

    # get latent heat of vapourisation and compute heat exchange due to evaporation
    latent_heat_vaporisation = enthalpy_of_vaporisation(air_temperature)
    evaporation_heat_flow = Unitful.uconvert(u"W", (m_eyes + m_cut) * latent_heat_vaporisation)
    # convert from kg/s to g/s
    m_eyes = uconvert(u"g/s", m_eyes)
    m_cut = uconvert(u"g/s", m_cut)
    return (; evaporation_heat_flow, m_cut, m_eyes)
end
