# Biophysics

"""
    conduction(; A_conduction, L, surface_temperature, substrate_temperature, k_substrate)
    conduction(A_conduction, L, surface_temperature, substrate_temperature, k_substrate)

Calculate conductive heat transfer between organism and substrate.

# Arguments / Keywords
- `A_conduction`: Contact area with substrate
- `L`: Conduction path length (depth into substrate)
- `surface_temperature`: Organism surface temperature
- `substrate_temperature`: Substrate temperature
- `k_substrate`: Thermal conductivity of substrate

# Returns
- Heat flux from organism to substrate (W). Positive = heat loss.
"""
function conduction(;
    A_conduction=0.0u"m^2",
    L=0.025u"m",
    surface_temperature=u"K"(25.0u"°C"),
    substrate_temperature=u"K"(10.0u"°C"),
    k_substrate=0.1u"W/m/K",
)
    return conduction(A_conduction, L, surface_temperature, substrate_temperature, k_substrate)
end
function conduction(A_conduction, L, surface_temperature, substrate_temperature, k_substrate)
    A_conduction * (k_substrate / L) * (surface_temperature - substrate_temperature)
end

"""
    solar(body, absorptivities, view_factors, solar_conditions; conduction_fraction=0.0)
    solar(body, absorptivities, view_factors, solar_conditions, A_silhouette, A_conduction)

Calculate solar energy balance.

# Arguments
- `body::AbstractBody`: Organism body with geometry
- `absorptivities::Absorptivities`: Body and ground absorptivities
- `view_factors::ViewFactors`: View factors to sky and ground
- `solar_conditions::SolarConditions`: Solar radiation conditions
- `conduction_fraction`: Fraction of body in contact with ground (default 0.0)

# Returns
NamedTuple with solar_flux, direct_flux, solar_sky_flux, solar_substrate_flux
"""
function solar(
    body::AbstractBody,
    absorptivities::Absorptivities,
    view_factors::ViewFactors,
    solar_conditions::SolarConditions;
    conduction_fraction=0.0,
)
    A_total = total_area(body)
    A_silhouette = silhouette_area(body, solar_conditions.zenith_angle)
    A_conduction = A_total * conduction_fraction
    return solar(body, absorptivities, view_factors, solar_conditions, A_silhouette, A_conduction)
end
function solar(
    body::AbstractBody,
    absorptivities::Absorptivities,
    view_factors::ViewFactors,
    solar_conditions::SolarConditions,
    A_silhouette,
    A_conduction,
)
    (; body_dorsal, body_ventral, ground) = absorptivities
    α_d, α_v, α_g = body_dorsal, body_ventral, ground
    F = view_factors
    (; zenith_angle, global_radiation, diffuse_fraction, shade) = solar_conditions

    A_total = total_area(body)

    direct_radiation = global_radiation * (1 - diffuse_fraction)
    diffuse_radiation = global_radiation * diffuse_fraction
    beam_radiation =
        zenith_angle < 90u"°" ? direct_radiation / cos(zenith_angle) : direct_radiation
    direct_flux = α_d * A_silhouette * beam_radiation * (1 - shade)
    solar_sky_flux = α_d * F.sky * A_total * diffuse_radiation * (1 - shade)
    solar_substrate_flux =
        α_v *
        F.ground *
        (A_total - A_conduction) *
        (1 - α_g) *
        global_radiation *
        (1 - shade)
    solar_flux = (direct_flux + solar_substrate_flux + solar_sky_flux)

    return (; solar_flux, direct_flux, solar_sky_flux, solar_substrate_flux)
end

"""
    radin(body, view_factors, emissivities, env_temps; conduction_fraction=0.0)

Calculate incoming longwave radiation from environment.

# Arguments
- `body::AbstractBody`: Organism body with geometry
- `view_factors::ViewFactors`: View factors to sky and ground
- `emissivities::Emissivities`: Body and environment emissivities
- `env_temps::EnvironmentTemperatures`: Environmental temperatures
- `conduction_fraction`: Fraction of body in contact with ground (default 0.0)

# Returns
NamedTuple with infrared_in_flux, infrared_sky_flux, infrared_substrate_flux
"""
function radin(
    body::AbstractBody,
    view_factors::ViewFactors,
    emissivities::Emissivities,
    env_temps::EnvironmentTemperatures;
    conduction_fraction=0.0,
)
    F = view_factors
    (; body_dorsal, body_ventral, ground, sky) = emissivities
    ϵ_d, ϵ_v, ϵ_g, ϵ_s = body_dorsal, body_ventral, ground, sky
    T = env_temps

    A_total = total_area(body)
    A_conduction = A_total * conduction_fraction

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    infrared_sky_flux = ϵ_d * F.sky * A_total * ϵ_s * σ * T.sky^4
    infrared_substrate_flux = ϵ_v * F.ground * (A_total - A_conduction) * ϵ_g * σ * T.ground^4
    infrared_in_flux = infrared_sky_flux + infrared_substrate_flux
    return (; infrared_in_flux, infrared_sky_flux, infrared_substrate_flux)
end

"""
    radout(body, view_factors, emissivities, conduction_fraction, dorsal_temperature, ventral_temperature)

Calculate outgoing longwave radiation from organism.

# Arguments
- `body::AbstractBody`: Organism body with geometry
- `view_factors::ViewFactors`: View factors to sky and ground
- `emissivities::Emissivities`: Body emissivities
- `conduction_fraction`: Fraction of body in contact with ground
- `dorsal_temperature`: Dorsal surface temperature
- `ventral_temperature`: Ventral surface temperature

# Returns
NamedTuple with infrared_out_flux, infrared_to_sky_flux, infrared_to_substrate_flux
"""
function radout(
    body::AbstractBody,
    view_factors::ViewFactors,
    emissivities::Emissivities,
    conduction_fraction,
    dorsal_temperature,
    ventral_temperature,
)
    F = view_factors
    (; body_dorsal, body_ventral) = emissivities
    ϵ_d, ϵ_v = body_dorsal, body_ventral

    A_total = total_area(body)
    A_conduction = A_total * conduction_fraction

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    infrared_to_sky_flux = A_total * F.sky * ϵ_d * σ * dorsal_temperature^4
    infrared_to_substrate_flux = (A_total - A_conduction) * F.ground * ϵ_v * σ * ventral_temperature^4
    infrared_out_flux = infrared_to_sky_flux + infrared_to_substrate_flux
    return (; infrared_out_flux, infrared_to_sky_flux, infrared_to_substrate_flux)
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
- `convection_flux`: Total convective heat loss
- `hc`: Combined heat transfer coefficient
- `hd`: Combined mass transfer coefficient
- `Sh`: Combined Sherwood number
- `free_convection_flux`, `forced_convection_flux`: Free and forced convective components
- `hc_free`, `hc_forc`: Free and forced heat transfer coefficients
- `Sh_free`, `Sh_forc`: Free and forced Sherwood numbers
- `hd_free`, `hd_forc`: Free and forced mass transfer coefficients
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
    β = 1 / air_temperature
    D = body.geometry.characteristic_dimension
    dry_air_out = dry_air_properties(air_temperature, atmospheric_pressure; gas_fractions)
    D_w = dry_air_out.vapour_diffusivity
    # checking to see if the fluid is water, not air
    if fluid == 1
        water_prop_out = water_properties(air_temperature)
        cp_fluid = water_prop_out.specific_heat
        ρ_air = water_prop_out.density
        k_fluid = water_prop_out.thermal_conductivity
        μ = water_prop_out.dynamic_viscosity
    else
        cp_fluid = 1.0057E+3u"J/K/kg"
        ρ_air = dry_air_out.density
        k_fluid = dry_air_out.thermal_conductivity
        μ = dry_air_out.dynamic_viscosity
    end

    # free convection
    Pr = cp_fluid * μ / u"J/s/K/m"(k_fluid)
    if fluid == 1
        #  water; no meaning
        Sc = 1
    else
        Sc = μ / (ρ_air * D_w)
    end
    δ_T = surface_temperature - air_temperature
    if δ_T <= 0.0u"K" # stability check - avoiding zero
        δ_T = δ_T + 0.00001u"K"
    end
    Gr = abs(((ρ_air^2) * β * Unitful.gn * (D^3) * δ_T) / (μ^2))
    Re = ρ_air * wind_speed * D / μ
    Nu_free = nusselt_free(body.shape, Gr, Pr)
    hc_free = (Nu_free * k_fluid) / D # heat transfer coefficient, free
    # calculating the Sherwood number from the Colburn analogy
    # Bird, Stewart & Lightfoot, 1960. Transport Phenomena. Wiley.
    Sh_free = Nu_free * (Sc / Pr)^(1 / 3) # Sherwood number, free
    # calculating the mass transfer coefficient from the Sherwood number
    hd_free = Sh_free * D_w / D # mass transfer coefficient, free
    free_convection_flux = hc_free * area * (surface_temperature - air_temperature) # free convective heat loss at surface
    # forced convection
    Nu = nusselt_forced(body.shape, Re) * convection_enhancement
    # forced convection for object
    hc_forc = Nu * k_fluid / D # heat transfer coefficient, forced
    Sh_forc = Nu * (Sc / Pr)^(1 / 3) # Sherwood number, forced
    hd_forc = Sh_forc * D_w / D # mass transfer coefficient
    forced_convection_flux = hd_forc * area * (surface_temperature - air_temperature) # forced convective heat transfer
    # combined free and forced convection
    # using Bird, Stewart & Lightfoot's mixed convection formula (p. 445, Transport Phenomena, 2002)
    Nu_comb = (Nu_free^3 + Nu^3)^(1 / 3)
    hc = Nu_comb * (k_fluid / D) # mixed convection heat transfer
    convection_flux = hc * area * (surface_temperature - air_temperature) # total convective heat loss
    Sh = Nu_comb * (Sc / Pr)^(1 / 3) # Sherwood number, combined
    hd_forc = Sh_forc * D_w / D  # mass transfer coefficient, forced
    hd = Sh * D_w / D # mass transfer coefficient, combined
    return (;
        convection_flux,
        hc,
        hd,
        Sh,
        free_convection_flux,
        forced_convection_flux,
        hc_free,
        hc_forc,
        Sh_free,
        Sh_forc,
        hd_free,
        hd_forc,
    )
end

"""
    nusselt_free(shape, Gr, Pr)

Calculate the Nusselt number for free convection based on body shape.

# Arguments
- `shape`: Body shape (Cylinder, Plate, Sphere, Ellipsoid, or animal-specific shapes)
- `Gr`: Grashof number (dimensionless)
- `Pr`: Prandtl number (dimensionless)

# Returns
- `Nu_free`: Free convection Nusselt number (dimensionless)
"""
function nusselt_free(shape::Union{Cylinder,DesertIguana,LeopardFrog}, Gr, Pr)
    #  free convection for a cylinder
    #  from p.334 Kreith (1965): Mc Adam's 1954 recommended coordinates
    Ra = Gr * Pr
    if Ra < 1.0E-05
        Nu_free = 0.4
    else
        if Ra < 0.1
            Nu_free = 0.976 * Ra^0.0784
        else
            if Ra <= 100
                Nu_free = 1.1173 * Ra^0.1344
            else
                if Ra < 10000.0
                    Nu_free = 0.7455 * Ra^0.2167
                else
                    if Ra < 1.0E+09
                        Nu_free = 0.5168 * Ra^0.2501
                    else
                        if Ra < 1.0E+12
                            Nu_free = 0.5168 * Ra^0.2501
                        end
                    end
                end
            end
        end
    end
    return Nu_free
end
function nusselt_free(shape::Plate, Gr, Pr)
    Ra = Gr * Pr
    #Nu_free = 0.55 * Ra ^ 0.25
    Nu_free = 0.13 * Ra ^ (1 / 3) # Gates 1980 eq. 9.77
    return Nu_free
end
function nusselt_free(shape::Union{Sphere,Ellipsoid}, Gr, Pr)
    #  sphere free convection
    #  from p.413 Bird et all (1960) Transport Phenomena
    Ra = (Gr ^ (1 / 4)) * (Pr ^ (1 / 3))
    Nu_free = 2.0 + 0.60 * Ra
    if Ra ≥ 200.0
        value = (Gr^(1/4)) * (Pr^(1/3))
        println("Ra = $Ra, (Gr^(1/4)) * (Pr^(1/3)) = $value")
    end
    return Nu_free
end

"""
    nusselt_forced(shape, Re)

Calculate the Nusselt number for forced convection based on body shape.

# Arguments
- `shape`: Body shape (Cylinder, Plate, Sphere, Ellipsoid, or animal-specific shapes)
- `Re`: Reynolds number (dimensionless)

# Returns
- `Nu`: Forced convection Nusselt number (dimensionless)
"""
function nusselt_forced(shape::Cylinder, Re)
    #  forced convection of a cylinder
    #  adjusting Nu - Re correlation for Re number (p. 260 McAdams, 1954)
    if Re < 4
        Nu = 0.891 * Re^0.33
    else
        if Re < 40
            Nu = 0.821 * Re^0.385
        else
            if Re < 4000.0
                Nu = 0.615 * Re^0.466
            else
                if Re < 40000.0
                    Nu = 0.174 * Re^0.618
                else
                    if Re < 400000.0
                        Nu = 0.0239 * Re^0.805
                    end
                end
            end
        end
    end
    return Nu
end
function nusselt_forced(shape::Plate, Re)
    # forced convection of a plate
    #0.102 * Re ^ 0.675 * Pr ^ (1 / 3)
    0.032 * Re ^ 0.8
end
function nusselt_forced(shape::Union{Ellipsoid,Sphere,DesertIguana,LeopardFrog}, Re)
    #  forced convection of a sphere
    0.35 * Re ^ 0.6 # from McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
end

"""
    evaporation(evap_pars, transfer, atmos, area, surface_temperature, air_temperature; kw...)

Compute surface evaporation based on mass transfer coefficient, wetness fractions,
and vapor density gradient between surface and air.

# Arguments
- `evap_pars::EvaporationParameters`: Evaporation parameters (wetness, eye_fraction, bare_skin_fraction)
- `transfer::TransferCoefficients`: Heat and mass transfer coefficients
- `atmos::AtmosphericConditions`: Atmospheric conditions (relative_humidity, atmospheric_pressure)
- `area`: Total surface area for evaporation
- `surface_temperature`: Surface temperature
- `air_temperature`: Air temperature

# Keywords
- `water_potential`: Body water potential (J/kg), default 0.0
- `gas_fractions::GasFractions`: Gas fractions for air properties

# Returns
NamedTuple with evaporation_flux, m_cut, m_eyes
"""
function evaporation(
    evap_pars::EvaporationParameters,
    transfer::TransferCoefficients,
    atmos::AtmosphericConditions,
    area,
    surface_temperature,
    air_temperature;
    water_potential=0.0u"J/kg",
    gas_fractions::GasFractions=GasFractions(),
)
    (; skin_wetness, eye_fraction, bare_skin_fraction) = evap_pars
    (; mass, mass_free) = transfer
    hd, hd_free = mass, mass_free
    (; relative_humidity, atmospheric_pressure) = atmos

    # effective areas for evaporation, partitioned into eye, insulated and bare
    effective_area_eye = area * eye_fraction
    effective_area_insulated = (area - effective_area_eye) * skin_wetness * (1 - bare_skin_fraction)
    effective_area_bare = (area - effective_area_eye) * skin_wetness * bare_skin_fraction

    # get vapour density at surface based on water potential of body
    M_w = (u"kg"(1u"molH₂O")) / 1u"mol" # molar mass of water
    rh_surf = exp(water_potential / (Unitful.R / M_w * surface_temperature))
    wet_air_out = wet_air_properties(surface_temperature, rh_surf, atmospheric_pressure; gas_fractions)
    ρ_vap_surf = wet_air_out.vapour_density

    # get air vapour density
    wet_air_out = wet_air_properties(air_temperature, relative_humidity, atmospheric_pressure; gas_fractions)
    ρ_vap_air = wet_air_out.vapour_density

    # mass of water lost
    m_eyes = hd * effective_area_eye * (ρ_vap_surf - ρ_vap_air) # forced + free
    m_cut_bare = hd * effective_area_bare * (ρ_vap_surf - ρ_vap_air) # forced + free
    m_cut_insulated = hd_free * effective_area_insulated * (ρ_vap_surf - ρ_vap_air) # free
    m_cut = m_cut_bare + m_cut_insulated

    # get latent heat of vapourisation and compute heat exchange due to evaporation
    L_v = enthalpy_of_vaporisation(air_temperature)
    evaporation_flux = Unitful.uconvert(u"W", (m_eyes + m_cut) * L_v)
    # convert from kg/s to g/s
    m_eyes = uconvert(u"g/s", m_eyes)
    m_cut = uconvert(u"g/s", m_cut)
    return (; evaporation_flux, m_cut, m_eyes)
end
