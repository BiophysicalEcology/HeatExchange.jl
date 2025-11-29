# Biophysics

"""
    conduction(; kw...)
    conduction(A_conduction, L, T_surface, T_substrate, k_sub)
"""
function conduction(;
    A_conduction=0.001325006u"m^2",
    L=0.025u"m",
    T_surface=u"K"(25u"°C"),
    T_substrate=u"K"(10u"°C"),
    k_sub=0.1u"W/m/K"
)
    return conduction(A_conduction, L, T_surface, T_substrate, k_sub)
end
conduction(A_conduction, L, T_surface, T_substrate, k_sub) = A_conduction * (k_sub / L) * (T_surface - T_substrate)

"""
    solar(; kw...)
    solar(α_body_dorsal, α_body_ventral, A_silhouette, A_total, A_conduction, F_ground, F_sky, α_ground, Q_solar, Q_direct, Q_diffuse)

Calculate solar energy balance.
"""
function solar(;
    α_body_dorsal,
    α_body_ventral,
    A_silhouette,
    A_total,
    A_conduction,
    F_ground,
    F_sky,
    α_ground,
    shade,
    Q_solar,
    Q_direct,
    Q_diffuse,
)
    return solar(α_body_dorsal, α_body_ventral, A_silhouette, A_total, A_conduction, F_ground, F_sky, α_ground, shade, Q_solar, Q_direct, Q_diffuse)
end
function solar(α_body_dorsal, α_body_ventral, A_silhouette, A_total, A_conduction, F_ground, F_sky, α_ground, shade, Q_solar, Q_direct, Q_diffuse)
    Q_direct = α_body_dorsal * A_silhouette * Q_direct * (1 - shade)
    Q_solar_sky = α_body_dorsal * F_sky * A_total * Q_diffuse * (1 - shade)
    Q_solar_substrate = α_body_ventral * F_ground * (A_total - A_conduction) * (1 - α_ground) * Q_solar * (1 - shade)
    Q_solar = (Q_direct + Q_solar_substrate + Q_solar_sky)
    return (; Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate)
end

"""
    radout(; kw...)
    radout(T_surface, A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral)

Calculate incoming radiation.
"""
function radin(;
    A_total,
    A_conduction,
    F_sky,
    F_ground,
    ϵ_body_dorsal,
    ϵ_body_ventral,
    ϵ_ground,
    ϵ_sky,
    T_sky,
    T_ground,
)
    return radin(A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral, ϵ_ground, ϵ_sky, T_sky, T_ground)
end
function radin(A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral, ϵ_ground, ϵ_sky, T_sky, T_ground)
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    Q_ir_sky = ϵ_body_dorsal * F_sky * A_total * ϵ_sky * σ * T_sky^4
    Q_ir_sub = ϵ_body_ventral * F_ground * (A_total - A_conduction) * ϵ_ground * σ * T_ground^4
    Q_ir_in = Q_ir_sky + Q_ir_sub
    return (; Q_ir_in, Q_ir_sky, Q_ir_sub)
end

"""
    radout(; kw...)
    radout(T_surface, A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral)

Calculate outgoing radiation.
"""
function radout(;
    T_surface,
    A_total,
    A_conduction,
    F_sky,
    F_ground,
    ϵ_body_dorsal,
    ϵ_body_ventral,
)
    radout(T_surface, A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral)
end
function radout(T_surface, A_total, A_conduction, F_sky, F_ground, ϵ_body_dorsal, ϵ_body_ventral)
    # computes longwave radiation lost
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    Q_ir_to_sky = A_total * F_sky * ϵ_body_dorsal * σ * T_surface^4
    Q_ir_to_sub = (A_total - A_conduction) * F_ground * ϵ_body_ventral * σ * T_surface^4
    Q_ir_out = Q_ir_to_sky + Q_ir_to_sub
    return (; Q_ir_out, Q_ir_to_sky, Q_ir_to_sub)
end

function convection(body, area, T_air, T_surface, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)
    G = Unitful.gn # acceleration due to gravity, m.s^2
    β = 1 / T_air
    D = body.geometry.characteristic_dimension
    dry_air_out = dry_air_properties(T_air, P_atmos; fO2, fCO2, fN2)
    D_w = dry_air_out.D_w
    # checking to see if the fluid is water, not air
    if fluid == 1
        water_prop_out = water_properties(T_air)
        cp_fluid = water_prop_out.c_p_H2O
        ρ_air = water_prop_out.ρ_H2O
        k_fluid = water_prop_out.k_H2O
        μ = water_prop_out.μ_H2O
    else
        cp_fluid = 1.0057E+3u"J/K/kg"
        ρ_air = dry_air_out.ρ_air
        k_fluid = dry_air_out.k_air
        μ = dry_air_out.μ
    end

    # free convection
    Pr = cp_fluid * μ / u"J/s/K/m"(k_fluid)
    if fluid == 1
        #  water; no meaning
        Sc = 1
    else
        Sc = μ / (ρ_air * D_w)
    end
    δ_T = max(1.0e-5u"K", T_surface - T_air)
    Gr = abs(((ρ_air^2) * β * G * (D^3) * δ_T) / (μ^2))
    Re = ρ_air * wind_speed * D / μ
    Nu_free = nusselt_free(body.shape, Gr, Pr)
    hc_free = (Nu_free * k_fluid) / D # heat transfer coefficient, free
    # calculating the Sherwood number from the Colburn analogy
    # Bird, Stewart & Lightfoot, 1960. Transport Phenomena. Wiley.
    Sh_free = Nu_free * (Sc / Pr)^(1 / 3) # Sherwood number, free
    # calculating the mass transfer coefficient from the Sherwood number
    hd_free = Sh_free * D_w / D # mass transfer coefficient, free
    Q_free = hc_free * area * (T_surface - T_air) # free convective heat loss at surface
    # forced convection
    Nu = nusselt_forced(body.shape, Re)
    # forced convection for object
    hc_forc = Nu * k_fluid / D # heat transfer coefficient, forced
    Sh_forc = Nu * (Sc / Pr)^(1 / 3) # Sherwood number, forced
    hd_forc = Sh_forc * D_w / D # mass transfer coefficient
    Q_forc = hd_forc * area * (T_surface - T_air) # forced convective heat transfer
    # combined free and forced convection
    # using Bird, Stewart & Lightfoot's mixed convection formula (p. 445, Transport Phenomena, 2002)
    Nu_comb = (Nu_free^3 + Nu^3)^(1 / 3)
    hc = Nu_comb * (k_fluid / D) # mixed convection heat transfer
    Q_conv = hc * area * (T_surface - T_air) # total convective heat loss
    Sh = Nu_comb * (Sc / Pr)^(1 / 3) # Sherwood number, combined
    hd = Sh * D_w / D # mass transfer coefficient, combined

    return (;Q_conv, hc, hd, Sh, Q_free, Q_forc, hc_free, hc_forc, Sh_free, Sh_forc, hd_free, hd_forc)
end

function nusselt_free(shape::Union{Cylinder, DesertIguana, LeopardFrog}, Gr, Pr)
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
   Nu_free = 0.55 * Ra ^ 0.25
   return Nu_free
end
function nusselt_free(shape::Ellipsoid, Gr, Pr)
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
    0.102 * Re ^ 0.675 * Pr ^ (1 / 3)
end
function nusselt_forced(shape::Union{Ellipsoid, Sphere, DesertIguana, LeopardFrog}, Re)
    #  forced convection of a sphere
    0.35 * Re ^ 0.6 # from McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
end

"""
    evaporation(; kw...)
    evaporation(T_core, T_surface, ψ_org, surface_wetness, A_total, hd, eye_fraction, bare_fraction, T_air, rh, elevation, P_atmos)
"""
function evaporation(;
    T_surface,
    ψ_org = 0.0u"J/kg",
    wetness,
    area,
    hd,
    hd_free = 0.0u"m/s", 
    eye_fraction,
    bare_fraction = 1.0, 
    T_air,
    rh,
    P_atmos,
    fO2 = 0.2095,
    fCO2 = 0.000412,
    fN2 = 0.7902,
)
    return evaporation(T_surface, ψ_org, wetness, area, hd, hd_free, eye_fraction, bare_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2)
end
function evaporation(T_surface, ψ_org, wetness, area, hd, hd_free, eye_fraction, bare_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2)
    # this subroutine computes surface evaporation based on the mass transfer
    # coefficient, fraction of surface of the skin acting as a free water surface
    # and exposed to the air, and the vapor density gradient between the
    # surface and the air, each at their own temperature.

	# effective areas for evaporation, partitioned into eye, insulated and bare
    effective_area_eye = area * eye_fraction
    effective_area_insulated = (area - effective_area_eye) * wetness * (1 - bare_fraction)
    effective_area_bare = (area - effective_area_eye) * wetness * bare_fraction

    # get vapour density at surface based on water potential of body
    M_w = (1u"molH₂O" |> u"kg")/1u"mol" # molar mass of water
    rh_surf = exp(ψ_org / (Unitful.R / M_w * T_surface))
    wet_air_out = wet_air_properties(T_surface, rh_surf, P_atmos; fO2, fCO2, fN2)
    ρ_vap_surf = wet_air_out.ρ_vap

    # get air vapour density
    wet_air_out = wet_air_properties(T_air, rh, P_atmos; fO2, fCO2, fN2)
    ρ_vap_air = wet_air_out.ρ_vap

    # mass of water lost
    m_eyes = hd * effective_area_eye * (ρ_vap_surf - ρ_vap_air) # forced + free
    m_cut_bare = hd * effective_area_bare * (ρ_vap_surf - ρ_vap_air) # forced + free
    m_cut_insulated = hd_free * effective_area_insulated * (ρ_vap_surf - ρ_vap_air) # free
    m_cut = m_cut_bare + m_cut_insulated

    # get latent heat of vapourisation and compute heat exchange due to evaporation
    L_v = enthalpy_of_vaporisation(T_air)
    Q_evap = Unitful.uconvert(u"W", (m_eyes + m_cut) * L_v)
    # convert from kg/s to g/s
    m_eyes = uconvert(u"g/s", m_eyes)
    m_cut = uconvert(u"g/s", m_cut)
    return (; Q_evap, m_cut, m_eyes)
end