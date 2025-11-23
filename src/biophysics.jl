# Biophysics

function get_nusselt_free(shape::Union{Cylinder, DesertIguana, LeopardFrog}, Gr, Pr)
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
function get_nusselt_free(shape::Plate, Gr, Pr)
   Ra = Gr * Pr
   Nu_free = 0.55 * Ra ^ 0.25
   return Nu_free
end
function get_nusselt_free(shape::Ellipsoid, Gr, Pr)
    #  sphere free convection
    #  from p.413 Bird et all (1960) Transport Phenomena
    Ra = (Gr ^ (1 /4)) * (Pr ^ (1 / 3))
    Nu_free = 2 + 0.60 * Ra
    # if Ra >= 200
    # print(Ra, '(Gr ^ 0.25) * (Pr ^ 0.333) is too large for correlation'))
    # end
    return Nu_free
end

function get_nusselt_forced(shape::Cylinder, Re)
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
function get_nusselt_forced(shape::Plate, Re)
    # forced convection of a plate
    0.102 * Re ^ 0.675 * Pr ^ (1 / 3)
end
function get_nusselt_forced(shape::Union{Ellipsoid, Sphere, DesertIguana, LeopardFrog}, Re)
    #  forced convection of a sphere
    0.35 * Re ^ 0.6 # from McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
end

function convection(body, A_convection, T_air, T_surface, wind_speed, P_atmos, fluid, fO2, fCO2, fN2)
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
    δ_T = T_surface - T_air
    Gr = abs(((ρ_air^2) * β * G * (D^3) * δ_T) / (μ^2))
    Re = ρ_air * wind_speed * D / μ
    Nu_free = get_nusselt_free(body.shape, Gr, Pr)
    hc_free = (Nu_free * k_fluid) / D # heat transfer coefficient, free
    # calculating the Sherwood number from the Colburn analogy
    # Bird, Stewart & Lightfoot, 1960. Transport Phenomena. Wiley.
    Sh_free = Nu_free * (Sc / Pr)^(1 / 3) # Sherwood number, free
    # calculating the mass transfer coefficient from the Sherwood number
    hd_free = Sh_free * D_w / D # mass transfer coefficient, free
    Q_free = hc_free * A_convection * (T_surface - T_air) # free convective heat loss at surface

    # forced convection
    Nu = get_nusselt_forced(body.shape, Re)
    # forced convection for object
    hc_forc = Nu * k_fluid / D # heat transfer coefficient, forced
    Sh_forc = Nu * (Sc / Pr)^(1 / 3) # Sherwood number, forced
    hd_forc = Sh_forc * D_w / D # mass transfer coefficient
    Q_forc = hd_forc * A_convection * (T_surface - T_air) # forced convective heat transfer

    # combined free and forced convection
    # using Bird, Stewart & Lightfoot's mixed convection formula (p. 445, Transport Phenomena, 2002)
    Nu_comb = (Nu_free^3 + Nu^3)^(1 / 3)
    hc = Nu_comb * (k_fluid / D) # mixed convection heat transfer
    Q_conv = hc * A_convection * (T_surface - T_air) # total convective heat loss
    Sh = Nu_comb * (Sc / Pr)^(1 / 3) # Sherwood number, combined
    hd = Sh * D_w / D # mass transfer coefficient, combined

    return (;Q_conv, hc, hd, Sh, Q_free, Q_forc, hc_free, hc_forc, Sh_free, Sh_forc, hd_free, hd_forc)
end


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
    solar(α_body_dorsal, α_body_ventral, A_silhouette, A_total, A_conduction, F_substrate, F_sky, α_substrate, Q_solar, Q_direct, Q_diffuse)

Calculate solar energy balance.
"""
function solar(;
    α_body_dorsal,
    α_body_ventral,
    A_silhouette,
    A_total,
    A_conduction,
    F_substrate,
    F_sky,
    α_substrate,
    shade,
    Q_solar,
    Q_direct,
    Q_diffuse,
)
    return solar(α_body_dorsal, α_body_ventral, A_silhouette, A_total, A_conduction, F_substrate, F_sky, α_substrate, shade, Q_solar, Q_direct, Q_diffuse)
end
function solar(α_body_dorsal, α_body_ventral, A_silhouette, A_total, A_conduction, F_substrate, F_sky, α_substrate, shade, Q_solar, Q_direct, Q_diffuse)
    Q_direct = α_body_dorsal * A_silhouette * Q_direct * (1 - shade)
    Q_solar_sky = α_body_dorsal * F_sky * A_total * Q_diffuse * (1 - shade)
    Q_solar_substrate = α_body_ventral * F_substrate * (A_total - A_conduction) * (1 - α_substrate) * Q_solar * (1 - shade)
    Q_solar = (Q_direct + Q_solar_substrate + Q_solar_sky)
    return (; Q_solar, Q_direct, Q_solar_sky, Q_solar_substrate)
end

"""
    radout(; kw...)
    radout(T_surface, A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral)

Calculate incoming radiation.
"""
function radin(;
    A_total,
    A_conduction,
    F_sky,
    F_substrate,
    ϵ_body_dorsal,
    ϵ_body_ventral,
    ϵ_substrate,
    ϵ_sky,
    T_sky,
    T_substrate,
)
    return radin(A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral, ϵ_substrate, ϵ_sky, T_sky, T_substrate)
end
function radin(A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral, ϵ_substrate, ϵ_sky, T_sky, T_substrate)
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    Q_ir_sky = ϵ_body_dorsal * F_sky * A_total * ϵ_sky * σ * T_sky^4
    Q_ir_sub = ϵ_body_ventral * F_substrate * (A_total - A_conduction) * ϵ_substrate * σ * T_substrate^4
    Q_ir_in = Q_ir_sky + Q_ir_sub
    return (; Q_ir_in, Q_ir_sky, Q_ir_sub)
end

"""
    radout(; kw...)
    radout(T_surface, A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral)

Calculate outgoing radiation.
"""
function radout(;
    T_surface,
    A_total,
    A_conduction,
    F_sky,
    F_substrate,
    ϵ_body_dorsal,
    ϵ_body_ventral,
)
    radout(T_surface, A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral)
end
function radout(T_surface, A_total, A_conduction, F_sky, F_substrate, ϵ_body_dorsal, ϵ_body_ventral)
    # computes longwave radiation lost
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    Q_ir_to_sky = A_total * F_sky * ϵ_body_dorsal * σ * T_surface^4
    Q_ir_to_sub = (A_total - A_conduction) * F_substrate * ϵ_body_ventral * σ * T_surface^4
    Q_ir_out = Q_ir_to_sky + Q_ir_to_sub
    return (; Q_ir_out, Q_ir_to_sky, Q_ir_to_sub)
end

"""
    evaporation(; kw...)
    evaporation(T_core, T_surface, m_resp, ψ_org, skin_wetness, A_total, hd, p_eyes, T_air, rh, elevation, P_atmos)
"""
function evaporation(;
    T_surface,
    m_resp,
    ψ_org,
    skin_wetness,
    A_total,
    hd,
    p_eyes,
    T_air,
    rh,
    P_atmos,
    fO2,
    fCO2,
    fN2,
)
    return evaporation(T_surface, m_resp, ψ_org, skin_wetness, A_total, hd, p_eyes, T_air, rh, P_atmos, fO2, fCO2, fN2)
end
function evaporation(T_surface, m_resp, ψ_org, skin_wetness, A_total, hd, p_eyes, T_air, rh, P_atmos, fO2, fCO2, fN2)
    # this subroutine computes surface evaporation based on the mass transfer
    # coefficient, fraction of surface of the skin acting as a free water surface
    # and exposed to the air, and the vapor density gradient between the
    # surface and the air, each at their own temperature.

    # get vapour density at surface based on water potential of body
    M_w = (1molH₂O |> u"kg")/1u"mol" # molar mass of water
    rh_surf = exp(ψ_org / (Unitful.R / M_w * T_surface))
    wet_air_out = wet_air_properties(T_surface, rh_surf, P_atmos; fO2, fCO2, fN2)
    ρ_vap_surf = wet_air_out.ρ_vap

    # get air vapour density
    wet_air_out = wet_air_properties(T_air, rh, P_atmos; fO2, fCO2, fN2)
    ρ_vap_air = wet_air_out.ρ_vap

    # water lost from eyes if present
    m_eyes = hd * p_eyes * A_total * (ρ_vap_surf - ρ_vap_air)
    if m_eyes > 0.0u"kg/s"
        m_cut = (A_total * skin_wetness - A_total * skin_wetness * p_eyes) * hd * (ρ_vap_surf - ρ_vap_air)
    else
        m_cut = A_total * skin_wetness * hd * (ρ_vap_surf - ρ_vap_air)
    end

    # total water lost
    m_evap = m_eyes + m_resp + m_cut

    # get latent heat of vapourisation and compute heat exchange due to evaporation
    L_v = enthalpy_of_vaporisation(T_air)
    Q_evap = Unitful.uconvert(u"W", (m_eyes + m_cut) * L_v)

    #onvert from kg/s to g/s
    m_eyes = uconvert(u"g/s", m_eyes)
    m_resp = uconvert(u"g/s", m_resp)
    m_cut = uconvert(u"g/s", m_cut)
    m_evap = uconvert(u"g/s", m_evap)

    return (; Q_evap, m_evap, m_resp, m_cut, m_eyes)
end

"""
    respiration(; kw...)
    respiration(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, elevation, P_atmos, fO2, fCO2, fN2) 

Computes respiratory heat and water loss via mass flow through the lungs 
given gas concentrations, pressure, respiration rate and humidity.
Note that there is no recovery of heat or moisture assumed in the nose.
If barometric preassure is known, elevation will be ignored. Otherwise, 
if atmospheric pressure is unknown, elevation will be used to estimate it.

# Keywords
- `T_x`: current core temperature guess, K
- `Q_metab`: metabolic rate, W
- `fO2_extract`: extraction efficiency, fractional
- `pant`: multiplier on breathing rate due to panting, -
- `rq`: respiratory quotient, (mol CO2 / mol O2)
- `T_air`: air temperature, K
- `rh`: relative humidity, fractional
- `elevation`: elevation, m
- `P_atmos`: barometric pressure, Pa
- `fO2`; fractional O2 concentration in atmosphere, -
- `fCO2`; fractional CO2 concentration in atmosphere, -
- `fN2`; fractional N2 concentration in atmosphere, -
"""
function respiration(;
    T_x,
    Q_metab,
    fO2_extract,
    pant,
    rq,
    T_air,
    rh,
    P_atmos,
    fO2,
    fCO2,
    fN2,
)
    return respiration(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, P_atmos, fO2, fCO2, fN2)
end
function respiration(T_x, Q_metab, fO2_extract, pant, rq, T_air, rh, P_atmos, fO2, fCO2, fN2)
    # adjust O2 to ensure sum to 1
    if fO2 + fCO2 + fN2 != 1
        fO2 = 1 - (fN2 + fCO2)
    end
    P_O2 = P_atmos * fO2
    Joule_m3_O2 = 20.1e6J / m^3 # joules of energy dissipated per m3 O2 consumed at STP (enthalpy of combustion)
    V_O2_STP = uconvert(u"m^3/s", Q_metab / Joule_m3_O2)

    # converting stp -> vol. of O2 at animal lung temperature, atm. press.
    T_lung = T_x
    V_O2 = (V_O2_STP * P_O2 / 273.15K) * (T_lung / P_O2)
    #n = PV/RT (ideal gas law: number of moles from press,vol,temp)
    J_O2 = uconvert(u"mol/s", P_atmos * V_O2 / (R * T_x)) # mol O2 consumed
    # moles/s of O2, N2, dry air at entrance [air flow = f(O2 consumption)]
    J_O2_in = J_O2 / fO2_extract # actual oxygen flow in (moles/s), accounting for efficiency of extraction
    J_N2_in = J_O2_in * (fN2 / fO2) #  actual nitrogen flow in (moles/s), accounting for efficiency of extraction
    V_air = V_O2 / fO2 # air flow
    V_CO2 = fCO2 * V_air #O2 flow
    J_CO2_in = P_atmos * V_CO2 / (R * T_lung)
    J_air_in = (J_O2_in + J_N2_in + J_CO2_in) * pant
    V_air = uconvert(u"m^3/s", (J_air_in * R * 273.15u"K" / 101325u"Pa")) # air volume @ stp (m3/s)
    # computing the vapor pressure at saturation for the subsequent calculation of 
    # actual moles of water based on actual relative humidity
    #wet_air_out = wet_air_properties(T_air, rh, P_atmos; fO2, fCO2, fN2)
    P_vap_sat = vapour_pressure(T_air)
    J_H2O_in = J_air_in * (P_vap_sat * rh) / (P_atmos - P_vap_sat * rh)
    # moles at exit
    J_O2_out = J_O2_in - J_O2 # remove consumed oxygen from the total
    J_N2_out = J_N2_in
    J_CO2_out = rq * J_O2 + J_CO2_in
    # total moles of air at exit will be approximately the same as at entrance, since 
    # the moles of O2 removed = approx. the # moles of co2 added
    J_air_out = (J_O2_out + J_N2_out + J_CO2_out) * pant
    # setting up call to wet_air_properties using temperature of exhaled air at body temperature, assuming saturated air
    rh_exit = 1.0
    #wet_air_out = wet_air_properties(T_x, rh_exit, P_atmos; fO2, fCO2, fN2)
    P_vap_sat = vapour_pressure(T_x)
    J_H2O_out = J_air_out * (P_vap_sat / (P_atmos - P_vap_sat))
    # enthalpy = U2-U1, internal energy only, i.e. lat. heat of vap. only involved, since assume 
    # P,T,V constant, so not significant flow energy, PV. (H = U + PV)

    # moles/s lost by breathing:
    J_evap = J_H2O_out - J_H2O_in
    # grams/s lost by breathing = moles lost * gram molecular weight of water:
    m_resp = J_evap * 18u"g/mol"
    # get latent heat of vapourisation and compute heat exchange due to respiration
    #L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(T_lung) - 273.15))J / kg # from wet_air_properties
    L_v = enthalpy_of_vaporisation(T_lung)
    # heat loss by breathing (J/s)=(J/kg)*(kg/s)
    Q_resp = uconvert(u"W", L_v * m_resp)

    return (;Q_resp, m_resp, J_air_in, J_air_out, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out, J_CO2_in, J_CO2_out)
end

metabolism(; mass, T_core, M1=0.013, M2=0.8, M3=0.038, M4=0.0) = metabolism(mass, T_core, M1, M2, M3, M4)
function metabolism(mass, T_core, M1, M2, M3, M4)
    mass_g = uconvert(u"g", mass)
    T_core = uconvert(u"°C", T_core)
    if T_core > 1u"°C"
        if T_core > 50u"°C"
            V_O2 = M1 * Unitful.ustrip(u"g", mass_g)^M2 * 10^(M3 * 50) * 10 ^ M4
         else
            V_O2 = M1 * Unitful.ustrip(u"g", mass_g)^M2 * 10^(M3 * Unitful.ustrip(u"°C", T_core)) * 10 ^ M4
        end
    else
        V_O2 = M1 * Unitful.ustrip(u"g", mass_g)^M2 * 10^(M3 * 1) * 10 ^ M4
    end
    Q_metab = (0.0056 * V_O2)u"W"

    V_O2 = (V_O2)u"ml/hr"

    return (; Q_metab, V_O2)
end

get_Tsurf_Tlung(body::AbstractBody, k_body, Q_gen_spec, T_core) = get_Tsurf_Tlung(shape(body), body, k_body, Q_gen_spec, T_core)
function get_Tsurf_Tlung(shape::Cylinder, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[2]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

    return (; T_surface, T_lung)  
end
function get_Tsurf_Tlung(shape::DesertIguana, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

    return (; T_surface, T_lung)  
end
function get_Tsurf_Tlung(shape::LeopardFrog, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

    return (; T_surface, T_lung)  
end
function get_Tsurf_Tlung(shape::Ellipsoid, body, k_body, Q_gen_spec, T_core)
    a = body.geometry.lengths[1] ^ 2
    b = body.geometry.lengths[2] ^ 2
    c = body.geometry.lengths[3] ^ 2
    x = ((a * b * c) / (a * b + a * c + b * c))
    T_surface = T_core - (Q_gen_spec / (2 * k_body)) * x
    T_lung = (Q_gen_spec / (4 * k_body)) * x + T_surface

    return (; T_surface, T_lung)  
end
#= function get_Tsurf_Tlung(shape::Sphere, body, k_body, Q_gen_spec, T_core)
    R_flesh = body.geometry.lengths[2]
    T_surface = T_core - (Q_gen_spec * R_flesh ^ 2) / (6 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (12 * k_body) + T_surface
    (T_surface = T_surface, T_lung = T_lung) 
end =#
# function get_Tsurf_Tlung(shape::DesertIguana, R_flesh, k_body, Q_gen_spec, T_core)
#     # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
#     T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
#     T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

#     return (; T_surface, T_lung)  
# end
# function get_Tsurf_Tlung(shape::LeopardFrog, R_flesh, k_body, Q_gen_spec, T_core)
#     # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
#     T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
#     T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

#     return (; T_surface, T_lung) 
# end

###############################################################################
#                          endotherm model code                               #
###############################################################################

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
    if isnothing(minimum_metabolic_rate) || minimum_metabolic_rate === missing
        mouse_elephant = (10^(-1.462 + 0.675 * log10(ustrip(u"g", mass))) * metabolic_multiplier)u"W" # TODO make a set of allometric models of metabolic rate
        Q_gen_min = mouse_elephant * Q10^((ustrip(u"°C", core_temperature) - 37) / 10)
    end

    # constants
    a_coef = 0.6
    b_coef = 0.5
    c_p_air = 1005.8u"J/kg/K"

    # body geometry
    V = mass / density # volume
    b = ((3 * V) / (4 * π * posture))^(1/3)
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
    A_o = 2 * π * b_o^2 + 2 * π * ((a_o * b_o) / e_o) * asin(e_o) # outer area of ellipsoid # TODO use general function
    R_ins = (b_o - b) / (k_fur * A_o) # insulation radius

    # air properties
    (; μ, k_air, ρ_air)  = dry_air_properties(T_f, P_atmos)

    V = (4 / 3) * π * a * b * c # recompute volume
    L_c = V^(1/3) # characteristic dimension
    e = sqrt(a^2 - c^2) / a # eccentricity
    A = 2 * π * b^2 + 2 * π * ((a * b) / e) * asin(e) # area of skin

    Re = ρ_air * v * L_c / μ # Reynold's number
    Pr = (μ * c_p_air) / k_air # Prandtl number

    q′′′_numerator = 2 * A * k_b * k_air * (2 + a_coef * (Re^b_coef) * Pr^(1/3)) * (T_c - T_f)
    q′′′_denominator = 2 * k_b * L_c * V + A * S2 * k_air * (2 + a_coef * Re^b_coef * Pr^(1/3))
    q′′′ = q′′′_numerator / q′′′_denominator
    T_s = T_c - (q′′′ * S2) / (2 * k_b) # skin temperature, Eq. 6

    Gr = abs((ρ_air^2) * (1 / T_f) * Unitful.gn * (L_c^3) * (T_s - T_f) / μ^2) # Grashof number
    Nu_free = 2 + 0.6 * (Gr^0.25) * (Pr^(1/3))
    Nu_forced = 0.37 * Re^0.6
    Nu_total = (Nu_free^3 + Nu_forced^3)^(1/3) #  Eq. 24
    h_cv = Nu_total * k_air / L_c # Eq. 25
    R_cv = 1 / (h_cv * A_o)
    R_rad = u"K/W"(1 / (4 * A_o * 0.95 * Unitful.σ * T_f^3)) # Eq. 39
    R_total = R_b + R_ins + (R_cv * R_rad) / (R_cv + R_rad)

    upper_critical_air_temperature = T_c - (Q_gen_min * stress_factor * R_total)
    lower_critical_air_temperature = T_c - Q_gen_min * R_total

    Q_gen_required = (T_c - T_f) / R_total
    Q_gen_final = max(Q_gen_required, Q_gen_min)

    O2_consumption_rate = Q_gen_final / 20.1u"J/mL" # TODO make general function to convert metabolic rate to O2
    ρ_vap_c = wet_air_properties(T_c, 1.0, P_atmos).ρ_vap
    ρ_vap_f = wet_air_properties(T_f, relative_humidity, P_atmos).ρ_vap

    respiratory_water_loss_rate = (O2_consumption_rate / f_O2 / oxygen_extraction_efficiency) * (ρ_vap_c - ρ_vap_f)
    latent_heat = enthalpy_of_vaporisation(T_f)
    Q_respiration = u"W"(respiratory_water_loss_rate * latent_heat)

    basal_metabolic_rate_fraction = Q_gen_final / Q_gen_min

    Q_evap = max(Q_respiration, (-Q_gen_final) + Q_gen_min)
    total_water_loss_rate = max((u"J/hr"(Q_evap) / latent_heat), 0.0u"kg/hr")
    fractional_mass_loss = total_water_loss_rate / mass

    return (;
        skin_temperature = u"°C"(T_s),
        upper_critical_air_temperature = u"°C"(upper_critical_air_temperature),
        lower_critical_air_temperature = u"°C"(lower_critical_air_temperature),
        Q_gen_required,
        Q_gen_final,
        Q_respiration,
        Q_evap,
        O2_consumption_rate = u"mL/hr"(O2_consumption_rate),
        respiratory_water_loss_rate = u"g/hr"(respiratory_water_loss_rate),
        total_water_loss_rate = u"g/hr"(total_water_loss_rate),
        basal_metabolic_rate_fraction,
        fractional_mass_loss,
    )
end
