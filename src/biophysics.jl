# Biophysics

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
    Ra = (Gr ^ (1 /4)) * (Pr ^ (1 / 3))
    Nu_free = 2 + 0.60 * Ra
    # if Ra >= 200
    # print(Ra, '(Gr ^ 0.25) * (Pr ^ 0.333) is too large for correlation'))
    # end
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
    Nu_free = nusselt_free(body.shape, Gr, Pr)
    hc_free = (Nu_free * k_fluid) / D # heat transfer coefficient, free
    # calculating the Sherwood number from the Colburn analogy
    # Bird, Stewart & Lightfoot, 1960. Transport Phenomena. Wiley.
    Sh_free = Nu_free * (Sc / Pr)^(1 / 3) # Sherwood number, free
    # calculating the mass transfer coefficient from the Sherwood number
    hd_free = Sh_free * D_w / D # mass transfer coefficient, free
    Q_free = hc_free * A_convection * (T_surface - T_air) # free convective heat loss at surface

    # forced convection
    Nu = nusselt_forced(body.shape, Re)
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
    evaporation(T_core, T_surface, ψ_org, skin_wetness, A_total, hd, eye_fraction, T_air, rh, elevation, P_atmos)
"""
function evaporation(;
    T_surface,
    ψ_org,
    skin_wetness,
    A_total,
    hd,
    eye_fraction,
    T_air,
    rh,
    P_atmos,
    fO2,
    fCO2,
    fN2,
)
    return evaporation(T_surface, ψ_org, skin_wetness, A_total, hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2)
end
function evaporation(T_surface, ψ_org, skin_wetness, A_total, hd, eye_fraction, T_air, rh, P_atmos, fO2, fCO2, fN2)
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
    m_eyes = hd * eye_fraction * A_total * (ρ_vap_surf - ρ_vap_air)
    if m_eyes > 0.0u"kg/s"
        m_cut = (A_total * skin_wetness - A_total * skin_wetness * eye_fraction) * hd * (ρ_vap_surf - ρ_vap_air)
    else
        m_cut = A_total * skin_wetness * hd * (ρ_vap_surf - ρ_vap_air)
    end

    # get latent heat of vapourisation and compute heat exchange due to evaporation
    L_v = enthalpy_of_vaporisation(T_air)
    Q_evap = Unitful.uconvert(u"W", (m_eyes + m_cut) * L_v)

    #onvert from kg/s to g/s
    m_eyes = uconvert(u"g/s", m_eyes)
    m_cut = uconvert(u"g/s", m_cut)

    return (; Q_evap, m_cut, m_eyes)
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

    return (; Q_resp, m_resp, J_air_in, J_air_out, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out, J_CO2_in, J_CO2_out)
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

Tsurf_and_Tlung(body::AbstractBody, k_body, Q_gen_spec, T_core) = Tsurf_and_Tlung(shape(body), body, k_body, Q_gen_spec, T_core)
function Tsurf_and_Tlung(shape::Cylinder, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[2]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

    return (; T_surface, T_lung)  
end
function Tsurf_and_Tlung(shape::DesertIguana, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

    return (; T_surface, T_lung)  
end
function Tsurf_and_Tlung(shape::LeopardFrog, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[1]
    T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

    return (; T_surface, T_lung)  
end
function Tsurf_and_Tlung(shape::Ellipsoid, body, k_body, Q_gen_spec, T_core)
    a = body.geometry.lengths[1] ^ 2
    b = body.geometry.lengths[2] ^ 2
    c = body.geometry.lengths[3] ^ 2
    x = ((a * b * c) / (a * b + a * c + b * c))
    T_surface = T_core - (Q_gen_spec / (2 * k_body)) * x
    T_lung = (Q_gen_spec / (4 * k_body)) * x + T_surface

    return (; T_surface, T_lung)  
end
#= function Tsurf_and_Tlung(shape::Sphere, body, k_body, Q_gen_spec, T_core)
    R_flesh = body.geometry.lengths[2]
    T_surface = T_core - (Q_gen_spec * R_flesh ^ 2) / (6 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (12 * k_body) + T_surface
    (T_surface = T_surface, T_lung = T_lung) 
end =#
# function Tsurf_and_Tlung(shape::DesertIguana, R_flesh, k_body, Q_gen_spec, T_core)
#     # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
#     T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
#     T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

#     return (; T_surface, T_lung)  
# end
# function Tsurf_and_Tlung(shape::LeopardFrog, R_flesh, k_body, Q_gen_spec, T_core)
#     # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
#     T_surface = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
#     T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surface 

#     return (; T_surface, T_lung) 
# end

radiant_temperature(body::AbstractBody, insulation, Q_evap_skin, T_core, T_skin, T_conduction, T_insulation, CD, longwave_depth_fraction, RRAD, i_ins) = radiant_temperature(shape(body), body, insulation, Q_evap_skin, T_core, T_skin, T_conduction, T_insulation, CD, longwave_depth_fraction, RRAD, i_ins)
function radiant_temperature(shape::Cylinder, body, insulation, Q_evap_skin, T_core, T_skin, T_conduction, T_insulation, CD, longwave_depth_fraction, RRAD, i_ins)

    k_compressed = insulation.insulation_conductivity_compressed
    r_skin = body.geometry.
    compression_fraction =
        (conduction_fraction * 2 * π * k_compressed * fibre_length) /
        log(r_compressed / r_skin)

    T_ins_compressed = conduction_fraction > 0 ?
        (compression_fraction * T_skin + CD * T_conduction) / (CD + compression_fraction) :
        0.0

    cd1 = (k_compressed / log(r_compressed / r_skin)) * conduction_fraction +
          (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)

    cd2 = (k_compressed / log(r_compressed / r_skin)) * conduction_fraction
    cd3 = (k_insulation / log(r_insulation / r_skin)) * (1 - conduction_fraction)

    dv1 = 1 +
        ((2 * π * fibre_length * r_flesh^2 * cd1) / (4 * k_flesh * volume)) +
        ((2 * π * fibre_length * r_flesh^2 * cd1) / (2 * k_fat * volume)) *
        log(r_skin / r_flesh)

    dv2 = Q_evap_skin * ((r_flesh^2 * cd1) / (4 * k_flesh * volume)) +
          Q_evap_skin * ((r_flesh^2 * cd1) / (2 * k_fat * volume)) *
          log(r_skin / r_flesh)

    dv3 = ((2 * π * fibre_length) / dv1) *
        (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) *
        r_flesh^2 / (2 * volume)

    dv4 = longwave_depth_fraction < 1 ?
        cd2 + (k_insulation / log(r_insulation / RRAD)) * (1 - conduction_fraction) :
        1.0

    return longwave_depth_fraction < 1 ?
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (T_insulation *
            ((k_insulation / log(r_insulation / RRAD)) *
            (1 - conduction_fraction))) / dv4 :
        T_insulation
end

function radiant_temperature(
    ::Sphere,
    conduction_fraction, k_compressed,
    r_compressed, r_skin, k_insulation, r_insulation,
    r_flesh, k_flesh, k_fat, volume,
    Q_evap_skin, 
    T_core, T_skin, T_conduction, T_insulation,
    CD, longwave_depth_fraction, RRAD
)

    compression_fraction =
        (conduction_fraction * 4 * π * k_compressed * r_compressed * r_skin) /
        (r_compressed - r_skin)

    T_ins_compressed = conduction_fraction > 0 ?
        (compression_fraction * T_skin + CD * T_conduction) / (CD + compression_fraction) :
        0.0

    cd1 = ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction +
          ((k_insulation * r_insulation) / (r_insulation - r_skin)) * (1 - conduction_fraction)

    cd2 = ((k_compressed * r_compressed) / (r_compressed - r_skin)) * conduction_fraction
    cd3 = ((k_insulation * r_insulation) / (r_insulation - r_skin)) * (1 - conduction_fraction)

    dv1 = 1 +
        ((4 * π * r_skin * r_flesh^2 * cd1) / (6 * k_flesh * volume)) +
        ((4 * π * r_skin * r_flesh^3 * cd1) / (3 * k_fat * volume)) *
            ((r_skin - r_flesh) / (r_flesh * r_skin))

    dv2 = Q_evap_skin * ((r_flesh^2 * cd1) / (6 * k_flesh * volume)) +
          Q_evap_skin * ((r_flesh^3 * cd1) / (3 * k_fat * volume)) *
          ((r_skin - r_flesh) / (r_flesh * r_skin))

    dv3 =
        ((4 * π * r_skin) / dv1) *
        (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) *
        r_flesh^3 / (3 * volume * RRAD)

    dv4 = longwave_depth_fraction < 1 ?
        cd2 + ((k_insulation * r_insulation) / (r_insulation - RRAD)) * (1 - conduction_fraction) :
        1.0

    return longwave_depth_fraction < 1 ?
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (T_insulation *
            ((k_insulation * r_insulation) / (r_insulation - RRAD) *
            (1 - conduction_fraction))) / dv4 :
        T_insulation
end

function radiant_temperature(
    ::Ellipsoid,
    conduction_fraction, k_compressed,
    ASEMAJ, BSEMIN, CSEMIN, FATTHK, SUBQFAT, ZL,
    BLCMP, r_flesh, k_flesh, k_fat, volume,
    Q_evap_skin, 
    T_core, T_skin, T_conduction, T_insulation,
    CD, longwave_depth_fraction
)

    FLSHASEMAJ = ASEMAJ - FATTHK
    FLSHBSEMIN = BSEMIN - FATTHK
    FLSHCSEMIN = CSEMIN - FATTHK

    ASQG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHASEMAJ^2 : ASEMAJ^2
    BSQG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHBSEMIN^2 : BSEMIN^2
    CSQG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHCSEMIN^2 : CSEMIN^2

    SSQG = (ASQG * BSQG * CSQG) / (ASQG * BSQG + ASQG * CSQG + BSQG * CSQG)

    BG = (Int(SUBQFAT) == 1 && FATTHK > 0) ? FLSHBSEMIN : BSEMIN
    BS = BSEMIN
    BL = BSEMIN + ZL
    BR = BS + longwave_depth_fraction * ZL

    compression_fraction =
        (conduction_fraction * 3 * k_compressed * volume * BLCMP * BS) /
        ((sqrt(3 * SSQG))^3 * (BLCMP - BS))

    T_ins_compressed = conduction_fraction > 0 ?
        (compression_fraction * T_skin + CD * T_conduction) / (CD + compression_fraction) :
        0.0

    cd1 = ((k_compressed * BLCMP) / (BLCMP - BS)) * conduction_fraction +
          ((k_insulation * BL) / (BL - BS)) * (1 - conduction_fraction)

    cd2 = ((k_compressed * BLCMP) / (BLCMP - BS)) * conduction_fraction
    cd3 = ((k_insulation * BL) / (BL - BS)) * (1 - conduction_fraction)

    dv1 = 1 +
        (3 * BS * SSQG * cd1) / (2 * k_flesh * (sqrt(3 * SSQG)^3)) +
        (BS * cd1) / k_fat * ((BS - BG) / (BS * BG))

    dv2 = Q_evap_skin *
        ((SSQG * cd1) / (2 * k_flesh * volume)) +
        Q_evap_skin *
        ((sqrt(3 * SSQG)^3 * cd1) / (3 * k_fat * volume)) *
        ((BS - BG) / (BS * BG))

    dv3 = (BS / dv1) * (T_core * cd1 - dv2 - T_ins_compressed * cd2 - T_insulation * cd3) / BR

    dv4 = longwave_depth_fraction < 1 ?
        cd2 + ((k_insulation * BL) / (BL - BR)) * (1 - conduction_fraction) :
        1.0

    return longwave_depth_fraction < 1 ?
        dv3 / dv4 +
        (T_ins_compressed * cd2) / dv4 +
        (T_insulation *
            ((k_insulation * BL) / (BL - BR) *
            (1 - conduction_fraction))) / dv4 :
        T_insulation
end


