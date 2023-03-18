# Biophysics
"""
    vapour_pressure

Calculates saturation vapour pressure (Pa) for a given air temperature.
    ...
    # Arguments
    - `T`: air temperature in K.
    ...    
"""
function vapour_pressure(T)
    T = Unitful.ustrip(T)# + 273.15
    logP_vap = T
    if T <= 273.15
        logP_vap = -9.09718 * (273.16 / T - 1) - 3.56654 * log10(273.16 / T) + 0.876793 * (1 - T / 273.16) + log10(6.1071)
    else
        logP_vap = -7.90298 * (373.16 / T - 1) + 5.02808 * log10(373.16 / T) - 1.3816E-07 * (10^(11.344 * (1 - T / 373.16)) - 1) + 8.1328E-03 * (10^(-3.49149 * (373.16 / T - 1)) - 1) + log10(1013.246)
    end
    (10^logP_vap) * 100Pa
end

"""
    wet_air

Calculates several properties of humid air as output variables below. The program
is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
Institution Press. Washington, DC. wet_air must be used in conjunction with function vapour_pressure.

Input variables are shown below. The user must supply known values for T_drybulb and P (P at one standard
atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
either (1) psychrometric data (T_wetbulb or rh), or (2) hygrometric data (T_dew)

(1) Psychrometric data:
If T_wetbulb is known but not rh, then set rh=-1 and dp=999
If rh is known but not T_wetbulb then set T_wetbulb=0 and dp=999

(2) Hygrometric data:
If T_dew is known then set T_wetublb = 0 and rh = 0.
...
# Arguments
- `T_drybulb`: Dry bulb temperature (K)
- `T_wetbulb`: Wet bulb temperature (K)
- `rh`: Relative humidity (%)
- `T_dew`: Dew point temperature (K)
- `P`: Barometric pressure (Pa)
# - `P_vap`: Vapour pressure (Pa)
# - `P_vap_sat`: Saturation vapour pressure (Pa)
# - `ρ_vap`: Vapour density (kg m-3)
# - `r_w Mixing`: ratio (kg kg-1)
# - `T_vir`: Virtual temperature (K)
# - `T_vinc`: Virtual temperature increment (K)
# - `ρ_air`: Density of the air (kg m-3)
# - `cp`: Specific heat of air at constant pressure (J kg-1 K-1)
# - `ψ`: Water potential (Pa)
# - `rh`: Relative humidity (%)
...    
"""
function wet_air(T_drybulb, T_wetbulb=T_drybulb, rh=0, T_dew=nothing, P_atmos=101325Pa)
    f_w = 1.0053 # (-) correction factor for the departure of the mixture of air and water vapour from ideal gas laws
    M_w = 0.018016kg/mol # molar mass of water
    M_a = 0.028965924869122257kg/mol # molar mass of air
    P_vap_sat = vapour_pressure(T_drybulb)
    if T_dew !== nothing
        P_vap = vapour_pressure(T_dew)
        rh = (P_vap / P_vap_sat) * 100
    else
        if rh !== nothing
            P_vap = P_vap_sat * rh / 100
        else
            δ_bulb = T_drybulb - T_wetbulb
            δ_P_vap = (0.000660 * (1 + 0.00115 * (Unitful.ustrip(T_wetbulb)-273.15)) * Unitful.ustrip(P) * Unitful.ustrip(δ_bulb))Pa
            P_vap = vapour_pressure(T_wetbulb) - δ_P_vap
            relhumid = (P_vap / P_vap_sat) * 100
        end
    end
    r_w = ((0.62197 * f_w * P_vap) / (P_atmos - f_w * P_vap))kg/kg
    ρ_vap = P_vap * M_w / (0.998 * Unitful.R * T_drybulb)
    ρ_vap = Unitful.uconvert(u"kg/m^3",ρ_vap) # simplify units
    T_vir = T_drybulb * ((1.0 + r_w / (18.016 / 28.966)) / (1 + r_w))
    T_vinc = T_vir - T_drybulb
    ρ_air = (M_a / Unitful.R) * P_atmos / (0.999 * T_vir)
    ρ_air = Unitful.uconvert(u"kg/m^3",ρ_air) # simplify units
    cp = ((1004.84 + (r_w * 1846.40)) / (1 + r_w))J/K/kg
    if min(rh) <= 0
        ψ = -999Pa
    else
        ψ = (4.615e+5 * Unitful.ustrip(T_drybulb) * log(rh / 100))Pa
    end
    (;P_vap, P_vap_sat, ρ_vap, r_w, T_vinc, ρ_air, cp, ψ, rh)
end

function dry_air(T_drybulb, P_atmos=nothing, elev=0m)
    σ = Unitful.uconvert(u"W/m^2/K^4",Unitful.σ) # Stefan-Boltzmann constant, W/m^2/K^4, extract σ when calling Unitful when units issue is fixed in Unitful
    M_a = 0.028965924869122257kg/mol # molar mass of air
    if P_atmos === nothing
     P_std = 101325Pa
     P_atmos = P_std * ((1 - (0.0065 * elev / 288m))^(1 / 0.190284))
    end
    ρ_air = (M_a / Unitful.R) * P_atmos / (T_drybulb)
    ρ_air = Unitful.uconvert(u"kg/m^3",ρ_air) # simplify units
    vis_not = 1.8325e-5kg/m/s
    T_not = 296.16K
    c = 120K
    μ = (vis_not * (T_not + c) / (T_drybulb + c)) * (T_drybulb / T_not)^1.5 # kg / m.s
    ν = μ / ρ_air # m2 / s or J.s/kg
    dif_vpr = 2.26e-5m^2/s * ((T_drybulb / 273.15K)^1.81) * (1.e5Pa / P_atmos) # m2 / s
    k_fluid = (0.02425 + (7.038e-5 * (Unitful.ustrip(T_drybulb) - 273.15)))W/m/K
    L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(T_drybulb) - 273.15))J/kg

    tcoeff = 1 / T_drybulb
    ggroup = 0.0980616m/s^2 * tcoeff / (ν^2) # 1 / m3.K
    bbemit = σ * ((T_drybulb)^4) # W/m2
    emtmax = 2.897e-3K*m / (T_drybulb) # m
    (;P_atmos, ρ_air, μ, ν, dif_vpr, k_fluid, L_v, tcoeff, ggroup, bbemit, emtmax)
end

function get_Nusselt_free(shape::Cylinder, Gr, Pr)
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
    (Nu_free)
end

function get_Nusselt_free(shape::Plate, Gr, Pr)
 Ra = Gr * Pr
 Nu_free = 0.55 * Ra ^ 0.25
 (Nu_free)
end

function get_Nusselt_free(shape::Ellipsoid, Gr, Pr)
    #  sphere free convection
    #  from p.413 Bird et all (1960) Transport Phenomena
    Ra = (Gr ^ (1 /4)) * (Pr ^ (1 / 3))
    Nu_free = 2 + 0.60 * Ra
    # if Ra >= 200
    # print(Ra, '(Gr ^ 0.25) * (Pr ^ 0.333) is too large for correlation'))
    # end
    (Nu_free)
end

function get_Nusselt_forced(shape::Cylinder, Re)
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
    (Nu)
end

function get_Nusselt_forced(shape::Plate, Re)
 # forced convection of a plate
 Nu = 0.102 * Re ^ 0.675 * Pr ^ (1 / 3)
 (Nu)
end

function get_Nusselt_forced(shape::Ellipsoid, Re)
    #  forced convection of a sphere
    Nu = 0.35 * Re ^ 0.6 # from McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
    (Nu)
end

function water_prop(T_water)
    # β = coefficient of expansion (1/K)
    # cp_fluid = specific heat (J/kg-K)
    # ρ_water = water density (kg/m3)
    # k_fluid = thermal conductivity (W/m-K)
    # μ = dynamic viscosity (kg/m-s)
    β = 0.21E-031 / K
    T_water = Unitful.ustrip(T_water) - 273.15
    cp_fluid = 4220.02 - 4.5531 * T_water + 0.182958 * T_water^2 - 0.00310614 * T_water^3 + 1.89399E-5 * T_water^4
    cp_fluid = (cp_fluid)J / kg / K
    if T_water < 30
        ρ_water = 1000.0
    else
        if T_water <= 60
            ρ_water = 1017 - 0.6 * T_water
        else
            T_water = 60
            ρ_water = 1017 - 0.6 * T_water
        end
    end
    k_fluid = (0.551666 + 0.00282144 * T_water - 2.02383E-5 * T_water^2)W / m / K
    μ = (0.0017515 - 4.31502E-5 * T_water + 3.71431E-7 * T_water^2)kg / m / s
    ρ_water = (ρ_water)kg / m^3
    (;β, cp_fluid, ρ_water, k_fluid, μ)
end

function convection(Body, A_conv, T_air, T_surf, vel, P_atmos, elev, fluid)
    shape = Body.geometry
    G = Unitful.gn # acceleration due to gravity, m.s^2
    β = 1 / T_air
    D = shape.characteristic_dimension
    dry_air_out = dry_air(T_air, P_atmos, elev)
    dif_vpr = dry_air_out.dif_vpr
    # checking to see if the fluid is water, not air
    if fluid == 1
        water_prop_out = water_prop(T_air)
        cp_fluid = water_prop_out.cp_fluid
        ρ_air = water_prop_out.ρ_water
        k_fluid = water_prop_out.k_fluid
        μ = water_prop_out.μ
    else
        cp_fluid = 1.0057E+3J/K/kg
        ρ_air = dry_air_out.ρ_air
        k_fluid = dry_air_out.k_fluid
        μ = dry_air_out.μ
    end

    # free convection
    Pr = cp_fluid * μ / k_fluid
    Pr = Unitful.ustrip(Pr) # is dimensionless, enforce this
    if fluid == 1
        #  water; no meaning
        Sc = 1
    else
        Sc = μ / (ρ_air * dif_vpr)
    end
    δ_T = T_surf - T_air
    Gr = abs(((ρ_air^2) * β * G * (D^3) * δ_T) / (μ^2))
    Re = ρ_air * vel * D / μ
    Nu_free = get_Nusselt_free(Body.shape, Gr, Pr)
    Hc_free = (Nu_free * k_fluid) / D # heat transfer coefficient, free
    # calculating the Sherwood number from the Colburn analogy
    # Bird, Stewart & Lightfoot, 1960. Transport Phenomena. Wiley.
    Sh_free = Nu_free * (Sc / Pr)^(1 / 3) # Sherwood number, free
    # calculating the mass transfer coefficient from the Sherwood number
    Hd_free = Sh_free * dif_vpr / D # mass transfer coefficient, free
    Q_free = Hc_free * A_conv * (T_surf - T_air) # free convective heat loss at surface

    # forced convection
    Nu = get_Nusselt_forced(Body.shape, Re)
    # forced convection for object
    Hc_forc = Nu * k_fluid / D # heat transfer coefficient, forced
    Sh_forc = Nu * (Sc / Pr)^(1 / 3) # Sherwood number, forced
    Hd_forc = Sh_forc * dif_vpr / D # mass transfer coefficient
    Q_forc = Hd_forc * A_conv * (T_surf - T_air) # forced convective heat transfer

    # combined free and forced convection
    # using Bird, Stewart & Lightfoot's mixed convection formula (p. 445, Transport Phenomena, 2002)
    Nu_comb = (Nu_free^3 + Nu^3)^(1 / 3)
    Hc = Nu_comb * (k_fluid / D) # mixed convection heat transfer
    Q_conv = Hc * A_conv * (T_surf - T_air) # total convective heat loss
    Sh = Nu_comb * (Sc / Pr)^(1 / 3) # Sherwood number, combined
    Hd = Sh * dif_vpr / D # mass transfer coefficient, combined
    (;Q_conv, Hc, Hd, Sh, Q_free, Q_forc, Hc_free, Hc_forc, Sh_free, Sh_forc, Hd_free, Hd_forc)
end

function conduction(A_cond = 0.001325006m^2,
    L=0.025m,
    T_surf=K(25°C),
    T_sub=K(10°C),
    k_sub=0.1W/m/K)
    A_cond * (k_sub / L) * (T_surf - T_sub)
end

function solar(α_org_dorsal=0.85,
    α_org_ventral=0.85,
    A_sil=0.004718043m^2,
    A_tot=0.01325006m^2,
    A_cond=0.001325006m^2,
    F_sub=0.4,
    F_sky=0.4,
    α_sub=0.8,
    Q_sol=1000W / m^2,
    Q_dir=964.177772475912W / m^2,
    Q_dif=100W / m^2)

    Q_direct = α_org_dorsal * A_sil * Q_dir
    Q_sol_sky = α_org_dorsal * F_sky * A_tot * Q_dif
    Q_sol_sub = α_org_ventral * F_sub * (A_tot - A_cond) * (1 - α_sub) * Q_sol
    Q_solar = (Q_direct + Q_sol_sub + Q_sol_sky)
    (;Q_solar, Q_direct, Q_sol_sky, Q_sol_sub)
end

function radin(A_tot=0.01325006m^2,
    A_cond=0.001325006m^2,
    F_sky=0.4,
    F_sub=0.4,
    ϵ_org_dorsal=0.95,
    ϵ_org_ventral=0.95,
    ϵ_sub=0.95,
    ϵ_sky=0.8,
    T_sky=(10 + 273.15)K,
    T_sub=(30 + 273.15)K)

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ) # Stefan-Boltzmann constant, W/m^2/K^4, extract σ when calling Unitful when units issue is fixed in Unitful
    Q_ir_sky = ϵ_org_dorsal * F_sky * A_tot * ϵ_sky * σ * T_sky^4
    Q_ir_sub = ϵ_org_ventral * F_sub * (A_tot - A_cond) * ϵ_sub * σ * T_sub^4
    Q_ir_in = Q_ir_sky + Q_ir_sub
    (;Q_ir_in, Q_ir_sky, Q_ir_sub)
end

function radout(
    T_surf=(25.1 + 273.15)K,
    A_tot=0.01325006m^2,
    A_cond=0.001325006m^2,
    F_sky=0.4,
    F_sub=0.4,
    ϵ_org_dorsal=0.95,
    ϵ_org_ventral=0.95)
    # computes longwave radiation lost
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ) # Stefan-Boltzmann constant, W/m^2/K^4, extract σ when calling Unitful when units issue is fixed in Unitful
    Q_ir_to_sky = A_tot * F_sky * ϵ_org_dorsal * σ * T_surf^4
    Q_ir_to_sub = (A_tot - A_cond) * F_sub * ϵ_org_ventral * σ * T_surf^4
    Q_ir_out = Q_ir_to_sky + Q_ir_to_sub
    (;Q_ir_out, Q_ir_to_sky, Q_ir_to_sub)
end

function evaporation(
    T_core=(25 + 273.15)K,
    T_surf=(25.1 + 273.15)K,
    m_resp=1.177235e-09kg / s,
    ψ_org=-7.07 * 100J / kg,
    p_wet=0.1 / 100,
    A_tot=0.01325006m^2,
    Hd=0.02522706m / s,
    p_eyes=0.03 / 100,
    T_air=(20 + 273.15)K,
    rh=50,
    elev=0m,
    P_atmos=101325Pa)

    # this subroutine computes surface evaporation based on the mass transfer
    # coefficient, fraction of surface of the skin acting as a free water surface
    # and exposed to the air, and the vapor density gradient between the
    # surface and the air, each at their own temperature.

    # get vapour density at surface based on water potential of body
    m_w = 0.018kg / mol #! molar mass of water, kg/mol
    rh_surf = exp(ψ_org / (Unitful.R / m_w * T_surf)) * 100 #
    wet_air_out = wet_air(T_surf, 0K, rh_surf, nothing, P_atmos)
    ρ_vap_surf = wet_air_out.ρ_vap

    # get air vapour density
    wet_air_out = wet_air(T_air, 0K, rh, nothing, P_atmos)
    ρ_vap_air = wet_air_out.ρ_vap

    # water lost from eyes if present
    m_eyes = Hd * p_eyes * A_tot * (ρ_vap_surf - ρ_vap_air)
    if m_eyes > 0kg / s
        m_cut = (A_tot * p_wet - A_tot * p_wet * p_eyes) * Hd * (ρ_vap_surf - ρ_vap_air)
    else
        m_cut = A_tot * p_wet * Hd * (ρ_vap_surf - ρ_vap_air)
    end

    # total water lost
    m_evap = m_eyes + m_resp + m_cut

    # get latent heat of vapourisation and compute heat exchange due to evaporation
    dry_air_out = dry_air(T_air, P_atmos, elev)
    L_v = dry_air_out.L_v
    Q_evap = Unitful.uconvert(u"W", (m_eyes + m_cut) * L_v)

    #onvert from kg/s to g/s
    m_eyes = uconvert(u"g/s", m_eyes)
    m_resp = uconvert(u"g/s", m_resp)
    m_cut = uconvert(u"g/s", m_cut)
    m_evap = uconvert(u"g/s", m_evap)

    (;Q_evap, m_evap, m_resp, m_cut, m_eyes)
end

"""
resp

Computes respiratory heat and water loss via mass flow through the lungs 
given gas concentrations, pressure, respiration rate and humidity.
Note that there is no recovery of heat or moisture assumed in the nose.
If barometric preassure is known, elevation will be ignored. Otherwise, 
if atmospheric pressure is unknown, elevation will be used to estimate it.
    ...
    # Arguments
    - `T_x`: current core temperature guess, K
    - `Q_metab`: metabolic rate, W
    - `fO2_extract`: extraction efficiency, fractional
    - `pant`: multiplier on breathing rate due to panting, -
    - `rq`: respiratory quotient, (mol CO2 / mol O2)
    - `T_air`: air temperature, K
    - `rh`: relative humidity, %
    - `elev`: elevation, m
    - `P_atmos`: barometric pressure, Pa
    - `fO2`; fractional O2 concentration in atmosphere, -
    - `fCO2`; fractional CO2 concentration in atmosphere, -
    - `fN2`; fractional N2 concentration in atmosphere, -
    ...    
"""
function respiration(
    T_x=296.15K,
    Q_metab=0.01241022W,
    fO2_extract=0.20,
    pant=1,
    rq=0.8,
    T_air=293.15K,
    rh=50,
    elev=nothing,
    P_atmos=101325Pa,
    fO2=0.2095,
    fCO2=0.0003,
    fN2=0.7902)

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
    V_air = uconvert(u"m^3/s", (J_air_in * R * 273.15K / 101325Pa)) # air volume @ stp (m3/s)
    # computing the vapor pressure at saturation for the subsequent calculation of 
    # actual moles of water based on actual relative humidity
    wet_air_out = wet_air(T_air, 273.15K, rh, nothing, P_atmos)
    P_vap_sat = wet_air_out.P_vap_sat
    J_H2O_in = J_air_in * (P_vap_sat * (rh / 100)) / (P_atmos - P_vap_sat * (rh / 100))
    # moles at exit
    J_O2_out = J_O2_in - J_O2 # remove consumed oxygen from the total
    J_N2_out = J_N2_in
    J_CO2_out = rq * J_O2 + J_CO2_in
    # total moles of air at exit will be approximately the same as at entrance, since 
    # the moles of O2 removed = approx. the # moles of co2 added
    J_air_out = (J_O2_out + J_N2_out + J_CO2_out) * pant
    # setting up call to wet_air using temperature of exhaled air at body temperature, assuming saturated air
    rh_exit = 100
    wet_air_out = wet_air(T_x, 0K, rh_exit, nothing, P_atmos)
    P_vap_sat = wet_air_out.P_vap_sat
    J_H2O_out = J_air_out * (P_vap_sat / (P_atmos - P_vap_sat))
    # enthalpy = U2-U1, internal energy only, i.e. lat. heat of vap. only involved, since assume 
    # P,T,V constant, so not significant flow energy, PV. (H = U + PV)

    # moles/s lost by breathing:
    J_evap = J_H2O_out - J_H2O_in
    # grams/s lost by breathing = moles lost * gram molecular weight of water:
    m_resp = J_evap * 18g / mol
    # get latent heat of vapourisation and compute heat exchange due to respiration
    L_v = (2.5012e6 - 2.3787e3 * (Unitful.ustrip(T_lung) - 273.15))J / kg # from dry_air
    # heat loss by breathing (J/s)=(J/kg)*(kg/s)
    Q_resp = uconvert(u"W", L_v * m_resp)
    (;Q_resp, m_resp, J_air_in, J_air_out, J_H2O_in, J_H2O_out, J_O2_in, J_O2_out, J_CO2_in, J_CO2_out)
end

function metabolism(mass=0.04kg, T_core=K(25°C), M1=0.013, M2=0.8, M3=0.038)
    mass_g = uconvert(u"g", mass)
    T_core = uconvert(u"°C", T_core)
    if T_core > 50°C
        V_O2 = M1 * Unitful.ustrip(mass_g)^M2 * 10^(M3 * 50)
    else
        V_O2 = M1 * Unitful.ustrip(mass_g)^M2 * 10^(M3 * Unitful.ustrip(T_core))
    end
    if T_core < 1°C
        Q_metab = 0.01W
    else
        Q_metab = (0.0056 * V_O2)W
    end
    V_O2 = (V_O2)ml / hr
    (;Q_metab, V_O2)
end

get_Tsurf_Tlung(body::AbstractBody, k_body, Q_gen_spec, T_core) = get_Tsurf_Tlung(shape(body), body, k_body, Q_gen_spec, T_core)

function get_Tsurf_Tlung(shape::Cylinder, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[2]
    T_surf = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surf 
    (;T_surf, T_lung)  
end

function get_Tsurf_Tlung(shape::Ellipsoid, body, k_body, Q_gen_spec, T_core)
    a = body.geometry.lengths[1] ^ 2
    b = body.geometry.lengths[2] ^ 2
    c = body.geometry.lengths[3] ^ 2
    T_surf = T_core - (Q_gen_spec / (2 * k_body)) * ((a * b * c) / (a * b + a * c + b * c))
    T_lung = (Q_gen_spec / (4 * k_body)) * ((a * b * c) / (a * b + a * c + b * c)) + T_surf
    (;T_surf, T_lung)  
end

#= function get_Tsurf_Tlung(shape::Sphere, body, k_body, Q_gen_spec, T_core)
    R_flesh = body.geometry.lengths[2]
    T_surf = T_core - (Q_gen_spec * R_flesh ^ 2) / (6 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (12 * k_body) + T_surf
    (T_surf = T_surf, T_lung = T_lung) 
end =#

function get_Tsurf_Tlung(shape::DesertIguana, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[2]
    T_surf = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surf 
    (;T_surf, T_lung)  
end

function get_Tsurf_Tlung(shape::LeopardFrog, body, k_body, Q_gen_spec, T_core)
    # cylinder: from P. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
    R_flesh = body.geometry.lengths[2]
    T_surf = T_core - Q_gen_spec * R_flesh ^ 2 / (4 * k_body)
    T_lung = (Q_gen_spec * R_flesh ^ 2) / (8 * k_body) + T_surf 
    (;T_surf, T_lung) 
end