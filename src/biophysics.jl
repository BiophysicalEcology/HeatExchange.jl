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
function wet_air(T_drybulb, T_wetbulb=T_drybulb, rh=0, T_dew=999K, P_atmos=101325Pa)
    #T = T_drybulb + 273.15°C
    f_w = 1.0053 # (-) correction factor for the departure of the mixture of air and water vapour from ideal gas laws
    M_w = 0.018016kg/mol # molar mass of water
    M_a = 0.028965924869122257kg/mol # molar mass of air
    P_vap_sat = vapour_pressure(T_drybulb)
    if T_dew < 999K
        P_vap = vapour_pressure(T_dew)
        rh = (P_vap / P_vap_sat) * 100
    else
        if min(rh) > -1
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
    (P_vap=P_vap, P_vap_sat, ρ_vap, r_w, T_vinc, ρ_air, cp, ψ, rh)
end

function dry_air(T_drybulb, P_atmos=101325Pa, elev=0m)
    σ = Unitful.k^4*π^2/(60*Unitful.ħ^3*Unitful.c0^2) # Stefan-Boltzmann constant, W/m^2/K^4, make Unitful.σ when error is fixed in Unitful
    M_a = 0.028965924869122257kg/mol # molar mass of air
    P_std = 101325Pa
    P_atmos = P_std * ((1 - (0.0065 * elev / 288m))^(1 / 0.190284))
    ρ_air = (M_a / Unitful.R) * P_atmos / (T_drybulb)
    ρ_air = Unitful.uconvert(u"kg/m^3",ρ_air) # simplify units
    vis_not = 1.8325e-5kg/m/s
    T_not = 296.16K
    c = 120K
    μ = (vis_not * (T_not + c) / (T_drybulb + c)) * (T_drybulb / T_not)^1.5 # kg / m.s
    ν = μ / ρ_air # m2 / s or J.s/kg
    dif_vpr = 2.26e-5m^2/s * ((T_drybulb / 273.15K)^1.81) * (1.e5Pa / P_atmos) # m2 / s
    k_fluid = (0.02425 + (7.038e-5 * (Unitful.ustrip(T_drybulb) - 273.15)))W/m/K
    L_v = (2.5012E6 - 2.3787e3 * Unitful.ustrip(T_drybulb) - 273.15)J/kg
    tcoeff = 1 / T_drybulb
    ggroup = 0.0980616m/s^2 * tcoeff / (ν^2) # 1 / m3.K
    bbemit = σ * ((T_drybulb)^4) # W/m2
    emtmax = 2.897e-3K*m / (T_drybulb) # m
    (P_atmos=P_atmos, ρ_air=ρ_air, μ=μ, ν=ν, dif_vpr=dif_vpr, k_fluid=k_fluid, L_v=L_v, tcoeff=tcoeff, ggroup=ggroup, bbemit=bbemit, emtmax=emtmax)
end

function get_Nusselt_free(shape::Cylinder, Gr, Pr)
    #C      FREE CONVECTION OF A CYLINDER
    #C      FROM P.334 KREITH (1965): MC ADAM'S 1954 RECOMMENDED COORDINATES
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

# function get_Nusselt_free(shape::Plate, Gr, Pr)
#     Ra = Gr * Pr
#     Nu_free = 0.55 * Ra ^ 0.25
#     (Nu_free)
# end

function get_Nusselt_free(shape::Ellipsoid, Gr, Pr)
    #C      SPHERE FREE CONVECTION
    #C      FROM P.413 BIRD ET AL (1960) TRANSPORT PHENOMENA)
    Ra = (Gr ^ (1 /4)) * (Pr ^ (1 / 3))
    Nu_free = 2 + 0.60 * Ra
    # if Ra >= 200
    #     print(Ra, '(Gr ^ 0.25) * (Pr ^ 0.333) IS TOO LARGE FOR CORREL.'))
    # end
    (Nu_free)
end

function get_Nusselt_forced(shape::Cylinder, Re)
    #C      FORCED CONVECTION OF A CYLINDER
    #C      ADJUSTING Nu - Re CORRELATION FOR Re NUMBER (P. 260 MCADAMS,1954)
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

# function get_Nusselt_forced(shape::Plate, Re)
#     #C      FORCED CONVECTION OF A PLATE
#     Nu = 0.102 * Re ^ 0.675 * Pr ^ (1 / 3)
#     (Nu)
# end

function get_Nusselt_forced(shape::Ellipsoid, Re)
    #C      FORCED CONVECTION OF A SPHERE
    Nu = 0.35 * Re ^ 0.6 # FROM McAdams, W.H. 1954. Heat Transmission. McGraw-Hill, New York, p.532
    (Nu)
end

function water_prop(T_water)
  # C     β = COEFFICIENT OF EXPANSION (1/C)
  # C     cp_fluid = SPECIFIC HEAT (J/KG-C)
  # C     ρ_air = WATER DENSITY (KG/M3)
  # C     k_fluid = THERMAL CONDUCTIVITY (W/M-C)
  # C     μ = DYNAMIC VISCOSITY (KG/M-S)
  β = 0.21E-031/K
  T_water = Unitful.ustrip(T_water)-273.15
  cp_fluid = 4220.02 - 4.5531 * T_water + 0.182958  * T_water ^ 2 - 0.00310614 *  T_water ^ 3 + 1.89399E-5 * T_water ^ 4
  cp_fluid = (cp_fluid)J/kg/K
  if T_water < 30
    ρ_water = 1000.
  else
    if T_water <= 60
        ρ_water = 1017 - 0.6 * T_water
    else
        T_water = 60
        ρ_water = 1017 - 0.6 * T_water
    end
end
  k_fluid = (0.551666+0.00282144 * T_water - 2.02383E-5 * T_water ^ 2)W/m/K
  μ = (0.0017515-4.31502E-5 * T_water + 3.71431E-7 * T_water ^ 2)kg/m/s
  ρ_water = (ρ_water)kg/m^3
  (β = β, cp_fluid = cp_fluid, ρ_water = ρ_water, k_fluid = k_fluid, μ = μ)
end

function convection(Body, area, T_air, T_surf, vel, P_atmos, elev, fluid)
    shape = Body.geometry
    G = Unitful.gn # acceleration due to gravity, m.s^2
    β = 1 / T_air
    A_conv = area # m2
    D = shape.characteristic_dimension
    dry_air_out = dry_air(T_air, P_atmos, elev)
    dif_vpr = dry_air_out.dif_vpr
    #C     CHECKING TO SEE IF THE FLUID IS WATER, NOT AIR
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
    Pr = cp_fluid * μ / k_fluid
    Pr = Unitful.ustrip(Pr) # is dimensionless, enforce this
    if fluid == 1
        #C      WATER; NO MEANING
        Sc = 1
    else
        Sc = μ / (ρ_air * dif_vpr)
    end
    δ_T = T_surf - T_air
    Gr = abs(((ρ_air^2) * β * G * (D^3) * δ_T) / (μ^2))
    Re = ρ_air * vel * D / μ
    Nu_free = get_Nusselt_free(Body.shape, Gr, Pr)
    Hc_free = (Nu_free * k_fluid) / D
    #C     CALCULATING THE SHERWOOD NUMBER FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.Nu_comb
    Sh_free = Nu_free * (Sc / Pr)^(1 / 3)
    #C     CALCULATING THE MASS TRANSFER COEFFICIENT FROM THE SHERWOOD NUMBER
    Hd_free = Sh_free * dif_vpr / D
    #C     CALCULATING THE CONVECTIVE HEAT LOSS AT THE SKIN
    Q_free = Hc_free * A_conv * (T_surf - T_air)
    Nu = get_Nusselt_forced(Body.shape, Re)
    #C     FORCED CONVECTION FOR ANIMAL
    Hc_forc = Nu * k_fluid / D # HEAT TRANFER COEFFICIENT
    Sh_forc = Nu * (Sc / Pr)^(1 / 3) # SHERWOOD NUMBER
    Hd_forc = Sh_forc * dif_vpr / D # MASS TRANSFER COEFFICIENT
    #C     USING BIRD, STEWART, & LIGHTFOOT'S MIXED CONVECTION FORMULA (P. 445, TRANSPORT PHENOMENA, 2002)
    Nu_comb = (Nu_free^3 + Nu^3)^(1 / 3)
    Hc = Nu_comb * (k_fluid / D) # MIXED CONVECTION HEAT TRANSFER
    Q_conv = Hc * A_conv * (T_surf - T_air) # TOTAL CONVECTIVE TRANSFER
    #C     CALCULATING THE SHERWOOD NUMBERs FROM THE COLBURN ANALOGY
    #C     (BIRD, STEWART & LIGHTFOOT, 1960. TRANSPORT PHENOMENA. WILEY.
    Sh = Nu_comb * (Sc / Pr)^(1 / 3) # SHERWOOD NUMBER
    Hd = Sh * dif_vpr / D # MASS TRANSFER COEFFICIENT
    (Q_conv=Q_conv, Hc=Hc, Hd=Hd, Sh=Sh, Q_free=Q_free, Hc_free=Hc_free, Hc_forc=Hc_forc, Sh_free=Sh_free, Sh_forc=Sh_forc, Hd_free=Hd_free, Hd_forc=Hd_forc)
end

function conduction(A, L, T_org, T_sub, k_sub)
    A * (k_sub / L) * (T_org - T_sub)
end

function solar(α_org_dorsal, α_org_ventral, A_sil, A_up, A_down, F_sub, F_sky, α_sub, Q_dir, Q_dif)
    Q_direct = α_org_dorsal * A_sil * Q_dir
    Q_sol_sky = α_org_dorsal * F_sky * A_up * Q_dif
    Q_sol_sub = α_org_ventral * F_sub * A_down * (1 - α_sub) * (Q_dir + Q_dif)
    (Q_direct + Q_sol_sub + Q_sol_sky)
end


function radin(A_tot = 0.01325006m^2,
    F_sky = 0.4,
    F_sub = 0.4,
    ϵ_org = 0.95,
    ϵ_sub = 0.95,
    ϵ_sky = 0.8,
    T_sky = (10+273.15)K,
    T_sub = (30+273.15)K)

σ = Unitful.k^4*π^2/(60*Unitful.ħ^3*Unitful.c0^2) # Stefan-Boltzmann constant, W/m^2/K^4, make Unitful.σ when error is fixed in Unitful
Q_ir_sky = ϵ_org * F_sky * A_tot * ϵ_sky * σ * T_sky ^ 4
Q_ir_sub = ϵ_org * F_sub * A_tot * ϵ_sub * σ * T_sub ^ 4
(Q_ir_sky + Q_ir_sub)
end

function radout(
    T_skin = (25.1+273.15)K,
    A_tot = 0.01325006m^2,
    F_sky = 0.4,
    F_sub = 0.4,
    ϵ_org = 0.95)
# C     COMPUTES LONGWAVE RADIATION LOST
σ = Unitful.k^4*π^2/(60*Unitful.ħ^3*Unitful.c0^2) # Stefan-Boltzmann constant, W/m^2/K^4, make Unitful.σ when error is fixed in Unitful
Q_ir_to_sky = A_tot * F_sky * ϵ_org * σ  * T_skin ^ 4
Q_ir_to_sub = A_tot * F_sub * ϵ_org * σ  * T_skin ^ 4
(Q_ir_to_sky + Q_ir_to_sub)
end

function evap(
    T_core = (25+273.15)K,
    T_skin = (25.1+273.15)K,
    J_resp = 1.177235e-09kg/s,
    ψ_org = -7.07 * 100J/kg,
    p_wet = 0.1,
    A_tot = 0.01325006m^2,
    Hd = 0.02522706m/s,
    p_eyes = 0.03 / 100,
    T_air = (20+273.15)K,
    rh = 50,
    P_atmos = 101325Pa)
    
  # C     THIS SUBROUTINE COMPUTES SURFACE EVAPORATION BASED ON THE MASS TRANSFER
  # C     COEFFICIENT, % OF SURFACE OF THE SKIN ACTING AS A FREE WATER SURFACE
  # C     AND EXPOSED TO THE AIR, AND THE VAPOR DENSITY GRADIENT BETWEEN THE
  # C     SURFACE AND THE AIR, EACH AT THEIR OWN TEMPERATURE.

  # get vapour density at surface based on water potential of body
  m_w = 0.018kg/mol #! molar mass of water, kg/mol
  RH = exp(ψ_org / (Unitful.R / m_w * T_skin)) * 100 #
  wet_air_out = wet_air(T_skin, 0K, RH, 999K, P_atmos)
  ρ_vap_surf = wet_air_out.ρ_vap

  # get air vapour density
  wet_air_out = wet_air(T_air, 0K, rh, 999K, P_atmos)
  ρ_vap_air = wet_air_out.ρ_vap

  # water lost from eyes if present
  J_eyes = Hd * p_eyes * A_tot * (ρ_vap_surf - ρ_vap_air)
  if J_eyes > 0kg/s
    J_cut = A_tot * p_wet * (1 - p_eyes) * Hd * (ρ_vap_surf - ρ_vap_air)
  else
    J_cut = A_tot * p_wet * Hd * (ρ_vap_surf - ρ_vap_air)
  end
  
  # total water lost
  J_evap = J_eyes + J_resp + J_cut

  # get latent heat of vapourisation and compute heat exchange due to evaporation
  dry_air_out = dry_air(T_air, P_atmos, elev)
  L_v = dry_air_out.L_v
  Q_evap = (J_eyes + J_cut) * L_v
  
  # convert from kg/s to g/s
  J_eyes = uconvert(u"g/s",J_eyes)
  J_resp = uconvert(u"g/s",J_resp)
  J_cut = uconvert(u"g/s",J_cut)
  J_evap = uconvert(u"g/s",J_evap)

  (Q_evap = Q_evap, J_evap = J_evap, J_resp = J_resp, J_cut = J_cut, J_eyes = J_eyes)
end