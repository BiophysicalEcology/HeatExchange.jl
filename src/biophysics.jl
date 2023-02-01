# Biopysics

function conduction(A, L, T_object, T_substrate, k_sub)
    A * (k_sub / L) * (T_object - T_substrate)
end

function solar(alpha_object, A_silhouette, A_up, A_down, alpha_substrate, Q_direct, Q_diffuse, Z)
    Q_norm = Q_direct / cos(Z)
    Q_direct = alpha_object * A_sillhouette * Q_norm
    Q_diffuse = alpha_object * A_up * Q_diffuse + alpha_object * A_down * (Q_direct + Q_diffuse) * (1 - alpha_substrate)
    Q_direct + Q_norm
end

function vapour_pressure(T)
    T = T+273.15
    loge=T
    if T <= 273.15
        loge=-9.09718*(273.16/T-1)-3.56654*log10(273.16/T)+0.876793*(1-T/273.16)+log10(6.1071)
    else
        loge=-7.90298*(373.16/T-1)+5.02808*log10(373.16/T)-1.3816E-07*(10^(11.344*(1-T/373.16))-1)+8.1328E-03*(10^(-3.49149*(373.16/T-1))-1)+log10(1013.246)
    end
    (10^loge)*100
end

function wetair(T_drybulb, T_wetbulb=T_drybulb, relhumid=0, dewpoint=999, P=101325)
    T = T_drybulb + 273.15
    esat = vapour_pressure(T)
    if dewpoint < 999
        e = vapour_pressure(dewpoint)
        relhumid = (e / esat) * 100
    else
        if min(relhumid) > -1
            e = esat * relhumid / 100
        else
            delta_bulb = T_drybulb - T_wetbulb
            wbsat = vapour_pressure(T_wetbulb)
            delta_e = 0.000660 * (1 + 0.00115 * T_wetbulb) * P * delta_bulb
            e = wbsat - delta_e
            relhumid = (e / esat) * 100
        end
    end
    rw = ((0.62197 * 1.0053 * e) / (P - 1.0053 * e))
    vd = e * 0.018016 / (0.998 * 8.31434 * T)
    tvir = T * ((1.0 + rw / (18.016 / 28.966)) / (1 + rw))
    tvinc = tvir - T
    denair = 0.0034838 * P / (0.999 * tvir)
    cp = (1004.84 + (rw * 1846.40)) / (1 + rw)
    if min(relhumid) <= 0
        theta = -999
    else
        theta = 4.615e+5 * T * log(relhumid / 100)
    end
    #const e2 = u"Pa"
    e2 = e
    (e=e2, esat, vd, rw, tvinc, denair, cp, theta, relhumid)
    #return(e, esat, vd, rw, tvinc, denair, cp, psi, relhumid)
end

function dryair(db, bp=101325, elev=0)
  pstd=101325
  T=db+273.15
  patmos=pstd*((1-(0.0065*elev/288))^(1/0.190284))
  densty=bp/(287.04*T)
  visnot=1.8325e-5
  tnot=296.16
  c=120
  visdyn=(visnot*(tnot+c)/(T+c))*(T/tnot)^1.5 # kg / m.s
  viskin=visdyn/densty # m2 / s or J.s/kg
  difvpr=2.26e-5*(((T)/273.15)^1.81)*(1.e5/bp) # m2 / s
  thcond=0.02425+(7.038e-5*db) # W / m.K
  htovpr=2.5012E6-2.3787e3*db # J/kg
  tcoeff=1/T
  ggroup=0.0980616*tcoeff/(viskin*viskin) # 1 / m3.K
  bbemit=5.670367e-8*((T)^4)
  emtmax=2.897e-3/(T)
  (patmos=patmos, densty=densty, visdyn=visdyn, viskin=viskin, difvpr=difvpr, thcond=thcond, htovpr=htovpr, tcoeff=tcoeff, ggroup=ggroup, bbemit=bbemit, emtmax=emtmax)
end
