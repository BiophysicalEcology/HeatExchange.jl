library(NicheMapR)

Ww_g = 40
alpha = 0.85
epsilon = 0.95
rho_body = 1000
fatosk = 0.4
fatosb = 0.4
shape = 2
shape_a = 1
shape_b = 3
shape_c = 2 / 3
conv_enhance = 1
#custom_shape = c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743)
pct_cond = 10
pct_touch = 0
postur = 1
orient = 0
k_flesh = 0.5
M_1 = 0.013
M_2 = 0.8
M_3 = 0.038
M_4 = 0
Q_act = 0
pct_wet = 0.1
pct_eyes = 0
pct_mouth = 0
psi_body = -707
pantmax = 1
F_O2 = 20
RQ = 0.8
delta_air = 0.1
leaf = 0
g_vs_ab = 0.3
g_vs_ad = 0
elev = 0
alpha_sub = 0.2
epsilon_sub = 1
epsilon_sky = 1
pres = 101325
fluid = 0
O2gas = 20.95
CO2gas = 0.03
N2gas = 79.02
K_sub = 0.5
PDIF = 0.1
SHADE = 0
QSOLR = 1000
Z = 20
TA = 20
TGRD = 30
TSUBST = 30
TSKY = -5
VEL = 1
RH = 5

ectoR_input <- c(Ww_g, alpha, epsilon, rho_body, fatosk, fatosb, shape, shape_a, shape_b, shape_c, conv_enhance, pct_cond, pct_touch, postur, orient, k_flesh, M_1, M_2, M_3, M_4, Q_act, pct_wet, pct_eyes, pct_mouth, psi_body, pantmax, F_O2, RQ, delta_air, leaf, g_vs_ab, g_vs_ad, elev, alpha_sub, epsilon_sub, epsilon_sky, pres, fluid, O2gas, CO2gas, N2gas, K_sub, PDIF, SHADE, QSOLR, Z, TA, TGRD, TSUBST, TSKY, VEL, RH)

ectoR.out <- ectoR_devel(
  Ww_g = Ww_g,
  alpha = alpha,
  epsilon = epsilon,
  rho_body = rho_body,
  fatosk = fatosk,
  fatosb = fatosb,
  shape = shape,
  shape_a = shape_a,
  shape_b = shape_b,
  shape_c = shape_c,
  conv_enhance = conv_enhance,
  custom_shape = custom_shape,
  pct_cond = pct_cond,
  pct_touch = pct_touch,
  postur = postur,
  orient = orient,
  k_flesh = k_flesh,
  M_1 = M_1,
  M_2 = M_2,
  M_3 = M_3,
  M_4 = M_4,
  Q_act = Q_act,
  pct_wet = pct_wet,
  pct_eyes = pct_eyes,
  pct_mouth = pct_mouth,
  psi_body = psi_body,
  pantmax = pantmax,
  F_O2 = F_O2,
  RQ = RQ,
  delta_air = delta_air,
  leaf = leaf,
  g_vs_ab = g_vs_ab,
  g_vs_ad = g_vs_ad,
  elev = elev,
  alpha_sub = alpha_sub,
  epsilon_sub = epsilon_sub,
  epsilon_sky = epsilon_sky,
  pres = pres,
  fluid = fluid,
  O2gas = O2gas,
  CO2gas = CO2gas,
  N2gas = N2gas,
  K_sub = K_sub,
  PDIF = PDIF,
  SHADE = SHADE,
  TA = TA,
  # air temperature at lizard height, deg C
  TGRD = TGRD,
  # ground temperature, deg C
  TSKY = TSKY,
  # sky temperature, deg C
  VEL = VEL,
  # wind speed, m/s
  RH = RH,
  # relative humidity, %
  QSOLR = QSOLR,
  # total horizontal plane solar radiation, W/m2
  Z = Z # solar zenith angle, degrees
)

temperatures <- c(ectoR.out$TC, ectoR.out$TSKIN, ectoR.out$TLUNG)
enbal <- unlist(ectoR.out$enbal)
masbal <- unlist(ectoR.out$masbal)
geom <- unlist(ectoR.out$GEOM.out)
solar <- unlist(ectoR.out$SOLAR.out)
resp <- unlist(ectoR.out$RESP.out)
conv <- unlist(ectoR.out$CONV.out)
evap <- unlist(ectoR.out$SEVAP.out)

write.csv(ectoR_input, '../data/ectoR_input.csv')
write.csv(temperatures, '../data/temperatures.csv')
write.csv(unlist(enbal), '../data/enbal.csv')
write.csv(unlist(masbal), '../data/masbal.csv')
write.csv(unlist(geom), '../data/geom.csv')
write.csv(unlist(solar), '../data/solar.csv')
write.csv(unlist(resp), '../data/resp.csv')
write.csv(unlist(conv), '../data/conv.csv')
write.csv(unlist(evap), '../data/evap.csv')
