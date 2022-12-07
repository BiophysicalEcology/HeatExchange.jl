using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using Plots
using UnitfulRecipes

"""
    Shape

Abstract type for the shape of the object being modelled.
"""
abstract type Shape end

abstract type Insulation end

abstract type AbstractBody end

struct Naked <: Insulation end

struct Fur{T} <: Insulation
    thickness::T
end

struct Geometry{V,C,L,A}
    volume::V
    characteristic_dimension::C
    lengths::L
    area::A
end

"""
    Body <: AbstractBody

Physical dimensions of a body or body part that may or may note be insulated.
"""
struct Body{S<:Shape,I<:Insulation,G} <: AbstractBody
    shape::S
    insulation::I
    geometry::G
end
# Add a Body constructor from shape, insulation, mass and density
function Body(shape::Shape, insulation::Insulation)
    Body(shape, insulation, geometry(shape, insulation))
end

# struct BodyContext{B<:Body,S,X}
#     body::B
#     silhouette:S
#     #other_time_location_dependent_thing::X
# end
# function BodyContext(body::Body, environment)
#     silhouette = calc_silhouette(body, environment)
#     #other_thing = calc_thing(body, environment)
#     BodyContext(body, silhouette)
# end


"""
    Cylinder <: Shape

A cylindrical organism.
"""
struct Cylinder{M,D,B} <: Shape
    mass::M
    density::D
    b::B
end

function geometry(shape::Cylinder, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    r = (volume / (shape.b * π * 2))^(1 / 3)
    length1 = shape.b * r * 2
    length2 = 2 * r
    area = area_of_cylinder(r, length1)
    return Geometry(volume, length, (length1, length2), area)
end

# function geometry(shape::Cylinder, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area, sil_area)
# end

"""
    Ellipsoid <: Shape

An ellipsoidal organism.
"""
struct Ellipsoid{M,D,B,C} <: Shape
    mass::M
    density::D
    b::B
    c::C
end

function geometry(shape::Ellipsoid, ::Naked)
    volume = shape.mass / shape.density
    length = volume^(1 / 3)
    a = ((3 / 4)* volume / (π * shape.b * shape.c)) ^ (1 / 3)
    b = a * shape.b
    c = a * shape.c
    p = 1.6075
    length1 = a * 2
    length2 = b * 2
    length3 = c * 2
    area = area_of_ellipsoid(a, b, c)
    return Geometry(volume, length, (length1, length2, length3), area)
end

function area_of_ellipsoid(a, b, c)
    e = ((a ^ 2 - c ^ 2) ^ 0.5 ) / a # eccentricity
    2 * π * b ^ 2 + 2 * π * (a * b / e) * asin(e)
end

"""
    sil_area_of_ellipsoid

Calculates the silhouette (projected) area of a prolate spheroid.
"""
function sil_area_of_ellipsoid(a, b, c, θ)
    a2 = cos(90) ^ 2 * (cos(θ) ^ 2 / a ^ 2 + sin(θ) ^ 2 / b ^ 2) + sin(90) ^ 2 / c ^ 2
    twohh = 2 * cos(90) * sin(90) * cos(θ) * (1 / b ^ 2 - 1 / a ^ 2)
    b2 = sin(θ) ^ 2 / a ^ 2 + cos(θ) ^ 2 / b ^ 2
    θ2 = 0.5 * atan(twohh, a2 - b2)
    sps = sin(θ2)
    cps = cos(θ2)
    a3 = cps * (a2 * cps + twohh * sps) + b2 * sps * sps
    b3 = sps * (a2 * sps - twohh * cps) + b2 * cps * cps
    semax1 = 1 / sqrt(a3)
    semax2 = 1 / sqrt(b3)
    π * semax1  * semax2
end

function area_of_cylinder(r, l)
    2 * π * r * l + 2 * π * r^2
end

"""
    sil_area_of_cylinder

Calculates the silhouette (projected) area of a cylinder.
Equation from Fig. 11.6 in Campbell, G. S., & Norman, J. M.
(1998). Environmental Biophysics. Springer.
"""
function sil_area_of_cylinder(r, l, θ)
    2 * r * l * sin(θ) + π * r^2 * cos(θ)
end
# function geometry(shape::Ellipsoid, fur::Fur)
#     # Something different with fur
#     # ...
#     return Geometry(volume, length, (length1, length2, length3), area)
# end


density = 1000#kg/m^3

trunkmass = 1#kg
trunkshapeb = 2
trunkshape = Cylinder(trunkmass, density, trunkshapeb)
trunk = Body(trunkshape, Naked())
θ = 90°
trunksilhouette = sil_area_of_cylinder(trunk.geometry.lengths[2]/2, trunk.geometry.lengths[1], θ)

zens = (0:90)°
sils = sil_area_of_cylinder.(trunk.geometry.lengths[2]/2, trunk.geometry.lengths[1], zens)
plot(zens, sils, xlabel = "zenith angle", ylabel = "silhouette area")

headmass = 0.5#kg
headshapeb = 2
headshapec = 2 / 3
headshape = Ellipsoid(headmass, density, headshapeb, headshapec)
head = Body(headshape, Naked())
headsilhouette = sil_area_of_ellipsoid(head.geometry.lengths[1], head.geometry.lengths[2], head.geometry.lengths[3], θ)

zens = (0:90)°
sils = sil_area_of_ellipsoid.(head.geometry.lengths[1], head.geometry.lengths[2], head.geometry.lengths[3], zens)
plot(zens, sils, xlabel = "zenith angle", ylabel = "silhouette area")

function conduction(A, L, T_object, T_substrate, k_sub)
    A * (k_sub / L) * (T_object - T_substrate)
end

function solar(alpha_object, A_silhouette, A_up, A_down, alpha_substrate, Q_direct, Q_diffuse, Z)
    Q_norm = Q_direct / cos(Z)
    Q_direct = alpha_object * A_sillhouette * Q_norm
    Q_diffuse = alpha_object * A_up * Q_diffuse + alpha_organism * A_down * (Q_direct + Q_diffuse) * (1 - alpha_substrate)
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

wetair = function(T_drybulb, T_wetbulb=T_drybulb, relhumid=0, dewpoint=999, P=101325)
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

wetair(20)

dryair = function(db, bp=101325, elev=0)
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
dryair(T)
