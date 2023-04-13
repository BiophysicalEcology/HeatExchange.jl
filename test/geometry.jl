using HeatExchange
using Unitful
using Test
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
# using Plots

density = 1000kg/m^3

trunkmass = 0.04kg
trunkshapeb = 2
trunkshape = Cylinder(trunkmass, density, trunkshapeb) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
trunk = Body(trunkshape, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct
θ = 90°
trunksilhouette = calc_silhouette_area(trunk, θ)

zens = (0:90)°
sils = calc_silhouette_area.(Ref(trunk), zens)
# plot(zens, sils, xlabel = "zenith angle", ylabel = "silhouette area")

headmass = 0.5#kg
headshapeb = 2
headshapec = 2 / 3
headshape = Ellipsoid(headmass, density, headshapeb, headshapec)
head = Body(headshape, Naked())
headsilhouette = calc_silhouette_area(head, θ)

zens = (0:90)°
sils = calc_silhouette_area.(Ref(head), zens)
# plot(zens, sils, xlabel = "zenith angle", ylabel = "silhouette area")

# frog
density = 1000kg/m^3
mass = 0.04kg

frogshape = LeopardFrog(mass, density) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
frog = Body(frogshape, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

mass_range = (1:100)u"g"
areas = []
for m in mass_range
    imass = Unitful.uconvert(u"kg", m)
    push!(areas, geometry(LeopardFrog(imass, density), Naked()).area)
end

# plot(mass_range, Unitful.uconvert.(u"cm^2",areas), xaxis=:log, yaxis=:log)
# xlims!(1e+0, 1e+2)
# ylims!(1e+0, 1e+2)

zens = (0:90)°
sils = calc_silhouette_area.(Ref(frog), zens)
# plot(zens, sils, xlabel = "zenith angle", ylabel = "silhouette area")

