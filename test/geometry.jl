

density = 1000kg/m^3

trunkmass = 0.04kg
trunkshapeb = 2
trunkshape = Cylinder(trunkmass, density, trunkshapeb) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
trunk = Body(trunkshape, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct
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

# frog
density = 1000kg/m^3
mass = 0.04kg

frogshape = Frog(mass, density) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
frog = Body(frogshape, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

mass_range = (1:100)u"g"
areas = []
for m in mass_range
    imass = Unitful.uconvert(u"kg", m)
    areas = push!(areas, geometry(Frog(imass, density), Naked()).area)
end

plot(mass_range, Unitful.uconvert.(u"cm^2",areas), xaxis=:log, yaxis=:log)
xlims!(1e+0, 1e+2)
ylims!(1e+0, 1e+2)

zens = (0:90)°
sils = sil_area_of_frog.(frog.geometry.area, zens)
plot(zens, sils, xlabel = "zenith angle", ylabel = "silhouette area")

