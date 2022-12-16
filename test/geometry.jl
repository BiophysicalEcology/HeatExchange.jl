

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
