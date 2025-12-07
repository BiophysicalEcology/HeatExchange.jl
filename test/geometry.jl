using HeatExchange
using Unitful
using Test
using ModelParameters
#using Plots

density = 1000.0u"kg/m^3"
mass = 2000.0u"kg"
shape_b = 2.0
shape_c = 1.0
fur_thickness = 10.0u"mm"
fibre_diameter = 30.0u"μm"
fibre_density = 3000u"cm^-2"
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"

fur = Fur(fur_thickness, fibre_diameter, fibre_density)
fat = Fat(fat_fraction, fat_density)

# plate
shape = Plate(mass, density, shape_b, shape_c)
body = Body(shape, Naked())

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fur)

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fat)

total_area(body)
skin_area(body)
evaporation_area(body)

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, CompositeInsulation(fur, fat))

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

# cylinder
shape = Cylinder(mass, density, shape_b)
body = Body(shape, Naked())

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fur)

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fat)

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, CompositeInsulation(fur, fat))

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

# sphere
shape = Sphere(mass, density)
body = Body(shape, Naked())

total_area(body)
skin_area(body)
evaporation_area(body)
silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fur)

total_area(body)
skin_area(body)
evaporation_area(body)
silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fat)

total_area(body)
skin_area(body)
evaporation_area(body)
silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, CompositeInsulation(fur, fat))

total_area(body)
skin_area(body)
evaporation_area(body)
silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

# ellipsoid
shape = Ellipsoid(mass, density, shape_b, shape_c)
body = Body(shape, Naked())

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fur)

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, fat)

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

body = Body(shape, CompositeInsulation(fur, fat))

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(body, 0.0u"°")
silhouette_area(body, NormalToSun())
silhouette_area(body, ParallelToSun())
silhouette_area(body, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

# leopard frog
shape = LeopardFrog(mass, density)
body = Body(shape, Naked())

total_area(body)
skin_area(body)
evaporation_area(body)

silhouette_area(shape, 45.0u"°")

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

# desert iguana
shape = DesertIguana(mass, density)
body = Body(shape, Naked())

total_area(body)
skin_area(body)
evaporation_area(body)
silhouette_area(shape, NormalToSun())
silhouette_area(shape, ParallelToSun())
silhouette_area(shape, Intermediate())

skin_radius(body)
insulation_radius(body)
flesh_radius(body)

body.geometry.area
body.geometry.length
body.geometry.volume

# trunkmass = 0.04u"kg"
# trunkshapeb = 2
# trunkshape = Cylinder(trunkmass, density, trunkshapeb) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
# trunk = Body(trunkshape, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct
# θ = 90u"°"
# trunksilhouette = silhouette_area(trunk, θ)

# zenith_angles = (0:90)u"°"
# sils = silhouette_area.(Ref(trunk), zenith_angles)
# plot(zenith_angles, sils, xlabel = "zenith angle", ylabel = "silhouette area")

# headmass = 0.5u"kg"
# headshapeb = 2
# headshapec = 2 / 3
# headshape = Ellipsoid(headmass, density, headshapeb, headshapec)
# head = Body(headshape, Naked())
# headsilhouette = silhouette_area(head, θ)

# zenith_angles = (0:90)u"°"
# sils = silhouette_area.(Ref(head), zenith_angles)
# plot(zenith_angles, sils, xlabel = "zenith angle", ylabel = "silhouette area")

# # frog
# density = 1000u"kg/m^3"
# mass = 0.04u"kg"

# frogshape = LeopardFrog(mass, density) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
# frog = Body(frogshape, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

# mass_range = (1:100)u"g"
# area = []
# for m in mass_range
#     imass = Unitful.uconvert(u"kg", m)
#     push!(area, geometry(LeopardFrog(imass, density), Naked()).area)
# end

# # plot(mass_range, Unitful.uconvert.(u"cm^2",area), xaxis=:log, yaxis=:log)
# # xlims!(1e+0, 1e+2)
# # ylims!(1e+0, 1e+2)

# zenith_angles = (0:90)u"°"
# sils = silhouette_area.(Ref(frog), zenith_angles)
# # plot(zenith_angles, sils, xlabel = "zenith angle", ylabel = "silhouette area")

