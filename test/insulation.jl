
insulation = InsulationParameters(;
    dorsal=FibreProperties(;
        diameter=Param(30.0u"μm"), # hair diameter, dorsal (m)
        length=Param(23.9u"mm"), # hair length, dorsal (m)
        density=Param(3000u"cm^-2"), # hair density, dorsal (1/m2)
        depth=Param(9.0u"mm"), # fur depth, dorsal (m)
        reflectance=Param(0.301), # fur reflectivity dorsal (fractional, 0-1)
        conductivity=Param(0.209u"W/m/K"), # hair thermal conductivity (W/m/K)
    ),
    ventral=FibreProperties(;
        diameter=Param(30.0u"μm"), # hair diameter, ventral (m)
        length=Param(23.9u"mm"), # hair length, ventral (m)
        density=Param(3000u"cm^-2"), # hair density, ventral (1/m2)
        depth=(9.0u"mm"), # fur depth, ventral (m)
        reflectance=Param(0.301), # fur reflectivity ventral (fractional, 0-1)
        conductivity=Param(0.209u"W/m/K"), # hair thermal conductivity (W/m/K)
    ),
    depth_compressed=Param(9.0u"mm"), # depth of compressed fur (for conduction) (m)
    longwave_depth_fraction=Param(1.0),
)

insulation_out = insulation_properties(insulation, (273.15 + 20.0)u"K", 0.3)
density = 1000.0u"kg/m^3"
mass = 0.04u"kg"
shapeb = 2
shape = Cylinder(mass, density, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(
    insulation_out.fibres.average.depth,
    insulation_out.fibres.average.diameter,
    insulation_out.fibres.average.density,
)
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))

body_geometry = geometry(shape, fur)

θ = 90u"°"
silhouette = silhouette_area(body, θ)

shape = Sphere(mass, density) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(
    insulation_out.fibres.average.depth,
    insulation_out.fibres.average.diameter,
    insulation_out.fibres.average.density,
)
composite_insulation = (fat, fur)

body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))

shape = Plate(mass, density, shapeb, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(
    insulation_out.fibres.average.depth,
    insulation_out.fibres.average.diameter,
    insulation_out.fibres.average.density,
)
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))

shape = Ellipsoid(mass, density, shapeb, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(
    insulation_out.fibres.average.depth,
    insulation_out.fibres.average.diameter,
    insulation_out.fibres.average.density,
)
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))
