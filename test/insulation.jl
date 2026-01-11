
insulation = InsulationParameters(;
    fibre_diameter_dorsal=Param(30.0u"μm"), # hair diameter, dorsal (m)
    fibre_diameter_ventral=Param(30.0u"μm"), # hair diameter, ventral (m)
    fibre_length_dorsal=Param(23.9u"mm"), # hair length, dorsal (m)
    fibre_length_ventral=Param(23.9u"mm"), # hair length, ventral (m)
    insulation_depth_dorsal=Param(9.0u"mm"), # fur depth, dorsal (m)
    insulation_depth_ventral=(9.0u"mm"), # fur depth, ventral (m)
    fibre_density_dorsal=Param(3000u"cm^-2"), # hair density, dorsal (1/m2)
    fibre_density_ventral=Param(3000u"cm^-2"), # hair density, ventral (1/m2)
    insulation_reflectance_dorsal=Param(0.301),  # fur reflectivity dorsal (fractional, 0-1)
    insulation_reflectance_ventral=Param(0.301),  # fur reflectivity ventral (fractional, 0-1)
    insulation_depth_compressed=Param(9.0u"mm"), # depth of compressed fur (for conduction) (m)
    fibre_conductivity=Param(0.209u"W/m/K"), # hair thermal conductivity (W/m°C)
    longwave_depth_fraction=Param(1.0),
)

insulation_out = insulation_properties(; insulation, insulation_temperature=(273.15 + 20.0)u"K", ventral_fraction=0.3)
insulation_properties(; insulation=insulation_pars,
            insulation_temperature=thermoregulation.T_insulation,
            ventral_fraction=radiation_pars.ventral_fraction)
density = 1000.0u"kg/m^3"
mass = 0.04u"kg"
shapeb = 2
shape = Cylinder(mass, density, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
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
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))


shape = Plate(mass, density, shapeb, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))

shape = Ellipsoid(mass, density, shapeb, shapeb) # define shape as a Cylinder struct of type 'Shape' and give it required values
fat_fraction = 0.2
fat_density = 901.0u"kg/m^3"
fat = Fat(fat_fraction, fat_density)
fur = Fur(insulation_out.insulation_depths[1], insulation_out.fibre_diameters[1], insulation_out.fibre_densities[1])
composite_insulation = (fat, fur)

body = Body(shape, Naked()) # construct a Body, which is furred
body = Body(shape, fur) # construct a Body, which is furred
body = Body(shape, fat) # construct a Body, which has a fat layer
body = Body(shape, CompositeInsulation(fur, fat))
