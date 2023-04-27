# using Microclim
using Plots
using Roots

using Unitful, 
      NCDatasets, 
      FieldMetadata, 
      FieldDefaults,
      ModelParameters

using Unitful: W, m, °C, hr, mol, K, s, J, Mg, g, kg, L, kPa, Pa, °, R

import FieldMetadata: default, units, limits, @units, @limits

include("/Users/usuario/Documents/GitHub/Microclim.jl//src/interpolation.jl")
include("/Users/usuario/Documents/GitHub/Microclim.jl/src/types.jl")
include("/Users/usuario/Documents/GitHub/Microclim.jl/src/netcdf.jl")
include("/Users/usuario/Documents/GitHub/Microclim.jl/src/files.jl")

include("../src/geometry.jl")
include("../src/organism.jl")
include("../src/environment.jl")
include("../src/biophysics.jl")
include("../src/heat_balance.jl")
include("../src/simulation.jl")

# basepath = "c:/Spatial_Data/microclimOz"
basepath = "/Volumes/igel/spatial/microclimOz"

years = 2000:2000
shade = 0
envgrid = load_grid(basepath, years, shade)
t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))

# plot(t1.airtemperature)
# plot(t1.soiltemperature)

# define the geometry
mass = 0.04kg
ρ_body = 1000kg/m^3
shapeb = 3
shapec = 2 / 3
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
geometric_traits = Body(shape_body, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

# construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
lizard = Model(Organism(geometric_traits, FunctionalTraits()))

# @time simulation(lizard, t1)

pred_tbs = simulation(lizard, t1)
plot(pred_tbs)

