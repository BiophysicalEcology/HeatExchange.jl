#=using Pkg=#
#=Pkg.add(path="/home/urtzi/Documents/GitHub/Microclim.jl");=#
using Microclim
using HeatExchange
using Plots
using Roots
using Test

using Unitful,
    NCDatasets,
    FieldMetadata,
    FieldDefaults,
    ModelParameters

using Unitful: W, m, °C, hr, mol, K, s, J, Mg, g, kg, L, kPa, Pa, °, R

import FieldMetadata: default, units, limits, @units, @limits

include("./src/simulation.jl")

# basepath = "c:/Spatial_Data/microclimOz"
basepath = "/media/urtzi/igel/spatial/microclimOz"

years = 2000:2000
shade = 0
envgrid = load_grid(basepath, years, shade)
t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))

#=plot(t1.airtemperature)=#
#=plot(t1.soiltemperature)=#

# define the geometry
mass = 0.04kg
ρ_body = 1000kg / m^3
shapeb = 3
shapec = 2 / 3
shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec)
geometric_traits = Body(shape_body, Naked())

# construct the Model 
lizard = Model(Organism(geometric_traits, MorphoPars(), PhysioPars()))

# @time simulation(lizard, t1)

pred_tbs = simulation(lizard, t1)
@test pred_tbs[1] == 13.420177854854956°C

plot(pred_tbs)

