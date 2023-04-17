using Microclim
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

# basepath = "c:/Spatial_Data/microclimOz"
basepath = "/Volumes/igel/spatial/microclimOz"

years = 2000:2000
shade = 0
envgrid = load_grid(basepath, years, shade)
t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))

plot(t1.airtemperature)
plot(t1.soiltemperature)

pred_tbs = Float32[]u"°C"
for i in eachindex(t1.radiation)
    # get the environmental parameters
    env_params = EnvironmentalParams()
    # get the variables for the environment
    env_vars = EnvironmentalVars(t1.airtemperature[i,1],
                                t1.skytemperature[i],
                                t1.soiltemperature[i,1],
                                t1.relhumidity[i,1],
                                t1.windspeed[i,1],
                                101325.0Pa,
                                t1.zenith[i],
                                0.5W/m/K,
                                t1.radiation[i],
                                t1.radiation[i] * 0.96,
                                t1.radiation[i] / 10)
    # define the geometry
    mass = 0.04kg
    ρ_body = 1000kg/m^3
    shapeb = 3
    shapec = 2 / 3
    shape_body = Ellipsoid(mass, ρ_body, shapeb, shapec) # define trunkshape as a Cylinder struct of type 'Shape' and give it required values
    geometric_traits = Body(shape_body, Naked()) # construct a Body, which is naked - this constructor will apply the 'geometry' function to the inputs and return a struct that has the struct for the 'Shape' type, as well as the insulation and the geometry struct

    # construct the Model which holds the parameters of the organism in the Organism concrete struct, of type AbstractOrganism
    lizard = Model(Organism(geometric_traits, FunctionalTraits()))

    # get the variables for the organism
    if i == 1 
        org_vars = OrganismalVars()
    else
        org_vars = OrganismalVars(heat_balance_out.T_core,
                                  heat_balance_out.T_surf,
                                  heat_balance_out.T_lung,
                                  -707J/kg)
    end
    
    variables = (organism = org_vars, environment = env_vars)
    
    T_air = env_vars.T_air
    try
        T_core_s = find_zero(t -> heat_balance(t, lizard, env_params, variables), (T_air - 40K, T_air + 100K), Bisection())
    catch e
        continue
    end
    T_core_C = (Unitful.ustrip(T_core_s) - 273.15)°C
    pred_tbs = push!(pred_tbs, T_core_C)
    heat_balance_out = heat_balance(T_core_s, lizard, env_params, variables)
end

plot(pred_tbs)

