# The endotherm, piece by piece

[`solve_metabolic_rate`](@ref) is a short function that calls others, and each can be called alone. Kearney et
al. (2021) presented the endotherm model of NicheMapR this way, as modules to be run alone or combined
differently for a particular animal, because endotherms differ too much in their responses for one arrangement
to suit all. This tutorial walks through the pieces in the order they are called, and puts them back together
by hand to get the answer of the [previous tutorial](endotherm.md). Each step names the subroutine of NicheMapR
it corresponds to.

```@setup components
using Main.FigureHelpers
using CairoMakie
```

## The animal and its surroundings

The default animal, at 0 °C:

```@example components
using HeatExchange, BiophysicalGeometry, Unitful
import HeatExchange: AtmosphericConditions, EnvironmentTemperatures, ViewFactors, GeometryVariables, MetabolicRates,
    Absorptivities, SolarConditions

shape = Ellipsoid(65.0u"kg", 1000.0u"kg/m^3", 1.1, 1.1)
insulation_pars = example_insulation_pars()
core_temperature = u"K"(37.0u"°C")
basal = metabolic_rate(Kleiber(), 65.0u"kg")

environment_vars = example_environment_vars(; air_temperature = u"K"(0.0u"°C"))
environment_pars = example_environment_pars()
air_temperature = environment_vars.air_temperature
nothing # hide
```

## Fur: `IRPROP`

First, the properties of the coat, which depend on the temperature of the air within it, taken as a weighted
mean of the guesses of the fur surface and skin temperatures:

```@example components
skin_guess, surface_guess = u"K"(34.0u"°C"), air_temperature
insulation = insulation_properties(insulation_pars, 0.7 * surface_guess + 0.3 * skin_guess, 0.5)
insulation.conductivities.dorsal, insulation.absorption_coefficients.dorsal
```

See [Insulation](../manual/insulation.md) for how the conductivity varies with the fibres.

## Geometry: `GEOM_ENDO`

The body is built from the shape, the fur and the fat, and gives the areas and lengths:

```@example components
fibres = insulation.fibres.dorsal
body = Body(shape, CompositeInsulation(FibrousLayer(fibres.depth, fibres.diameter, fibres.density), FatLayer(0.0, 901.0u"kg/m^3")))
total_area(body), skin_area(body), evaporation_area(body), characteristic_dimension(VolumeCubeRoot(), body)
```

## Sunlight: `SOLAR_ENDO`

There is no sun in this example, so the call returns zeros, but it is made the same way when there is, see
[Radiation](../manual/radiation.md):

```@example components
absorbed = solar(body, Absorptivities(example_radiation_pars(), environment_pars), ViewFactors(0.5, 0.5, 0.0, 0.0),
                 SolarConditions(environment_vars), silhouette(body, Intermediate()), 0.0u"m^2")
absorbed.solar_flow
```

## Convection: `CONV_ENDO`

Convection is from the outer surface of the fur. At the guessed surface temperature, the air temperature, there
is none yet, but the coefficients are found, see
[Convection and conduction](../manual/convection_conduction.md):

```@example components
air = convection(; body, area = total_area(body), air_temperature, surface_temperature = surface_guess + 10.0u"K",
                   wind_speed = environment_vars.wind_speed, atmospheric_pressure = environment_vars.atmospheric_pressure,
                   fluid = environment_pars.fluid)
air.heat_transfer_coefficient.combined, air.mass_transfer_coefficient.combined
```

## Evaporation from the skin: `SEVAP_ENDO`

The skin under fur is sheltered, and loses water by free convection only, see
[Evaporation and respiration](../manual/evaporation_respiration.md):

```@example components
atmosphere = AtmosphericConditions(environment_vars)
skin = AnimalEvaporationParameters(; skin_wetness = 0.005, eye_fraction = 0.0, bare_skin_fraction = 0.0)
evaporation(skin, air.mass_transfer_coefficient, atmosphere, evaporation_area(body), skin_guess, air_temperature).evaporation_heat_flow
```

## Skin and fur temperatures: `SIMULSOL`

The core of the model is the simultaneous solution for the skin and fur surface temperatures of one side of
the body, see [Solving a heat balance](../manual/heat_balance.md#With-insulation).
[`solve_temperatures`](@ref) takes the body, the coat, the surroundings of that side and the traits. The
convection and evaporation above are computed again inside it, at each trial pair of temperatures.

The dorsal side sees the sky and the ventral side the ground. Each is solved as a whole animal with that coat
and view, so the view factor of each is doubled:

```@example components
traits = (; core_temperature, flesh_conductivity = 0.9u"W/m/K", fat_conductivity = 0.23u"W/m/K", ϵ_body = 0.99,
            skin_wetness = 0.005, insulation_wetness = 0.0, bare_skin_fraction = 0.0, eye_fraction = 0.0)

function side(body_side, view_factors)
    solve_temperatures(;
        body, insulation_pars, insulation,
        geometry_vars = GeometryVariables(; side = body_side, conductance_coefficient = 0.0u"W/K", ventral_fraction = 0.5,
                                            conduction_fraction = 0.0, longwave_depth_fraction = 1.0),
        environment_vars = (; temperature = EnvironmentTemperatures(environment_vars), view_factors, atmos = atmosphere,
                              fluid = environment_pars.fluid, solar_flow = absorbed.solar_flow,
                              gas_fractions = environment_pars.gas_fractions,
                              convection_enhancement = environment_pars.convection_enhancement),
        traits, temperature_tolerance = 1e-3u"K", skin_temperature = skin_guess, insulation_temperature = surface_guess)
end

dorsal = side(Dorsal(), ViewFactors(1.0, 0.0, 0.0, 0.0))     # sky, ground, bush, vegetation
ventral = side(Ventral(), ViewFactors(0.0, 1.0, 0.0, 0.0))
u"°C"(dorsal.skin_temperature), u"°C"(dorsal.insulation_temperature), dorsal.flows.net_metabolic
```

`flows` is a [`HeatFlows`](@ref HeatExchange.HeatFlows) with each heat flow of that side. `net_metabolic` is
``Q_{gen,net}``, the heat that must be conducted from core to skin to hold those temperatures:

```@example components
f = dorsal.flows
(; f.convection, f.longwave, f.skin_evaporation, f.net_metabolic)
```

Convection, longwave radiation and evaporation from the skin add up to it.

## The two sides together

The heat required is the mean of the two sides, weighted by the view of the sky and of the ground. With both at
air temperature and the same coat on each, the two sides are equal here:

```@example components
net_metabolic = 0.5 * dorsal.flows.net_metabolic + 0.5 * ventral.flows.net_metabolic
skin_temperature = (dorsal.skin_temperature + ventral.skin_temperature) / 2
lung_temperature = (core_temperature + skin_temperature) / 2
net_metabolic, u"°C"(lung_temperature)
```

## Respiration: `ZBRENT_ENDO` and `RESPFUN`

The animal must generate that heat and the heat it loses in breathing, which depends on how much it generates.
[`respiration`](@ref) returns the heat lost for a trial metabolic rate, and the residual `balance`. The
metabolic rate is where the balance is zero:

```@example components
resp_pars = example_respiration_pars()
breath(metabolic) = respiration(MetabolicRates(; metabolic, sum = net_metabolic, minimum = basal), resp_pars, atmosphere,
                                65.0u"kg", lung_temperature, air_temperature; O2conversion = Kleiber1961())

trial = 60.0:5.0:160.0   # W
fig, ax = figure_axis("Trial metabolic rate (W)", "Respiration balance (W)")
lines!(ax, trial, [ustrip(u"W", breath(q * u"W").balance) for q in trial]; linewidth = 2)
hlines!(ax, [0.0]; color = :black, linewidth = 1)
fig
```

```@example components
metabolic_heat_flow = zbrent(q -> ustrip(u"W", breath(q * u"W").balance), -100.0, 1000.0, 1e-6) * u"W"
```

Below the basal rate the animal still breathes at its basal rate, which is the bend in the line.

## The answer

```@example components
direct = solve_metabolic_rate(Organism(body, example_heat_exchange_traits()),
                              (; environment_pars, environment_vars), skin_guess, surface_guess)
metabolic_heat_flow, direct.energy_flows.metabolic_heat_flow
```

The pieces assembled by hand give the answer of the solver. The breath carries the water lost with it:

```@example components
final = breath(metabolic_heat_flow)
final.respiration_heat_flow, final.respiration_mass_flow
```

## Arranging the pieces differently

This is the default arrangement. Others are made of the same pieces:

- **One side only**, for a part of a body that faces one way: [`solve_part_surface`](@ref).
- **Several parts**, each solved with [`solve_part_surface`](@ref) and added before the one respiration
  balance: [`solve_coupled_metabolic_rate`](@ref), see [A human of many parts](human.md).
- **A core temperature that is not given**, with the same surface solve inside a search for it:
  [`solve_temperature`](@ref).
- **No search at all**, with the temperatures and metabolic rate handed in and the residuals handed back, for
  an optimiser: [`solve_part_heat_balance`](@ref), see
  [Differentiability and the NLP interface](../manual/autodiff.md).
- **A loop around the whole**, changing a trait between solves, for thermoregulation: see
  [Endotherm thermoregulation by rules](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/endotherm_rules)
  in the documentation of BiophysicalBehaviour.jl.
