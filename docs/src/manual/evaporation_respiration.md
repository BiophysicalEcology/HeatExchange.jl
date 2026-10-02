# Evaporation and respiration

Water that evaporates from an organism takes heat with it, about 2.4 kJ for each gram. This couples the heat
budget to the water budget (Kearney and Porter 2020), see [Flows of mass](gradients.md#Flows-of-mass). Water
evaporates from the skin, the eyes, wet fur, the stomata of a leaf, and the lungs.

```@setup evaporation
using Main.FigureHelpers
using CairoMakie
```

## From the surface

The rate at which water leaves a wet surface is set by the difference in vapour density between the air at the
surface and the air beyond it, the area that is wet, and the mass transfer coefficient ``h_d`` from
[`convection`](@ref):

```math
\dot{m} = h_d \, A_{wet} \, (\rho_{v,s} - \rho_{v,a}), \qquad Q_{evap} = \lambda \, \dot{m}
```

where ``\lambda`` is the latent heat of vaporisation. The air at the surface is saturated at the surface
temperature or, for an organism not fully hydrated, at the humidity in equilibrium with its water potential
``\psi``:

```math
h_s = \exp\left(\frac{\psi \, M_w}{R \, T_s}\right)
```

This is 1 at a water potential of zero, and 0.995 at the −707 J/kg of the default ectotherm. Vapour densities
are from [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl).

The wet area is described by [`AnimalEvaporationParameters`](@ref), after Tracy (1976) and Kearney and Porter
(2004):

| Field | Meaning |
|:--|:--|
| `skin_wetness` | the fraction of the skin that acts as a free water surface. About 0.001 for a desert lizard, near 1 for a frog, and rising with sweating in a mammal |
| `eye_fraction` | the fraction of the surface that is open eye, always wet |
| `bare_skin_fraction` | the fraction of the wet skin that is bare. Bare skin evaporates by free and forced convection. Skin under fur is sheltered from wind and evaporates by free convection only |
| `insulation_wetness` | the fraction of the outer surface of the fur that is wet, from rain or licking |

```@example evaporation
using HeatExchange, BiophysicalGeometry, Unitful
import HeatExchange: AtmosphericConditions

body = Body(Cylinder(40.0u"g", 1000.0u"kg/m^3", 6.0), Naked())
air_temperature, skin_temperature = u"K"(25.0u"°C"), u"K"(30.0u"°C")
air = convection(; body, area = total_area(body), air_temperature, surface_temperature = skin_temperature,
                   wind_speed = 1.0u"m/s", atmospheric_pressure = 101325.0u"Pa", fluid = Air())
atmosphere = AtmosphericConditions(0.3, 1.0u"m/s", 101325.0u"Pa")   # relative humidity, wind speed, pressure

lizard_skin = AnimalEvaporationParameters(; skin_wetness = 0.001, eye_fraction = 0.0003, bare_skin_fraction = 1.0)
frog_skin = AnimalEvaporationParameters(; skin_wetness = 1.0, eye_fraction = 0.0003, bare_skin_fraction = 1.0)
dry, wet = map((lizard_skin, frog_skin)) do skin
    evaporation(skin, air.mass_transfer_coefficient, atmosphere, total_area(body), skin_temperature, air_temperature)
end
dry.evaporation_heat_flow, wet.evaporation_heat_flow
```

A wet-skinned animal the size of the lizard in [Get started](../get_started.md) would lose more heat by
evaporation than that lizard gained from the sun. The water lost is returned too, from the skin and the eyes:

```@example evaporation
uconvert(u"g/hr", wet.cutaneous_mass_flow), uconvert(u"g/hr", wet.eyes_mass_flow)
```

## From a leaf

A leaf loses water through its stomata, and how far they are open is a conductance, not a wetted fraction.
[`LeafEvaporationParameters`](@ref) holds the vapour conductance of the lower (abaxial) and upper (adaxial)
surface, and a cuticular conductance that remains when the stomata are closed, in mol m⁻² s⁻¹. Each surface is
in series with the boundary layer, the same ``h_d`` from [`convection`](@ref), and the two are in parallel:

```math
h_{leaf} = \frac{1}{2} \frac{h_{ab} \, h_d}{h_{ab} + h_d} + \frac{1}{2} \frac{h_{ad} \, h_d}{h_{ad} + h_d}
```

with each stomatal conductance converted to a velocity by ``h = g \, R T / P``. The method of
[`evaporation`](@ref) is chosen by the type of the parameters, so a leaf is an [`Organism`](@ref) with
[`LeafEvaporationParameters`](@ref), and the rest of the heat budget is shared with animals. See
[A leaf](../tutorials/leaf.md).

## From the lungs

Air is breathed in at the temperature and humidity of the surroundings and out warm and saturated. ``Q_{resp}``
is the heat used to evaporate the water added, less the heat given up by the air if it leaves cooler than it
came. [`respiration`](@ref) computes it by a balance of moles through the lungs:

1. The metabolic rate is converted to a rate of oxygen consumption, see [Metabolism](metabolism.md).
2. The air needed to supply that oxygen follows from the fraction of oxygen in the air and the fraction the
   lungs extract, `oxygen_extraction_efficiency`. Panting multiplies the air flow by `pant`.
3. The moles of oxygen, carbon dioxide, nitrogen and water in and out follow, with carbon dioxide produced in
   the ratio `respiratory_quotient` to the oxygen consumed.
4. The water evaporated is the vapour leaving, saturated to `exhaled_relative_humidity` at lung temperature,
   less that entering.

```@example evaporation
import HeatExchange: MetabolicRates

metabolic_heat_flow = metabolic_rate(Kleiber(), 65.0u"kg")
breath = respiration(MetabolicRates(; metabolic = metabolic_heat_flow), example_respiration_pars(),
                     AtmosphericConditions(0.3, 1.0u"m/s", 101325.0u"Pa"), 65.0u"kg",
                     u"K"(35.0u"°C"), u"K"(20.0u"°C"))
breath.respiration_heat_flow, breath.respiration_mass_flow
```

The parameters are in [`RespirationParameters`](@ref), and the composition of the air in `gas_fractions` of
[`EnvironmentalPars`](@ref), which can be changed for a burrow. The molar flows are returned as
[`MolarFluxes`](@ref HeatExchange.MolarFluxes) in `molar_fluxes_in` and `molar_fluxes_out`:

```@example evaporation
breath.molar_fluxes_in.oxygen, breath.molar_fluxes_out.oxygen
```

The flows of gas are computed from demand, not from gradients of partial pressure, see
[Flows of mass](gradients.md#Flows-of-mass).

The lung temperature is between that of the core and the skin. For a bare body it comes from
[`surface_and_lung_temperature`](@ref), and for an insulated body it is the mean of core and skin temperature.

!!! note "Temperature of exhaled air"
    Air is exhaled at the lung temperature in this version. NicheMapR exhales it at the air temperature plus an
    offset, `DELTAR` or `delta_air`, where that is lower, for the recovery of heat and water in the nasal
    passages. The parameter exists here, `exhaled_temperature_offset`, but is not yet used. For an endotherm in
    the cold this gives a higher respiratory heat and water loss than NicheMapR with its default `DELTAR = 0`,
    see [An endotherm: metabolic rate](../tutorials/endotherm.md).

### The respiration balance

When metabolic rate is the unknown, respiration and metabolism depend on each other: the heat generated must
cover the heat conducted to the skin and the heat lost in the breath, which depends on the heat generated.
[`respiration`](@ref) returns the residual,

```math
\mathrm{balance} = Q_{gen} - Q_{resp}(Q_{gen}) - Q_{gen,net}
```

where ``Q_{gen,net}`` is passed in as `sum` in `MetabolicRates`. [`solve_metabolic_rate`](@ref) finds the
``Q_{gen}`` at which it is zero, see [Solving a heat balance](heat_balance.md#Closing-the-budget). The metabolic
rate used for breathing is not allowed below `minimum`, so an animal under a heat load still breathes at its
basal rate.

## With panting and sweating

Evaporation is the only way to lose heat to an environment hotter than the body. Panting raises `pant`, and
sweating or licking raises `skin_wetness`. In a bare-skinned animal, panting also adds `mouth_fraction` to the
skin wetness, for the open mouth. The amounts are decided by
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), see
[Endotherm thermoregulation by rules](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/endotherm_rules).
This package gives the heat and water that result.
