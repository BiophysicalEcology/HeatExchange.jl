# Evaporation and respiration

Water that evaporates from an organism takes heat with it, about 2.4 kJ for each gram. This couples the heat
budget to the water budget: the evaporation term of one is the loss term of the other (Kearney and Porter 2020).
Water evaporates from the skin, from the eyes, from wet fur, from the stomata of a leaf, and from the lungs.

```@setup evaporation
using Main.FigureHelpers
using CairoMakie
```

## From the surface

The rate at which water leaves a wet surface is set by the difference in vapour density between the air at the
surface and the air beyond it, by the area that is wet, and by the mass transfer coefficient ``h_d`` from
[`convection`](@ref):

```math
\dot{m} = h_d \, A_{wet} \, (\rho_{v,s} - \rho_{v,a}), \qquad Q_{evap} = \lambda \, \dot{m}
```

where ``\lambda`` is the latent heat of vaporisation. The air at the surface is taken to be saturated at the
temperature of the surface, or, for an organism that is not fully hydrated, at the humidity in equilibrium with its
water potential ``\psi``:

```math
h_s = \exp\left(\frac{\psi \, M_w}{R \, T_s}\right)
```

This is 1 at a water potential of zero, and 0.995 at the −707 J/kg of the default ectotherm. Vapour densities
are from [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl).

The wet area is described by [`AnimalEvaporationParameters`](@ref), after Tracy (1976) and Kearney and Porter
(2004):

| Field | Meaning |
|:--|:--|
| `skin_wetness` | the fraction of the skin that acts as a free water surface. It is about 0.001 for a desert lizard, near 1 for a frog, and rises with sweating in a mammal |
| `eye_fraction` | the fraction of the surface that is open eye, always wet |
| `bare_skin_fraction` | the fraction of the wet skin that is bare. Bare skin evaporates by free and forced convection. Skin under fur is sheltered from the wind and evaporates by free convection only |
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

A wet-skinned animal of the size of the lizard in [Get started](../get_started.md) would lose more heat by
evaporation than that lizard gained from the sun. The water lost is returned as well, from the skin and from the
eyes:

```@example evaporation
uconvert(u"g/hr", wet.cutaneous_mass_flow), uconvert(u"g/hr", wet.eyes_mass_flow)
```

## From a leaf

A leaf loses water through its stomata, and how far they are open is a conductance, not a wetted fraction.
[`LeafEvaporationParameters`](@ref) holds the vapour conductance of the lower (abaxial) and upper (adaxial)
surface, and a cuticular conductance that remains when the stomata are closed, in mol m⁻² s⁻¹. Each surface is in
series with the boundary layer of the leaf, which is the same mass transfer coefficient from [`convection`](@ref),
and the two surfaces are in parallel:

```math
h_{leaf} = \frac{1}{2} \frac{h_{ab} \, h_d}{h_{ab} + h_d} + \frac{1}{2} \frac{h_{ad} \, h_d}{h_{ad} + h_d}
```

with each stomatal conductance converted to a velocity by ``h = g \, R T / P``. The method of
[`evaporation`](@ref) is chosen by the type of the parameters, so a leaf is an [`Organism`](@ref) with
[`LeafEvaporationParameters`](@ref) as its `evaporation_pars`, and everything else in the heat budget is shared
with animals. See the tutorial [A leaf](../tutorials/leaf.md).

## From the lungs

Air is breathed in at the temperature and humidity of the surroundings and out warm and saturated. The heat lost,
``Q_{resp}``, is the heat used to evaporate the water added, less the heat given up by the air if it leaves
cooler than it came in. [`respiration`](@ref) computes it by a balance of moles through the lungs:

1. The metabolic rate is converted to a rate of oxygen consumption, see [Metabolism](metabolism.md).
2. The air that must be breathed to supply that oxygen follows from the fraction of oxygen in the air and the
   fraction of it that the lungs extract, `oxygen_extraction_efficiency`. Panting multiplies the air flow by
   `pant`.
3. The moles of oxygen, carbon dioxide, nitrogen and water in and out follow, with carbon dioxide produced in the
   ratio `respiratory_quotient` to the oxygen consumed.
4. The water evaporated is the difference between the water vapour leaving, saturated to
   `exhaled_relative_humidity` at the lung temperature, and that entering.

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

The lung temperature is between that of the core and the skin. For a bare body it comes from
[`surface_and_lung_temperature`](@ref), and for an insulated body it is the mean of the core and skin temperatures.

!!! note "Temperature of exhaled air"
    Air is exhaled at the lung temperature in this version. NicheMapR exhales it at the air temperature plus an
    offset, `DELTAR` or `delta_air`, where that is lower, to represent the recovery of heat and water in the nasal
    passages. The parameter for that offset exists here, `exhaled_temperature_offset`, but is not yet used. For an
    endotherm in the cold this gives a higher respiratory heat and water loss than NicheMapR does with its default
    `DELTAR = 0`, see the tutorial [An endotherm: metabolic rate](../tutorials/endotherm.md).

### The respiration balance

When the metabolic rate is the unknown, respiration and metabolism depend on each other: the heat generated must
cover the heat conducted to the skin and the heat lost in the breath, and the heat lost in the breath depends on
the heat generated. [`respiration`](@ref) returns the residual of this,

```math
\mathrm{balance} = Q_{gen} - Q_{resp}(Q_{gen}) - Q_{gen,net}
```

where ``Q_{gen,net}`` is passed in as `sum` in `MetabolicRates`. [`solve_metabolic_rate`](@ref) finds the ``Q_{gen}``
at which it is zero, see [Solving a heat balance](heat_balance.md). The metabolic rate used for breathing is not
allowed to fall below `minimum`, so that an animal under a heat load still breathes at its basal rate.

## With panting and sweating

Evaporation is the only way to lose heat to an environment that is hotter than the body. Panting raises `pant`,
and sweating or licking raises `skin_wetness`. In a bare-skinned animal, panting also adds `mouth_fraction` to the
skin wetness, for the open mouth. The amounts are decided by
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), and this package gives
the heat and water that result.
