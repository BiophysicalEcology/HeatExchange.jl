# For NicheMapR users

HeatExchange.jl descends from the heat budget code of [NicheMapR](https://github.com/mrke/NicheMapR): the
ectotherm model (Kearney and Porter 2020), the endotherm model (Kearney et al. 2021) and its human extension,
HomoTherm (Kearney et al. 2026). The physics is the same. The structure is not, and this page says what moved
where.

## Two programs become one question with two answers

NicheMapR is organised by kind of organism. The ectotherm model finds the body temperature at which the heat
budget balances, for bare skin and small metabolic heat. The endotherm model holds the core temperature and
finds the metabolic rate, for an animal that may have fur and fat. Each has its own Fortran library and R
interface.

Kearney et al. (2021) noted that this split reflects a real biological difference, the far greater energetic
intensity of endotherms, but it is not a physical one. Both programs solve the same steady-state heat budget.
They differ in two things only:

- **which quantity is unknown**: the core temperature, or the rate of metabolic heat production;
- **what the surface is**: bare skin, or skin under fur or feathers.

HeatExchange.jl is organised by those two things, and any combination can be solved:

| | Bare skin | Insulated |
|:--|:--|:--|
| **Solve for temperature** | the NicheMapR ectotherm model, and leaves | *new*: a torpid or dead endotherm, a furred insect, a chick before it can thermoregulate |
| **Solve for metabolic rate** | *new*: a naked endotherm, a thermogenic flower | the NicheMapR endotherm model |

The first row is [`solve_temperature`](@ref) and the second [`solve_metabolic_rate`](@ref). The column is
decided by the body: `Naked()` takes the bare-skin path and a fibrous layer the insulated path, see
[Temperature or metabolic rate](solvers.md). The heat-flow functions underneath are shared, so there is one
`convection`, one `evaporation` and one `respiration` where NicheMapR has two of each.

## What moved where

| NicheMapR | Here |
|:--|:--|
| `ectotherm` / `ectoR_devel`: `FUN` solved for `TC` by `ZBRENT` (`uniroot` in R) | [`heat_balance`](@ref) solved by [`solve_temperature`](@ref) with [`zbrent`](@ref) |
| `SOLAR`, `RADIN`, `RADOUT`, `CONV`, `COND`, `SEVAP`, `RESP`, `MET` and their `_ecto` and `_ENDO` versions | [`solar`](@ref), [`radiation_in`](@ref), [`radiation_out`](@ref), [`convection`](@ref), [`conduction`](@ref), [`evaporation`](@ref), [`respiration`](@ref), [`metabolic_rate`](@ref) |
| `endoR` / `endoR_devel` / `SOLVENDO` | [`solve_metabolic_rate`](@ref) |
| `IRPROP`, `GETKFUR` | [`insulation_properties`](@ref), [`insulation_thermal_conductivity`](@ref) |
| `SIMULSOL`, called for the dorsal and the ventral side | [`solve_temperatures`](@ref) for each side, on the residuals of [`solve_part_heat_balance`](@ref) |
| `ZBRENT_ENDO` on `RESPFUN` | [`zbrent`](@ref) on the `balance` returned by [`respiration`](@ref) |
| the closed-form flesh, fat and fur conduction of each shape | a stack of radial layers, see [Layers as a radial graph](radial_layers.md) |
| the allometric metabolic rate equations in `MET` and `QBASAL` | [`MetabolicRateEquation`](@ref) types here for now, moving to [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl), see [Metabolism](metabolism.md) |
| `ellipsoid_endo` | [`ellipsoid_endotherm`](@ref) |
| `leaf = 1` in `ectotherm` | [`LeafEvaporationParameters`](@ref) in place of [`AnimalEvaporationParameters`](@ref) |
| `HomoTherm`: `endoR` once per body part, summed | [`solve_coupled_metabolic_rate`](@ref) over the parts of a body, see [Bodies of many parts](multipart.md) |
| `GEOM`, `GEOM_ENDO`, the geometric shape codes | [BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl): `Cylinder`, `Sphere`, `Plate`, `Ellipsoid` |
| the empirical surface areas in `GEOM` and `GEOM_ENDO`: the desert iguana and leopard frog, and the bird and mammal options of `SAMODE` | [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl), which already has the bird and mammal surface areas (`surface_area(EutherianMammal(), mass)`). These are allometric relations, not geometry. `DesertIguana` and `LeopardFrog` are still in BiophysicalGeometry.jl in this version, and are to move |
| `DRYAIR`, `WETAIR`, `VAPPRS`, `WATER` | [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl) |
| the thermoregulation loops of `ectotherm` (shade, posture, burrow, climb) and of `endoR` (uncurl, vasodilate, raise core temperature, pant, sweat) | [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), see [For NicheMapR users](https://biophysicalecology.github.io/BiophysicalBehaviour.jl/dev/manual/nichemapr) in its documentation |
| transient heat budgets: `onelump`, `onelump_var`, `twolump`, `trans_behav`, and the transient option of `ectotherm` | to be added to this package: the residual of [`heat_balance`](@ref) tracked through time, in a varying environment, and turned into body temperature by the heat capacity. `trans_behav` and the behaviour of the transient `ectotherm` are for BiophysicalBehaviour.jl |
| the Dynamic Energy Budget model of `ectotherm` | not in this package: to come through AnimalMapper.jl, by way of DEBtool_J.jl, see [Flows of mass](gradients.md#Flows-of-mass) |
| `micro_global` and the other microclimate functions | [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) and [SolarRadiation.jl](https://github.com/BiophysicalEcology/SolarRadiation.jl) |

The last four rows are the largest change in scope. `ectotherm` is a whole simulation: it reads a year of
microclimate, chooses where the animal is each hour, and can run an energy budget alongside. `endoR` includes a
thermoregulatory sequence. This package does neither. It solves the heat budget of one organism in one
environment, and is the layer that the behaviour and life-cycle packages call, see
[Environments and the ecosystem](ecosystem.md). It corresponds most closely to `ectoR_devel`, and to
`endoR_devel` with `THERMOREG = 0`.

## There is no separate transient program

The transient functions of NicheMapR are separate solutions for a lump of tissue with heat storage, written
apart from the steady-state code. Here a transient calculation needs no new physics. [`heat_balance`](@ref)
returns the net heat flow at a given temperature, which is the rate of heat storage. A steady state is where
that is zero. A transient is its integral through time, divided by the heat capacity. The transient calculation
of this package is to be exactly that, in an environment that may vary, as in `onelump_var`. It is not yet in
the version documented here, see [Solving a heat balance](heat_balance.md#Steady-state-and-storage).

## Dorsal and ventral

`endoR` calls `SIMULSOL` twice, for a whole animal in dorsal fur facing the sky and a whole animal in ventral
fur facing the ground, and takes a weighted mean of the two heat flows. This is kept, as the
[`MultiSided`](@ref) strategy, and reproduces `endoR`.

It is also now a special case. A body can be two real halves joined at a shared core, each with its own coat,
area and view, see [Back and belly: two halves](../tutorials/two_parts.md).

## Layers

`endoR` has one layer of flesh, an optional layer of fat set by `FATPCT`, and one of fur, with the heat flow
written out for each shape. Here the core-to-skin conduction is an ordered list of layers, each with a
resistance that depends on its shape, added in series. One flesh layer and one fat layer reproduce `endoR` to
machine precision. See [Layers as a radial graph](radial_layers.md).

## Many parts

`HomoTherm` builds an animal from several parts: run the endotherm model for each part without respiration, add
up the heat each must be supplied with, and close the respiration balance once. Each part gives up a fixed
fraction of its area, `PJOINs`, to its joins. That scheme is [`solve_coupled_metabolic_rate`](@ref). Here the
parts are those of a `CompositeBody`, so the hidden areas and the view each part has of the sky, the ground and
the other parts come from where the parts are. Parts can also hold different core temperatures and conduct heat
to each other. [A human of many parts](../tutorials/human.md) rebuilds the `HomoTherm` human.

The body is built with BiophysicalGeometry.jl: each part a `Body`, joined at named surfaces by a `Join` that
says where on each part the join is and how large. Its documentation has tutorials that build a dog, a cow and
this human, and an interactive page,
[Build an animal](https://biophysicalecology.github.io/BiophysicalGeometry.jl/dev/builder), where an animal is
assembled with sliders and the code can be copied.

## Parameter names

The parameters are grouped into structs by process, with names in words, see [Parameters](parameters.md).

### Ectotherm

| `ectoR_devel` | Here | In |
|:--|:--|:--|
| `Ww_g`, `rho_body`, `shape`, `shape_b`, `shape_c` | mass, density, shape type and axis ratios | the shape, from BiophysicalGeometry.jl |
| `alpha` | `body_absorptivity_dorsal`, `body_absorptivity_ventral` | [`RadiationParameters`](@ref) |
| `epsilon` | `body_emissivity_dorsal`, `body_emissivity_ventral` | [`RadiationParameters`](@ref) |
| `fatosk`, `fatosb` | `sky_view_factor`, `ground_view_factor` | [`RadiationParameters`](@ref) |
| `postur` | `solar_orientation`: `Intermediate()`, `NormalToSun()`, `ParallelToSun()` | [`RadiationParameters`](@ref) |
| `pct_cond` (%) | `conduction_fraction` (0 to 1) | [`ExternalConductionParameters`](@ref) |
| `k_flesh` | `flesh_conductivity` | [`InternalConductionParameters`](@ref) |
| `pct_wet`, `pct_eyes` (%) | `skin_wetness`, `eye_fraction` (0 to 1) | [`AnimalEvaporationParameters`](@ref) |
| `pct_mouth` (%) | `mouth_fraction` (0 to 1) | [`RespirationParameters`](@ref) |
| `psi_body` | `water_potential` | [`HydraulicParameters`](@ref) |
| `F_O2` (%), `RQ`, `pantmax`, `delta_air` | `oxygen_extraction_efficiency`, `respiratory_quotient`, `pant`, `exhaled_temperature_offset` | [`RespirationParameters`](@ref) |
| `M_1`, `M_2`, `M_3`, `M_4` | `mass_normalisation`, `mass_exponent`, `thermal_sensitivity`, `metabolic_state` | [`AndrewsPough2`](@ref) |
| `g_vs_ab`, `g_vs_ad` | `abaxial_vapour_conductance`, `adaxial_vapour_conductance` | [`LeafEvaporationParameters`](@ref) |
| `1 - alpha_sub`, `epsilon_sub`, `epsilon_sky`, `elev`, `fluid`, `conv_enhance` | `ground_albedo`, `ground_emissivity`, `sky_emissivity`, `elevation`, `fluid`, `convection_enhancement` | [`EnvironmentalPars`](@ref) |
| `O2gas`, `CO2gas`, `N2gas` (%) | `GasFractions(oxygen, carbon_dioxide, nitrogen)` | [`EnvironmentalPars`](@ref) |
| `TA`, `TSKY`, `TGRD`, `TSUBST`, `VEL`, `RH` (%), `pres` | `air_temperature`, `sky_temperature`, `ground_temperature`, `substrate_temperature`, `wind_speed`, `relative_humidity` (0 to 1), `atmospheric_pressure` | [`EnvironmentalVars`](@ref) |
| `QSOLR`, `Z`, `PDIF`, `SHADE` (%), `K_sub` | `global_radiation`, `zenith_angle`, `diffuse_fraction`, `shade` (0 to 1), `substrate_conductivity` | [`EnvironmentalVars`](@ref) |

### Endotherm

| `endoR_devel` | Here | In |
|:--|:--|:--|
| `AMASS`, `ANDENS`, `SHAPE`, `SHAPE_B`, `SHAPE_C` | mass, density, shape type and axis ratios | the shape |
| `FATPCT` (%), `FATDEN` | `FatLayer(fraction, density)`; `fat_fraction`, `fat_density` | the body; [`InternalConductionParameters`](@ref) |
| `AK1`, `AK2` | `flesh_conductivity`, `fat_conductivity` | [`InternalConductionParameters`](@ref) |
| `DHAIRD`, `LHAIRD`, `ZFURD`, `RHOD`, `REFLD`, `KHAIR` (and `…V`) | `diameter`, `length`, `depth`, `density`, `reflectance`, `conductivity` of `dorsal` (and `ventral`) | [`FibreProperties`](@ref) in [`InsulationParameters`](@ref) |
| `ZFURCOMP`, `XR` | `depth_compressed`, `longwave_depth_fraction` | [`InsulationParameters`](@ref) |
| `PVEN`, `PCOND` | `ventral_fraction`, `conduction_fraction` | [`RadiationParameters`](@ref), [`ExternalConductionParameters`](@ref) |
| `EMISAN`, `FSKREF`, `FGDREF`, `FABUSH`, `ORIENT` | `body_emissivity_…`, `sky_view_factor`, `ground_view_factor`, `bush_view_factor`, `solar_orientation` | [`RadiationParameters`](@ref) |
| `1 - REFLD`, `1 - REFLV` | `body_absorptivity_dorsal`, `body_absorptivity_ventral` | [`RadiationParameters`](@ref) |
| `PCTWET`, `FURWET`, `PCTEYES`, `PCTBAREVAP` (%) | `skin_wetness`, `insulation_wetness`, `eye_fraction`, `bare_skin_fraction` (0 to 1) | [`AnimalEvaporationParameters`](@ref) |
| `EXTREF` (%), `RQ`, `PANT`, `DELTAR`, `RELXIT` (%) | `oxygen_extraction_efficiency`, `respiratory_quotient`, `pant`, `exhaled_temperature_offset`, `exhaled_relative_humidity` | [`RespirationParameters`](@ref) |
| `TC`, `QBASAL`, `Q10` | `core_temperature`, `metabolic_heat_flow`, `q10` | [`MetabolismParameters`](@ref) |
| `RESPIRE`, `DIFTOL`, `BRENTOL` | `respire`, `temperature_error_tolerance`, `resp_tolerance` | [`SolveMetabolicRateOptions`](@ref) |
| `TS`, `TFA` | the starting `skin_temperature` and `insulation_temperature` | arguments of [`solve_metabolic_rate`](@ref) |
| `1 - ABSSB` | `ground_albedo` | [`EnvironmentalPars`](@ref) |
| `TAREF`, `TBUSH`, `TCONDSB`, `KSUB` | `reference_air_temperature`, `bush_temperature`, `substrate_temperature`, `substrate_conductivity` | [`EnvironmentalVars`](@ref) |
| `UNCURL`, `SHAPE_B_MAX`, `AK1_INC`, `AK1_MAX`, `TC_INC`, `TC_MAX`, `PANT_INC`, `PANT_MAX`, `PCTWET_INC`, `PCTWET_MAX`, `TREGMODE` | limits and steps of thermoregulation | BiophysicalBehaviour.jl |

Percentages in NicheMapR are fractions here, and temperatures in °C are
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities in any unit of temperature.

### Outputs

| NicheMapR table | Here |
|:--|:--|
| `enbal` | `energy_balance` from [`solve_temperature`](@ref); `energy_flows` from [`solve_metabolic_rate`](@ref) |
| `treg` | `thermoregulation` |
| `masbal` | `mass_balance`; `mass_flows` |
| `morph` | `morphology`, and the body itself |

## How closely the numbers agree

The tests compare the package with output saved from NicheMapR, and the tutorials show the comparisons, see
[An ectotherm: body temperature](../tutorials/ectotherm.md) and
[An endotherm: metabolic rate](../tutorials/endotherm.md).

**Ectotherm.** Body temperature and every heat and mass flow agree with `ectoR_devel` to one part in ten
thousand, and the geometry to rounding.

**Endotherm.** Geometry agrees with `endoR_devel` to rounding, and core and lung temperatures to 0.1 %. Skin
and fur temperatures, fur conductivity and surface heat flows agree to about 0.3 %, and metabolic rate to about
1 %. The cause is known. `SIMULSOL` updates the fur surface temperature with a linearised surface balance and
stops when successive guesses agree within `DIFTOL`. Here the two surface temperatures are found by a Newton
iteration that balances the surface energy budget to 10⁻⁹. Metabolic rate depends on the small difference
between core and skin temperatures, so a small change in skin temperature is a larger relative change in
metabolic rate.
