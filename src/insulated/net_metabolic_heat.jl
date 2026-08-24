"""
    net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature)

Calculate net metabolic heat conducted core→skin to maintain the core-to-skin gradient.

Thin wrapper over the radial-layer conduction stack (`radial_net_metabolic_heat`,
radial_layers.jl): the flesh core plus fat shell as an ordered series of shape-dispatched
resistances. See that file for the per-shape resistances and their derivation.

# Keywords
- `body::AbstractBody`: Body geometry
- `conductivities::ThermalConductivities`: Thermal conductivities (flesh, fat, insulation)
- `core_temperature`: Core body temperature
- `skin_temperature`: Skin temperature

# Returns
- `net_metabolic_heat_production`: Net metabolic heat generation (W)
"""
net_metabolic_heat(; body::AbstractBody, conductivities::ThermalConductivities,
                     core_temperature, skin_temperature, smoothing::SmoothingStrategy=HardBound()) =
    radial_net_metabolic_heat(body, conductivities, core_temperature, skin_temperature; smoothing)
