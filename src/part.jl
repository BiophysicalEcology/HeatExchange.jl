"""
    Part(; body, insulation, traits, environment, ventral_fraction, side = Dorsal(),
        conduction_fraction = 0.0, conductance_coefficient = 0.0u"W/K", covered_area,
        characteristic_dimension, weight = 1.0, neighbours = nothing)

A surface of an organism that comes to its own skin and insulation temperatures around a
core: its geometry `body`, the `side` of the `InsulationParameters` `insulation` it wears,
the `traits` of its tissue and surface, and the `environment` it sees.

`traits` holds `flesh_conductivity`, `fat_conductivity`, `ϵ_body`, `skin_wetness`,
`insulation_wetness`, `bare_skin_fraction` and `eye_fraction`. `environment` holds the
`temperature`s and `view_factors` of the part's surroundings, the `atmos`, `fluid`,
`gas_fractions` and `convection_enhancement`, and the `solar_flow` the part absorbs.

`conduction_fraction` of the part's surface touches the ground through
`conductance_coefficient`. `covered_area` is the surface joined to other parts, which
exchanges nothing with the environment; it defaults to none. The `characteristic_dimension`
for convection defaults to the part's own [`VolumeCubeRoot`](@ref). `weight` is the share of
the organism's heat the part carries. `neighbours`, keyed like the parts of its body, are
the shares of its view that the other parts take.
"""
struct Part{B,I,S,T,E,F,C,V,A,AE,D,W,N}
    body::B
    insulation::I
    side::S
    traits::T
    environment::E
    conduction_fraction::F
    conductance_coefficient::C
    ventral_fraction::V
    area::A
    evaporation_area::AE
    characteristic_dimension::D
    weight::W
    neighbours::N
end
function Part(; body, insulation, traits, environment, ventral_fraction, side = Dorsal(),
    conduction_fraction = 0.0, conductance_coefficient = 0.0u"W/K",
    covered_area = zero(BiophysicalGeometry.total_area(body)),
    characteristic_dimension = characteristic_dimension(VolumeCubeRoot(), body), weight = 1.0,
    neighbours = nothing,
)
    area = BiophysicalGeometry.total_area(body) - covered_area
    # The evaporating area is a share of the skin, and the covered skin takes the same share.
    exposed_evaporation_area = evaporation_area(body) * (1 - covered_area / skin_area(body))
    Part(body, insulation, side, traits, environment, conduction_fraction, conductance_coefficient,
        ventral_fraction, area, exposed_evaporation_area, characteristic_dimension, weight, neighbours)
end

"""
    parts(body::CompositeBody, setups::NamedTuple, views::NamedTuple; shade = 0.0, understory_cover = 0.0)

The [`Part`](@ref)s of `body`, from `setups` keyed like `body.parts`, each holding the
keywords of `Part` other than `body`, `covered_area` and `neighbours`, and an `environment`
without `view_factors`.

`views`, keyed like `body.parts`, splits each part's view into `sky`, `ground` and
`neighbours`, as [`silhouette_factors`](@ref) does. The canopy hides `shade` of the sky and
the understory `understory_cover` of the ground.

Halves joined over their whole cut faces are one body to the air around them, so by
default each takes the [`VolumeCubeRoot`](@ref) of all the halves it is joined to as its
`characteristic_dimension`.
"""
function parts(body::CompositeBody{<:NamedTuple{K}}, setups::NamedTuple, views::NamedTuple;
    shade = 0.0, understory_cover = 0.0,
) where {K}
    covered = BiophysicalGeometry.covered_areas(body.parts, body.joins)
    group = _groups(body, map(_joins_halves, body.joins))
    volumes = map(part -> part.geometry.volume, values(body.parts))
    dimensions = NamedTuple{K}(map(values(body.parts), group) do part, g
        cbrt(sum(map((volume, h) -> h == g ? volume : zero(volume), volumes, group))) +
            _outer_insulation_thickness(part.insulation)
    end)
    map(body.parts, setups, views, covered, dimensions) do part, setup, view, covered_area, characteristic_dimension
        split = _split_view(view.sky, view.ground; shade, understory_cover)
        view_factors = ViewFactors(split.sky, split.ground, split.understory, split.canopy)
        Part(; characteristic_dimension, setup..., body = part, covered_area,
            environment = (; setup.environment..., view_factors), neighbours = view.neighbours)
    end
end

_joins_halves(::Join{<:Any,<:Any,<:Attachment{<:Flat,FullCover},<:Attachment{<:Flat,FullCover}}) = true
_joins_halves(::Join) = false

# The open sky and ground, and the canopy and understory that hide `shade` of the one and
# `understory_cover` of the other.
function _split_view(sky, ground; shade, understory_cover)
    canopy = sky * shade
    understory = ground * understory_cover
    (; sky = sky - canopy, canopy, ground = ground - understory, understory)
end

"""
    heat_balance(temperatures::OrganismTemperatures, metabolic_heat_flow, part::Part;
        respiration_pars, neighbours = (), minimum_metabolic_heat = 0 W,
        k_flesh, pant, skin_wetness, smoothing = HardBound())

The heat budget of `part` at the given temperatures and metabolic heat flow, with the
residuals that vanish at its steady state:

- `residual_energy_balance`: metabolic + solar − respiration − evaporation − convection −
  radiation − conduction (W)
- `residual_internal_conduction`: metabolic − respiration − the heat conducted from core
  to skin (W)
- `surface_balance`: the difference of those two, in which metabolism and respiration
  cancel (W)
- `residual_skin_temperature`: the skin temperature less the one the insulation side of
  the budget implies (K)

`neighbours` are the parts this one sees, as `(; fraction, temperature)` of its view.
`k_flesh`, `pant` and `skin_wetness` override the part's flesh conductivity, the panting
of `respiration_pars` and the part's skin wetness.
"""
function heat_balance(temperatures::OrganismTemperatures, metabolic_heat_flow, part::Part;
    respiration_pars,
    neighbours = (),
    minimum_metabolic_heat = zero(metabolic_heat_flow),
    k_flesh = part.traits.flesh_conductivity,
    pant = respiration_pars.pant,
    skin_wetness = part.traits.skin_wetness,
    smoothing::SmoothingStrategy = HardBound(),
)
    (; core_temperature, skin_temperature, insulation_temperature) = temperatures
    (; body, side, conduction_fraction, conductance_coefficient) = part
    insulation_pars = part.insulation
    (; temperature, view_factors, atmos, fluid, solar_flow, gas_fractions, convection_enhancement) = part.environment
    T = temperature
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; fat_conductivity, ϵ_body, insulation_wetness, bare_skin_fraction, eye_fraction) = part.traits
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    area_convection = part.area * (1 - conduction_fraction)

    insulation = insulation_properties(insulation_pars, insulation_temperature * 0.7 + skin_temperature * 0.3,
        part.ventral_fraction; smoothing)
    (; insulation_conductivity, effective_conductivity) = _insulation_conductivity(
        insulation, insulation_pars, side, insulation_temperature, skin_temperature, σ; smoothing)
    conductivities = ThermalConductivities(k_flesh, fat_conductivity, insulation_conductivity)
    insulation = setproperties(insulation, (; conductivity_compressed = effective_conductivity))

    conv = convection(;
        body,
        area = area_convection,
        air_temperature = T.air,
        surface_temperature = insulation_temperature,
        wind_speed,
        atmospheric_pressure,
        fluid,
        gas_fractions,
        convection_enhancement,
        characteristic_dim = part.characteristic_dimension,
        smoothing,
    )
    atmos = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
    skin_evaporation_heat_flow = evaporation(
        AnimalEvaporationParameters(; skin_wetness, eye_fraction, bare_skin_fraction),
        conv.mass_transfer_coefficient,
        atmos,
        part.evaporation_area,
        skin_temperature,
        T.air;
        gas_fractions,
    ).evaporation_heat_flow

    radiant_out = radiant_temperature(;
        body,
        insulation,
        insulation_pars,
        org_temps = temperatures,
        conductivities,
        side,
        conductance_coefficient,
        conduction_fraction,
        evaporation_flow = skin_evaporation_heat_flow,
        substrate_temperature = T.substrate,
        smoothing,
    )
    (; compressed_insulation_temperature, conductances) = radiant_out
    radiant = radiant_out.radiant_temperature

    insulation_evaporation_heat_flow = _insulation_evaporation(
        conv, atmos, area_convection, insulation_temperature, T.air,
        insulation_wetness, insulation.insulation_test; gas_fractions, smoothing)

    coefficients = _radiation_coefficients(area_convection, view_factors, ϵ_body, σ, radiant, T)
    sky_radiation_flow = coefficients.sky * (radiant - T.sky)
    understory_radiation_flow = coefficients.understory * (radiant - T.understory)
    canopy_radiation_flow = coefficients.canopy * (radiant - T.canopy)
    ground_radiation_flow = coefficients.ground * (radiant - T.ground)
    radiation_heat_flow = _add_neighbour_radiation(
        sky_radiation_flow + understory_radiation_flow + canopy_radiation_flow + ground_radiation_flow,
        area_convection, neighbours, ϵ_body, σ, radiant)

    convection_heat_flow = conv.heat_transfer_coefficient.combined * area_convection * (insulation_temperature - T.air)
    conduction_heat_flow = u"W"(conductance_coefficient * (compressed_insulation_temperature - T.substrate))
    net_metabolic_heat_internal = net_metabolic_heat(; body, conductivities, core_temperature, skin_temperature, smoothing)

    environment_flow = radiation_heat_flow + convection_heat_flow + conduction_heat_flow +
        insulation_evaporation_heat_flow - solar_flow
    skin_temperature_from_heat_balance = mean_skin_temperature(;
        body,
        insulation,
        insulation_pars,
        conductivities,
        conductances,
        conduction_fraction,
        environment_flow,
        skin_evaporation_flow = skin_evaporation_heat_flow,
        core_temperature,
        calculated_insulation_temperature = insulation_temperature,
        compressed_insulation_temperature,
        smoothing,
    ).mean_skin_temperature

    respiration_out = respiration(
        MetabolicRates(; metabolic = metabolic_heat_flow, sum = metabolic_heat_flow, minimum = minimum_metabolic_heat),
        setproperties(respiration_pars, (; pant)),
        atmos,
        mass(body.shape),
        (core_temperature + skin_temperature) / 2,
        T.air;
        gas_fractions,
        O2conversion = Kleiber1961(),
        smoothing,
    )
    respiration_heat_flow = uconvert(u"W", respiration_out.respiration_heat_flow)

    residual_energy_balance = metabolic_heat_flow + solar_flow - respiration_heat_flow -
        skin_evaporation_heat_flow - insulation_evaporation_heat_flow -
        convection_heat_flow - radiation_heat_flow - conduction_heat_flow
    residual_internal_conduction = metabolic_heat_flow - respiration_heat_flow - net_metabolic_heat_internal

    return (;
        solar_heat_flow = solar_flow,
        convection_heat_flow,
        radiation_heat_flow,
        conduction_heat_flow,
        skin_evaporation_heat_flow,
        insulation_evaporation_heat_flow,
        respiration_heat_flow,
        net_metabolic_heat_internal,
        sky_radiation_flow,
        understory_radiation_flow,
        canopy_radiation_flow,
        ground_radiation_flow,
        residual_energy_balance,
        residual_internal_conduction,
        surface_balance = residual_energy_balance - residual_internal_conduction,
        residual_skin_temperature = skin_temperature - skin_temperature_from_heat_balance,
    )
end

"""
    solve_surface(part::Part, core_temperature; skin_temperature, insulation_temperature,
        tolerance, neighbours = (), smoothing = HardBound())

The skin and insulation temperatures of `part` around `core_temperature`, starting from the
given ones, with `neighbours` as for [`heat_balance`](@ref).

Returns `(; skin_temperature, insulation_temperature, net_metabolic, flesh_conductance, flows,
insulation_conductivity, tolerance, success, ntry)`, where `net_metabolic` is the heat
conducted from core to skin and `flesh_conductance` that heat over the temperature drop.
"""
function solve_surface(part::Part, core_temperature; skin_temperature, insulation_temperature, tolerance,
    neighbours = (), smoothing::SmoothingStrategy = HardBound(),
)
    insulation_test = _insulation_test(part.insulation)
    result = if insulation_test > zero(insulation_test)
        _solve_insulated(part, core_temperature, skin_temperature, insulation_temperature, tolerance, neighbours, smoothing)
    else
        _solve_bare(part, core_temperature, skin_temperature, insulation_temperature, tolerance, neighbours, smoothing)
    end
    net_metabolic = result.flows.net_metabolic
    flesh_conductance = _flesh_conductance(net_metabolic, core_temperature, result.skin_temperature)
    return (; result..., net_metabolic, flesh_conductance)
end

# The heat conducted to the skin over the drop in temperature to it, with the drop kept
# from vanishing.
function _flesh_conductance(net_metabolic, core_temperature, skin_temperature)
    gradient = core_temperature - skin_temperature
    floor_gradient = 1e-6u"K"
    return net_metabolic / (abs(gradient) < floor_gradient ? oftype(gradient, floor_gradient) : gradient)
end

# Metabolism and respiration cancel out of the surface balance, so any will do here.
function _solve_insulated(part, core_temperature, skin_temperature, insulation_temperature, tolerance, neighbours, smoothing)
    respiration_pars = stripparams(RespirationParameters())
    surface_residuals(skin, insulation) =
        heat_balance(OrganismTemperatures(core_temperature, skin, insulation), 1.0u"W", part;
            respiration_pars, neighbours, smoothing)
    solution = _newton_surface(surface_residuals, skin_temperature, insulation_temperature)
    b = solution.balance

    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)
    insulation = insulation_properties(part.insulation,
        insulation_temperature * 0.7 + skin_temperature * 0.3, part.ventral_fraction; smoothing)
    (; insulation_conductivity) = _insulation_conductivity(
        insulation, part.insulation, part.side, solution.insulation, solution.skin, σ; smoothing)

    flows = HeatFlows(
        b.convection_heat_flow,
        b.conduction_heat_flow,
        b.net_metabolic_heat_internal,
        b.skin_evaporation_heat_flow,
        b.insulation_evaporation_heat_flow,
        b.radiation_heat_flow,
        b.solar_heat_flow,
        b.sky_radiation_flow,
        b.understory_radiation_flow,
        b.canopy_radiation_flow,
        b.ground_radiation_flow,
    )
    return (;
        skin_temperature = solution.skin,
        insulation_temperature = solution.insulation,
        flows,
        insulation_conductivity,
        tolerance,
        success = solution.success,
        ntry = solution.iters,
    )
end

# Damped Newton on the skin and insulation temperatures, with residuals `surface_balance`
# (W) and `residual_skin_temperature` (K) taken without units. The Jacobian is a finite
# difference. Each step is halved until the residual shrinks and stays finite, and a
# singular Jacobian falls back to a bounded fixed-point step.
function _newton_surface(f, skin0, insulation0; tol = 1.0e-9, maxiter = 200)
    Tunit = oneunit(skin0)
    scale1(x) = x / oneunit(x)
    resnorm(r) = hypot(scale1(r.surface_balance), scale1(r.residual_skin_temperature))
    s = skin0 / Tunit
    n = insulation0 / Tunit
    r = f(s * Tunit, n * Tunit)
    nr = resnorm(r)
    success = false
    iters = 0
    for iter in 1:maxiter
        iters = iter
        if nr < tol
            success = true
            break
        end
        R1 = scale1(r.surface_balance)
        R2 = scale1(r.residual_skin_temperature)
        h = 1.0e-4
        rs = f((s + h) * Tunit, n * Tunit)
        rn = f(s * Tunit, (n + h) * Tunit)
        J11 = (scale1(rs.surface_balance) - R1) / h
        J21 = (scale1(rs.residual_skin_temperature) - R2) / h
        J12 = (scale1(rn.surface_balance) - R1) / h
        J22 = (scale1(rn.residual_skin_temperature) - R2) / h
        det = J11 * J22 - J12 * J21
        if !isfinite(det) || abs(det) < 1.0e-12
            ds = clamp(R2, -1.0, 1.0)
            dn = clamp(R1, -1.0, 1.0)
        else
            ds = (J22 * R1 - J12 * R2) / det
            dn = (-J21 * R1 + J11 * R2) / det
        end
        α = 1.0
        snew, nnew, rnew, nrnew = s, n, r, nr
        for _ in 1:40
            snew = s - α * ds
            nnew = n - α * dn
            rnew = f(snew * Tunit, nnew * Tunit)
            nrnew = resnorm(rnew)
            (isfinite(nrnew) && nrnew < nr) && break
            α *= 0.5
        end
        s, n, r, nr = snew, nnew, rnew, nrnew
    end
    return (; skin = s * Tunit, insulation = n * Tunit, balance = r, success, iters)
end

# Without insulation the skin is the surface, and its temperature is a fixed point of the
# linearised balance. All of the skin is bare.
function _solve_bare(part, core_temperature, skin_temperature, insulation_temperature, tolerance, neighbours, smoothing)
    (; body) = part
    (; temperature, view_factors, atmos, fluid, solar_flow, gas_fractions, convection_enhancement) = part.environment
    T = temperature
    (; relative_humidity, wind_speed, atmospheric_pressure) = atmos
    (; flesh_conductivity, ϵ_body, skin_wetness, eye_fraction) = part.traits
    σ = Unitful.uconvert(u"W/m^2/K^4", Unitful.σ)

    ntry = 0
    volume = flesh_volume(body)
    area = part.area
    r_skin = skin_radius(body)
    evaporation_pars = AnimalEvaporationParameters(; skin_wetness, eye_fraction, bare_skin_fraction = 1.0)
    atmos = AtmosphericConditions(relative_humidity, wind_speed, atmospheric_pressure)
    flesh = 4 * flesh_conductivity * volume / r_skin^2

    while true
        ntry += 1
        for i in 1:20
            conv = convection(;
                body,
                area,
                air_temperature = T.air,
                surface_temperature = skin_temperature,
                wind_speed,
                atmospheric_pressure,
                fluid,
                gas_fractions,
                convection_enhancement,
                characteristic_dim = part.characteristic_dimension,
                smoothing,
            )
            heat_transfer_coefficient = conv.heat_transfer_coefficient.combined
            skin_evaporation_flow = evaporation(
                evaporation_pars, conv.mass_transfer_coefficient, atmos, area, skin_temperature, T.air; gas_fractions,
            ).evaporation_heat_flow

            c = _radiation_coefficients(area, view_factors, ϵ_body, σ, skin_temperature, T)
            nb = _neighbour_coefficients(area, neighbours, ϵ_body, σ, skin_temperature)
            gain = flesh * core_temperature - skin_evaporation_flow + heat_transfer_coefficient * area * T.air + solar_flow +
                c.sky * T.sky + c.understory * T.understory + c.canopy * T.canopy + c.ground * T.ground + nb.cTsum
            conductance = flesh + heat_transfer_coefficient * area + c.sky + c.understory + c.canopy + c.ground + nb.csum
            skin_temperature_calc = gain / conductance

            sky_radiation_flow = c.sky * (skin_temperature_calc - T.sky)
            understory_radiation_flow = c.understory * (skin_temperature_calc - T.understory)
            canopy_radiation_flow = c.canopy * (skin_temperature_calc - T.canopy)
            ground_radiation_flow = c.ground * (skin_temperature_calc - T.ground)
            longwave_flow = sky_radiation_flow + understory_radiation_flow + canopy_radiation_flow + ground_radiation_flow +
                nb.csum * skin_temperature_calc - nb.cTsum
            convection_flow = heat_transfer_coefficient * area * (skin_temperature_calc - T.air)
            net_metabolic = flesh * (core_temperature - skin_temperature_calc)
            flows = HeatFlows(
                convection_flow,
                0.0u"W",
                net_metabolic,
                skin_evaporation_flow,
                0.0u"W",
                longwave_flow,
                solar_flow,
                sky_radiation_flow,
                understory_radiation_flow,
                canopy_radiation_flow,
                ground_radiation_flow,
            )

            if abs(skin_temperature - skin_temperature_calc) < tolerance
                return (; skin_temperature = skin_temperature_calc, insulation_temperature, flows,
                    insulation_conductivity = nothing, tolerance, success = true, ntry)
            end
            skin_temperature = skin_temperature_calc
            insulation_temperature = skin_temperature_calc
            ntry += 1
            if ntry == 101
                if tolerance <= 0.001u"K"
                    tolerance = 0.01u"K"
                    ntry = 0
                else
                    return (; skin_temperature = skin_temperature_calc, insulation_temperature, flows,
                        insulation_conductivity = nothing, tolerance, success = false, ntry)
                end
            end
        end
    end
end
