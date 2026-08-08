using HeatExchange
using BiophysicalGeometry
using Unitful
using Test

const ρ_rl = 1000.0u"kg/m^3"

# Fur (only affects skin→insulation, not the core→skin stack under test) + genuine fat
# so that flesh_radius < skin_radius and the fat shell term is non-trivial.
_ins_pars = example_insulation_pars()
_fibre    = _ins_pars.dorsal
_fur      = FibrousLayer(_fibre.depth, _fibre.diameter, _fibre.density)
_fat      = FatLayer(0.1, 901.0u"kg/m^3")   # 10% fat fraction → real fat shell

furred(shape) = Body(shape, CompositeInsulation(_fur, _fat))

# Distinct flesh / fat conductivities so both stack terms matter.
const _cond = ThermalConductivities(0.5u"W/m/K", 0.2u"W/m/K", 0.05u"W/m/K")
const _core = 310.15u"K"
const _skin = 305.0u"K"

# Reference values captured from the original shape-specific closed forms (the code the
# radial stack replaced), so this pins the numbers independently rather than comparing the
# stack to its own wrapper. core = 310.15 K, skin = 305.0 K.
@testset "radial stack pins the original closed-form values" begin
    fur = FibrousLayer(5.8u"mm", 30.0u"μm", 5e7u"1/m^2")
    ins_nofat = CompositeInsulation(fur, FatLayer(0.0, 901.0u"kg/m^3"))
    ins_fat   = CompositeInsulation(fur, FatLayer(0.15, 901.0u"kg/m^3"))
    c_nofat = ThermalConductivities(0.9u"W/m/K", 0.23u"W/m/K", 0.05u"W/m/K")
    c_fat   = ThermalConductivities(0.5u"W/m/K", 0.2u"W/m/K",  0.05u"W/m/K")
    cases = (
        ("cyl_nofat", Body(Cylinder(1.0u"kg", ρ_rl, 3.0), ins_nofat), c_nofat, 13.131383303665624),
        ("cyl_fat",   Body(Cylinder(1.0u"kg", ρ_rl, 3.0), ins_fat),   c_fat,    5.267026928959943),
        ("sph_nofat", Body(Sphere(1.0u"kg", ρ_rl), ins_nofat),        c_nofat,  7.226478724342244),
        ("sph_fat",   Body(Sphere(1.0u"kg", ρ_rl), ins_fat),          c_fat,    2.918771069163008),
        ("ell_nofat", Body(example_ellipsoid_shape_pars(mass=33.7u"g", axis_ratio_b=1.1, axis_ratio_c=1.1), ins_nofat), c_nofat, 2.3434263008914344),
        ("ell_fat",   Body(example_ellipsoid_shape_pars(mass=1.0u"kg", axis_ratio_b=1.3, axis_ratio_c=1.5), ins_fat),   c_fat,   2.8978347465563004),
        # Gram-mass ellipsoid *with* fat: flesh/skin b-semi-minor axes come out in different
        # composite unit representations, so the fat shell's inner/outer radii are differently
        # typed. Guards ConductiveShell against requiring identical inner/outer types.
        ("ell_g_fat", Body(example_ellipsoid_shape_pars(mass=33.7u"g", axis_ratio_b=1.1, axis_ratio_c=1.5), ins_fat),   c_fat,   0.9336494822479774),
    )
    for (name, body, c, expected) in cases
        @testset "$name" begin
            got = radial_net_metabolic_heat(body, c, _core, _skin)
            @test ustrip(u"W", got) ≈ expected rtol = 1e-10
            # The public wrapper delegates to the same stack.
            @test net_metabolic_heat(; body, conductivities=c,
                core_temperature=_core, skin_temperature=_skin) == got
        end
    end
end

@testset "radial stack: series resistance adds, generating core is not a plain shell" begin
    for build in (Cylinder(1.0u"kg", ρ_rl, 3.0),
                  Sphere(1.0u"kg", ρ_rl),
                  example_ellipsoid_shape_pars(mass=1.0u"kg", axis_ratio_b=1.3, axis_ratio_c=1.5))
        body = furred(build)
        stack = core_to_skin_stack(body, _cond)
        @test length(stack) == 2
        @test stack[1] isa GeneratingCore
        @test stack[2] isa ConductiveShell
        # Series resistance is the sum of the two layer resistances.
        r_total = stack_resistance(stack, body)
        r_core  = HeatExchange._layer_resistance(stack[1], HeatExchange.shape(body), body)
        r_shell = HeatExchange._layer_resistance(stack[2], HeatExchange.shape(body), body)
        @test r_total ≈ r_core + r_shell
        @test r_core > zero(r_core) && r_shell > zero(r_shell)
    end
end
