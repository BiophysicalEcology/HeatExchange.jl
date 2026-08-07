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

@testset "radial stack reproduces net_metabolic_heat (core→skin)" begin
    shapes = (
        "cylinder" => (m -> Cylinder(m, ρ_rl, 3.0)),
        "sphere"   => (m -> Sphere(m, ρ_rl)),
    )
    for (name, build) in shapes
        @testset "$name" begin
            for mass in (0.1u"kg", 1.0u"kg", 10.0u"kg", 100.0u"kg")
                body = furred(build(mass))
                # Fat is genuinely present, so the shell term is exercised.
                @test flesh_radius(body) < skin_radius(body)

                ref = net_metabolic_heat(; body, conductivities = _cond,
                                         core_temperature = _core, skin_temperature = _skin)
                got = radial_net_metabolic_heat(body, _cond, _core, _skin)
                @test got ≈ ref rtol = 1e-10
            end
        end
    end
end

@testset "radial stack: series resistance adds, generating core is not a plain shell" begin
    body = furred(Cylinder(1.0u"kg", ρ_rl, 3.0))
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
