using HeatExchange
using Unitful
using Test
using Unitful: °, rad, °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using CSV, DataFrames

testdir = realpath(joinpath(dirname(pathof(HeatExchange)), "../test"))

names_in = [
    :GEOMETRY, :shape_b, :shape_c, :AMASS, :ANDENS, :SKINW, :SKINT, :RINSUL, :PTCOND, :PMOUTH, :PANT
]
names_out = [
    :AREA, :AV, :AT, :AL, :VOL, :R, :R1, :ASEMAJR, :BSEMINR, :CSEMINR, :ASILN, :ASILP, :AEFF
]

# plate test
plate_in_NMR_vec = (DataFrame(CSV.File("$testdir/data/plate_in.csv")))[:, 2]
plate_in_NMR = (; zip(names_in, plate_in_NMR_vec)...)
plate_out_NMR_vec = (DataFrame(CSV.File("$testdir/data/plate_out.csv")))[:, 2]
plate_out_NMR = (; zip(names_out, plate_out_NMR_vec)...)

density = (plate_in_NMR.ANDENS)u"kg/m^3"
mass = (plate_in_NMR.AMASS)u"kg"
shapeb = plate_in_NMR.shape_b
shapec = plate_in_NMR.shape_c
shape = Plate(mass, density, shapeb, shapec)
plate_out = Body(shape, Naked())

@testset "plate shape comparisons" begin
    @test plate_out.geometry.area ≈ (plate_out_NMR.AREA)u"m^2" atol=1e-6u"m^2"
    # TODO work out best strategy for characteristic dimension of plate for NicheMapR and HeatExchange
    @test plate_out.geometry.characteristic_dimension ≈ (plate_out_NMR.VOL ^ (1/3))u"m" atol=1e-6u"m"
    @test plate_out.geometry.volume ≈ (plate_out_NMR.VOL)u"m^3" atol=1e-6u"m^3"
end

# cylinder test
cylinder_in_NMR_vec = (DataFrame(CSV.File("$testdir/data/cylinder_in.csv")))[:, 2]
cylinder_in_NMR = (; zip(names_in, cylinder_in_NMR_vec)...)
cylinder_out_NMR_vec = (DataFrame(CSV.File("$testdir/data/cylinder_out.csv")))[:, 2]
cylinder_out_NMR = (; zip(names_out, cylinder_out_NMR_vec)...)

density = (cylinder_in_NMR.ANDENS)u"kg/m^3"
mass = (cylinder_in_NMR.AMASS)u"kg"
shapeb = cylinder_in_NMR.shape_b
shape = Cylinder(mass, density, shapeb)
cylinder_out = Body(shape, Naked())

@testset "cylinder shape comparisons" begin
    @test cylinder_out.geometry.area ≈ (cylinder_out_NMR.AREA)u"m^2" atol=1e-6u"m^2"
    @test cylinder_out.geometry.characteristic_dimension ≈ (cylinder_out_NMR.VOL ^ (1/3))u"m" atol=1e-6u"m"
    @test cylinder_out.geometry.volume ≈ (cylinder_out_NMR.VOL)u"m^3" atol=1e-6u"m^3"
end 

# ellipsoid test
ellipsoid_in_NMR_vec = (DataFrame(CSV.File("$testdir/data/ellipsoid_in.csv")))[:, 2]
ellipsoid_in_NMR = (; zip(names_in, ellipsoid_in_NMR_vec)...)
ellipsoid_out_NMR_vec = (DataFrame(CSV.File("$testdir/data/ellipsoid_out.csv")))[:, 2]
ellipsoid_out_NMR = (; zip(names_out, ellipsoid_out_NMR_vec)...)

density = (ellipsoid_in_NMR.ANDENS)u"kg/m^3"
mass = (ellipsoid_in_NMR.AMASS)u"kg"
shapeb = ellipsoid_in_NMR.shape_b
shapec = ellipsoid_in_NMR.shape_c
shape = Ellipsoid(mass, density, shapeb, shapec)
ellipsoid_out = Body(shape, Naked())

@testset "ellipsoid shape comparisons" begin
    @test ellipsoid_out.geometry.area ≈ (ellipsoid_out_NMR.AREA)u"m^2" atol=1e-6u"m^2"
    @test ellipsoid_out.geometry.characteristic_dimension ≈ (ellipsoid_out_NMR.VOL ^ (1/3))u"m" atol=1e-6u"m"
    @test ellipsoid_out.geometry.volume ≈ (ellipsoid_out_NMR.VOL)u"m^3" atol=1e-6u"m^3"
end

# desert iguana test
iguana_in_NMR_vec = (DataFrame(CSV.File("$testdir/data/iguana_in.csv")))[:, 2]
iguana_in_NMR = (; zip(names_in, iguana_in_NMR_vec)...)
iguana_out_NMR_vec = (DataFrame(CSV.File("$testdir/data/iguana_out.csv")))[:, 2]
iguana_out_NMR = (; zip(names_out, iguana_out_NMR_vec)...)

density = (iguana_in_NMR.ANDENS)u"kg/m^3"
mass = (iguana_in_NMR.AMASS)u"kg"
shape = DesertIguana(mass, density)
iguana_out = Body(shape, Naked())

@testset "iguana shape comparisons" begin
    @test iguana_out.geometry.area ≈ (iguana_out_NMR.AREA)u"m^2" atol=1e-6u"m^2"
    @test iguana_out.geometry.characteristic_dimension ≈ (iguana_out_NMR.VOL ^ (1/3))u"m" atol=1e-6u"m"
    @test iguana_out.geometry.volume ≈ (iguana_out_NMR.VOL)u"m^3" atol=1e-6u"m^3"
end

# leopard frog test
frog_in_NMR_vec = (DataFrame(CSV.File("$testdir/data/frog_in.csv")))[:, 2]
frog_in_NMR = (; zip(names_in, frog_in_NMR_vec)...)
frog_out_NMR_vec = (DataFrame(CSV.File("$testdir/data/frog_out.csv")))[:, 2]
frog_out_NMR = (; zip(names_out, frog_out_NMR_vec)...)

density = (frog_in_NMR.ANDENS)u"kg/m^3"
mass = (frog_in_NMR.AMASS)u"kg"
shape = LeopardFrog(mass, density)
frog_out = Body(shape, Naked())

@testset "frog shape comparisons" begin
    @test frog_out.geometry.area ≈ (frog_out_NMR.AREA)u"m^2" atol=1e-6u"m^2"
    @test frog_out.geometry.characteristic_dimension ≈ (frog_out_NMR.VOL ^ (1/3))u"m" atol=1e-6u"m"
    @test frog_out.geometry.volume ≈ (frog_out_NMR.VOL)u"m^3" atol=1e-6u"m^3"
end