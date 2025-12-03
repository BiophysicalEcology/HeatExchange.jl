shape(body::AbstractBody) = body.shape
insulation(body::AbstractBody) = body.insulation
geometry(body::AbstractBody) = body.geometry
surface_area(body::AbstractBody) = surface_area(shape(body), body)

# functions to extract appropriate surface areas from different objects

get_total_area(body::AbstractBody) = get_total_area(shape(body), insulation(body), body)
get_skin_area(body::AbstractBody) = get_skin_area(shape(body), insulation(body), body)
get_evaporation_area(body::AbstractBody) = get_evaporation_area(shape(body), insulation(body), body)

# for composite insulation cases (fat and fur/feathers)
outer_insulation(ins::Insulation) = ins
outer_insulation(ins::CompositeInsulation) = begin
    # find fur layer if present
    fur_layer = findlast(i -> i isa Fur, ins.layers)
    if fur_layer !== nothing
        ins.layers[fur_layer]
    else
        # otherwise the last layer
        ins.layers[end]
    end
end

get_total_area(shape, ins::CompositeInsulation, body) =
    get_total_area(shape, outer_insulation(ins), body)
get_skin_area(shape, ins::CompositeInsulation, body) =
    get_skin_area(shape, outer_insulation(ins), body)
get_evaporation_area(shape, ins::CompositeInsulation, body) =
    get_evaporation_area(shape, outer_insulation(ins), body)

# functions to get the appropriate radii
get_r_skin(body::AbstractBody) = get_r_skin(shape(body), insulation(body), body)
get_r_insulation(body::AbstractBody) = get_r_insulation(shape(body), insulation(body), body)
get_r_flesh(body::AbstractBody) = get_r_flesh(shape(body), insulation(body), body)

# functions to compute silhouette area as a function of zenith angle (passive) or 
# orientation (thermoregulating)

"""
    silhouette_area(body::AbstractBody, θ)

Calculates the silhouette (projected) area of an object given a solar zenith angle, θ.

"""
silhouette_area(body::AbstractBody, θ) = silhouette_area(shape(body), insulation(body), body, θ)

"""
    silhouette_area(body::AbstractBody, orientation)

Calculates the silhouette (projected) area of an object given an orientation towards the sun (normal, parallel or inbetween).

"""
silhouette_area(body::AbstractBody) = silhouette_area(shape(body), insulation(body), body)

# Orientation-specific implementations
silhouette_area(body::AbstractBody, ::NormalToSun) =
    silhouette_area(body).normal

silhouette_area(body::AbstractBody, ::ParallelToSun) =
    silhouette_area(body).parallel

silhouette_area(body::AbstractBody, ::Intermediate) =
    (silhouette_area(body).normal + silhouette_area(body).parallel) * 0.5

function calc_area_hair(fibre_diameter, fibre_density, skin)
    π * (fibre_diameter / 2) ^ 2 * (fibre_density * skin)
end