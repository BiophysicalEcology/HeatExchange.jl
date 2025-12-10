"""
    mammal_skin_area

Skin area for a bird based on Walsberg & King (1978)
    Journal of Experimental Biology, 76, 185-189.
"""
function mammal_skin_area(shape)
    return Unitful.uconvert(u"m^2", (1110.0 * Unitful.ustrip(u"kg", shape.mass) ^ 0.65)u"cm^2")
end

"""
    mammal_fur_area

Fur area for a mammal based on Stahl (1967)
    Journal of Applied Physiology, 22, 453-460.
"""
function mammal_fur_area(body::AbstractBody)
    fur_multiplier = total_area(body) / skin_area(body)
    return Unitful.uconvert(u"m^2", (1110.0 * Unitful.ustrip(u"kg", body.shape.mass) ^ 0.65 * fur_multiplier)u"cm^2")
end