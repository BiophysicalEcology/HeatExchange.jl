"""
    bird_skin_area

Skin area for a bird based on Walsberg & King (1978)
    Journal of Experimental Biology, 76, 185-189.
Note that skin area is greater than plumage area.
"""
function bird_skin_area(shape)
    return Unitful.uconvert(u"m^2", (10.0 * Unitful.ustrip(u"g", shape.mass) ^ 0.667)u"cm^2")
end

"""
    bird_plumage_area

Plumage area for a bird based on Walsberg & King (1978)
    Journal of Experimental Biology, 76, 185-189.
Note that plumage area is less than skin area.
"""
function bird_plumage_area(shape)
    return Unitful.uconvert(u"m^2", (8.11 * Unitful.ustrip(u"g", shape.mass) ^ 0.667)u"cm^2")
end