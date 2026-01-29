"""
    interpolate(a::FibreProperties, b::FibreProperties, t::Number)
    interpolate(a::Number, b::Number, t)

Linearly interpolate between two `FibreProperties` objects.

Returns `a * (1 - t) + b * t` for each field.
"""
function interpolate(a::FibreProperties, b::FibreProperties, t::Number)
    interpolated = map(getproperties(a), getproperties(b)) do a, b
        interpolate(a, b, t)
    end
    return setproperties(a; interpolated...)
end
interpolate(a::Number, b::Number, t::Number) = a * (1 - t) + b * t
