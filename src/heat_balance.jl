
# Heat balance


function heat_balance() end

heat_balance(b::AbstractBody, parameters, environment) = begin
Qrad_in = radin()
Qsolar = solar()
end