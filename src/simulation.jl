
function simulation end

function simulation(o::Model, t::MicroclimPoint)
    
    pred_tbs = Vector{Float64}(undef,length(t.radiation))u"°C"
    local T_core_s
    local heat_balance_out
    
    for i in eachindex(t.radiation)
        # get the environmental parameters
        env_params = EnvironmentalParams()
        # get the variables for the environment
        env_vars = EnvironmentalVars(t.airtemperature[i,1],
                                    t.skytemperature[i],
                                    t.soiltemperature[i,1],
                                    t.relhumidity[i,1],
                                    t.windspeed[i,1],
                                    101325.0Pa,
                                    t.zenith[i],
                                    0.5W/m/K,
                                    t.radiation[i],
                                    t.radiation[i] * 0.96,
                                    t.radiation[i] / 10)
        # get the variables for the organism
        org_vars = if i == 1 
                        OrganismalVars()
                    elseif T_core_s === nothing
                        OrganismalVars()
                    else
                        OrganismalVars(heat_balance_out.T_core,
                           heat_balance_out.T_surf,
                           heat_balance_out.T_lung,
                           -707J/kg)
                    end
        
        # combine variables for both the environment and the organism
        variables = (organism = org_vars, environment = env_vars)
        
        T_air = env_vars.T_air
        try
            T_core_s = find_zero(t -> heat_balance(t, o, env_params, variables), (T_air - 40K, T_air + 100K), Bisection())
        catch e
            T_core_s = nothing
            continue
        end
        pred_tbs[i] = (Unitful.ustrip(T_core_s) - 273.15)°C
        heat_balance_out = heat_balance(T_core_s, o, env_params, variables)
    end

    return pred_tbs
end

