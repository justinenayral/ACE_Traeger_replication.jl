#  Joos (2013) Carbon impuls response
#  Best fit to CMIP 5 specification

# Fraction staying forever
a0 = 0.2173

# Fractions moving into other boxes
a = [0.2240, 0.2824, 0.2763]

# Half life of different reservoirs
tau_life = [394.4, 36.54, 4.304]

half_life = tau_life .* log(2)
gamma_annual = exp.(-1 ./ tau_life)
gamma_5year = exp.(-5 ./ tau_life)
gamma_decadal = exp.(-10 ./ tau_life)

if timestep in [1, 5, 10]
    if timestep == 1
        CarbonMatrix = diagm([1; gamma_annual])
    elseif timestep == 5
        CarbonMatrix = diagm([1; gamma_5year])
    elseif timestep == 10
        CarbonMatrix = diagm([1; gamma_decadal])
    end
else
    # default is decadal
    CarbonMatrix = diagm([1; gamma_decadal])
end

CarbonWeights = [a0; a]