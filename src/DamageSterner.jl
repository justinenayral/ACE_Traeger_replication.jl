#This function define the Sterner damage function: add into a package

#For sterner: get the xi0
function xi0Sterner(initial_guess::Vector, a, matchtemp, T=2.5)
    function xi0solve_Sterner!(F, xi0)
        F[1] = a_Sterner * matchtemp^2 - (1 - exp(-xi0[1] * exp(xi1 * matchtemp) + xi0[1]))
    end
    result = nlsolve(xi0solve_Sterner!, initial_guess) # Solve the equation using NLsolve
    x = result.zero[1] # Extract the solution
    return x
end

#For DICE: get the xi0
function xi0DICE(initial_guess::Vector, a, matchtemp, T=2.5)
    function xi0solve_DICE!(F, xi0)
        F[1] = 1 / (1 + a * matchtemp^2) - exp(-xi0[1] * exp(xi1 * matchtemp) + xi0[1])
    end
    result = nlsolve(xi0solve_DICE!, initial_guess) # Solve the equation using NLsolve
    x = result.zero[1] # Extract the solution
    return x
end

#Define the damage functions
#For Dice:
function dam_DICE(a, T)
    x=1 .- 1 ./ (1 .+ a.* T.^2)
    return x
end

#For Howard and Sterner(2017): with the default value a, similar as in Sterner
function dam_Sterner(T, a=0.01145)
    x= a .* T.^2
    return x
end

#For Dam_Weitz
function dam_Weitz(T)
    x=1 .- 1 ./ (1 .+ (T ./ 20.46).^2 .+ (T ./ 6.081).^6.754)
    return x
end

#For ACE:
function dam_ACE(xi0, xi1, T)
    x=1 .- exp.(-xi0 * exp.(xi1 .* T) .+ xi0)
    return x
end

##Damage functions:

function squeeze(A::Vector{Float64})
    singleton_dims = tuple((d for d in 1:ndims(vec(A)) if size(vec(A), d) == 1)...)
    return squeeze(vec(A), singleton_dims)
end