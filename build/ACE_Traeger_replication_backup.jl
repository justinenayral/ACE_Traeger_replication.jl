module ACE_Traeger_replication

import MAT #To open .mat datasets
import Printf #For labelling graphs
import LinearAlgebra
import NLsolve
import DataFrames
import Plots
import CSV
import XLSX

using MAT #To open .mat datasets
using Printf #For labelling graphs
using LinearAlgebra
using NLsolve
using DataFrames
using Plots
using CSV
using XLSX

greet() = print("Hello ACE!")

function figure2() 
    
    #Need to run the julia file containing all the functions: DamageSterner.jl : need to create a package with documentation
    
    
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
    
    #Plot specification:
    range1_inC=[0, 10]
    range2_inC=[0, 10]
    st=0.01
    
    T = range(range2_inC[1], stop=range2_inC[2], step=st)
    length(T) #Same length as in Matlab code
    range1=[1:1:(range1_inC[2]-range1_inC[1])/st+1]
    range2 = [((i - range1_inC[1]) / st + 1) for i in range(range2_inC[1], stop=range2_inC[2], step=st)]
    
    # DICE 2013R by Nordhaus:
    # Climate damage parameters
    # a1: Damage intercept                                 /0       /
    # a2: Damage quadratic term                            /0.00267 /
    # a3: Damage exponent                                  /2.00    /
    a_Nord_13 = 0.00267
    
    # Damage Parameters Nordhaus DICE 2016
    a_Nord_16 = 0.00236
    
    # Damage Parameters Nordhaus DICE 2007
    a_Nord = 0.0028
    
    # ACE model parameters
    cs = 3
    xi1 = 0.23104906018 * 3 / cs
    
    ##Matching Nordhaus and HSP, damages parameters
    #matchtemp=2.5
    # Extract the solution
    xi0 = xi0DICE([0.0], a_Nord, 2.5)
    print("For DICE: cs=", round(Int, cs), " find xi0=", xi0, " and xi1=", xi1 )
    
    
    #Matching Howard Sterner 2017 Damages parameters
    #Damage Parameters Howard Sterner:
    a_Sterner=0.01145
    
    #matchtemp = 3  # 3 for paper's match
    
    # Solve for Sterner's xi0
    xi0_Sterner=xi0Sterner([0.0], a_Sterner, 3)
    
    print("For HSP1: cs=", round(Int, cs), " find xi0=", xi0_Sterner, " and xi1=", xi1 )
    
    # Define the function to Solve for Sterner's second xi0
    xi0_Sterner2=xi0DICE([0.0], a_Sterner, 3)
    
    # Display the result for Sterner's second xi0
    print("For HSP1: cs=", round(Int, cs), " find xi0=", xi0_Sterner2, " and xi1=", xi1 )
    
    
    ##Damage functions: 
    Dam_Nord = dam_DICE(a_Nord, T)
    Dam_ACE = dam_ACE(xi0, xi1, T)
    Dam_Weitz = dam_Weitz(T)
    Dam_Sterner = dam_Sterner(T)
    Dam_Sterner2 = dam_DICE(a_Sterner, T)
    Dam_Nord_13 = dam_DICE(a_Nord_13, T)
    Dam_Nord_16 = dam_DICE(a_Nord_16, T)
    
    #Define +/- 50% variations:
    Dam_ACE_V1 = dam_ACE(xi0 * 0.5, xi1, T)
    Dam_ACE_V2 = dam_ACE(xi0 * 1.5, xi1, T)
    
    #Sterner matched at original function:
    Dam_ACE_Sterner = dam_ACE(xi0_Sterner, xi1, T)
    Dam_ACE_S_V1 = dam_ACE(0.5*xi0_Sterner, xi1, T)
    Dam_ACE_S_V2 = dam_ACE(1.5*xi0_Sterner, xi1, T)
    
    # Sterner matched at 1/1+aT^2 function:
    Dam_ACE_Sterner2 = dam_ACE(xi0_Sterner2, xi1, T)
    Dam_ACE_S2_V1 =  dam_ACE(0.5*xi0_Sterner2, xi1, T)
    Dam_ACE_S2_V2 = dam_ACE(1.5*xi0_Sterner2, xi1, T)
    
    #Draw the graph
    
    plotrange = 0.71
    h1 = plot(title="Damage Function Calibration", xlabel="Degree C above Preindustrial", ylabel="Damages in Percent of Output", ylim=(0, plotrange), xlim=(minimum(T), maximum(T)), legend=:topleft, size=(400,300))
    
    plot!(T, Dam_Nord, color=:olivedrab2, linestyle=:solid, label="DICE 2007", linewidth=2.0)
    plot!(T, Dam_Nord_13, color=:forestgreen, linestyle=:solid, label="DICE 2013", linewidth=2.0) # Assuming Nordall == 1
    plot!(T, Dam_Nord_16, color=:darkgreen, linestyle=:solid, label="DICE 2016", linewidth=2.0) # Assuming Nordall == 1
    plot!(T, Dam_ACE, color=:olivedrab2, linestyle=:dash, label="ACE base", linewidth=2.0)
    plot!(T, Dam_ACE_V2, color=:turquoise3, linestyle=:dash, label="ACE base ± 50%", linewidth=2.0)
    plot!(T, Dam_Sterner2, color=:red, linestyle=:solid, label="HSP-norm", linewidth=2.0)
    plot!(T, Dam_Sterner, color=:gray, linestyle=:solid, label="HSP-nn", linewidth=2.0)
    plot!(T, Dam_ACE_Sterner2, color=:red, linestyle=:dash, label="ACE HSP", linewidth=2.0)
    plot!(T, Dam_ACE_Sterner, color=:black, linestyle=:dash, label="ACE HSP*", linewidth=2.0)
    plot!(T, Dam_ACE_V1, color=:turquoise3, linestyle=:dash, label="", linewidth=2.0)
    #plot!(T, Dam_ACE_S_V1, color=:magenta, linestyle=:dot, label="ACE HSP ± 50%") # Assuming HSPpm == 1
    #plot!(T, Dam_ACE_S_V2, color=:magenta, linestyle=:dot, label="ACE HSP ± 50%") # Assuming HSPpm == 1
    
    yticks!([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], ["0", "10", "20", "30", "40", "50", "60", "70"])
    
    #Plot for larger range 2: 
    plotrange = 0.31
    h2 = plot(title="Damage Function Calibration", xlabel="Degree C above Preindustrial", ylabel="Damages in Percent of Output", ylim=(0, plotrange), xlim=(minimum(T), 5), legend=:topleft, size=(400,300))
    
    plot!(T, Dam_Nord, color=:olivedrab2, linestyle=:solid, label="", linewidth=2.0)
    plot!(T, Dam_Nord_13, color=:forestgreen, linestyle=:solid, label="", linewidth=2.0) # Assuming Nordall == 1
    plot!(T, Dam_Nord_16, color=:darkgreen, linestyle=:solid, label="", linewidth=2.0) # Assuming Nordall == 1
    plot!(T, Dam_ACE, color=:olivedrab2, linestyle=:dash, label="", linewidth=2.0)
    plot!(T, Dam_ACE_V2, color=:turquoise3, linestyle=:dash, label="", linewidth=2.0)
    plot!(T, Dam_Sterner2, color=:red, linestyle=:solid, label="", linewidth=2.0)
    plot!(T, Dam_Sterner, color=:gray, linestyle=:solid, label="", linewidth=2.0)
    plot!(T, Dam_ACE_Sterner2, color=:red, linestyle=:dash, label="", linewidth=2.0)
    plot!(T, Dam_ACE_Sterner, color=:black, linestyle=:dash, label="", linewidth=2.0)
    plot!(T, Dam_ACE_V1, color=:turquoise3, linestyle=:dash, label="", linewidth=2.0)
    #plot!(T, Dam_ACE_S_V1, color=:magenta, linestyle=:dot, label="ACE HSP ± 50%") # Assuming HSPpm == 1
    #plot!(T, Dam_ACE_S_V2, color=:magenta, linestyle=:dot, label="ACE HSP ± 50%") # Assuming HSPpm == 1
    
    yticks!([0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3], ["0", "5", "10", "15", "20", "25", "30"])
    
    h = plot(h1, h2, layout=(1,2), size=(800,550))

end


end # module ACE_Traeger_replication
