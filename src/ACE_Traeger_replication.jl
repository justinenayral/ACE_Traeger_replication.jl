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

path = "C:/Users/Guada/Dropbox/Traeger/Julia code/"

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

function figure3()
##First, the authors use the MAGICC6.0 model by Meinshausen, Raper and Wigley(2011) to simulate
#the RCP scenario over a time horizon of 400 years.

#The calibration of ACE uses 2 ocean layers(upper and deep) compared to MAGGIC's 50 layers and DICE's ocean layers

#Setup:
Startdate=2015
Enddate=2400

#Fixed settings of final Calibration
cs=3 #climate sensitivity
η=3.8 #forcing parameter used to calculate CO2equivalents from RF
timestep=10 #Number of time step
Oceanlayers=2 #Number of ocean layers used in calibration
Selectlayers=[5, 40] # Ocean layer's of Magicc used in calibration for shallow ocean (layer 5) and deep ocean (layer 40)
weight=[4, 2, 1] # sets weights for error returned, weight on atmosphere, shallow, and deep ocean
σ=[0.538194273157518,	0.0765198358638229,	0.00399267991987196,	0.461805699284704,	0.0274865322002226] #Corresponding sigma (heat flows between adjacent layers), same as in TMDynamics.m

#importing Magicc 6.0 data generated by CT using downloaded version of Magicc 6.
data=matread(path*"//MagiccOcean.mat")


#Get the name of all the variables: MagiccOcean is a Dictionnary with Dictionnaries inside
scenario = collect(keys(data["MagiccOcean"]))
S=length(scenario) #Number of scenarios

#define the first and the end dates
startpos = findfirst(data["MagiccOcean"][scenario[1]]["Year"].==Startdate)[1]
endpos = findfirst(data["MagiccOcean"][scenario[1]]["Year"].==Enddate)[1]

#CHeck that the index for the first and the last dates in the data are consistent
if endpos-startpos!=Enddate-Startdate
    error("Years Missing, interrrupt")
end

#Read the Layer's names and initialize the Temp matrix (layers (atm and ocean) , scenarios , horizon )
Layer=data["MagiccOcean"][scenario[1]]["Layer"] #Read in layer names
Temp = zeros(Oceanlayers + 1, S, floor(Int, (endpos - startpos) / timestep + 1))
RF = zeros(S, floor(Int, (endpos - startpos) / timestep + 1)) #Radiative forcing (greenhouse effect)

for s=1:S
    RF[s, :] = data["MagiccOcean"][scenario[s]]["RFtot"][startpos:timestep:endpos]'
    Temp[1, s, :] = data["MagiccOcean"][scenario[s]]["Tempatm"][startpos:timestep:endpos]' 
    for l = 1:Oceanlayers
        Temp[l+1, s, :] = data["MagiccOcean"][scenario[s]]["Temp"][startpos:timestep:endpos, Selectlayers[l]]' 
    end
end 

#forcing
forcing=exp.(log(2)/η*RF)

#We define the parameter xi1=log(2)/s:
xi_base=zeros(Oceanlayers+1)

for i in 1:length(xi_base)
xi_base[i]=log(2)/cs
end

#Get the size of the three vectors
(layers, scenarios, horizon) = size(Temp)

#Temp: 3-D temperature field: layers (atmosphere and ocean) , scenarios , horizon 

logic = true

# Call the function TempSimulation_ACE to run the simulation
Simulation=TempSimulation_ACE(σ, Temp, forcing, xi_base, weight, logic)


##Evaluating results for print out on screen and saving, using Magicc data
# Assign sigma_up (first layers) and sigma_down (upper layers):
optim=Dict()
optim["σ_up"]=σ[1:layers]
optim["σ_down"]=σ[layers+1:2*layers-1]

# Construct TempHelpDown and TempHelpUp matrices
TempHelpDown = [zeros(layers-1, 1) diagm(optim["σ_down"][1:layers-1]); zeros(1, layers)]
TempHelpUp = [zeros(1, layers); diagm(optim["σ_up"][2:layers]) zeros(layers-1, 1)] #OK

# Assign taubar
taubar = 1
optim["TempMatrix"]=diagm(1 .- optim["σ_up"] - [optim["σ_down"]; 0]) + TempHelpDown + TempHelpUp
optim["xi"]=Simulation["xi"]

pos_2015 = findfirst(data["MagiccOcean"][scenario[1]]["Year"].==2015)[1] 

RF_init=zeros(S,1)
Temp_init=zeros(Oceanlayers+1, S)
for s = 1:S
    RF_init[s, :] .= data["MagiccOcean"][scenario[s]]["RFtot"][pos_2015]'
    Temp_init[1, s, :] .= data["MagiccOcean"][scenario[s]]["Tempatm"][pos_2015]' # Note: Use square brackets for indexing dictionaries
    for l = 1:Oceanlayers
        Temp_init[l + 1, s, :] .= data["MagiccOcean"][scenario[s]]["Temp"][pos_2015]'
    end
end

optim["τ_initial"] = exp.(diagm(optim["xi"]) * Temp_init)
optim["Temp_initial"] = Temp_init # [layers scenarios] in 2015


####Plotting
#Export extra datasets to get other specifications
d=matread(path*"//TempDataCalel.mat")

#Necessary variables for graph
d["data"]["year_ACE"] = Startdate:timestep:Enddate
d["data"]["year_Magicc"] = Startdate:timestep:Enddate

#CHeck the minimum and maximum year of each specifications
maxyearplot = maximum(d["data"]["year_ACE"]) + 5
minyearplot = minimum(d["data"]["year_ACE"]) - 5

maxyearplot_DICE = findlast(x -> x <= min(maximum(d["data"]["DICE13"][:, 2]), maxyearplot), d["data"]["DICE13"][:, 2])
minyearplot_DICE = findfirst(x -> x >= max(minimum(d["data"]["DICE13"][:, 2]), minyearplot), d["data"]["DICE13"][:, 2])

maxyearplot_PAGE = findlast(x -> x <= min(maximum(d["data"]["PAGE"][:, 2]), maxyearplot), d["data"]["PAGE"][:, 2])
minyearplot_PAGE = findfirst(x -> x >= max(minimum(d["data"]["PAGE"][:, 2]), minyearplot), d["data"]["PAGE"][:, 2])

maxyearplot_FUND = findlast(x -> x <= min(maximum(d["data"]["FUND"][:, 2]), maxyearplot), d["data"]["FUND"][:, 2])
minyearplot_FUND = findfirst(x -> x >= max(minimum(d["data"]["FUND"][:, 2]), minyearplot), d["data"]["FUND"][:, 2])

#Define the legend
legend=["RCP 6 to 4.5", "RCP 3", "RCP 4.5", "RCP 8.5", "RECP 4.5 to 3", "RCP 6"]
#Initialization of the variables for building the plot
Temp_plot=zeros(39,3)
TempSim_plot=zeros(39,3)
Temp_combined_plot_atm=zeros(39,6)
Temp_combined_plot_atm_sim=zeros(39,6)

#
for i=1:scenarios # This loop over scenarios closes only for final all scenario plots

    for l = 1:Oceanlayers + 1
        Temp_plot[:, l] = Temp[l, i, :]
        TempSim_plot[:, l] .= Simulation["Tempsim"][l, i, :]
    
        if l == 1
            # Saving for plot containing only atm layer but all scenarios
            Temp_combined_plot_atm[:, i] = Temp[l, i, :]
            Temp_combined_plot_atm_sim[:, i] = Simulation["Tempsim"][l, i, :]
        end
    end

end

# Joint scenario plots
g = plot(title="Temp Fit all scen", xlabel="year", ylabel="Temperature")
legend=["RCP 6 to 4.5", "RCP 3", "RCP 4.5", "RCP 8.5", "RECP 4.5 to 3", "RCP 6"]
colors = ["deepskyblue2",  "midnightblue",   "red2",   "olivedrab3",   "darkgreen",   "indigo"]

order=[2, 5, 3, 1, 6, 4]

for i in order
    plot!(Startdate:timestep:Enddate, Temp_combined_plot_atm_sim[:,i], linestyle=:dash, color=colors[i], linewidth=2, label="")
    plot!(Startdate:timestep:Enddate, Temp_combined_plot_atm[:,i], color=colors[i], linewidth=2, label=legend[i])
end

hline!([-Inf], color=:black, linestyle=:solid, label="MAGICC")
hline!([-Inf], color=:black, linestyle=:dash, label="ACE", linewidth=2)

# Plot DICE13, PAGE, and FUND
plot!(d["data"]["DICE13"][minyearplot_DICE:maxyearplot_DICE, 2], d["data"]["DICE13"][minyearplot_DICE:maxyearplot_DICE, 3], linestyle=:dash, linewidth=1, label="DICE", color=:pink)
plot!(d["data"]["PAGE"][minyearplot_PAGE:maxyearplot_PAGE, 2], d["data"]["PAGE"][minyearplot_PAGE:maxyearplot_PAGE, 3], linestyle=:dash, linewidth=1, label="PAGE", color=:lightskyblue2)
plot!(d["data"]["FUND"][minyearplot_FUND:maxyearplot_FUND, 2], d["data"]["FUND"][minyearplot_FUND:maxyearplot_FUND, 3], linestyle=:dash, linewidth=1, label="FUND", color=:grey)

for i in 2:4
    plot!(d["data"]["DICE13"][minyearplot_DICE:maxyearplot_DICE, 2], d["data"]["DICE13"][minyearplot_DICE:maxyearplot_DICE, 2+i], linestyle=:dash, linewidth=1, label="", color=:pink)
    plot!(d["data"]["PAGE"][minyearplot_PAGE:maxyearplot_PAGE, 2], d["data"]["PAGE"][minyearplot_PAGE:maxyearplot_PAGE, 2+i], linestyle=:dash, linewidth=1, label="", color=:lightskyblue2)
    plot!(d["data"]["FUND"][minyearplot_FUND:maxyearplot_FUND, 2], d["data"]["FUND"][minyearplot_FUND:maxyearplot_FUND, 2+i], linestyle=:dash, linewidth=1, label="", color=:grey)
end
plot!(legendfontsize=6)

# Set xlabel and ylabel
xlabel!("year")
ylabel!("Degrees (Celsius) above preindustrial")
title!("Temperature dynamics calibration")

end


function TempSimulation_ACE(σ, Temp, forcing, xi, weight, logic)

#Inputs:
#Temp: 3-D temparature field: layers (atm and ocean) , scenarios , horizon 
#Forcing: 2-D scenarios , time 
#sigma: 
#first "1:layer" entries: sigma_up
#next  "1:layer-1" entries: sigma_down
#xi: vector of length layers
#weight: vector of length layers, weights on squared differences

#Output: 
#squared and weighted sum of errors, Simulated T 

    (layers, scenarios, horizon) = size(Temp)

    ## Extracing different endogneously optimized parameters from sigma
    xione=xi[1]
    σ_up=σ[1:layers] 
    σ_down=vcat(σ[layers+1:2*layers-1],0) 
    τ_bar=1
    #Standard version: no heat exchange with lower boundary layer



    ##Initializing and Simulating
    τ=zeros(layers, scenarios, horizon) 
    τ[:,:,1]=exp.(diagm(xi)*Temp[:,:,1])
    TempSim=zeros(layers, scenarios, horizon)
    TempSim[:,:,1]=Temp[:,:,1]

    for y = 1:horizon-1
        # Need time in last, layers in rows of the remaining matrix, so scenarios in columns
        τ[:, :, y+1] = diagm(1 .- σ_up .- σ_down) * τ[:, :, y] +
                        diagm(σ_up) * [forcing[:, y]'; τ[1:layers-1, :, y]] +
                        diagm(σ_down) * [τ[2:layers, :, y]; fill(τ_bar, 1, scenarios)]
        TempSim[:, :, y+1] = diagm(1 ./ xi) * log.(τ[:, :, y+1])
    end #OK


    ##Evaluating
    diff_sum=0

    for l = 1:layers
        diff = weight[l] .* (TempSim[l, :, :] .- Temp[l, :, :]) .^ 2
        diff_sum += sum(diff)
    end #OK

    # Set err to diff_sum
    err = Float64(diff_sum)

    # Display the result
    println(" Residual ", err, " σ-vector ", σ)
    println()  # Print an empty line

    return Dict("err"=>err, "Tempsim"=>TempSim, "xi"=>xi)

end

end # module ACE_Traeger_replication
