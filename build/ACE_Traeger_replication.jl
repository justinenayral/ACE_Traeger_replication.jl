module ACE_Traeger_replication()

import MAT
import Printf
import LinearAlgebra
import NLsolve
import DataFrames
import Plots
import CSV
import XLSX

using MAT
using Printf
using LinearAlgebra
using NLsolve
using DataFrames
using Plots
using CSV
using XLSX

timestep = 10

############################################### 

############################################### DAMAGE FUNCTIONS

    
#Damage function
#Define the damage functions
#For Dice:
"""
DICE damage function replace the quadratic damage term ``aT^2``, to limit damages to 100 percent of production. 
``D(T)=1-(1/(1+aT^2))``
    
*Input:*
- ``a`` : Damage coefficients (``a=0.0028`` for 2007, ``a=0.00267`` for 2013, ``a=0.00236`` for 2016)
- ``T``: temperature
    
*Syntax:*
```julia
ACE_Traeger_replication.dam_DICE(a, T)
```
*Output:*
- The calculated DICE damage estimate.

*Example*
```julia-repl
julia> ACE_Traeger_replication.dam_Sterner(0.0028,3)
0.02458056964494726
```
"""
function dam_DICE(a, T)
    x=1 .- 1 ./ (1 .+ a.* T.^2)
    return x
end
    
#For Howard and Sterner(2017): with the default value a, similar as in Sterner
"""
This damage function is based on Howard and Sterner (2017):
``D(T)=aT^2``
    
*Input:*
- ``a`` : Damage coefficient. By default, it is set to 0.01145. The user must write the value of a explicitely, if he wants to use the dam_Sterner function with another ``a``.
- ``T``: temperature
    
*Syntax:*
```julia
ACE_Traeger_replication.dam_Sterner(T, a)
```
    
*Output:*
- The calculated DICE-Howard-Sterner damage estimate.

*Example*
```julia-repl
julia> ACE_Traeger_replication.dam_Sterner(3)
0.10305
```
"""
function dam_Sterner(T, a=0.01145)
    x= a .* T.^2
    return x
end
    
#For Dam_Weitz (in the initial code but not used to generate results)
function dam_Weitz(T)
    x=1 .- 1 ./ (1 .+ (T ./ 20.46).^2 .+ (T ./ 6.081).^6.754)
    return x
end

#For ACE:
"""
The damage function in ACE takes an exponential functional form: 
- ``D(T)=1-exp(-ξ0*exp(ξ1*T)+ξ0)``
with ``ξ0``, a free damage parameters and the paramteter ``ξ1=log(2)/s`` pinned down by climate sensitivity s.
    
*Input:*
- ``ξ0``: xi0 (Calibrated to 0.022 in the paper. The base calibration is an exact match of the two calibration points 0° and 2.5° in the 2007 model)
- ``ξ1``: xi1 (Estimated to ``0.25`` in the paper)
- ``T``: temperature
    
*Syntax:*
```julia
ACE_Traeger_replication.dam_ACE(xi0, xi1, T)
```
*Output:*
- The calculated ACE damage estimate.

*Example:*
```julia-repl
julia> ACE_Traeger_replication.dam_ACE(0.022, 0.25, 3)
0.024274517795511485
```
"""
function dam_ACE(xi0, xi1, T)
    x=1 .- exp.(-xi0 * exp.(xi1 .* T) .+ xi0)
    return x
end
    
##Damage functions: 
#function squeeze(A::Vector{Float64})
#    singleton_dims = tuple((d for d in 1:ndims(vec(A)) if size(vec(A), d) == 1)...)
#    return squeeze(vec(A), singleton_dims)
#end


############################################### FIGURE 2

"""
The function below replicates original figure II in the paper. This figure shows the predicted damages of different models depending on the temperature degrees above the preindustrial level.

The different models are:
- DICE 2007, 2013, 2016
- ACE base damage calibration. It matches the two calibration points 0 and 2.5° in the 2007 model. 
- HSP-norm: Howard and Sterner's damage function using DICE's approach to limit damages to 100 percents
- HSP-nn: Howard and Sterner's damage function using DICE's approach, not normalized to limit damages to 100 percents. Damages exceed production at 9.5°. 
    
The left side of the plot refers to the temperature range of the IPCC secenarios in Figure 3 and the right side is focuses on lower degrees of warming. 
    
*Input:*
- path
    
*Syntax:*

```julia
ACE_Traeger_replication.Damage_function_plot(path)
```
    
*Output:*
- Figure 2 saved into path
"""
function Damage_function_plot(path=pwd())
    path=path   
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
    savefig(h, path*"/Damage_function_plot.pdf")
end

###############################################

############################################### FIGURE 3
"""
This function replicates Figure III of the paper. The authors use the MAGICC6.0 model by Meinshausen, Raper and Wigley(2011) to simulate
the Representative Concentration Pathway (RCP) scenario over a time horizon of 400 years. The calibration of ACE uses 2 ocean layers(upper and deep) compared to MAGGIC's 50 layers and DICE's ocean layers.
The figure shows ACE's temperature response from present times to 2400 compared to MAGICC6.0 using the color-coded radiative forcing scenarios of the latest IPCC assessment report.
    
Note that the dataframe MagiccOcean.mat and TempDataCalel.mat are required to run this function. Please save these dataframes to your datapath of choice (default is current working directory).        

Make sure you have the following dataframes:
- MagiccOcean.mat
- TempDataCalel.mat
        
*Input:*
- path: The user should provide their preferred path to save the figure. Otherwise the function takes the working directory by default.
- datapath: The user should provide the location of the data needed for the figure (found in replication package). Otherwise the function takes the current working directory as default. 

*Necessary dataframes are:*
- MagiccOcean.mat
- TempDataCalel.mat
        
*Syntax:*
```julia
ACE_Traeger_replication.TempFitSim(path, datapath)
```
        
*Output:*
- Figure 3 saved into path
"""
function TempFitSim(path = pwd(), datapath = pwd())
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
    data=matread(datapath*"//MagiccOcean.mat")


    #Get the name of all the variables: MagiccOcean is a Dictionnary with Dictionnaries inside
    scenario = collect(keys(data["MagiccOcean"]))
    S=length(scenario) #Number of scenarios

    #define the first and the end dates
    startpos = findfirst(data["MagiccOcean"][scenario[1]]["Year"].==Startdate)[1]
    endpos = findfirst(data["MagiccOcean"][scenario[1]]["Year"].==Enddate)[1]

    #Check that the index for the first and the last dates in the data are consistent
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
    d=matread(datapath*"//TempDataCalel.mat")

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

    savefig(g, path*"/TempFitSim.pdf")

end

###############################################

############################################### TEMPSIMULATION

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
    end


    ##Evaluating
    diff_sum=0

    for l = 1:layers
        diff = weight[l] .* (TempSim[l, :, :] .- Temp[l, :, :]) .^ 2
        diff_sum += sum(diff)
    end 

    # Set err to diff_sum
    err = Float64(diff_sum)

    # Display the result
    println(" Residual ", err, " σ-vector ", σ)
    println()  # Print an empty line

    return Dict("err"=>err, "Tempsim"=>TempSim, "xi"=>xi)

end

###############################################

############################################### SUB-FIGURE 4
"""
The function below replicates a part of figure 4 in the paper. This figure shows temperature increases (in degree celsius) resulting from 
    a one-time release of certain amount of gigatons of carbon across different models. Users can specify the amount of carbon released into the atmosphere. 
    The default is 100 GtC as in the paper.

The different models are:
- ACE-DICE ACE using DICE's carbon cycle
- ACE-Joos ACE using Joos et al. (2013) impulse-response model
- Dynamic Integrated Climate-Economy (DICE) 2013, 2016
    
*Input:*
- impulse: The user can provide an impulse of their choice (in gigatons of carbon-equivalent)
- path: The user should provide their preferred path to save the figure. Otherwise the function takes the working directory by default.
    
*Syntax:*

```julia
ACE_Traeger_replication.Impulse_response(pulse, path)
```
    
*Output:*
- figure saved into path
"""
function Impulse_response(pulse = 100, path=pwd())

    logic_plot=1 # turns plotting on

    logic_timing_pulse = 0;
       # =0: adds pulse to M stock at beginning of present period (default)
    logic_emission_decay = 0;  
       # =0: Emissions added straight to atmosphere (default)
    logic_forcing_delay = 1;
       # =1: T_t+1(M_t+1) DICE and ACE preferred timing   (default)


    #Parameters
    horizon=210   # years
    pulse = pulse  # Carbon pulse in GtC
       
    # Select versions to calculate, save, and plot (ACE-DICE always on)
    logic_ACE_Joos=1 # Includes Joosetal 2013
    logic_DICE=1   # Includes DICE 2013
    logic_DICE16=1 # Includes DICE 2016

        # Joos Decay boxes
            if logic_ACE_Joos == 1
                #  Joos (2013) Carbon impulse response
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

                CarbonMatrixJoos = copy(CarbonMatrix)
            end
    
    param_timestep = timestep  # needed in TMDynamics for scaling of carbon cycle dynamics

    # TMDynamics
    cs=3
    Mpre=596.4
    param = Dict()
    param["xi"] = [0.231049060186648, 0.231049060186648, 0.231049060186648]
    sigma_temp_up = [0.538194275028702]

    TempMatrix = [0.000000027559817 0.461805697411481 0;
              0.076519839162246 0.895993627259103 0.027486533578651;
              0 0.003992679633153 0.996007320366847]

    param_xi1= log(2) / cs #Approximately 0.23

    if param_xi1-param["xi"][1]>1e-10
        print("Inconsistency with xi1 definition, check underlying cs")
    end

    #RCP 6 2015 values for possible initialization (not required for paper):
    Initial = Dict("Temperature" => [0.98901461, 0.73635941, 0.73635941],
               "tau" => [1.256727218547292, 1.185465055669960, 1.185465055669960],
               "M" => [818.985, 1527, 10010])

    #Equilibrium temperature increase of given specification:
    tau_eq_vec =(Matrix(1.0I, length(Initial["Temperature"]), length(Initial["Temperature"]))-TempMatrix)^(-1)*2*[sigma_temp_up[1]; zeros(length(Initial["Temperature"])-1)]
    #and translating back into temperature in degree Celsius:
    T_eq_vec=log.(tau_eq_vec)./param["xi"]
    print("Test of long-run equilibrium Temp at 2xCO2", round.(Int, T_eq_vec))

    #Parameters Carbon Cycle Matrix
    sigma_carb = [0.088, 0.0025, 0.03832888889, 3.3750e-04]  # DICE
    sigma_carb = sigma_carb * timestep / 5  # rescales to 10 year time step. Extremely similar to 5 year

    # Carbon Cycle Matrix
    CarbonMatrix = [1 - sigma_carb[1]  sigma_carb[3]  0; 
                sigma_carb[1]  1 - sigma_carb[2] - sigma_carb[3]  sigma_carb[4];
                0  sigma_carb[2]  1 - sigma_carb[4]]


    xi = param["xi"]
    sigma_forc = sigma_temp_up[1]
    
    # Path Initializations
    hori = ceil(Int, horizon / timestep)
    horizon_legend = (0:hori-1) * timestep
    M = zeros(3, hori)
    tau = zeros(3, hori)
    tau_pulse = zeros(3, hori)
    tau_Joos = zeros(3, hori)
    Temp = zeros(3, hori)
    Temp_pulse = zeros(3, hori)
    Temp_Joos = zeros(3, hori)
    M_pulse = zeros(3, hori)
    M[:, 1] = Initial["M"]
    tau[:, 1] = Initial["tau"]
    tau_pulse[:, 1] = Initial["tau"]  # initial temp for pulse experiment using DICE 2013 carbon cycle
    tau_Joos[:, 1] = Initial["tau"]   # initial temp for pulse experiment using Joos et al 2013
    Temp[:, 1] = diagm(1 ./ xi) * log.(tau[:, 1])
    Temp_pulse[:, 1] = diagm(1 ./ xi) * log.(tau[:, 1])
    Temp_Joos[:, 1] = Temp[:, 1]
    M_pulse[:, 1] = Initial["M"] .+ [pulse, 0, 0]

    # DICE 2013 Dynamic Variables
    ParaDice = Dict()
    ParaDice["Mpre"] = 588
    ParaDice["climate_sens"] = 3 
    println("Climate sensitivity is ", ParaDice["climate_sens"])
    println("DICE13 itself would use 3.2")
    ParaDice["eta"] = 3.8
    ParaDice["Initial_M"] = [830.4 1527 10010]  # from GAMS. In EXCEL [818.985 1527 10010].

    # Parameters Carbon Cycle Matrix:
    sigma_carb = [0.088, 0.0025, 0.03832888889, 3.3750e-04]
    sigma_carb .= sigma_carb .* timestep / 5  # rescales time step based on DICE's 5 years.
    ParaDice["CarbonMatrix"] = [1 - sigma_carb[1] sigma_carb[3] 0;
                         sigma_carb[1] 1 - sigma_carb[2] - sigma_carb[3] sigma_carb[4];
                         0 sigma_carb[2] 1 - sigma_carb[4]]

    # Temperature Dynamics
    ParaDice["Initial_T"] = [0.80 0.0068]  # GAMS
    ParaDice["c1"] = 0.098
    ParaDice["c3"] = 0.088
    ParaDice["c4"] = 0.025

    sigma_tempDICE = zeros(3)
    sigma_tempDICE[1] = ParaDice["c1"] * ParaDice["eta"] / ParaDice["climate_sens"]
    ParaDice["sigma_tempDICE"] = copy(sigma_tempDICE)
    sigma_tempDICE[2] = ParaDice["c1"] * ParaDice["c3"]
    sigma_tempDICE[3] = ParaDice["c4"]

    if timestep != 5
        println("Using timestep not equal to 5, DICE results will be off")
    end

    sigma_tempDICE .= 1 .- (1 .- sigma_tempDICE) .^ (timestep / 5)  # rescales time step based on DICE's 5 years.
    ParaDice["TempDICEMatrix"] = [1 - sigma_tempDICE[1] - sigma_tempDICE[2] sigma_tempDICE[2];
                            sigma_tempDICE[3] 1 - sigma_tempDICE[3]]


    # DICE 2016 Dynamic Variables
    if logic_DICE16 == 1

    ParaDice16 = Dict()    
    # Main parameters
    ParaDice16["eta"] = 3.6813
    ParaDice16["Mpre"] = 588
    ParaDice16["climate_sens"] = 3.1
    
    # Carbon Dynamics DICE 2016
    ParaDice16["Initial_M"] = [851, 460, 1740]  # GtC from GAMS 2016
    # Parameters for Carbon Cycle Matrix:
    ParaDice16["Eq16"] = [588, 360, 720]  # Equilibrium concentration atmosphere (GtC) / upper strata (GtC) / lower strata (GtC)
    ParaDice16["b12"] = 0.12
    ParaDice16["b23"] = 0.007
    ParaDice16["b11"] = 1 - ParaDice16["b12"]
    ParaDice16["b21"] = ParaDice16["b12"] * ParaDice16["Eq16"][1] / ParaDice16["Eq16"][2]
    ParaDice16["b22"] = 1 - ParaDice16["b21"] - ParaDice16["b23"]
    ParaDice16["b32"] = ParaDice16["b23"] * ParaDice16["Eq16"][2] / ParaDice16["Eq16"][3]
    ParaDice16["b33"] = 1 - ParaDice16["b32"]
    # transformations for time scaling
    sigma_carb16 = [1 - ParaDice16["b11"], ParaDice16["b23"], ParaDice16["b21"], ParaDice16["b32"]]  # DICE
    sigma_carb16 .= sigma_carb16 .* timestep / 5  # rescales time step based on DICE's 5 years.
    ParaDice16["CarbonMatrix"] = [1 - sigma_carb16[1] sigma_carb16[3] 0;
                               sigma_carb16[1] 1 - sigma_carb16[2] - sigma_carb16[3] sigma_carb16[4];
                               0 sigma_carb16[2] 1 - sigma_carb16[4]]
    # Tempdynamics DICE 16 GAMS
    ParaDice16["Initial_T"] = [0.85 0.0068]
    ParaDice16["c1"] = 0.1005
    ParaDice16["c3"] = 0.088
    ParaDice16["c4"] = 0.025
    ParaDice16["climate_sens"] = 3.1
    ParaDice16["eta"] = 3.6813
    sigma_tempDICE16 = zeros(3)
    sigma_tempDICE16[1] = ParaDice16["c1"] * ParaDice16["eta"] / ParaDice16["climate_sens"]
    ParaDice16["sigma_tempDICE"] = copy(sigma_tempDICE16)
    sigma_tempDICE16[2] = ParaDice16["c1"] * ParaDice16["c3"]
    sigma_tempDICE16[3] = ParaDice16["c4"]
    # Timestep-scaling does not work for temperature in DICE (here alternative version to above - if you see how to do it right let me know)
    sigma_tempDICE16 .= 1 .- (1 .- sigma_tempDICE16) .^ (timestep / 5)
    ParaDice16["TempDICEMatrix"] = [1 - sigma_tempDICE16[1] - sigma_tempDICE16[2] sigma_tempDICE16[2];
                                 sigma_tempDICE16[3] 1 - sigma_tempDICE16[3]]
    end

    if logic_DICE == 1
        Tvec = zeros(2, hori)
        Tvec_pulse = zeros(2, hori)
        Temp_DICE = zeros(1, hori)
        Temp_DICE_pulse = zeros(1, hori)

        Tvec[:, 1] = ParaDice["Initial_T"]
        Tvec_pulse[:, 1] = ParaDice["Initial_T"]
        Temp_DICE[1] = Tvec[1, 1]
        Temp_DICE_pulse[1] = Tvec_pulse[1, 1]
    end

    if logic_DICE16 == 1
        M16 = zeros(3, hori)
        M_pulse16 = zeros(3, hori)
        Tvec16 = zeros(2, hori)
        Tvec_pulse16 = zeros(2, hori)
        Temp_DICE16 = zeros(1, hori)
        Temp_DICE_pulse16 = zeros(1, hori)

        M16[:, 1] = ParaDice16["Initial_M"]
        M_pulse16[:, 1] = ParaDice16["Initial_M"] .+ [pulse, 0, 0]
        Tvec16[:, 1] = ParaDice16["Initial_T"]
        Tvec_pulse16[:, 1] = ParaDice16["Initial_T"]
        Temp_DICE16[1] = Tvec16[1, 1]
        Temp_DICE_pulse16[1] = Tvec_pulse16[1, 1]
    end

    for t = 1:hori - 1
        # Carbon cycle:
        # Calculate emissions that keep concentration constant. Equation equivalent to CarbonMatrix*M+E=M.
        if logic_emission_decay == 0
            Evec = (diagm(ones(3)) - CarbonMatrix) * M[:, t]
        elseif logic_emission_decay == 1  # carbon cycle already hits current emissions.
            Evec = (diagm(ones(3)) - CarbonMatrix) * M[:, t] / CarbonMatrix[1, 1]
        end

        # Picking only atmospheric part keeps atm concentration fix, but not ocean (saturates)
        E = zeros(hori)
        E[t] = Evec[1]
        if logic_timing_pulse == 0
            M[:, t + 1] = CarbonMatrix * M[:, t] + [E[t], 0, 0]
            M_pulse[:, t + 1] = CarbonMatrix * M_pulse[:, t] + [E[t], 0, 0]
            elseif logic_timing_pulse == 1  # reset pulse stock in present to baseline stock
                if logic_emission_decay == 0  # all of pulse appears next period
                M[:, t + 1] = CarbonMatrix * M[:, t] + [E[t], 0, 0]
                    if t == 1
                        M_pulse[:, t] = M[:, t]
                        M_pulse[:, t + 1] = CarbonMatrix * M[:, t] + [pulse + E[t], 0, 0]
                    else
                    M_pulse[:, t + 1] = CarbonMatrix * M_pulse[:, t] + [E[t], 0, 0]
                end
            elseif logic_emission_decay == 1  # carbon cycle already hits current emissions.
                M[:, t + 1] = CarbonMatrix * (M[:, t] + [E[t], 0, 0])
                    if t == 1
                        M_pulse[:, t] = M[:, t]
                        M_pulse[:, t + 1] = CarbonMatrix * (M[:, t] + [pulse + E[t], 0, 0])
                    else
                M_pulse[:, t + 1] = CarbonMatrix * (M_pulse[:, t] + [E[t], 0, 0])
            end
        end
    else
        error("option for timing not available")
    end

    # Temperature evolution:
    forc = zeros(hori)
    forc_pulse = zeros(hori)
    if logic_forcing_delay == 1
        forc[t] = sigma_forc * M[1, t + 1] / Mpre
    else
        forc[t] = sigma_forc * M[1, t] / Mpre
    end
    tau[:, t + 1] = TempMatrix * tau[:, t] + [forc[t], 0, 0]
    Temp[:, t + 1] = diagm(1.0 ./ xi) * log.(tau[:, t + 1])

    if logic_forcing_delay == 1
        forc_pulse[t] = sigma_forc * M_pulse[1, t + 1] / Mpre
    else
        forc_pulse[t] = sigma_forc * M_pulse[1, t] / Mpre
    end
    tau_pulse[:, t + 1] = TempMatrix * tau_pulse[:, t] + [forc_pulse[t], 0, 0]
    Temp_pulse[:, t + 1] = diagm(1.0 ./ xi) * log.(tau_pulse[:, t + 1])

    if logic_ACE_Joos == 1
        Imp = zeros(hori)
        forc_Joos = zeros(hori)
        # Joos Impulse response
        Imp[t + 1] = sum(CarbonMatrixJoos^(t) * CarbonWeights) * pulse
        forc_Joos[t] = sigma_forc * (M[1, t] + Imp[t + 1]) / Mpre
        tau_Joos[:, t + 1] = TempMatrix * tau_Joos[:, t] + [forc_Joos[t], 0, 0]
        Temp_Joos[:, t + 1] = diagm(1.0 ./ xi) * log.(tau_Joos[:, t + 1])
    end

    if logic_DICE == 1
        Fback = zeros(hori)
        Fpul = zeros(hori)
        Fback[t] = ParaDice["eta"] * (log(M[1, t + 1]) - log(ParaDice["Mpre"])) / log(2)  
        Tvec[:, t + 1] = [1 0]' * ParaDice["climate_sens"] * ParaDice["sigma_tempDICE"][1] / ParaDice["eta"] * (ParaDice["eta"] * (log(M[1, t + 1]) - log(ParaDice["Mpre"])) / log(2)) + ParaDice["TempDICEMatrix"] * Tvec[:, t]
        Temp_DICE[t + 1] = Tvec[1, t + 1]
        Fpul[t] = ParaDice["eta"] * (log(M_pulse[1, t + 1]) - log(ParaDice["Mpre"])) / log(2) 
        Tvec_pulse[:, t + 1] = [1 0]' * ParaDice["climate_sens"] * ParaDice["sigma_tempDICE"][1] / ParaDice["eta"] * (ParaDice["eta"] * (log(M_pulse[1, t + 1]) - log(ParaDice["Mpre"])) / log(2)) + ParaDice["TempDICEMatrix"] * Tvec_pulse[:, t]
        Temp_DICE_pulse[t + 1] = Tvec_pulse[1, t + 1]
    end
    
    if logic_DICE16 == 1 || logic_DICE == 1
        # Carbon cycle:
        # Calculate emissions that keep concentration constant.
        if logic_emission_decay == 0
            Evec16 = (diagm(ones(3)) - ParaDice16["CarbonMatrix"]) * M16[:, t]
        elseif logic_emission_decay == 1  # carbon cycle already hits current emissions.
            Evec16 = (diagm(ones(3)) - ParaDice16["CarbonMatrix"]) * (M16[:, t] / ParaDice16["CarbonMatrix"][1, 1])
        end
        # Picking only atmospheric part keeps atm concentration fix, but not ocean (saturates)
        E16 = zeros(hori)
        E16[t] = Evec16[1]
        if logic_timing_pulse == 0
            M16[:, t + 1] = ParaDice16["CarbonMatrix"] * M16[:, t] + [E16[t], 0, 0]
            M_pulse16[:, t + 1] = ParaDice16["CarbonMatrix"] * M_pulse16[:, t] + [E16[t], 0, 0]
        elseif logic_timing_pulse == 1     # reset pulse stock in present to baseline stock
            if logic_emission_decay == 0    # all of pulse appears next period
                M16[:, t + 1] = ParaDice16["CarbonMatrix"] * M16[:, t] + [E16[t], 0, 0]
                if t == 1
                    M_pulse16[:, t] = M16[:, t]
                    M_pulse16[:, t + 1] = ParaDice16["CarbonMatrix"] * M16[:, t] + [pulse + E16[t], 0, 0]
                else
                    M_pulse16[:, t + 1] = ParaDice16["CarbonMatrix"] * M_pulse16[:, t] + [E16[t], 0, 0]
                end
            elseif logic_emission_decay == 1  # carbon cycle already hits current emissions.
                M16[:, t + 1] = ParaDice16["CarbonMatrix"] * (M16[:, t] + [E16[t], 0, 0])
                if t == 1
                    M_pulse16[:, t] = M16[:, t]
                    M_pulse16[:, t + 1] = ParaDice16["CarbonMatrix"] * (M16[:, t] + [pulse + E16[t], 0, 0])
                else
                    M_pulse16[:, t + 1] = ParaDice16["CarbonMatrix"] * (M_pulse16[:, t] + [E16[t], 0, 0])
                end
            end
        else
            error("option for timing not available")
        end
    end

    # Temperature evolution
    if logic_DICE16 == 1
        Fback16 = zeros(hori)
        Fpul16 = zeros(hori)
        Fback16[t] = ParaDice16["eta"] * (log(M16[1, t + 1]) - log(ParaDice16["Mpre"])) / log(2)
        Tvec16[:, t + 1] = [1 0]' * ParaDice16["climate_sens"] * ParaDice16["sigma_tempDICE"][1] / ParaDice16["eta"] * (ParaDice16["eta"] * (log(M16[1, t + 1]) - log(ParaDice16["Mpre"])) / log(2)) + ParaDice16["TempDICEMatrix"] * Tvec16[:, t]
        Temp_DICE16[t + 1] = Tvec16[1, t + 1]
        Fpul16[t] = ParaDice16["eta"] * (log(M_pulse16[1, t + 1]) - log(ParaDice16["Mpre"])) / log(2)
        Tvec_pulse16[:, t + 1] = [1 0]' * ParaDice16["climate_sens"] * ParaDice16["sigma_tempDICE"][1] / ParaDice16["eta"] * (ParaDice16["eta"] * (log(M_pulse16[1, t + 1]) - log(ParaDice16["Mpre"])) / log(2)) + ParaDice16["TempDICEMatrix"] * Tvec_pulse16[:, t]
        Temp_DICE_pulse16[t + 1] = Tvec_pulse16[1, t + 1]
    end
    

    end

    if logic_plot == 1
    
        linewidth = 1.75
        titlefontsize = 18

        h1 = plot(horizon_legend, Temp_pulse[1, 1:hori] .- Temp[1, 1:hori], linewidth=linewidth, label = "ACE-DICE")
        if logic_ACE_Joos == 1
            plot!(horizon_legend, Temp_Joos[1, 1:hori] .- Temp[1, 1:hori], linewidth=linewidth, label = "ACE-Joos")
        end
        if logic_DICE == 1
            plot!(horizon_legend, Temp_DICE_pulse[1, 1:hori] .- Temp_DICE[1, 1:hori], linewidth=linewidth, label = "DICE 2013")
        end
        if logic_DICE16 == 1
            plot!(horizon_legend, Temp_DICE_pulse16[1, 1:hori] .- Temp_DICE16[1, 1:hori], linewidth=linewidth, label = "DICE 2016")
        end
        title!("Temperature Impulse Response to $(pulse)GtC", fontsize=titlefontsize)
        xlabel!("Year")
        ylabel!("Temperature Diff in C")
        xticks!(0:20:horizon)
        savefig(h1, path*"/ImpulseResponse_Temperature_$(pulse)GtC_timestep_$(timestep).pdf")
    
    end

end    

############################################### Figure 4
"""
The function below replicates figure 4 in the paper. This figure shows temperature increases (in degree celsius) resulting from 
    a one-time release of 100 of gigatons of carbon (GtC) across different models.

The different models are:
- ACE-DICE ACE using DICE's carbon cycle
- ACE-Joos ACE using Joos et al. (2013) impulse-response model
- Dynamic Integrated Climate-Economy (DICE) 2013, 2016
- Climate Framework for Uncertainty, Negotiation and Distribution (FUND)
- Policy Analysis of the Greenhouse Effect (PAGE)
- Coupled Model Intercomparison Project Phase 5 (CMIP5) 3.0 and 3.25 (DPRV)
- Bern present day (PD)
- Bern pre-industrial (PI)
- Transient climate response to cumulative carbon emissions model (TCRE)
    
Make sure you have the following dataframes:
    - impulse_timestep_5_logi_001.mat
    - impulse_timestep_1_logi_001.mat    
    - Venmans_CS3_lam_1p06.mat
    - Impulse_Response_Bern2p5_Joos.mat
    - Impulse_Response_Bern2p5_PD_Joos.mat
    - Venmans_CS3p1.mat

*Input:*
- path: The user should provide their preferred path to save the figure. Otherwise the function takes the current working directory by default.
- datapath: The user should provide the location of the data needed for the figure (found in replication package). Otherwise the function takes the current working directory as default. 
    Necessary dataframes are:
        - impulse_timestep_5_logi_001.mat
        - impulse_timestep_1_logi_001.mat    
        - Venmans_CS3_lam_1p06.mat
        - Impulse_Response_Bern2p5_Joos.mat
        - Impulse_Response_Bern2p5_PD_Joos.mat
        - Venmans_CS3p1.mat
*Syntax:*

```julia
ACE_Traeger_replication.Impulse_response_combined(path, datapath)
```
    
*Output:*
- Figure 4 saved into path
"""
function Impulse_response_combined(path = pwd(), datapath = pwd())

    horizon=210
    pulse = 100
    logic_plot = 1

    logi = "001" # 001 is standard timing: pulse added to M0, no emission decay during period, T_t+1(M_t+1), gives same result as 111
    
    horizon_cutoff = 160
    logic_ACE_DICE=1 
    logic_ACE_Joos=1 # Includes ACE-Joos based on Joos et al 2013
    logic_DICE=1 # Includes DICE 2013 
    logic_DICE16=1 # Includes DICE 2016 
    logic_VenBest=1 # Includes Venman's best ensemble fit (FairGeoffrey)
    logic_VenBest_alt=1 # takes CMIP5 climate sensitivity average and inverts to get lambda. Then reduces F2CO2 to get CS 0f 3. 
    logic_VenFUND=1 # Includes Venman-Dietz's FUND
    logic_VenPAGE=1 # Includes Venman-Dietz's PAGE
    
    # Timesteps for ACE with 1 and 5 years:
    logic_ACE_Joos1=1 
    logic_ACE_DICE1=0  
    logic_Bern2p5_Joos_PI = 1
    logic_Bern2p5_Joos_PD = 1
    logic_TCRE = 1
    
    # Take DICE 2013 and 16 from 5 year time step
    
    data1 = matread(datapath*"/impulse_timestep_5_logi_$(logi).mat")
    
    M13 = data1["M"]
    M_pulse13 = data1["M_pulse"]
    DICE_horizon = data1["horizon_legend"]
    DICE_Impulse = data1["Temp_DICE_pulse"][1,1:trunc(Int, data1["hori"])]-data1["Temp_DICE"][1,1:trunc(Int, data1["hori"])]
    DICE_Impulse16 = data1["Temp_DICE_pulse16"][1,1:trunc(Int, data1["hori"])]-data1["Temp_DICE16"][1,1:trunc(Int, data1["hori"])]
    ACE_horizon5 = data1["horizon_legend"]
    ACE_Impulse_DICE5 = data1["Temp_pulse"][1,1:trunc(Int, data1["hori"])]-reshape(data1["Temp"][1,:,:], 6*39)[1:trunc(Int, data1["hori"])]
    ACE_Impulse_Joos5 = data1["Temp_Joos"][1,1:trunc(Int, data1["hori"])]-reshape(data1["Temp"][1,:,:], 6*39)[1:trunc(Int, data1["hori"])]
    pulse5y = data1["pulse"]
    
    # Take 1 year time step version of ACE
    data2 = matread(datapath*"/impulse_timestep_1_logi_$(logi).mat")
    ACE_horizon1 = data2["horizon_legend"]
    ACE_Impulse_DICE1 = data2["Temp_pulse"][1,1:trunc(Int, data2["hori"])]-reshape(data2["Temp"][1,:,:], 6*39)[1:trunc(Int, data2["hori"])]
    ACE_Impulse_Joos1 = data2["Temp_Joos"][1,1:trunc(Int, data2["hori"])]-reshape(data2["Temp"][1,:,:], 6*39)[1:trunc(Int, data2["hori"])]
    pulse1y = data2["pulse"]
    if pulse1y != pulse5y
        error("loading models where 1 year pulse differs from 5 year pulse")
    end
    
    if logic_VenBest_alt==1
        Ven = matread(datapath*"/Venmans_CS3_lam_1p06.mat")
        JGbest_alt = Ven["Ven"]["JGbest"]
    end 
    
    data3 = matread(datapath*"/impulse_timestep_1_logi_$(logi).mat")
    ACE_horizon = data3["horizon_legend"]
    ACE_Impulse_DICE = data3["Temp_pulse"][1,1:trunc(Int, data3["hori"])]-reshape(data3["Temp"][1,:,:], 6*39)[1:trunc(Int, data3["hori"])]
    ACE_Impulse_Joos = data3["Temp_Joos"][1,1:trunc(Int, data3["hori"])]-reshape(data3["Temp"][1,:,:], 6*39)[1:trunc(Int, data3["hori"])]
    if data3["pulse"] != pulse5y
        error("loading 10 year models where pulse differs from 5 year pulse")
    end
    
    # FROM AR6: 
    TCRE_best = 1.65*data3["pulse"]/1000
    # likely interval (by guidance notes AR5) 
    TCRE_high = 2.3*data3["pulse"]/1000
    TCRE_low = 1*data3["pulse"]/1000
    
    if logic_Bern2p5_Joos_PI==1
        data4 = matread(datapath*"/Impulse_Response_Bern2p5_Joos.mat")   # preindustrial, by courtesy of Fortunat Joos
        Impulse_Response_Bern2p5 = data4["Impulse_Response_Bern2p5"][1:201]
        horizon_Bern = 0.5:200.5
    end
    
    if logic_Bern2p5_Joos_PD==1
        data5 = matread(datapath*"/Impulse_Response_Bern2p5_PD_Joos.mat")  # present day, by courtesy of Fortunat Joos
        Impulse_Response_Bern2p5_PD = data5["Impulse_Response_Bern2p5_Joos"][1:201]
        horizon_Bern = 0.5:200.5
    end
    
    data6 = matread(datapath*"/Venmans_CS3p1.mat")  # Loads time paths VenBest courtesy of Frank Venmans.
    

    if logic_plot == 1
        
        linewidth = 1.75
        titlefontsize = 18
    
        h3 = plot()
        if logic_ACE_Joos == 1
            plot!(0:maximum(ACE_horizon), ACE_Impulse_Joos, linewidth=linewidth, label = "ACE-Joos")
        end
        if logic_ACE_DICE == 1
            plot!(0:maximum(ACE_horizon), ACE_Impulse_DICE, linewidth=linewidth, label = "ACE-DICE")
        end
        if logic_ACE_Joos1 == 1
            plot!(0:maximum(ACE_horizon1), ACE_Impulse_Joos1, linewidth=linewidth, label = "ACE-Joos-1y")
        end
        if logic_ACE_DICE1 ==1
            plot!(0:maximum(ACE_horizon1), ACE_Impulse_DICE1, linewidth=linewidth, label = "ACE-DICE-1")
        end
        if logic_DICE == 1
            plot!(0:5:maximum(DICE_horizon), DICE_Impulse, linewidth=linewidth, label = "DICE 2013")
        end
        if logic_DICE16 == 1
           plot!(0:5:maximum(DICE_horizon), DICE_Impulse16, linewidth=linewidth, label = "DICE 2016")
        end
        if logic_VenBest ==1
            plot!(0:maximum(vec(Ven["Ven"]["year"])), vec(Ven["Ven"]["JGbest"]), linewidth=linewidth, label = "CMIP5 DPRV")
        end
        if logic_VenBest_alt ==1
            plot!(0:maximum(vec(Ven["Ven"]["year"])), vec(JGbest_alt), linewidth=linewidth, label = "CMIP5 3.0*")
        end
        if logic_Bern2p5_Joos_PD==1
            plot!(horizon_Bern, Impulse_Response_Bern2p5_PD, linewidth=linewidth, label = "Bern PD")
        end
        if logic_Bern2p5_Joos_PI==1
            plot!(horizon_Bern, Impulse_Response_Bern2p5, linewidth=linewidth, label = "Bern PI")
        end
        if logic_VenFUND ==1
            plot!(0:maximum(vec(Ven["Ven"]["year"])), vec(Ven["Ven"]["FUND"]), linewidth=linewidth, label = "FUND")
        end
        if logic_VenPAGE ==1
            plot!(0:maximum(vec(Ven["Ven"]["year"])), vec(Ven["Ven"]["PAGE"]), linewidth=linewidth, label = "PAGE")
        end
        if logic_TCRE == 1
            hline!((TCRE_best, TCRE_best), linewidth=linewidth, linestyle=:dot, label = "TCRE-best", color = "black")
            hline!((TCRE_low, TCRE_low), linewidth=linewidth, linestyle=:dot, label = "TCRE-likely", color = "grey")
            hline!((TCRE_high, TCRE_high), linewidth=linewidth, linestyle=:dot, color = "grey", label = "")
        end
        title!("Temperature Impulse Response to 100GtC", fontsize=titlefontsize)
        xlabel!("Year")
        ylabel!("Temperature Difference in C")
        xticks!(0:20:horizon)
        savefig(h3, path*"/ImpulseResponse_Temp_Complot_$(pulse)GtC.pdf")
    
    end
    

end

###############################################

############################################### TABLE 1 (SCC)


"""
This function replicates Table I of the paper. The table shows the social cost of carbon (SCC) estimates based on a variety of assumptions. 
It also provides information on the cent per gallon and euro per liter cost of gasoline.
These assumptions are included in the scenarios and can take the following values:
 - various annual discount factors (rho).
 - damages can be calibrated to the DICE model (1) or damages are  calibrated to Howard-Sterner (2017) and Pyndick (2020) models (2)
 - the scenario can be based on the use of the carbon cycle matrix given by DICE or the carbon cycle matrix given by Joos et al. (2013)
 - the scenario can include population weighting or not
 - the scenario can include a capital share in production (kappa) of value 0.3 (stylized fact) or 0.4 (more recent estimates) as well as calibrated ones.

Note that population recalibration is not implemented in this replication, therefore 3 scenarios are missing.

*Input:*
- path

*Syntax:*
```julia
ACE_Traeger_replication.SCC(path)
```

*Output:*
- Table 1 saved into path
"""
function SCC(path = pwd())

# Initialize globals
global replicate_table = 1
global tab = Dict()
global param = Dict()


# Initialize table entries
tab["scen"] = zeros(40)
tab["carb_mult"] = zeros(40)
tab["wo_delay"] = zeros(40)
tab["SCC"] = zeros(40)
tab["cent_gallon"] = zeros(40)
tab["cent_liter"] = zeros(40)
tab["rho_discount_annual_exo"] = zeros(40)
tab["damages"] = zeros(40)
tab["boxmodel"] = zeros(40)
tab["population"] = zeros(40)
tab["pop_recalibrate"] = zeros(40)
tab["kappa"] = zeros(40)
tab["kappa_calib"] = zeros(40)

# List of scenarios
global list = Dict()

    list["rho_discount_annual_exo"] = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01, 0.01,0.01,0.01,0.01, 
                                      0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.001,0.001,
                                      0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001]

    list["damages"] = [1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2]
    
    list["boxmodel"] = [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]
    list["population"] = [0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1]
    list["pop_recalibrate"] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    list["kappa"] = [0.3,  0.3,  0.3,  0.4,  0.4,  0.3,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.3,  0.3,  0.3,  0.4,  0.3,   0.3,  0.4,  0.3,  0.4,  0.3,  0.4,  0.4,  0.3,  0.3,  0.3,  0.4,  0.3,  0.3,  0.4,  0.3,  0.4,  0.3,  0.4,  0.4]
    list["kappa_calib"] = [0.3,  0.3,  0.3,  0.4,  0.3,  0.3,  0.4,  0.4,  0.3,  0.4,  0.3,  0.4,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3]


# loop over the different scenarios
for k in 1:length(list["damages"])
    println("Preparing scerario $k")
    global boxmodel = list["boxmodel"][k]
    global population = list["population"][k]
    global jj = k

    # If scenario uses standard carbon matrix and no population weighting
    if boxmodel == 0 && population == 0
        include(joinpath(@__DIR__, "ACE_SCC_Deterministic_boxmodel0.jl"))
    end
    # If scenario uses standard carbon matrix and population weighting
    if boxmodel == 0 && population == 1
        include(joinpath(@__DIR__, "ACE_SCC_Deterministic_population1_boxmodel0.jl"))
    end
    # If scenario uses Joos carbon matrix and no population weighting
    if boxmodel == 1 && population == 0
        include(joinpath(@__DIR__, "ACE_SCC_Deterministic_boxmodel1.jl"))
    end
    # If scenario uses Joos carbon matrix and population weighting
    if boxmodel == 1 && population == 1
        include(joinpath(@__DIR__, "ACE_SCC_Deterministic_population1_boxmodel1.jl"))
    end

end

# Create data frame in csv and in excel format
df = DataFrame(Scenario = tab["scen"],
                rho_discount_annual_exo = tab["rho_discount_annual_exo"],
                damages = tab["damages"],
                boxmodel = tab["boxmodel"],
                population = tab["population"],
                pop_recalibrate = tab["pop_recalibrate"],
                kappa = tab["kappa"],
                kappa_calib = tab["kappa_calib"],
                Carb_Mult = tab["carb_mult"], 
                w_o_TD = tab["wo_delay"], 
                SCC = tab["SCC"], 
                Gallon = tab["cent_gallon"], 
                liter = tab["cent_liter"])
display(df)
CSV.write(path*"/output.csv", df)
XLSX.writetable(path*"/output.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

end

end # module ACE_Traeger_replication
