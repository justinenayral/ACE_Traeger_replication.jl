#This file
#####is the master file for the replication of Traeger (2023)

#Path
path = "C:/Users/justi/Dropbox/Traeger/Julia code/"
#path = "C:/Users/Guada/Dropbox/Traeger/Julia code/"
graphpath=path*"figure/"

#Packages
using MAT #To open .mat datasets
using Printf #For labelling graphs
using LinearAlgebra
using DataFrames
using CSV
using XLSX
using NLsolve
using Plots

include("DamageSterner.jl")
include("DamageFctPlot_ACE.jl") # Plot Figure 2
include("TempSimulation_ACE.jl")
include("TempFitSim_ACE.jl") # Plot Figure 3


# Plot Figure 4
list_impulse_timestep = [1, 1, 5, 5, 5, 10, 10, 10]
list_impulse_logic_timing_pulse = [0, 1, 0, 1, 1, 0, 1, 1]
list_impulse_logic_emission_decay = [0, 1, 0, 0, 1, 0, 0, 1]
list_impulse_logic_forcing_delay = [1, 1, 1, 1, 1, 1, 1, 1]

for kk in 1:length(list_impulse_timestep)
    include("Impulse_Response_ACE.jl")
end

include("Impulse_Response_Complot.jl")

replicate_table = 1
global tab = Dict()

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

list = Dict()

    list["rho_discount_annual_exo"] = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01, 0.01,0.01,0.01,0.01,0.01,0.01, 
                                      0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.001,0.001,
                                      0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001]

    list["damages"] = [1, 1, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 2]
    
    list["boxmodel"] = [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]
    list["population"] = [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1]
    list["pop_recalibrate"] = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    list["kappa"] = [0.3,  0.3,  0.3,  0.4,  0.4,  0.3,  0.3,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.4,  0.3,  0.3,  0.3,  0.4,  0.3,   0.3,  0.4,  0.3,  0.4,  0.3,  0.4,  0.4,  0.3,  0.3,  0.3,  0.4,  0.3,  0.3,  0.4,  0.3,  0.4,  0.3,  0.4,  0.4]
    list["kappa_calib"] = [0.3,  0.3,  0.3,  0.4,  0.3,  0.3,  0.3,  0.4,  0.4,  0.3,  0.4,  0.4,  0.3,  0.4,  0.4,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3]

for k in 1:length(list["damages"])


    global boxmodel = list["boxmodel"][k]
    global population = list["population"][k]
    global jj = k

    if boxmodel == 0 && population == 0
    include("ACE_SCC_Deterministic_boxmodel0.jl")
    end

    if boxmodel == 0 && population == 1
    include("ACE_SCC_Deterministic_population1_boxmodel0")
    end

    if boxmodel == 1 && population == 0
    include("ACE_SCC_Deterministic_boxmodel1.jl")
    end

    if boxmodel == 1 && population == 1
    include("ACE_SCC_Deterministic_population1_boxmodel1")
    end

end

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
CSV.write("output.csv", df)
XLSX.writetable("output.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))