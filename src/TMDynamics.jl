##This julia file:
####Initializes ACE CLimate Dynamics

#Package
using LinearAlgebra

cs=3
Mpre=596.4
param = Dict()
param["xi"] = [0.231049060186648, 0.231049060186648, 0.231049060186648]
sigma_temp_up = [0.538194275028702]

TempMatrix = [0.000000027559817 0.461805697411481 0;
              0.076519839162246 0.895993627259103 0.027486533578651;
              0 0.003992679633153 0.996007320366847]

#Testing for possible mistakes/incompatibility
param_xi1= log(2) / cs #Approximately 0.23, temperature trafo parameter

if param_xi1-param["xi"][1]>1e-10
    print("Inconsistency with xi1 definition, check underlying cs")
end

#RCP 6 2015 values for possible initialization (not required for paper):
Initial = Dict("Temperature" => [0.98901461, 0.73635941, 0.73635941],
               "tau" => [1.256727218547292, 1.185465055669960, 1.185465055669960],
               "M" => [818.985, 1527, 10010])

#Testing equilibrium temperature increase of given specification:
tau_eq_vec =(Matrix(1.0I, length(Initial["Temperature"]), length(Initial["Temperature"]))-TempMatrix)^(-1)*2*[sigma_temp_up[1]; zeros(length(Initial["Temperature"])-1)]
#and translating back into temperature in degree Celsius:
T_eq_vec=log.(tau_eq_vec)./param["xi"]
print("Test of long-run equilibrium Temp at 2xCO2", round.(Int, T_eq_vec))

##Carbon Dynamics ##########Redo with parameters !!!!!

#Parameters Carbon Cycle Matrix
sigma_carb = [0.088, 0.0025, 0.03832888889, 3.3750e-04]  # DICE
sigma_carb = sigma_carb * timestep / 5  # rescales to 10 year time step. Extremely similar to 5 year
CarbonMatrix = [1 - sigma_carb[1]  sigma_carb[3]  0; 
                sigma_carb[1]  1 - sigma_carb[2] - sigma_carb[3]  sigma_carb[4];
                0  sigma_carb[2]  1 - sigma_carb[4]]

