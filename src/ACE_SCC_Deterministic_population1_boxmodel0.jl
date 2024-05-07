# ACE_SCC_Deterministic.jl
# Calculates deterministic shadow values and SCC for ACE model


# MAIN SETTINGS

cs=3
xi1 = 0.23104906018 * 3 / cs
Mpre=596.4
param = Dict()
param["xi"] = [0.231049060186648, 0.231049060186648, 0.231049060186648]
sigma_temp_up = [0.538194275028702]

TempMatrix = [0.000000027559817 0.461805697411481 0;
              0.076519839162246 0.895993627259103 0.027486533578651;
              0 0.003992679633153 0.996007320366847]

param_xi1= log(2) / cs #Approximately 0.23, temperature trafo parameter

if param_xi1-param["xi"][1]>1e-10
    print("Inconsistency with xi1 definition, check underlying cs")
end

Initial = Dict("Temperature" => [0.98901461, 0.73635941, 0.73635941],
               "tau" => [1.256727218547292, 1.185465055669960, 1.185465055669960],
               "M" => [818.985, 1527, 10010])

tau_eq_vec =(Matrix(1.0I, length(Initial["Temperature"]), length(Initial["Temperature"]))-TempMatrix)^(-1)*2*[sigma_temp_up[1]; zeros(length(Initial["Temperature"])-1)]
T_eq_vec=log.(tau_eq_vec)./param["xi"]

##Carbon Dynamics

#Parameters Carbon Cycle Matrix
sigma_carb = [0.088, 0.0025, 0.03832888889, 3.3750e-04]  # DICE
sigma_carb = sigma_carb * timestep / 5  # rescales to 10 year time step. Extremely similar to 5 year
CarbonMatrix = [1 - sigma_carb[1]  sigma_carb[3]  0; 
                sigma_carb[1]  1 - sigma_carb[2] - sigma_carb[3]  sigma_carb[4];
                0  sigma_carb[2]  1 - sigma_carb[4]]


param["rho_discount_annual_exo"] = list["rho_discount_annual_exo"][jj]  
damages = list["damages"][jj] 
use_beta_invrate_parity = false
param["kappa"] = list["kappa"][jj]         
param["kappa_calib"] = list["kappa_calib"][jj]    
EURperUSD = 1 / 1.1348  
switch_damvar = 0  
   
    if damages == 1
        matchtemp = 2.5  # DICE std: 2.5 - defines temperature at which matching DICE damages
        a_Nord = 0.0028  # DICE 2007  % BASE SCENARIO together with matchtemp = 2.5
        #a_Nord = 0.00267  # DICE 2013R
        #a_Nord = 0.00236  # DICE 2016
    elseif damages == 2
        # Calibration damages to HSP scenario (uncomment next two lines):
        matchtemp = 3   # for Howard & Sterner matched at 3C  % HSP Variation
        a_Nord = 0.01145  # Howard & Sterner (2017)             % HSP Variation
    else
        error("You have to select damages = 1 (Nordhaus) or damages = 2 (HSP)")
    end
       
    # Set timestep (has to corresponds to TMDynamics, here 10 years)
    param["timestep"] = 10 # in years

    # IMF forecast data 
    param["Y20"]=130.186703e12; # IMF Oct 2020 forecast for 2020 in PPP (trillion=10^12)
    # Investment rate: 
    param["InvRate20"]=0.26108	; # IMF 2020 investment rate estimate, forecasted in Oct 2020
    # Implied factors:
    param["beta_exo"]=exp(-param["rho_discount_annual_exo"]*param["timestep"]);

    function xi0DICE(initial_guess::Vector, a, matchtemp, T=2.5)
        function xi0solve_DICE!(F, xi0)
            F[1] = 1 / (1 + a * matchtemp^2) - exp(-xi0[1] * exp(xi1 * matchtemp) + xi0[1])
        end
        result = nlsolve(xi0solve_DICE!, initial_guess) # Solve the equation using NLsolve
        x = result.zero[1] # Extract the solution
        return x
    end

    # Solve for damage coefficient
    param["xi0"] = xi0DICE([0.0], a_Nord, matchtemp)

    if switch_damvar == 1
        param["xi0"] *= dam_fact
        println("!! Changing damages by factor $(dam_fact)")
    elseif switch_damvar == 2
        param["xi0"] = param["xi0_exo"]
        println("!! Changing damages parameter $(param["xi0"]) to exogenously given $(param["xi0"]) probably for deterministic base value in stochastic run...")
    end

# Population input
    if population == 1 # load UN population data, 2020-2100 in 5 year steps
        pop_5year = [7794799, 8184437,	8548487, 8887524, 9198847, 9481803, 9735034, 9958099,	10151470, 10317879, 10459240, 10577288, 10673904, 10750662, 10809892, 10851860, 10875394]
        pop_decadal = pop_5year[1:2:17] # first 2020, last 2100
        pop = pop_decadal
        if param["timestep"] != 10
            error("population uses decadal time step, different from model time step")
        end 
    end

# Prepare for calculating deterministic shadow values
# Calculate rho and beta for endogenous time preference calibration
    beta_implied20 = param["InvRate20"] / param["kappa_calib"]
    rho_implied20 = -1 / param["timestep"] * log(beta_implied20)

# Arranging discount parameters into a list for looping
# Eliminated outdated world output calibration as list entry 2 (thus 2 is skipped in for loop below)
    param["beta"] = zeros(3)
    param["beta"][1] = param["beta_exo"]
    param["beta"][3] = beta_implied20
    param["rho"] = zeros(3)
    for i in 1:3
        param["rho"][i] = -1 / param["timestep"] * log(param["beta"][i]) # putting rho's into list for loop
    end


    # Initialization
        phi = Dict()
        phi_NoTempLag = Dict()
        phi["k"] = zeros(3)
        phi["M"] = zeros(3,3)
        phi["M_edec"] = zeros(3,3)
        phi_NoTempLag["M"] = zeros(3,3)
        PsiTau_scen = zeros(3,3,3)
        PsiM_scen = zeros(3,3,3)
        phi["tau"] = zeros(3,3)
        phi["M_immediate_forc"] = zeros(3,3)
        phi["M_immediate_forc_tot"] = zeros(3,3)
        phi["k_pop"] = zeros(3,9)
        phi["tau_pop"] = zeros(3,9,3)
        phi["M_pop"] = zeros(3,9,3)
        phi["M_pop_help"] = zeros(3,9,3)
        phi_NoTempLag["tau_pop"] = zeros(3,9,3)
        phi_NoTempLag["M_pop"] = zeros(3,9,3)
        phi_NoTempLag["M_pop_help"] = zeros(3,9,3)
        x = zeros(3)
        x_pop = zeros(3,8)
        C = zeros(3)
        SCC_inC = zeros(3,3)
        SCC_inC_immediate_forc = zeros(3,3)
        SCC_inCO2_immediate_forc = zeros(3,3)
        SCC_inCO2_immediate_forc_tot = zeros(3,3)
        SCC_inC_immediate_forc_tot = zeros(3,3)
        SCC_gall = zeros(3)
        SCC_inCO2_15 = zeros(3,3)
        SCC_lit = zeros(3)
        SCC_inC_NoTempLag = zeros(3,3)
        SCC_inCO2_NoTempLag = zeros(3,3)
        SCC_inC_edec = zeros(3,3) 
        SCC_inCO2_15_edec = zeros(3,3) 
        SCC_inC_direct = []
        SCC_inC_pop = zeros(3,3) 
        SCC_inCO2_pop = zeros(3,3) 
        SCC_inC_NoTempLag_pop = zeros(3,3) 
        SCC_inCO2_NoTempLag_pop = zeros(3,3) 
        SCC_gall_pop = zeros(3) 
        SCC_lit_pop = zeros(3) 

    for i in (1,3)
    # STATIONARY CASE (no population weighting scenario)
    # Calculate phi["M"]
    phi["k"][i] = param["kappa"] / (1 - param["kappa"] * param["beta"][i])
    # Below "Psi_sigma" matrix: (1-beta sigma)^-1
    PsiTau = inv(diagm(ones(3)) - param["beta"][i] * TempMatrix)
    PsiTau_scen[i, :, :] = PsiTau
    phi["tau"][i, :] = -param["xi0"] * (1 + param["beta"][i] * phi["k"][i]) * PsiTau[1, :]
    # Below "Psi_Phi" matrix: (1-beta Phi)^-1
    PsiM = inv(diagm(ones(3)) - param["beta"][i] * CarbonMatrix) 
    PsiM_scen[i, :, :] = PsiM


    phi["M"][i, :] = param["beta"][i] * sigma_temp_up[1] * phi["tau"][i, 1] / Mpre * PsiM[1, :]
    phi["M_edec"][i, :] = param["beta"][i] * sigma_temp_up[1] * phi["tau"][i, 1] / Mpre * PsiM * CarbonMatrix[:, 1]
    # w/o temp delay replace sigma_forc*phi.tau by - xi0 (1_beta*phi.k):
    phi_NoTempLag["M"][i, :] = param["beta"][i] * 1 * (-param["xi0"]) * (1 + param["beta"][i] * phi["k"][i]) * 1 / Mpre * PsiM[1, :]
    # Timing where T_t+1 is forced by M_t+1: 
    phi["M_immediate_forc"][i, :] = param["beta"][i] * sigma_temp_up[1] * phi["tau"][i, 1] / Mpre * CarbonMatrix[1, :]' * PsiM[:,:]
    phi["M_immediate_forc_tot"][i, :] = phi["M_immediate_forc"][i, :] .+ sigma_temp_up[1] / Mpre * phi["tau"][i, 1]
    
    # SCC calculation based on population weighting
    if population == 1
        growth_pop_factor = pop[2:end] ./ pop[1:end-1]
        growth_pop_rate = growth_pop_factor .- 1
        
        # Initialize year 2100
        phi["k_pop"][i, length(growth_pop_factor) + 1, :] .= phi["k"][i]
        phi["tau_pop"][i, length(growth_pop_factor) + 1, :] .= phi["tau"][i, :]
        phi_NoTempLag["tau_pop"][i, length(growth_pop_factor) + 1, :] .= phi["tau"][i, :]
        

            # Initialize shadow value with 2100 stationary values
            phi["M_pop_help"][i, length(growth_pop_factor) + 1, :] .= phi["M"][i, :]
            phi["M_pop"][i, length(growth_pop_factor) + 1, :] .= phi["M_pop_help"][i, length(growth_pop_factor) + 1, :]
            phi_NoTempLag["M_pop_help"][i, length(growth_pop_factor) + 1, :] .= phi_NoTempLag["M"][i, :]
            phi_NoTempLag["M_pop"][i, length(growth_pop_factor) + 1, :] .= phi_NoTempLag["M_pop_help"][i, length(growth_pop_factor) + 1, :]
            Vec = [1; zeros(2)]
        
    
        
        for lind = length(growth_pop_factor):-1:1
            phi["k_pop"][i, lind, 1] = param["kappa"] + param["beta"][i] * growth_pop_factor[lind] * param["kappa"] * phi["k_pop"][i, lind+1, 1]
            phi["tau_pop"][i, lind, :] = param["beta"][i] * growth_pop_factor[lind] * transpose(phi["tau_pop"][i, lind+1, :]) * TempMatrix -
                                        reshape([(1 + param["beta"][i] * growth_pop_factor[lind] * phi["k_pop"][i, lind+1, 1]) * param["xi0"];
                                        zeros(2)], (3,1))'
            phi_NoTempLag["tau_pop"][i, lind, :] .= -[(1 + param["beta"][i] * growth_pop_factor[lind] * phi["k_pop"][i, lind+1, 1]) * param["xi0"];
                                                    zeros(2)]
            phi["M_pop_help"][i, lind, :] .= param["beta"][i] * growth_pop_factor[lind] * phi["tau_pop"][i, lind+1, 1] * sigma_temp_up[1] / Mpre * Vec +
                                           reshape(param["beta"][i] * growth_pop_factor[lind] * transpose(phi["M_pop_help"][i, lind+1, :]) * CarbonMatrix, (3,1))
            phi_NoTempLag["M_pop_help"][i, lind, :] .= param["beta"][i] * growth_pop_factor[lind] * phi_NoTempLag["tau_pop"][i, lind+1, 1] / Mpre * Vec +
                                                     reshape(param["beta"][i] * growth_pop_factor[lind] * transpose(phi_NoTempLag["M_pop_help"][i, lind+1, :]) * CarbonMatrix, (3,1))
            x_pop[i, lind] = 1 / (1 + param["beta"][i] * growth_pop_factor[lind] * phi["k_pop"][i, lind, 1])
            
            phi["M_pop"][i, lind, :] .= phi["M_pop_help"][i, lind, :]
            phi_NoTempLag["M_pop"][i, lind, :] .= phi_NoTempLag["M_pop_help"][i, lind, :]
        end
    end
    


        if use_beta_invrate_parity == 1
            x[i] = 1 - param["beta"][i] * param["kappa"]  # Consumption rate is based on discount rate calibration
        else
            x[i] = 1 - param["InvRate20"]
        end
        

        C[i] = param["Y20"] * x[i] * param["timestep"]  # Consumption level is rate times actual consumption
        param["Y"]=param["Y20"]
        SCC_inC[i, :] = -phi["M"][i, :] * C[i] / 10^9  # /10^9 transforms per Gt to per ton
        SCC_inCO2_15[i, :] = SCC_inC[i, :] / (44 / 12)
    

    # Added in August 21 for timing with M_t+1 in T_t+1
        SCC_inC_immediate_forc[i, :] = -phi["M_immediate_forc"][i, :] * C[i] / 10^9  # /10^9 transforms per Gt to per ton
        SCC_inCO2_immediate_forc[i, :] = SCC_inC_immediate_forc[i, :] / (44 / 12)
    
        SCC_inC_immediate_forc_tot[i, :] = -phi["M_immediate_forc_tot"][i, :] * C[i] / 10^9  # /10^9 transforms per Gt to per ton
        SCC_inCO2_immediate_forc_tot[i, :] = SCC_inC_immediate_forc_tot[i, :] / (44 / 12)
    
        SCC_inC_edec[i, :] = -phi["M_edec"][i, :] * C[i] / 1e9  # edec = emission decay already over the course of the period
        SCC_inCO2_15_edec[i, :] = SCC_inC_edec[i, :] / (44 / 12)
        SCC_gall[i, 1] = SCC_inCO2_15[i, 1] * 8.78 / 1000 * 100  # 8.78 kgCO2/gallon (and ton->kg, and USD->cent)
        SCC_lit[i, 1] = SCC_gall[i] / 3.78541 * EURperUSD  # 1 gallon = 3.78541 liter

    # Cutting out temperature delay would result in SCC
        SCC_inC_NoTempLag[i, :] = -phi_NoTempLag["M"][i, :] * C[i] / 10^9
        SCC_inCO2_NoTempLag[i, :] .= SCC_inC_NoTempLag[i] / (44 / 12)


        if population == 1
            # SCC calculation based on population weighting
            SCC_inC_pop[i,:] = -phi["M_pop"][i,1,:] * C[i] / 10^9  # where second position is year 2020
            SCC_inCO2_pop[i,:] = SCC_inC_pop[i,:] / (44/12)
        
            # Cutting out temperature delay would result in SCC
            SCC_inC_NoTempLag_pop[i,:] = -phi_NoTempLag["M_pop"][i,1,:] * C[i] / 10^9  # where second position is year 2020
            SCC_inCO2_NoTempLag_pop[i,:] = SCC_inC_NoTempLag_pop[i,:] / (44/12)
        
            # Convert SCC to gallons and liters
            SCC_gall_pop[i,1] = SCC_inCO2_pop[i,1] * 8.78 / 1000 * 100  # 8.78 kgCO2/gallon (and ton->kg, and USD->cent)
            SCC_lit_pop[i,1] = SCC_gall_pop[i,1] / 3.78541 * EURperUSD  # 1 gallon = 3.78541 liter

        end
        

# Save information for table (if in table replication scenario)
if replicate_table == 1
    if (param["rho_discount_annual_exo"] < 0.006 && i == 1) || (param["rho_discount_annual_exo"] > 0.005 && i == 3)
        tab["scen"][jj] = jj
        
        tab["rho_discount_annual_exo"][jj] = list["rho_discount_annual_exo"][jj]
        tab["damages"][jj] = list["damages"][jj]
        tab["boxmodel"][jj] = list["boxmodel"][jj]
        tab["population"][jj] = list["population"][jj]
        tab["pop_recalibrate"][jj] = list["pop_recalibrate"][jj]
        tab["kappa"][jj] = list["kappa"][jj]
        tab["kappa_calib"][jj] = list["kappa_calib"][jj]

        if population == 0
            if boxmodel != 1
                tab["carb_mult"][jj] = PsiM[1, 1]
            elseif boxmodel == 1
                tab["carb_mult"][jj] = CarbBoxMultiplier
            end
        else
            tab["carb_mult"][jj] = 0
        end

        if population == 0
            tab["wo_delay"][jj] = SCC_inC_NoTempLag[i] / (44 / 12)
            tab["SCC"][jj] = SCC_inCO2_15[i, 1]
            tab["cent_gallon"][jj] = SCC_gall[i]
            tab["cent_liter"][jj] = SCC_lit[i]
        else
            tab["wo_delay"][jj] = SCC_inC_NoTempLag_pop[i] / (44 / 12)
            tab["SCC"][jj] = SCC_inCO2_pop[i, 1]
            tab["cent_gallon"][jj] = SCC_gall_pop[i]
            tab["cent_liter"][jj] = SCC_lit_pop[i]
        end
    end
end

end