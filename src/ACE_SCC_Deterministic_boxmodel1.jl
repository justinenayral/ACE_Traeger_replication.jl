# ACE_SCC_Deterministic.jl
# Last edited CT 2021 (created CT 2014 merging earlier files)
#
# Calculates deterministic shadow values and SCC for ACE model
# - calls TMDynamics.jl to initialize Climate Dynamics
# - if using boxmodel calls Joos_decay_boxes.jl 
# - if using population weighting with time preference recalibration calls
#      CalibrateBetaPopulation.jl
#

    #For DICE: get the xi0
    function xi0DICE(initial_guess::Vector, a, matchtemp, T=2.5)
        function xi0solve_DICE!(F, xi0)
            F[1] = 1 / (1 + a * matchtemp^2) - exp(-xi0[1] * exp(xi1 * matchtemp) + xi0[1])
        end
        result = nlsolve(xi0solve_DICE!, initial_guess) # Solve the equation using NLsolve
        x = result.zero[1] # Extract the solution
        return x
    end

# MAIN SETTINGS
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


    # Solve for damage coefficient
    param["xi0"] = xi0DICE([0.0], a_Nord, matchtemp)
    println("Using cs=$(cs), xi1=$(param["xi"][1]), matching xi0=$(param["xi0"])")

    if switch_damvar == 1
        param["xi0"] *= dam_fact
        println("!! Changing damages by factor $(dam_fact)")
    elseif switch_damvar == 2
        param["xi0"] = param["xi0_exo"]
        println("!! Changing damages parameter $(param["xi0"]) to exogenously given $(param["xi0"]) probably for deterministic base value in stochastic run...")
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

    println("Matched quadratic damage coefficient of $a_Nord at $matchtemp degree Celsius")
    println("DETERMINISTIC MODEL OUTPUT (for different calibrations): ")
    

    phi = Dict()
    phi_NoTempLag = Dict()
    phi["k"] = zeros(3)
    phi["M"] = zeros(3,5)
    phi["M_edec"] = zeros(3,5)
    phi["Ebox_contr"] = zeros(3,4)
    phi["Ebox_contr_perc"] = zeros(3,4)
    phi_NoTempLag["Ebox_contr"] = zeros(3,4)
    phi_NoTempLag["Ebox_contr_perc"] = zeros(3,4)
    phi_NoTempLag["M"] = zeros(3,5)
    PsiTau_scen = zeros(3,3,3)
    PsiM_scen = zeros(3,4,4)
    phi["tau"] = zeros(3,3)
    x = zeros(3)
    C = zeros(3)
    SCC_inC = zeros(3,5)
    SCC_inCO2_15 = zeros(3,5)
    SCC_inC_edec = zeros(3,5)
    SCC_inCO2_15_edec = zeros(3,5)
    SCC_inCO2_NoTempLag = zeros(3)
    SCC_gall = zeros(3)
    SCC_lit = zeros(3)
    SCC_inC_NoTempLag = zeros(3,5)
    SCC_inC_edec = zeros(3,5) 
    SCC_inC_direct = []


    for i in (1,3)
    # STATIONARY CASE (no population weighting scenario)
    # Calculate phi["M"]
    phi["k"][i] = param["kappa"] / (1 - param["kappa"] * param["beta"][i])
    # Below "Psi_sigma" matrix: (1-beta sigma)^-1
    PsiTau = inv(diagm(ones(3)) - param["beta"][i] * TempMatrix)
    PsiTau_scen[i, :, :] = PsiTau
    phi["tau"][i, :] = -param["xi0"] * (1 + param["beta"][i] * phi["k"][i]) * PsiTau[1, :]
    # Below "Psi_Phi" matrix: (1-beta Phi)^-1
    PsiM = inv(diagm(ones(4)) - param["beta"][i] * CarbonMatrix) 
    PsiM_scen[i, :, :] = PsiM


    # Here first box shadow value does not correspond to emissions shadow value
    # Code makes SCC first entry of phi.M and then adds box shadow values as additional entries to phi.M
    # Formula change: first two lines use diagonals diag(PsiM) rather than first row of PsiM to get at shadow value of individual boxes
    # Shadow value by unweighted(!) box ("virtual reservoir"):
    temp = param["beta"][i] * sigma_temp_up[1] * phi["tau"][i, 1] / Mpre * diag(PsiM)  # for ppm
    temp_edec = param["beta"][i] * sigma_temp_up[1] * phi["tau"][i, 1] / Mpre * PsiM * CarbonMatrix  # for ppm
    temp_NoTempLag = param["beta"][i] * 1 * (-param["xi0"]) * (1 + param["beta"][i] * phi["k"][i]) * 1 / Mpre * diag(PsiM)  # for ppm
    # First entry: Shadow value of a unit of emission spreading into the different boxes
    # then: temp' also adds unweighted entries back into phi.M
    phi["M"][i, :] = [dot(CarbonWeights,temp); temp]  # Note: temp contains individual box values, first entry is total
    phi["M_edec"][i, :] = [dot(CarbonWeights, ones(1, 4) * temp_edec); diag(temp_edec)]
    # Individual box's contribution to SCC:
    phi["Ebox_contr"][i, :] = CarbonWeights .* temp
    phi["Ebox_contr_perc"][i, :] = CarbonWeights .* temp ./ dot(CarbonWeights, temp) * 100
    CarbBoxMultiplier = CarbonWeights' * diag(PsiM)
    Emultiplier = CarbonWeights .* diag(PsiM)
    # Same for case w/o temperature lag:
    phi_NoTempLag["M"][i, :] = [dot(CarbonWeights, temp_NoTempLag); temp_NoTempLag]
    phi_NoTempLag["Ebox_contr"][i, :] = CarbonWeights .* temp_NoTempLag
    phi_NoTempLag["Ebox_contr_perc"][i, :] = CarbonWeights .* temp_NoTempLag ./ dot(CarbonWeights, temp_NoTempLag) * 100
    # others remain same as with temp-lag: NoTempLag.CarbBoxMultiplier and NoTempLag.Emultiplier
    

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
        SCC_inC_edec[i, :] = -phi["M_edec"][i, :] * C[i] / 1e9  # edec = emission decay already over the course of the period
        SCC_inCO2_15_edec[i, :] = SCC_inC_edec[i, :] / (44 / 12)
        SCC_gall[i, 1] = SCC_inCO2_15[i, 1] * 8.78 / 1000 * 100  # 8.78 kgCO2/gallon (and ton->kg, and USD->cent)
        SCC_lit[i, 1] = SCC_gall[i] / 3.78541 * EURperUSD  # 1 gallon = 3.78541 liter

    # Cutting out temperature delay would result in SCC
        SCC_inC_NoTempLag[i, :] = -phi_NoTempLag["M"][i, :] * C[i] / 10^9
        SCC_inCO2_NoTempLag[i, :] .= SCC_inC_NoTempLag[i] / (44 / 12)


    # Decomposition terms:
    println(" Case: ", jj, " - implied annual rho ", round(param["rho"][i] * 100, digits=2), "%")
    println("  Varphi_k (utility shadow value of capital) ", round(phi["k"][i], digits=5))
    println("  Varphi_tau (utility shadow value generalized temperature) ", round.(phi["tau"][i, :], digits=5))
    println("  Varphi_M (utility shadow value carbon) ", round.(phi["M"][i, :], digits=5))
    println(" SCC year 2020: ", round(SCC_inC[i], digits=5), " & in CO2: ", round(SCC_inC[i] / (44 / 12), digits=5), ". Decomposing: ")
    println("    beta 10Y/Mpre= ", round(param["beta"][i] * 10 * param["Y"] / (Mpre * 10^9), digits=5), " & in CO2: ", round(param["beta"][i] * param["Y"] / (Mpre * 10^9) / (44 / 12) * 10, digits=5), ".")
    println("    damage parameter xi0 ", round(param["xi0"], digits=5))
    println("    beta 10Y/Mpre*xi0= ", round(param["beta"][i] * 10 * param["Y"] / (Mpre * 10^9) * param["xi0"], digits=5), " & in CO2: ", round(param["beta"][i] * param["Y"] / (Mpre * 10^9) * param["xi0"] / (44 / 12) * 10, digits=5))
    println("   Base: beta 10C/(1-beta*kappa)/Mpre*xi0= ", round(param["beta"][i] * C[i] / (1 - param["beta"][i] * param["kappa"]) / (Mpre * 10^9) * param["xi0"], digits=5), " & in CO2: ", round(param["beta"][i] * C[i] / (1 - param["beta"][i] * param["kappa"]) / (Mpre * 10^9) * param["xi0"] / (44 / 12), digits=5))
    println("    sigma^forc ", round(sigma_temp_up[1], digits=5))
    println("    CarbBoxMultiplier contribution =", round(CarbBoxMultiplier, digits=5))
    println("      with $(round.(Emultiplier, digits=5)) coming from each of the boxes")
    println("      and individual SCC contributions of $(round.(SCC_inCO2_15[i] / CarbBoxMultiplier * Emultiplier, digits=5)) coming from each of the boxes")
    println("    PsiTau_{1,1} contribution ", round(PsiTau[1, 1], digits=5))
    println(" SCC w/o temp delay (cut sigma Matrix and sig^forc): ", round(SCC_inC_NoTempLag[i], digits=5), " & in CO2 : ", round(SCC_inC_NoTempLag[i] / (44 / 12), digits=5))
    println(" PERCENT of orignal tax if ignoring carbon dynamics ", round(100 / (PsiM[1, 1]), digits=5))
    println("")

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
