# ACE_SCC_Deterministic.jl
# Calculates deterministic shadow values and SCC for ACE model


# MAIN SETTINGS
# Fraction staying forever
a0 = 0.2173
xi1 = 0.23104906018 * 3 / cs

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
    phi["M"] = zeros(3,5)
    phi["M_edec"] = zeros(3,5)
    phi["Ebox_contr"] = zeros(3,4)
    phi["Ebox_contr_perc"] = zeros(3,4)
    phi["k_pop"] = zeros(3,9)
    phi["tau_pop"] = zeros(3,9,3)
    phi["M_pop"] = zeros(3,9,5)
    phi["M_pop_help"] = zeros(3,9,4)
    phi_NoTempLag["tau_pop"] = zeros(3,9,3)
    phi_NoTempLag["M_pop"] = zeros(3,9,5)
    phi_NoTempLag["M_pop_help"] = zeros(3,9,4)
    x = zeros(3)
    x_pop = zeros(3,8)
    phi_NoTempLag["Ebox_contr"] = zeros(3,4)
    phi_NoTempLag["Ebox_contr_perc"] = zeros(3,4)
    phi_NoTempLag["M"] = zeros(3,5)
    PsiTau_scen = zeros(3,3,3)
    PsiM_scen = zeros(3,4,4)
    phi["tau"] = zeros(3,3)
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
    SCC_inC_pop = zeros(3,5) 
    SCC_inCO2_pop = zeros(3,5) 
    SCC_inC_NoTempLag_pop = zeros(3,5) 
    SCC_inCO2_NoTempLag_pop = zeros(3,5) 
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
    PsiM = inv(diagm(ones(4)) - param["beta"][i] * CarbonMatrix) 
    PsiM_scen[i, :, :] = PsiM


    temp = param["beta"][i] * sigma_temp_up[1] * phi["tau"][i, 1] / Mpre * diag(PsiM)  # for ppm
    temp_edec = param["beta"][i] * sigma_temp_up[1] * phi["tau"][i, 1] / Mpre * PsiM * CarbonMatrix  # for ppm
    temp_NoTempLag = param["beta"][i] * 1 * (-param["xi0"]) * (1 + param["beta"][i] * phi["k"][i]) * 1 / Mpre * diag(PsiM)  # for ppm

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
    # others remain same as with temp-lag
    
    # SCC calculation based on population weighting
    if population == 1
        growth_pop_factor = pop[2:end] ./ pop[1:end-1]
        growth_pop_rate = growth_pop_factor .- 1
        
        # Initialize year 2100
        phi["k_pop"][i, length(growth_pop_factor) + 1, :] .= phi["k"][i]
        phi["tau_pop"][i, length(growth_pop_factor) + 1, :] .= phi["tau"][i, :]
        phi_NoTempLag["tau_pop"][i, length(growth_pop_factor) + 1, :] .= phi["tau"][i, :] # Initialization doesn't matter for phi.tau (overwrite)
        

            phi["M_pop_help"][i, length(growth_pop_factor)+1,:] = phi["M"][i,2:length(temp)+1]
            phi["M_pop"][i,length(growth_pop_factor)+1,:] = phi["M"][i,1:length(temp)+1]  # one more element for SCC (2100, rest assigned ex-post)
            phi_NoTempLag["M_pop_help"][i,length(growth_pop_factor)+1,:] = phi_NoTempLag["M"][i,2:length(temp_NoTempLag)+1]
            phi_NoTempLag["M_pop"][i,length(growth_pop_factor)+1,:] = phi_NoTempLag["M"][i,1:length(temp_NoTempLag)+1] # one more element for SCC (2100, rest assigned ex-post)
            Vec = ones(4) # all boxes contribute to forcing
        
    
        # Loop works backwards taking 2100 as eql values calculated above w/o population growth
        for lind = length(growth_pop_factor):-1:1
            phi["k_pop"][i, lind, 1] = param["kappa"] + param["beta"][i] * growth_pop_factor[lind] * param["kappa"] * phi["k_pop"][i, lind+1, 1]
            phi["tau_pop"][i, lind, :] = param["beta"][i] * growth_pop_factor[lind] * transpose(phi["tau_pop"][i, lind+1, :]) * TempMatrix -
                                        reshape([(1 + param["beta"][i] * growth_pop_factor[lind] * phi["k_pop"][i, lind+1, 1]) * param["xi0"];
                                        zeros(2)], (3,1))'
            phi_NoTempLag["tau_pop"][i, lind, :] .= -[(1 + param["beta"][i] * growth_pop_factor[lind] * phi["k_pop"][i, lind+1, 1]) * param["xi0"];
                                                    zeros(2)]
            phi["M_pop_help"][i, lind, :] .= param["beta"][i] * growth_pop_factor[lind] * phi["tau_pop"][i, lind+1, 1] * sigma_temp_up[1] / Mpre * Vec +
                                           reshape(param["beta"][i] * growth_pop_factor[lind] * transpose(phi["M_pop_help"][i, lind+1, :]) * CarbonMatrix, (4,1))
            phi_NoTempLag["M_pop_help"][i, lind, :] .= param["beta"][i] * growth_pop_factor[lind] * phi_NoTempLag["tau_pop"][i, lind+1, 1] / Mpre * Vec +
                                                     reshape(param["beta"][i] * growth_pop_factor[lind] * transpose(phi_NoTempLag["M_pop_help"][i, lind+1, :]) * CarbonMatrix, (4,1))
            x_pop[i, lind] = 1 / (1 + param["beta"][i] * growth_pop_factor[lind] * phi["k_pop"][i, lind, 1])
            
            phi["M_pop"][i,lind,:]=[dot(CarbonWeights, phi["M_pop_help"][i,lind,:]); phi["M_pop_help"][i,lind,:]]
            phi_NoTempLag["M_pop"][i,lind,:]=[dot(CarbonWeights, phi_NoTempLag["M_pop_help"][i,lind,:]); phi_NoTempLag["M_pop_help"][i,lind,:]]
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
    

    
        SCC_inC_edec[i, :] = -phi["M_edec"][i, :] * C[i] / 1e9  # edec = emission decay already over the course of the period
        SCC_inCO2_15_edec[i, :] = SCC_inC_edec[i, :] / (44 / 12)
        SCC_gall[i, 1] = SCC_inCO2_15[i, 1] * 8.78 / 1000 * 100  # 8.78 kgCO2/gallon (and ton->kg, and USD->cent)
        SCC_lit[i, 1] = SCC_gall[i] / 3.78541 * EURperUSD  # 1 gallon = 3.78541 liter

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

