using Plots
using LinearAlgebra

timestep = 10

# timestep can be chosen arbitrarily and will generate (wrong) results. Only reasonable time steps:
# timestep = 5 required for DICE models 
# timestep = 10 is main ACE specification
# timestep = 1 generates ACE specification calibrated to annual temperature dynamics
    
logic_plot=1 # turns plotting on


logic_timing_pulse = 0;
       # =0: adds pulse to M stock at beginning of present period (default for simple)
       # =1: adds pulse to E in current period  (default for advanced timing)
logic_emission_decay = 0;  
       # =0: Emissions added straight to atmosphere (default for simple)
       # =1: Emissions already partly decay on way to atmosphere - carbon cycle is applied as CarbonMatrix*(M_t+E_t)  (default for advanced timing)
       # setting 1 requires also logic_timing_pulse = 1 embedding pulse in emissions.
logic_forcing_delay = 1;       
       # =0: T_t+1(M_t) "common timing with delay"  
       # =1: T_t+1(M_t+1) DICE and ACE preferred timing   (default for simple & advanced timing)


#Parameters
    horizon=210;   # years
    pulse = 100 ;  # Carbon pulse in GtC
       
# Select versions to calculate, save, and plot (ACE-DICE always on)
    logic_ACE_Joos=1; # Includes Joosetal 2013
    logic_DICE=1;   # Includes DICE 2013
    logic_DICE16=1; # Includes DICE 2016


    if logic_ACE_Joos == 1
        include("Joos_decay_boxes.jl")
        CarbonMatrixJoos = copy(CarbonMatrix)
    end
    
param_timestep = timestep  # needed in TMDynamics for scaling of carbon cycle dynamics
include("TMDynamics.jl") 
xi = param["xi"]
    
    if timestep == 1
        # overwrite temperature dynamics with system calibrated to 1 year time step
        param.xi = [0.231049060186648, 0.231049060186648, 0.231049060186648]  # same as above
        sigma_temp_up = [0.146170047477258]
        TempMatrix = [0.711637064433489  0.142192888089253  0;
                      0.007035452978769  0.991128650190495  0.001835896830736;
                      0                  0.000378982121102  0.999621017878898]
    end
    
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

if logic_timing_pulse == 0 && logic_emission_decay == 1
    error("invalid option, if logic_emission_decay=1, then also need logic_timing_pulse=1")
end

# DICE 2013 Dynamic Variables
ParaDice = Dict()
ParaDice["Mpre"] = 588
ParaDice["climate_sens"] = 3  # 3.2;  3;
println("Climate sensitivity is ", ParaDice["climate_sens"])
println("DICE13 itself would use 3.2")
ParaDice["eta"] = 3.8
ParaDice["Initial_M"] = [830.4 1527 10010]  # from GAMS. In EXCEL [818.985 1527 10010]. But initialization doesn't matter for linear model's impulse response

# Parameters Carbon Cycle Matrix:
sigma_carb = [0.088, 0.0025, 0.03832888889, 3.3750e-04]
sigma_carb .= sigma_carb .* timestep / 5  # rescales time step based on DICE's 5 years.
ParaDice["CarbonMatrix"] = [1 - sigma_carb[1] sigma_carb[3] 0;
                         sigma_carb[1] 1 - sigma_carb[2] - sigma_carb[3] sigma_carb[4];
                         0 sigma_carb[2] 1 - sigma_carb[4]]

# Temperature Dynamics - unfortunately timestep scaling not working
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


# DICE 2016 Climate Dynamics
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
        # Joos Impulse response: Note that it's already an impulse response model for carbon
        Imp[t + 1] = sum(CarbonMatrixJoos^(t) * CarbonWeights) * pulse
        # (result is independent of base-concentration by nature of impulse response model)
        forc_Joos[t] = sigma_forc * (M[1, t] + Imp[t + 1]) / Mpre
        tau_Joos[:, t + 1] = TempMatrix * tau_Joos[:, t] + [forc_Joos[t], 0, 0]
        Temp_Joos[:, t + 1] = diagm(1.0 ./ xi) * log.(tau_Joos[:, t + 1])
    end

    if logic_DICE == 1
        Fback = zeros(hori)
        Fpul = zeros(hori)
        # Using same M as in ACE-DICE 
        Fback[t] = ParaDice["eta"] * (log(M[1, t + 1]) - log(ParaDice["Mpre"])) / log(2)  # only needed if using them in Dietz et al for comparison purposes
        Tvec[:, t + 1] = [1 0]' * ParaDice["climate_sens"] * ParaDice["sigma_tempDICE"][1] / ParaDice["eta"] * (ParaDice["eta"] * (log(M[1, t + 1]) - log(ParaDice["Mpre"])) / log(2)) + ParaDice["TempDICEMatrix"] * Tvec[:, t]
        Temp_DICE[t + 1] = Tvec[1, t + 1]
        Fpul[t] = ParaDice["eta"] * (log(M_pulse[1, t + 1]) - log(ParaDice["Mpre"])) / log(2)  # only needed if using them in Dietz et al for comparison purposes
        Tvec_pulse[:, t + 1] = [1 0]' * ParaDice["climate_sens"] * ParaDice["sigma_tempDICE"][1] / ParaDice["eta"] * (ParaDice["eta"] * (log(M_pulse[1, t + 1]) - log(ParaDice["Mpre"])) / log(2)) + ParaDice["TempDICEMatrix"] * Tvec_pulse[:, t]
        Temp_DICE_pulse[t + 1] = Tvec_pulse[1, t + 1]
    end
    
    if logic_DICE16 == 1 || logic_DICE == 1
        # Carbon cycle:
        # Calculate emissions that keep concentration constant. Equation equivalent to CarbonMatrix*M+E=M.
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
    title!("Temperature Imp Response to 100GtC", fontsize=titlefontsize)
    xlabel!("Year")
    ylabel!("Temperature Diff in C")
    xticks!(0:20:horizon)
    savefig(h1, "ImpulseResponse_Temperature_$(pulse)GtC_timestep_$(timestep).png")

end
