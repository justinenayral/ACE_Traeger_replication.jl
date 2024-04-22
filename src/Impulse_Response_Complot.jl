using MAT

include("Impulse_Response_ACE.jl")

logi = "001" # 001 is standard timing: pulse added to M0, no emission decay during period, T_t+1(M_t+1), gives same result as 111
# logi = 111 # 111 is advanced timing: adds pulse to current emissions and has emissions decay already, T_t+1(M_t+1) 

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

# Take DICE 2013 and 16 from 5 year time step (scaling of temperature model not working)

data1 = matread(path*"impulse_timestep_5_logi_$(logi).mat")

M13 = data1["M"]  # 5 rather than 10 year time step
M_pulse13 = data1["M_pulse"] # 5 rather than 10 year time step"]
DICE_horizon = data1["horizon_legend"]
DICE_Impulse = data1["Temp_DICE_pulse"][1,1:trunc(Int, data1["hori"])]-data1["Temp_DICE"][1,1:trunc(Int, data1["hori"])]
DICE_Impulse16 = data1["Temp_DICE_pulse16"][1,1:trunc(Int, data1["hori"])]-data1["Temp_DICE16"][1,1:trunc(Int, data1["hori"])]
ACE_horizon5 = data1["horizon_legend"]
ACE_Impulse_DICE5 = data1["Temp_pulse"][1,1:trunc(Int, data1["hori"])]-reshape(data1["Temp"][1,:,:], 6*39)[1:trunc(Int, data1["hori"])]
ACE_Impulse_Joos5 = data1["Temp_Joos"][1,1:trunc(Int, data1["hori"])]-reshape(data1["Temp"][1,:,:], 6*39)[1:trunc(Int, data1["hori"])]
pulse5y = data1["pulse"]

# Take 1 year time step version of ACE

data2 = matread(path*"impulse_timestep_1_logi_$(logi).mat")
ACE_horizon1 = data2["horizon_legend"]
ACE_Impulse_DICE1 = data2["Temp_pulse"][1,1:trunc(Int, data2["hori"])]-reshape(data2["Temp"][1,:,:], 6*39)[1:trunc(Int, data2["hori"])]
ACE_Impulse_Joos1 = data2["Temp_Joos"][1,1:trunc(Int, data2["hori"])]-reshape(data2["Temp"][1,:,:], 6*39)[1:trunc(Int, data2["hori"])]
pulse1y = data2["pulse"]
if pulse1y != pulse5y
    error("loading models where 1 year pulse differs from 5 year pulse")
end

if logic_VenBest_alt==1
    Ven = matread(path*"Venmans_CS3_lam_1p06.mat")
    JGbest_alt = Ven["Ven"]["JGbest"]
end 

data3 = matread(path*"impulse_timestep_1_logi_$(logi).mat")
ACE_horizon = data3["horizon_legend"]
ACE_Impulse_DICE = data3["Temp_pulse"][1,1:trunc(Int, data3["hori"])]-reshape(data3["Temp"][1,:,:], 6*39)[1:trunc(Int, data3["hori"])]
ACE_Impulse_Joos = data3["Temp_Joos"][1,1:trunc(Int, data3["hori"])]-reshape(data3["Temp"][1,:,:], 6*39)[1:trunc(Int, data3["hori"])]
if data3["pulse"] != pulse5y
    error("loading 10 year models where pulse differs from 5 year pulse")
end

# FROM AR6: 
TCRE_best = 1.65*data3["pulse"]/1000;
# likely interval (by guidance notes AR5) 
TCRE_high = 2.3*data3["pulse"]/1000;
TCRE_low = 1*data3["pulse"]/1000;

if logic_Bern2p5_Joos_PI==1;
    data4 = matread(path*"Impulse_Response_Bern2p5_Joos.mat")   # preindustrial, by courtesy of Fortunat Joos
    Impulse_Response_Bern2p5 = data4["Impulse_Response_Bern2p5"][1:201]
    horizon_Bern = 0.5:200.5
end

if logic_Bern2p5_Joos_PD==1;
    data5 = matread(path*"Impulse_Response_Bern2p5_PD_Joos.mat")  # present day, by courtesy of Fortunat Joos
    Impulse_Response_Bern2p5_PD = data5["Impulse_Response_Bern2p5_Joos"][1:201]
    horizon_Bern = 0.5:200.5
end

data6 = matread(path*"Venmans_CS3p1.mat")  # Loads time paths VenBest courtesy of Frank Venmans.


if logic_plot == 1
    
    linewidth = 1.75
    titlefontsize = 18

    #h1 = plot()
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
    title!("Temperature Imp Response to 100GtC", fontsize=titlefontsize)
    xlabel!("Year")
    ylabel!("Temperature Diff in C")
    xticks!(0:20:horizon)
    savefig(h3, "ImpulseResponse_Temp_Complot_$(pulse)GtC.png")

end
