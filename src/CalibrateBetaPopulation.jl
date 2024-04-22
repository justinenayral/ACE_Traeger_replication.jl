function CalibrateBetaPopulation(param, phi, TempMatrix, sigma_temp_up, CarbonMatrix, Mpre, growth_pop_factor, beta_cali)
    i = 3  # only running recalibration for third scenario

    if param["popsim"] == 1  # get shadow values re-running routine after optimization
         global phi 
         global phi_NoTempLag 
         global x_pop_recal 
         global CarbBoxMultiplier_cali 
         global Emultiplier_cali
    end

    # Recalculate equilibrium shadow values for 2100
    phi["k_cali"] = param["kappa"] / (1 - param["kappa"] * beta_cali)
    PsiTau = inv(diagm(ones(3)) - beta_cali * TempMatrix)
    phi["tau_cali"][:] = -param["xi0"] * (1 + beta_cali * phi["k_cali"]) * PsiTau[1, :]
    PsiM = inv(diag(ones(3)) - beta_cali * CarbonMatrix)
    
    if param["boxmodel"] != 1
        phi["M_cali"][:] = beta_cali * sigma_temp_up[1] * phi["tau_cali"][1] / Mpre * PsiM[1, :]  # for ppm
        phi_NoTempLag["M_cali"][:] = beta_cali * 1 * (-param["xi0"]) * (1 + beta_cali * phi["k_cali"]) * 1 / Mpre * PsiM[1, :]  # for ppm
        Vec = [1; zeros(2)]  # only atm carbon contributes forcing
    elseif param["boxmodel"] == 1  # boxmodel
        temp = beta_cali * sigma_temp_up[1] * phi.tau_cali[1] / Mpre * diag(PsiM)  # for ppm
        temp_NoTempLag = beta_cali * 1 * (-param["xi0"]) * (1 + beta_cali * phi["k_cali"]) * 1 / Mpre * diag(PsiM)  # for ppm
        phi["M_cali"][:] = temp'  # [param.CarbonWeights * temp temp']
        phi_NoTempLag["M_cali"][:] = temp_NoTempLag'  # [param.CarbonWeights * temp_NoTempLag temp_NoTempLag']
        Vec = ones(length(param["CarbonWeights"]))  # all boxes contribute to forcing

        # Loop back from 2100 under population growth
        phi["k_pop_cali"][length(growth_pop_factor) + 1] = phi["k_cali"]
        phi["tau_pop_cali"][length(growth_pop_factor) + 1, :] = phi["tau_cali"]
        phi_NoTempLag["tau_pop_cali"][length(growth_pop_factor) + 1, :] = phi["tau_cali"]
        phi["M_pop_cali"][length(growth_pop_factor) + 1, :] = phi["M_cali"]
        phi_NoTempLag["M_pop_cali"][length(growth_pop_factor) + 1, :] = phi_NoTempLag["M_cali"]
        
        for lind = length(growth_pop_factor):-1:1
            phi["k_pop_cali"][lind] = param["kappa"] + beta_cali * growth_pop_factor[lind] * param["kappa"] * phi["k_pop_cali"][lind + 1]
            phi["tau_pop_cali"][lind, :] = beta_cali * growth_pop_factor[lind] * phi["tau_pop_cali"][lind + 1, :] * TempMatrix - [(1 + beta_cali * growth_pop_factor[lind] * phi["k_pop_cali"][lind + 1]) * param["xi0"] zeros(1, length(TempMatrix) - 1)]
            phi_NoTempLag["tau_pop_cali"][lind, :] = -[(1 + beta_cali * growth_pop_factor[lind] * phi["k_pop_cali"][lind + 1]) * param["xi0"] zeros(2)]
            phi["M_pop_cali"][lind, :] = beta_cali * growth_pop_factor[lind] * squeeze(phi["tau_pop_cali"][lind + 1, 1]) * sigma_temp_up[1] / Mpre * Vec + beta_cali * growth_pop_factor[lind] * phi["M_pop_cali"][lind + 1, :] * CarbonMatrix
            phi_NoTempLag["M_pop_cali"][lind, :] = beta_cali * growth_pop_factor[lind] * squeeze(phi_NoTempLag["tau_pop_cali"][lind + 1, 1]) / Mpre * Vec + beta_cali * growth_pop_factor[lind] * phi_NoTempLag["M_pop_cali"][lind + 1, :] * CarbonMatrix
            x_pop_recal[lind] = 1 / (1 + beta_cali * growth_pop_factor[lind] * phi["k_pop_cali"][lind])
        end

        # Rename M_pop_cali to save SCC in addition to shadow values
        if param["boxmodel"] == 1
            # Append SCC to beginning of phi_M which currently only contains box-shadow values
            phiMcal_tempo = phi["M_pop_cali"][lind, :]  # save temporary
            phi["M_pop_cali"] = [phiMcal_tempo * param["CarbonWeights"]'; phiMcal_tempo]
            phiMcal_tempo_NoTempLag = phi_NoTempLag["M_pop_cali"][lind, :]  # save temporary
            phi_NoTempLag["M_pop_cali"] = [phiMcal_tempo_NoTempLag * param["CarbonWeights"]'; phiMcal_tempo_NoTempLag]
        end

        consrate = x_pop_recal[1]
    end
end
