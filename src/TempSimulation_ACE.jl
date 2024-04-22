#Need to create a package for this function !!!!
using Printf

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
function TempSimulation_ACE(σ, Temp, forcing, xi, weight, logic)
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