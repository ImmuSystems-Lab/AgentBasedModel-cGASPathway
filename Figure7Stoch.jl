using Distributed
addprocs(4)

@everywhere include("ProblemGenerator.jl")
@everywhere include("VirusCallBacks.jl")

########################################################################
#Virus simulation (Figure 7)
########################################################################

#Define the percentages of cells producing IFN
#Each cell has its own unique parameter set

@everywhere function VirusStoch(percentIFN,prob)
    #Make a copy of the parameter container
    θ = deepcopy(prob.p)
    #The 11th parameter (kcat7) controls IFN production, if zero then no IFN
    if percentIFN == 0.0
        θ.par[11] .= 0.0
    else
        θ.par[11][rand(N,N) .> percentIFN] .= 0.0
    end

    #Assign the new parameters to the model
    probStoch = remake(prob; p=θ)
    sol = @time solve(probStoch,CVODE_BDF(linear_solver=:GMRES),callback=cb)

    #Calculate the cell state dynamics over time
    saveTimePoints = range(prob.tspan[1],prob.tspan[2],step=0.1)
    allStates = zeros(Int,length(saveTimePoints),3) #Healthy, Infected or Dead

    for (i,t) in enumerate(saveTimePoints)
        allStates[i,:] = cellStates(t,θ)
    end

    return allStates
end


#Define a new problem
@everywhere prob = ModelSetup(:Virus,:Stochastic)

percentIFN = 0.0:0.1:1.0
#Number of repititions for each IFN percentage
reps = 10
states = Vector(undef,reps)

for i = 1:reps
    states[i] = pmap(x -> VirusStoch(x,prob),percentIFN)
end


plot(range(prob.tspan[1],prob.tspan[2],step=0.1),states ./ nCells)
