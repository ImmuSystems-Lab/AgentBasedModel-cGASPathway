using Distributed
addprocs(50)

@everywhere include("ProblemGenerator.jl")
@everywhere include("VirusCallBacks.jl")

########################################################################
#Virus simulation (Figure 7)
########################################################################


#Define the percentages of cells producing IFN
#Each cell has its own unique parameter set
@everywhere function VirusStoch(rep,percentIFN,prob)
    #Make a copy of the parameter container and the initial condition
    θ = deepcopy(prob.p)
	u₀ = deepcopy(prob.u0)

	#Assign cells to be initially infected
	probDistInfected = Poisson(moi)
	u₀[:,:,2] = @. m2c(1e3*rand(probDistInfected,N,N))

	#Reset what cells are initially infected and dead
	θ.cellsInfected = fill(Inf,N,N)
	θ.cellsInfected[findall(u₀[:,:,2] .> 0.0), 1] .= 0.0
	θ.cellsDead = fill(Inf,N,N)
	θ.infectFirstAttempt = trues(N,N)

    #The 11th parameter (kcat7) controls IFN production, if zero then no IFN
    if percentIFN == 0.0
        θ.par[11] .= 0.0 #Set all parameters to zero
    else
		θ.par[11] .= 47639.70295 .* rand(Bernoulli(percentIFN),N,N)
    end

    #Assign the new parameters to the model
	#println(θ)
	@show rep, percentIFN
    probStoch = remake(prob; p=θ,u0=u₀)
    sol = @time solve(probStoch,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)

    #Calculate the cell state dynamics over time
	#Percent, Sample, Time, Healthy, Infected, Dead
	numTime = length(saveTimePoints)
    allStates = zeros(numTime,6)

	allStates[:,1] = fill(percentIFN,numTime)
	allStates[:,2] = fill(rep,numTime)
	allStates[:,3] = saveTimePoints

    for (i,t) in enumerate(saveTimePoints)
        allStates[i,4:6] = cellStates(t,probStoch.p)
    end

    return allStates
end


#Define a new problem
@everywhere prob = ModelSetup(:Virus,:Stochastic,:Hetero)

#Time points to save during simulation
@everywhere const saveTimePoints = range(prob.tspan[1],prob.tspan[2],step=0.1)

#What percent of the cell population will be capable of producing IFNb?
percentIFN = 0.0:0.1:1.0
#Number of repititions for each IFN percentage
reps = 10
#Vector to save the simulations
#states = Vector(undef,reps)

iter = Iterators.product(1:reps,percentIFN)

#Simulate all of the different percentages in parallel and write data
states = pmap(x -> VirusStoch(x[1],x[2],prob),iter)

#Saving the Data
#Need to create a DataFrame with the following columns:
#Percent, Sample, Time, Healthy, Infected, Dead
statesLong = vcat(states...)
statesLong[:,4:6] = statesLong[:,4:6]./nCells #Get cell Percents

VirusSim = convert(DataFrame,statesLong)
rename!(VirusSim,[:Percent, :Sample, :Time, :Healthy, :Infected, :Dead])

#Save the simulation to a CSV to plot
CSV.write("./ServerSimulations/VirusFigure7DataHetero.csv",VirusSim)
