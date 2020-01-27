using Distributed
addprocs(11)

@everywhere include("ProblemGenerator.jl")
@everywhere include("VirusCallBacks.jl")

########################################################################
#Virus simulation (Figure 7)
########################################################################


#Define the percentages of cells producing IFN
#Each cell has its own unique parameter set
@everywhere function VirusStoch(percentIFN,prob)
    #Make a copy of the parameter container and the initial condition
    θ = deepcopy(prob.p)
	u₀ = deepcopy(prob.u0)

	#Assign cells to be initially infected
	probDistInfected = Poisson(moi)
	u₀[:,:,2] = @. m2c(1e3*rand(probDistInfected,N,N))

	#Reset what cells are initially infected
	θ.cellsInfected = fill(Inf,N,N)
	θ.cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0

    #The 11th parameter (kcat7) controls IFN production, if zero then no IFN
    if percentIFN == 0.0
        θ.par[11] .= 0.0 #Set all parameters to zero
    else
		θ.par[11] .= 47639.70295 #Reset the parameters all to nonzero
        θ.par[11][rand(N,N) .> percentIFN] .= 0.0 #Set some to zero
    end

    #Assign the new parameters to the model
	println(θ)
    probStoch = remake(prob; p=θ,u0=u₀)
    sol = @time solve(probStoch,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)

    #Calculate the cell state dynamics over time
    allStates = zeros(Int,length(saveTimePoints),3) #Healthy, Infected or Dead

    for (i,t) in enumerate(saveTimePoints)
        allStates[i,:] = cellStates(t,θ)
    end

    return allStates
end


#Define a new problem
@everywhere prob = ModelSetup(:Virus,:Stochastic,:Homo)

#Time points to save during simulation
@everywhere const saveTimePoints = range(prob.tspan[1],prob.tspan[2],step=0.1)

#What percent of the cell population will be capable of producing IFNb?
percentIFN = 0.0:0.1:1.0
#Number of repititions for each IFN percentage
reps = 10
#Vector to save the simulations
states = Vector(undef,reps)

#Loop through all the replicates
for i = 1:reps
	@show i
	#Simulate all of the different percentages in parallel and write data
    states[i] = pmap(x -> VirusStoch(x,prob),percentIFN)
end


#Saving the Data
#Need to create a DataFrame with the following columns:
#Percent, Sample, Time, Healthy, Infected, Dead
numPercent = length(percentIFN)
numTime = length(saveTimePoints)
numSample = reps
statesLong = vcat(vcat(states...)...) #There has to be a better way to do this
statesLong /= nCells

VirusSim = DataFrame()
	VirusSim.Percent = repeat(percentIFN,inner=numTime,outer=numSample)
	VirusSim.Sample = repeat(1:reps,inner=numTime*numPercent)
	VirusSim.Time = repeat(saveTimePoints,outer=numPercent*numSample)
	VirusSim.Healthy = statesLong[:,1]
	VirusSim.Infected = statesLong[:,2]
	VirusSim.Dead = statesLong[:,3]


#Save the simulation to a CSV to plot
CSV.write("./ServerSimulations/VirusFigure7DataHomo.csv",VirusSim)
