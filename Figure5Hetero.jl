using RCall
using Distributed
addprocs(2)
@everywhere include("ProblemGenerator.jl")


#Need a function that creates distributions of heterogenous states
@everywhere AddNoise2States(σ) = [TruncatedNormal(μ,σ*μ,0,Inf) for μ in nonZeroSpeciesValues]

#Create a function that calculate the peak IFN distribution given:
 #The number of cells initially infected
 #The amount of noise desired for the initial protien concentrations
@everywhere function HeteroSimulation(primeCellPercent,varIC,prob)
    #Make a copy of the initial conditions
    u0 = copy(prob.u0)
    #Get number of primary cells wanted from percent
    primeCellNum = round(Int,nCells * primeCellPercent)
    secondaryCellNum = nCells-primeCellNum
    #Randomly place them in the ABM
    cellLocations = shuffle([ones(primeCellNum);zeros(secondaryCellNum)])
    #Update the number of initially infected cells (DNA conc.)
    u0[:,:,2] = m2c.(1e3.*reshape(cellLocations,N,N))

    #Add noise to the nonzero initial conditions
    noiseDistributions = AddNoise2States(varIC)
    for (i,index) in enumerate(nonZeroSpeciesIdx)
        u0[:,:,index] = rand(noiseDistributions[i],N,N)
    end

    #Make a copy of the parameter container
    θ = deepcopy(prob.p)
    #Reset the number of cells infected initially
    θ.cellsInfected = fill(Inf,N,N)
    θ.cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0
    #Reset the mass balances
    θ.mass = [u0[:,:,i] for i in nonZeroSpeciesIdx]

    #Define the new problem and solve
    probHetero = remake(prob; u0=u0, p=θ)
    solHetero = solve(probHetero,CVODE_BDF(linear_solver=:GMRES))

    #Seperate by primary and secondary cells
    isPrimary = u0[:,:,2] .> 0.0
    maxIFN = [maximum(solHetero[coord,7,:]) for coord in cellIndicies]

    return [maxIFN[isPrimary],maxIFN[.!isPrimary]]
end

@everywhere prob = ModelSetup(:ISD,:notStochastic,:Hetero)

varICRange = 0.0:0.1:1.0
varIClength = length(varICRange)
primeCellPercent = [0.01, 0.10, 0.63]
maxIFNDis = Vector(undef,length(primeCellPercent))

for i = 1:length(primeCellPercent)
    maxIFNDis[i] = pmap(x -> HeteroSimulation(primeCellPercent[i],x,prob),varICRange)
end



primary = [maxIFNDis[3][i][1] for i=1:length(varICRange)]
secondary = [maxIFNDis[3][i][2] for i=1:length(varICRange)]

primary = DataFrame(percent=repeat(1:varIClength,inner= round(Int,nCells * primeCellPercent[3])), maxIFN = vcat(primary...))
secondary = DataFrame(percent=repeat(1:varIClength,inner= round(Int,nCells * (1.0-primeCellPercent[3]))), maxIFN = vcat(secondary...))


@df primary violin(:percent, :maxIFN,side=:left,frame=:box,label=:primary,legend=:topleft,show_median=true,trim=false)
@df secondary violin!(:percent, :maxIFN,side=:right,label=:secondary,show_median=true,trim=false)
xticks!(1:varIClength,string.(varICRange))
xlabel!("Parameter Variability")
ylabel!("Maximum IFN Production")
