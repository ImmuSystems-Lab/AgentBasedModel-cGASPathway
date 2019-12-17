#Run the Chronic simulations in parallel
#change the number of cores used as needed
using RCall
using Distributed
addprocs(56)
#The @everywhere macro runs the command/defines the thing on all the processors
@everywhere include("ProblemGenerator.jl")

#Given a new parameter set, determine average IFN at end of simulation
@everywhere function ChronicIFN(kcat8,τ7,prob)
    #Make a copy of the parameter container
    θ = deepcopy(prob.p)

    #Need to check if each cell has a unique parameter set
    if isa(prob.p.par[1],Array)
      #Set the parameters to desired input
      θ.par[13] .= kcat8
      θ.par[27] .= τ7
    else
      #Set the parameters to desired input
      θ.par[13] = kcat8
      θ.par[27] = τ7
    end


    @show (kcat8,τ7)
    #Generate the new problem (run for longer to ensure steady state)
    probNew = remake(prob; p=θ,tspan=(0.0,168.0))
    sol = solve(probNew,CVODE_BDF(linear_solver=:GMRES),saveat=1.0)

    return mean(sol[:,:,7,end])
end


#Generate a problem
@everywhere prob = ModelSetup(:ISD,:Stochastic)

#We need to know if every cell has been given a unique parameter set
if isa(prob.p.par[1],Array) #Is the first parameter a single number or array?
  kcat8True = prob.p.par[13][1]
  τ7True = prob.p.par[27][1]
else
  kcat8True,τ7True = prob.p.par[[13,27]]
end

#Vary STAT production (kcat8) and IFN Degradation (τ7)
parRange = 20
kcat8Vals = range(0.5*kcat8True,2.0*kcat8True,length=parRange)
τ7Vals = range(0.5*τ7True,2.0*τ7True,length=parRange)

iter = Iterators.product(kcat8Vals,τ7Vals)
meanIFNβ = pmap(x -> ChronicIFN(x[1],x[2],prob),iter)

#plt = heatmap(kcat8Vals,τ7Vals,meanIFNβ,color=:viridis)

kvals = zeros(parRange^2,2)
for ((kcat8,τ7),i) in zip(iter,1:parRange^2)
    kvals[i,:] = [kcat8,τ7]
end


chronicIFN = DataFrame(kcat8 = kvals[:,1],tau7 = kvals[:,2], IFN = reverse(meanIFNβ[:]))

#Save the simulation to plot
CSV.write("StochFigure6Data.csv",chronicIFN)
