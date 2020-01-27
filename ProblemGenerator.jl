#Import Packages needed for simulation
using DifferentialEquations, Sundials #For implementing Differential equations
using LinearAlgebra, SparseArrays, Distributions, Statistics, Random #Linear Algebra and Statistics
using CSV, DataFrames #Data handling
using StatsPlots #For graphing
using JLD2, FileIO #Saving simulations

###############################################################
# 1. Define all the constants that will not change
###############################################################

#Constants for all cells
const N=200 #number of grid points along one dimensions
const nCells = N^2 #number of cells in the simulation
const cellVol = 3e-12 #Cell Volume (liters)
const Na = 6.02e23 #Avagadro's number
const species = 14 #Number of states within each cell (including virus)
const moi = 1.0e-3 #Multicity of infection
const cellIndicies = CartesianIndices(zeros(N,N)) #set of indices for looping through every cell

#Function that converts molecules to nM
m2c(molecule) = @. 1e9*molecule/(cellVol*Na)

#Constants for all simulations
const tspan = (0.0,48.0) #Time span for simulation
const tstop = sort(rand(Uniform(tspan[1],tspan[2]),1000)) #Times where the simulation stops and check for virus movement
const ΔIFNβ = zeros(N,N) #Define memory space to hold the Laplacian
const statesNames = ["cGAS","DNA","Sting","cGAMP","IRF3","IFNbm","IFNb","STAT",
                     "SOCSm","IRF7m","TREX1m","IRF7","TREX1","Virus"] #for plotting
const θNames = [:k1f, :k1r, :k3f, :k3r, :k4f, :kcat5, :Km5, :k5r, :kcat6, :Km6, :kcat7,
 :Km7, :kcat8, :Km8, :k8f, :k9f, :k10f1, :k10f2, :k11f, :k12f, :k13f, :k6f, :kcat2,
 :Km2, :τ4, :τ6, :τ7, :τ8, :τ9, :τ10, :τ11, :τ12, :τ13, :k14f,:τ14] #Parameter names
 #These species are not starting at zero
 const nonZeroSpeciesIdx = [1,3,5] #cGAS, Sting, IRF3
 const nonZeroSpeciesValues = m2c([1e3, 1e3, 1e4]) #convert to concentration

#THis function modifies initial values by adding guassian noise
AddNoise2States(σ) = [TruncatedNormal(μ,σ*μ,0,Inf) for μ in nonZeroSpeciesValues]


 #Often it is useful to pass parameters between functions during the ODE solve,
 #This struct hold all the parameters that we need to keep track of
 mutable struct ParContainer{T}
   par::T #Rate Constants
   mass::Vector{Array{Float64,2}} #Mass balances
   deathParameter::Array{Float64,2} # 0 or 1 indicating cell is dead
   DNAReplicate::Int64 #0 if ISD, 1 if Virus
   cellsInfected::Array{Float64,2} #Time cell was infected
   cellsDead::Array{Float64,2} #Time cell was killed
   infectFirstAttempt::BitArray{2} #Has cell tried to infect neighbors?
 end


 ###############################################################
 # 2. Helper function that discretizes the Laplacian
 ###############################################################
 #∇u is updated with Laplacian estimate from u
 function ∇²(Δu,u)
   #Get dimensions of the input and define some constants
     n1, n2 = size(u)
     Δx = 32.0 #Grid spacing (diameter of cell in μm)
     D=97.5*3600.0 #Diffusion coefficient (μm^2/hr)
     h = D/Δx^2

     # internal nodes
     for j = 2:n2-1
         for i = 2:n1-1
             @inbounds  Δu[i,j] = h*(u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - 4*u[i,j])
         end
     end

     # left/right edges
     for i = 2:n1-1
         @inbounds Δu[i,1] = h*(u[i+1,1] + u[i-1,1] + 2*u[i,2] - 4*u[i,1])
         @inbounds Δu[i,n2] = h*(u[i+1,n2] + u[i-1,n2] + 2*u[i,n2-1] - 4*u[i,n2])
     end

     # top/bottom edges
     for j = 2:n2-1
         @inbounds Δu[1,j] = h*(u[1,j+1] + u[1,j-1] + 2*u[2,j] - 4*u[1,j])
         @inbounds Δu[n1,j] = h*(u[n1,j+1] + u[n1,j-1] + 2*u[n1-1,j] - 4*u[n1,j])
     end

     # corners
     @inbounds Δu[1,1] = h*(2*(u[2,1] + u[1,2]) - 4*u[1,1])
     @inbounds Δu[n1,1] = h*(2*(u[n1-1,1] + u[n1,2]) - 4*u[n1,1])
     @inbounds Δu[1,n2] = h*(2*(u[2,n2] + u[1,n2-1]) - 4*u[1,n2])
     @inbounds Δu[n1,n2] = h*(2*(u[n1-1,n2] + u[n1,n2-1]) - 4*u[n1,n2])
 end

 ###############################################################
 # 3. This is the ODE model for the cGAS pathway
 ###############################################################

 # Define the discretized PDE as an ODE function
 function Model!(du,u,p,t)

   #Species
   cGAS = @view u[:,:,1]
   DNA = @view u[:,:,2]
   Sting = @view u[:,:,3]
   cGAMP = @view u[:,:,4]
   IRF3 = @view u[:,:,5]
   IFNβm = @view u[:,:,6]
   IFNβ = @view u[:,:,7]
   STAT = @view u[:,:,8]
   SOCSm = @view u[:,:,9]
   IRF7m = @view u[:,:,10]
   TREX1m = @view u[:,:,11]
   IRF7 = @view u[:,:,12]
   TREX1 = @view u[:,:,13]
   Virus = @view u[:,:,14]

   #Derivatives
   d_cGAS = @view du[:,:,1]
   d_DNA = @view du[:,:,2]
   d_Sting = @view du[:,:,3]
   d_cGAMP = @view du[:,:,4]
   d_IRF3 = @view du[:,:,5]
   d_IFNβm = @view du[:,:,6]
   d_IFNβ = @view du[:,:,7]
   d_STAT = @view du[:,:,8]
   d_SOCSm = @view du[:,:,9]
   d_IRF7m = @view du[:,:,10]
   d_TREX1m = @view du[:,:,11]
   d_IRF7 = @view du[:,:,12]
   d_TREX1 = @view du[:,:,13]
   d_Virus = @view du[:,:,14]

   #Parameters
   k1f, k1r, k3f, k3r, k4f, kcat5, Km5, k5r, kcat6, Km6, kcat7, Km7, kcat8, Km8, k8f, k9f, k10f1, k10f2, k11f, k12f, k13f, k6f, kcat2, Km2, τ4, τ6, τ7, τ8, τ9, τ10, τ11, τ12, τ13, k14f, τ14 = p.par
   #Constants from the mass balances
   cGAStot, Stingtot, IRF3tot = p.mass
   #Which cells are dead?
   💀 = p.deathParameter
   #Should DNA be allowed to replicate (only with virus, not with ISD)?
   🔁 = p.DNAReplicate

   #Calculate the diffusion of IFNβ
   ∇²(ΔIFNβ,IFNβ)

   #Update derivatives for each species according to model
   @. d_cGAS = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS)
   @. d_DNA = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS) - kcat2*TREX1*DNA / (Km2 + DNA) + 🔁*💀*DNA*(0.55-DNA)/0.55
   @. d_Sting = -k3f*cGAMP*Sting + k3r*(Stingtot - Sting)
   @. d_cGAMP = k4f*(cGAStot - cGAS) - k3f*cGAMP*Sting + k3f*(Stingtot - Sting) - τ4*cGAMP
   @. d_IRF3 = -kcat5*IRF3*(Stingtot - Sting) / (Km5 +IRF3) + k5r*(IRF3tot - IRF3)
   @. d_IFNβm = 💀*kcat6*(IRF3tot - IRF3) / (Km6 + (IRF3tot - IRF3)) + 💀*k6f*IRF7 - τ6*IFNβm
   @. d_IFNβ = 💀*kcat7*IFNβm / (Km7 + IFNβm) - τ7*IFNβ + ΔIFNβ #Add the diffusion in here
   @. d_STAT = 💀*kcat8*IFNβ / (Km8 + IFNβ) * 1.0/(1.0+k8f*SOCSm) - τ8*STAT
   @. d_SOCSm = 💀*k9f*STAT - τ9*SOCSm
   @. d_IRF7m = 💀*k10f1*STAT + 💀*k10f2*IRF7 - τ10*IRF7m
   @. d_TREX1m = 💀*k11f*STAT - τ11*TREX1m
   @. d_IRF7 = 💀*k12f*IRF7m - τ12*IRF7
   @. d_TREX1 = 💀*k13f*TREX1m - τ13*TREX1
   @. d_Virus = 💀*k14f*DNA - τ14*Virus
 end

 ###############################################################
# 4. Set up function that will return an ODE problem to solve
 ###############################################################

function ModelSetup(infectionMethod,IFNStoch,Hetero)
    #Paramter values for the ODEs
    θVals = [2.6899, 4.8505, 0.0356, 7.487, 517.4056, 22328.3852, 11226.3682,0.9341,
             206.9446, 10305.461, 47639.70295,3.8474, 13.006, 78.2048, 0.0209,
             0.0059, 0.001, 0.0112, 0.001, 99.9466, 15.1436,0.0276, 237539.3249,
             61688.259, 0.96, 0.347, 12.2428736,1.2399, 1.5101, 0.347, 0.165, 6.9295,
             0.0178]
    θVirus = [1.0, 1.0] # k14f τ14 (Virus Parameters)
    append!(θVals,θVirus) #Append the virus parameters to the orginal parameters

    #Should IFN signaling be stochastic?
    if IFNStoch == :Stochastic
        #Keep most parameters the same
        θ = fill.(θVals,N,N)
        #kcat8 produces IFN, make it nonzero ~20% of the time (can be changed later)
        θ[11] .= rand([zeros(4)...,θVals[11]],N,N)
    else
        θ = θVals #Just keep the parameters as is (same for each cell)
    end

#Now we need to deal with the initial condition
  #Define the initial conditions
  u0 = zeros(N,N,species)

if Hetero == :Hetero
  noiseDistributions = AddNoise2States(0.5)
  for (i,index) in enumerate(nonZeroSpeciesIdx)
      u0[:,:,index] = rand(noiseDistributions[i],N,N)
  end

else
  #Loop through non-zero species and update their concentrations
  for (idx,val) in zip(nonZeroSpeciesIdx,nonZeroSpeciesValues)
      u0[:,:,idx] .=  val
  end
end

  #Finally need to set the DNA initial condition
  if infectionMethod == :ISD
    #Define a region on the domain where cells will be infected
    circleOrigin = [0,0] #Where is the center of the drop?
    circleRadiusSquared = 200^2 #How big is the drop?
    #Calculate squared distances
    sqDist(x,c) = reduce(+, @. (x-c)^2)
    #Loop though cells and check if they are infected
    for currentCell in cellIndicies
        #Are the cells inside the infected region?
        if sqDist([currentCell[1],currentCell[2]],circleOrigin) <= circleRadiusSquared
            u0[currentCell,2] = m2c(1e3)
        end
    end

  elseif infectionMethod == :Virus
    #Assume a poisson ditribution to randomly choose each cell's level of infection
    probDistInfected = Poisson(moi)
    u0[:,:,2] = @. m2c(1e3*rand(probDistInfected,N,N))
  end

#Need to wrap everything up into the ParContainer struct
  mass = [u0[:,:,i] for i in nonZeroSpeciesIdx]
  deathParameter = ones(N,N)
  DNAReplicate = infectionMethod==:ISD ? 0 : 1
  #Keep track of infected cells (save time when infected, Inf means not infected)
  cellsInfected = fill(Inf,N,N) #Make constant when not testing
  cellsInfected[findall(u0[:,:,2] .> 0.0), 1] .= 0.0
  #Keep track of time of death (TOD)
  cellsDead = fill(Inf,N,N) #Inf implies alive
  #Create an array that keeps track of whether or not a cell has tried to infect neighbors
  infectFirstAttempt = trues(N,N)

#Create an instance of the structures
θ = ParContainer(
  θ, #Rate constants
  mass, #Mass balances
  deathParameter, # cell is dead? (1==alive, 0==dead)
  DNAReplicate, #can DNA replicate? (1 if virus)
  cellsInfected, #Time cell was infected
  cellsDead, #Time cell was killed
  infectFirstAttempt) #Has cell tried to infect neighbors?)


return ODEProblem(Model!,u0,tspan,θ)
end


###############################################################
# 5. Function to count the cell states (healthy,infected,dead)
###############################################################
function cellStates(t,θ)
    #Number of healthy cells at time t
    totaHealthy = sum(θ.cellsInfected .> t)
    #Number of dead cells at time t
    totalDead = sum(θ.cellsDead .< t)
    #Number of infected cells at time t
    totalInfected = nCells - totaHealthy - totalDead
    return [totaHealthy,totalInfected,totalDead]
end
