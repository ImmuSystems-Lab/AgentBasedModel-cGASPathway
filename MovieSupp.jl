using RCall, CSV, DataFrames
include("ProblemGenerator.jl")
include("VirusCallBacks.jl")


function DataRangle(sol)
  timeLength = length(sol.t)
  df = DataFrame()
  df.x = repeat(repeat(1:N,N),timeLength)
  df.y = repeat(repeat(1:N,inner=N),timeLength)
  df.t = repeat(sol.t,inner=N^2)
  df.IFN = vec(sol[:,:,7,:])

  return df
end

function CellStateGrid(θ::ParContainer,t::Float64)
  grid = zeros(Int64,size(θ.cellsInfected))
  grid[θ.cellsInfected .<= t] .= 1
  grid[θ.cellsDead .<= t] .= 2

  return vec(grid)
end

#-----------ISD, Not Stochastic, Homo-----------
probISD_nS_Ho = ModelSetup(:ISD,:notStochastic,:Homo)
    solISD_nS_Ho = @time solve(probISD_nS_Ho,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

ISD_nS_Ho = DataRangle(solISD_nS_Ho)
CSV.write("/home/robertgregg/R/Scripts/ABMMovie/ISD_nS_Ho.csv", ISD_nS_Ho)

#-----------Virus, Not Stochastic, Homo-----------
probVirus_nS_Ho = ModelSetup(:Virus,:notStochastic,:Homo)
    solVirus_nS_Ho = @time solve(probVirus_nS_Ho,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)

Virus_nS_Ho = DataRangle(solVirus_nS_Ho)
CSV.write("/home/robertgregg/R/Scripts/ABMMovie/Virus_nS_Ho.csv", Virus_nS_Ho)

#-----------ISD, Stochastic, Hetero-----------
probISD_S_He = ModelSetup(:ISD,:Stochastic,:Hetero)
    solISD_S_He = @time solve(probISD_S_He,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

ISD_S_He = DataRangle(solISD_S_He)
CSV.write("/home/robertgregg/R/Scripts/ABMMovie/ISD_S_He.csv", ISD_S_He)

#-----------Virus, Stochastic, Hetero-----------
probVirus_S_He = ModelSetup(:Virus,:Stochastic,:Hetero)
    solVirus_S_He = @time solve(probVirus_S_He,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)

Virus_S_He = DataRangle(solVirus_S_He)
CSV.write("/home/robertgregg/R/Scripts/ABMMovie/Virus_S_He.csv", Virus_S_He)


#-----------Cell States----------
timeLength = length(solVirus_nS_Ho.t)
CellStateVirus_nS_Ho = DataFrame()
CellStateVirus_nS_Ho.x = repeat(repeat(1:N,N),timeLength)
CellStateVirus_nS_Ho.y = repeat(repeat(1:N,inner=N),timeLength)
CellStateVirus_nS_Ho.t = repeat(solVirus_S_He.t,inner=N^2)
CellStateVirus_nS_Ho.Cell = vcat([CellStateGrid(probVirus_nS_Ho.p,t) for t in solVirus_nS_Ho.t]...)
CSV.write("/home/robertgregg/R/Scripts/ABMMovie/CellStateVirus_nS_Ho.csv", CellStateVirus_nS_Ho)



CellStateVirus_S_He = DataFrame()
CellStateVirus_S_He.x = repeat(repeat(1:N,N),timeLength)
CellStateVirus_S_He.y = repeat(repeat(1:N,inner=N),timeLength)
CellStateVirus_S_He.t = repeat(solVirus_S_He.t,inner=N^2)
CellStateVirus_S_He.Cell = vcat([CellStateGrid(probVirus_S_He.p,t) for t in solVirus_S_He.t]...)
CSV.write("/home/robertgregg/R/Scripts/ABMMovie/CellStateVirus_S_He.csv", CellStateVirus_S_He)
