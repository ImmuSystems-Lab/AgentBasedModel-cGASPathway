#Shows what a typical ISD and Viral Infections looks like
using RCall
include("ProblemGenerator.jl")
include("VirusCallBacks.jl")
########################################################################
#Heatmaps of IFN dynamics (Figure 2)
########################################################################

#-----------ISD, Not Stochastic, Homo-----------
probISD_nS_Ho = ModelSetup(:ISD,:notStochastic,:Homo)
    solISD_nS_Ho = @time solve(probISD_nS_Ho,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
ISD_nS_Ho = DataFrame()
ISD_nS_Ho.x = repeat(1:N,N)
ISD_nS_Ho.y = repeat(1:N,inner=N)
ISD_nS_Ho.IFN = vec(solISD_nS_Ho(10.0)[:,:,7])

#-----------Virus, Not Stochastic, Homo-----------
probVirus_nS_Ho = ModelSetup(:Virus,:notStochastic,:Homo)
    solVirus_nS_Ho = @time solve(probVirus_nS_Ho,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)
Virus_nS_Ho = DataFrame()
Virus_nS_Ho.x = repeat(1:N,N)
Virus_nS_Ho.y = repeat(1:N,inner=N)
Virus_nS_Ho.IFN = vec(solVirus_nS_Ho(10.0)[:,:,7])

#-----------Virus, Not Stochastic, Homo, Cell States-----------
function CellStateGrid(θ::ParContainer,t::Float64)
  grid = zeros(Int64,size(θ.cellsInfected))
  grid[θ.cellsInfected .<= t] .= 1
  grid[θ.cellsDead .<= t] .= 2

  return grid
end

CellStateVirus_nS_Ho_raw = CellStateGrid(probVirus_nS_Ho.p,10.0)
CellStateVirus_nS_Ho = DataFrame()
CellStateVirus_nS_Ho.x = repeat(1:N,N)
CellStateVirus_nS_Ho.y = repeat(1:N,inner=N)
CellStateVirus_nS_Ho.Cell = vec(CellStateVirus_nS_Ho_raw)



#-----------ISD, Stochastic, Hetero-----------
probISD_S_He = ModelSetup(:ISD,:Stochastic,:Hetero)
    solISD_S_He = @time solve(probISD_S_He,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
ISD_S_He = DataFrame()
ISD_S_He.x = repeat(1:N,N)
ISD_S_He.y = repeat(1:N,inner=N)
ISD_S_He.IFN = vec(solISD_S_He(10.0)[:,:,7])

#-----------Virus, Stochastic, Hetero-----------
probVirus_S_He = ModelSetup(:Virus,:Stochastic,:Hetero)
    solVirus_S_He = @time solve(probVirus_S_He,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)
Virus_S_He = DataFrame()
Virus_S_He.x = repeat(1:N,N)
Virus_S_He.y = repeat(1:N,inner=N)
Virus_S_He.IFN = vec(solVirus_S_He(10.0)[:,:,7])

#-----------Virus, Not Stochastic, Homo, Cell States-----------
function CellStateGrid(θ::ParContainer,t::Float64)
  grid = zeros(Int64,size(θ.cellsInfected))
  grid[θ.cellsInfected .<= t] .= 1
  grid[θ.cellsDead .<= t] .= 2

  return grid
end

CellStateVirus_S_He_raw = CellStateGrid(probVirus_S_He.p,10.0)
CellStateVirus_S_He = DataFrame()
CellStateVirus_S_He.x = repeat(1:N,N)
CellStateVirus_S_He.y = repeat(1:N,inner=N)
CellStateVirus_S_He.Cell = vec(CellStateVirus_S_He_raw)




#-----------R Plot-----------
@rput ISD_nS_Ho
@rput Virus_nS_Ho
@rput CellStateVirus_nS_Ho
@rput ISD_S_He
@rput Virus_S_He
@rput CellStateVirus_S_He
R"""
  library(ggplot2)
  library(ggpubr)

  #To draw an arc on the heatmap to show initial condition
  arc <- data.frame(
  x0 = seq(0, 100, length.out = 1000),
  y0 = sqrt(100^2 - (seq(0, 100, length.out = 1000)^2))
  )

  #Color Range
  clim = c(min(ISD_S_He$IFN), max(ISD_nS_Ho$IFN))
  commonFigureOptions <- list(scale_fill_distiller(palette = "Spectral",guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),limits = clim),
    scale_x_continuous(expand=c(0,0)),
    scale_y_continuous(expand=c(0,0)),
    theme_bw(base_size = 12),
    ylab("Cell"),
    xlab("Cell"),
    labs(fill="IFN (nM)"),
    theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

  p1 <- ggplot() +
    geom_raster(ISD_nS_Ho, mapping =aes(x, y, fill=IFN)) +
    geom_line(data=arc, aes(x=x0, y=y0),size=0.2,linetype="dashed") +
    labs(fill="IFN (nM)") +
    ggtitle("ISD Transfection") +
    commonFigureOptions

  p2 <- ggplot(Virus_nS_Ho, aes(x, y, fill=IFN)) +
    geom_raster(aes(fill=IFN)) +
    ggtitle("Virus Infection") +
    commonFigureOptions

  p3 <- ggplot(CellStateVirus_nS_Ho, aes(x, y, fill=factor(Cell))) +
    geom_raster(aes(fill=factor(Cell))) +
    scale_fill_manual(values=c("#377EB8","#4DAF4A","#E41A1C"),labels = c("Healthy", "Infected","Dead")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    labs(title ="Cell States", x="Cell", y="Cell", fill = "Cell State") +
    theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1)

  p4 <- ggplot() +
    geom_raster(ISD_S_He, mapping =aes(x, y, fill=IFN)) +
    geom_line(data=arc, aes(x=x0, y=y0),size=0.2,linetype="dashed") +
    labs(fill="IFN (nM)") +
    commonFigureOptions

  p5 <- ggplot(Virus_S_He, aes(x, y, fill=IFN)) +
    geom_raster(aes(fill=IFN)) +
    commonFigureOptions

  p6 <- ggplot(CellStateVirus_S_He, aes(x, y, fill=factor(Cell))) +
    geom_raster(aes(fill=factor(Cell))) +
    scale_fill_manual(values=c("#377EB8","#4DAF4A","#E41A1C"),labels = c("Healthy", "Infected","Dead")) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    labs( x="Cell", y="Cell", fill = "Cell State") +
    theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1)

  figure1 <- ggarrange(p1, p2, p3, p4, p5, p6,
                      labels = "AUTO",
                      align = "h",
                      ncol = 3, nrow = 2)

  ggsave("./Figures/Figure2.pdf",width=10,height=5,units="in")
"""


stateToPlot = 7
plotState=[solVirussh[coord,stateToPlot,:] for coord in cellIndicies]
plot(solVirussh.t,plotState[:],leg=false,framestyle=:box,xtickfontsize=16,ytickfontsize=16)
title!(statesNames[stateToPlot])
xlabel!("Time (hrs)")
ylabel!("(nM)")


#Get the max value for IFNβ concentration for whole simulation
Maxval = maximum(solVirussh[:,:,7,:])
#Loop through time points to make an animation
anim = @animate for i = tspan[1]:0.1:tspan[2]
    heatmap(solVirussh(i)[:,:,7],clims=(0.0,Maxval),
            title="Time = " * string(i) * " hrs")
end

gif(anim,"./Figures/IFNAnimation.gif")
