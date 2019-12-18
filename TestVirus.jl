using RCall

include("ProblemGenerator.jl")
include("VirusCallBacks.jl")

prob = ModelSetup(:Virus,:Stochastic)

sol = solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)

#-----------IFN Animation-----------

#Get the max value for IFNβ concentration for whole simulation
Maxval = maximum(sol[:,:,7,:])
#Loop through time points to make an animation
anim = @animate for i = tspan[1]:0.1:tspan[2]
    heatmap(sol(i)[:,:,7],clims=(0.0,Maxval),
            title="Time = " * string(i) * " hrs")
end

gif(anim,"./Figures/IFNAnimation.gif")



#-----------Cell States Plots-----------
saveTimePoints = range(tspan[1],tspan[2],step=0.1)
allCellStates = zeros(Int64,length(saveTimePoints),3)

for k=1:length(saveTimePoints)
    allCellStates[k,:] = cellStates(saveTimePoints[k],prob.p)
end

plot(saveTimePoints,allCellStates./nCells,framestyle = :box,linewidth=2)


#-----------IRNb State Plots-----------

allStateData = DataFrame()
for (i,st) in enumerate(statesNames)
    allStateData[!,Symbol(st)] = vec(sol[:,:,i,:])
end

timeLength = length(sol.t)
allStateData.Cell = repeat(1:nCells,timeLength)
allStateData.Time = repeat(sol.t,inner=nCells)
allStateData.Infected = repeat(vec(prob.u0[:,:,2] .> 0.0),timeLength)

@rput allStateData
R"""
library(ggplot2)
library(ggpubr)

stateAve = aggregate(allStateData[,1:14], list(allStateData$Infected,allStateData$Time), mean)
colnames(stateAve)[1] <- "Cell.State"
colnames(stateAve)[2] <- "Time"

logic<- unlist(lapply(stateAve$Cell.State, function(x) x == TRUE))
stateAve$Cell.State[which(logic == TRUE)] <- "Primary"
stateAve$Cell.State[which(logic == FALSE)] <- "Secondary"
low = aggregate(allStateData[,1:14], list(allStateData$Infected,allStateData$Time), FUN = 'quantile',probs=0.05)
high = aggregate(allStateData[,1:14], list(allStateData$Infected,allStateData$Time), FUN = 'quantile',probs=0.95)

commonFigureOptions <- list( scale_x_continuous(breaks=seq(0, 48, 12)),
  theme_pubr(border=TRUE),
  ylab("nM"),
  xlab("Time (hours)")
  theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

p1 <- ggplot(stateAve) + geom_line(aes(y=DNA, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$DNA, ymax=high$DNA, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("DNA") +
  commonFigureOptions

p2 <- ggplot(stateAve) + geom_line(aes(y=IFNb, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$IFNb, ymax=high$IFNb, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IFNb") +
  commonFigureOptions

p3 <- ggplot(stateAve) + geom_line(aes(y=IRF7, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$IRF7, ymax=high$IRF7, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IRF7") +
  commonFigureOptions

p4 <- ggplot(stateAve) + geom_line(aes(y=TREX1, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$TREX1, ymax=high$TREX1, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("TREX1") +
  commonFigureOptions

figure <- ggarrange(p1, p2, p3, p4,
                    labels = c("A", "B", "C", "D"),
                    common.legend = TRUE, legend = "right",
                    align = "hv",
                    ncol = 2, nrow = 2)

ggsave("./Figures/Figure2Stoch.pdf")
"""

#-----------HeatMap Cell States -----------

function CellStateGrid(θ::ParContainer,t::Float64)
  grid = zeros(Int64,size(θ.cellsInfected))
  grid[θ.cellsInfected .<= t] .= 1
  grid[θ.cellsDead .<= t] .= 2

  return grid
end


anim = @animate for i = tspan[1]:0.1:tspan[2]
    heatmap(CellStateGrid(prob.p,i),clims=(0.0,2.0),
            title="Time = " * string(i) * " hrs")
end

gif(anim,"./Figures/CellStateGrid.gif")
