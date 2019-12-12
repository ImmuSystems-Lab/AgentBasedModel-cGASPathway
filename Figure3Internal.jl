using RCall
include("ProblemGenerator.jl")


########################################################################
#Plotting States (Figure 3)
########################################################################
probFigure3 = ModelSetup(:ISD,:notStochastic)
  solFigure3 = @time solve(probFigure3,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

allStateData = DataFrame()
for (i,st) in enumerate(statesNames)
    allStateData[!,Symbol(st)] = vec(solFigure3[:,:,i,:])
end

timeLength = length(solFigure3.t)
allStateData.Cell = repeat(1:nCells,timeLength)
allStateData.Time = repeat(solFigure3.t,inner=nCells)
allStateData.Infected = repeat(vec(probFigure3.u0[:,:,2] .> 0.0),timeLength)

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

ggsave("./Figures/Figure3.pdf")
"""
