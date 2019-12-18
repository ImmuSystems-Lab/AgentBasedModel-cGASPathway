using RCall
include("ProblemGenerator.jl")
include("VirusCallBacks.jl")

########################################################################
#Plotting States (Figure 3)
########################################################################

#--------------ISD--------------
ISDprob = ModelSetup(:ISD,:notStochastic)
  ISDsol = @time solve(ISDprob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

allStateDataISD = DataFrame()
for (i,st) in enumerate(statesNames)
    allStateDataISD[!,Symbol(st)] = vec(ISDsol[:,:,i,:])
end

timeLength = length(ISDsol.t)
allStateDataISD.Cell = repeat(1:nCells,timeLength)
allStateDataISD.Time = repeat(ISDsol.t,inner=nCells)
allStateDataISD.Infected = repeat(vec(ISDprob.u0[:,:,2] .> 0.0),timeLength)

#--------------Virus--------------
Virusprob = ModelSetup(:Virus,:Stochastic)
  Virussol = @time solve(Virusprob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1,callback=cb)

allStateDataVirus = DataFrame()
for (i,st) in enumerate(statesNames)
    allStateDataVirus[!,Symbol(st)] = vec(Virussol[:,:,i,:])
end

timeLength = length(Virussol.t)
allStateDataVirus.Cell = repeat(1:nCells,timeLength)
allStateDataVirus.Time = repeat(Virussol.t,inner=nCells)
allStateDataVirus.Infected = repeat(vec(Virusprob.u0[:,:,2] .> 0.0),timeLength)



@rput allStateDataISD
@rput allStateDataVirus
R"""
library(ggplot2)
library(ggpubr)

#--------------ISD--------------
stateAveISD = aggregate(allStateDataISD[,1:14], list(allStateDataISD$Infected,allStateDataISD$Time), mean)
colnames(stateAveISD)[1] <- "Cell.State"
colnames(stateAveISD)[2] <- "Time"

logic<- unlist(lapply(stateAveISD$Cell.State, function(x) x == TRUE))
stateAveISD$Cell.State[which(logic == TRUE)] <- "Primary"
stateAveISD$Cell.State[which(logic == FALSE)] <- "Secondary"
lowISD = aggregate(allStateDataISD[,1:14], list(allStateDataISD$Infected,allStateDataISD$Time), FUN = 'quantile',probs=0.05)
highISD = aggregate(allStateDataISD[,1:14], list(allStateDataISD$Infected,allStateDataISD$Time), FUN = 'quantile',probs=0.95)

#--------------Virus--------------
stateAveVirus = aggregate(allStateDataVirus[,1:14], list(allStateDataVirus$Infected,allStateDataVirus$Time), mean)
colnames(stateAveVirus)[1] <- "Cell.State"
colnames(stateAveVirus)[2] <- "Time"

logic<- unlist(lapply(stateAveVirus$Cell.State, function(x) x == TRUE))
stateAveVirus$Cell.State[which(logic == TRUE)] <- "Primary"
stateAveVirus$Cell.State[which(logic == FALSE)] <- "Secondary"
lowVirus = aggregate(allStateDataVirus[,1:14], list(allStateDataVirus$Infected,allStateDataVirus$Time), FUN = 'quantile',probs=0.05)
highVirus = aggregate(allStateDataVirus[,1:14], list(allStateDataVirus$Infected,allStateDataVirus$Time), FUN = 'quantile',probs=0.95)



commonFigureOptions <- list( scale_x_continuous(breaks=seq(0, 48, 12)),
  theme_pubr(border=TRUE),
  ylab("nM"),
  xlab("Time (hours)"),
  theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))


#--------------ISD Plots--------------
p1 <- ggplot(stateAveISD) + geom_line(aes(y=DNA, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISD$DNA, ymax=highISD$DNA, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("DNA") +
  commonFigureOptions

p2 <- ggplot(stateAveISD) + geom_line(aes(y=IFNb, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISD$IFNb, ymax=highISD$IFNb, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IFNb") +
  commonFigureOptions

p3 <- ggplot(stateAveISD) + geom_line(aes(y=IRF7, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISD$IRF7, ymax=highISD$IRF7, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IRF7") +
  commonFigureOptions

p4 <- ggplot(stateAveISD) + geom_line(aes(y=TREX1, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISD$TREX1, ymax=highISD$TREX1, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("TREX1") +
  commonFigureOptions

#--------------Virus Plots--------------
p5 <- ggplot(stateAveVirus) + geom_line(aes(y=DNA, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$DNA, ymax=highVirus$DNA, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("DNA") +
  commonFigureOptions

p6 <- ggplot(stateAveVirus) + geom_line(aes(y=IFNb, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$IFNb, ymax=highVirus$IFNb, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IFNb") +
  commonFigureOptions

p7 <- ggplot(stateAveVirus) + geom_line(aes(y=IRF7, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$IRF7, ymax=highVirus$IRF7, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IRF7") +
  commonFigureOptions

p8 <- ggplot(stateAveVirus) + geom_line(aes(y=TREX1, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$TREX1, ymax=highVirus$TREX1, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("TREX1") +
  commonFigureOptions

figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
                    labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                    common.legend = TRUE, legend = "right",
                    align = "hv",
                    ncol = 2, nrow = 4)

ggsave("./Figures/Figure3ISD_Virus.pdf",width = 12,height=6,units="in")
"""
