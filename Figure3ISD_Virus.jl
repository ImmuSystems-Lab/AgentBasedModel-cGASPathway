using RCall
include("ProblemGenerator.jl")
include("VirusCallBacks.jl")

########################################################################
#Plotting States (Figure 3)
########################################################################

#--------------ISD (Det)--------------
ISDprob = ModelSetup(:ISD, :notStochastic, :Homo)
ISDsol = @time solve(ISDprob, CVODE_BDF(linear_solver = :GMRES), saveat = 0.1)

allStateDataISDDet = DataFrame()
for (i, st) in enumerate(statesNames)
  allStateDataISDDet[!, Symbol(st)] = vec(ISDsol[:, :, i, :])
end

timeLength = length(ISDsol.t)
allStateDataISDDet.Cell = repeat(1:nCells, timeLength)
allStateDataISDDet.Time = repeat(ISDsol.t, inner = nCells)
allStateDataISDDet.Infected = repeat(vec(ISDprob.u0[:, :, 2] .> 0.0), timeLength)

#--------------ISD (Stoch)--------------
ISDprob = ModelSetup(:ISD, :Stochastic, :Hetero)
ISDsol = @time solve(ISDprob, CVODE_BDF(linear_solver = :GMRES), saveat = 0.1)

allStateDataISDStoch = DataFrame()
for (i, st) in enumerate(statesNames)
  allStateDataISDStoch[!, Symbol(st)] = vec(ISDsol[:, :, i, :])
end

timeLength = length(ISDsol.t)
allStateDataISDStoch.Cell = repeat(1:nCells, timeLength)
allStateDataISDStoch.Time = repeat(ISDsol.t, inner = nCells)
allStateDataISDStoch.Infected = repeat(vec(ISDprob.u0[:, :, 2] .> 0.0), timeLength)


#--------------Virus--------------
Virusprob = ModelSetup(:Virus, :Stochastic, :Hetero)
Virussol = @time solve(
  Virusprob,
  CVODE_BDF(linear_solver = :GMRES),
  saveat = 0.1,
  callback = cb,
)

allStateDataVirus = DataFrame()
for (i, st) in enumerate(statesNames)
  allStateDataVirus[!, Symbol(st)] = vec(Virussol[:, :, i, :])
end

timeLength = length(Virussol.t)
allStateDataVirus.Cell = repeat(1:nCells, timeLength)
allStateDataVirus.Time = repeat(Virussol.t, inner = nCells)
allStateDataVirus.Infected = repeat(vec(Virusprob.u0[:, :, 2] .> 0.0), timeLength)



@rput allStateDataISDDet
@rput allStateDataISDStoch
@rput allStateDataVirus
R"""
library(ggplot2)
library(ggpubr)

#--------------ISD (Det)--------------
stateAveISDDet = aggregate(allStateDataISDDet[,1:14], list(allStateDataISDDet$Infected,allStateDataISDDet$Time), mean)
colnames(stateAveISDDet)[1] <- "Cell.State"
colnames(stateAveISDDet)[2] <- "Time"

logic<- unlist(lapply(stateAveISDDet$Cell.State, function(x) x == TRUE))
stateAveISDDet$Cell.State[which(logic == TRUE)] <- "Primary"
stateAveISDDet$Cell.State[which(logic == FALSE)] <- "Secondary"
lowISDDet = aggregate(allStateDataISDDet[,1:14], list(allStateDataISDDet$Infected,allStateDataISDDet$Time), FUN = 'quantile',probs=0.05)
highISDDet = aggregate(allStateDataISDDet[,1:14], list(allStateDataISDDet$Infected,allStateDataISDDet$Time), FUN = 'quantile',probs=0.95)

#--------------ISD (Stoch)--------------
stateAveISDStoch = aggregate(allStateDataISDStoch[,1:14], list(allStateDataISDStoch$Infected,allStateDataISDStoch$Time), mean)
colnames(stateAveISDStoch)[1] <- "Cell.State"
colnames(stateAveISDStoch)[2] <- "Time"

logic<- unlist(lapply(stateAveISDStoch$Cell.State, function(x) x == TRUE))
stateAveISDStoch$Cell.State[which(logic == TRUE)] <- "Primary"
stateAveISDStoch$Cell.State[which(logic == FALSE)] <- "Secondary"
lowISDStoch = aggregate(allStateDataISDStoch[,1:14], list(allStateDataISDStoch$Infected,allStateDataISDStoch$Time), FUN = 'quantile',probs=0.05)
highISDStoch = aggregate(allStateDataISDStoch[,1:14], list(allStateDataISDStoch$Infected,allStateDataISDStoch$Time), FUN = 'quantile',probs=0.95)

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
p1 <- ggplot(stateAveISDDet) + geom_line(aes(y=DNA, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDDet$DNA, ymax=highISDDet$DNA, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("DNA") +
  commonFigureOptions

p2 <- ggplot(stateAveISDDet) + geom_line(aes(y=IFNb, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDDet$IFNb, ymax=highISDDet$IFNb, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IFNb") +
  commonFigureOptions

p3 <- ggplot(stateAveISDDet) + geom_line(aes(y=IRF7, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDDet$IRF7, ymax=highISDDet$IRF7, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IRF7") +
  commonFigureOptions

p4 <- ggplot(stateAveISDDet) + geom_line(aes(y=TREX1, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDDet$TREX1, ymax=highISDDet$TREX1, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("TREX1") +
  commonFigureOptions

p5 <- ggplot(stateAveISDStoch) + geom_line(aes(y=DNA, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDStoch$DNA, ymax=highISDStoch$DNA, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("DNA") +
  commonFigureOptions

p6 <- ggplot(stateAveISDStoch) + geom_line(aes(y=IFNb, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDStoch$IFNb, ymax=highISDStoch$IFNb, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IFNb") +
  commonFigureOptions

p7 <- ggplot(stateAveISDStoch) + geom_line(aes(y=IRF7, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDStoch$IRF7, ymax=highISDStoch$IRF7, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IRF7") +
  commonFigureOptions

p8 <- ggplot(stateAveISDStoch) + geom_line(aes(y=TREX1, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowISDStoch$TREX1, ymax=highISDStoch$TREX1, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("TREX1") +
  commonFigureOptions

#--------------Virus Plots--------------
p9 <- ggplot(stateAveVirus) + geom_line(aes(y=DNA, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$DNA, ymax=highVirus$DNA, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("DNA") +
  commonFigureOptions

p10 <- ggplot(stateAveVirus) + geom_line(aes(y=IFNb, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$IFNb, ymax=highVirus$IFNb, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IFNb") +
  commonFigureOptions

p11 <- ggplot(stateAveVirus) + geom_line(aes(y=IRF7, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$IRF7, ymax=highVirus$IRF7, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IRF7") +
  commonFigureOptions

p12 <- ggplot(stateAveVirus) + geom_line(aes(y=TREX1, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=lowVirus$TREX1, ymax=highVirus$TREX1, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("TREX1") +
  commonFigureOptions

figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
                    labels = "AUTO",
                    common.legend = TRUE, legend = "right",
                    align = "hv",
                    ncol = 4, nrow = 3)

ggsave("./Figures/Figure3Stoch.pdf",width = 15,height=8,units="in")
"""
