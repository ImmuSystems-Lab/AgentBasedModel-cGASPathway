using RCall

#Run the KD simulations in parallel
#change the number of cores used as needed
using Distributed
addprocs(2)
#The @everywhere macro runs the command/defines the thing on all the processors
@everywhere include("ProblemGenerator.jl")

########################################################################
#Knockdown simulations of IRF7 and TREX1 (Figure 4)
########################################################################

#Putting the calculations in a function will allow for easy parallelization
@everywhere function KnockDownSim(par2Change,prob)
    #Percentages to run the KD simulations
    KnockDownVals=1:-0.25:0
    percentLabels = [string(i)*"%" for i=0:25:100]
    kdSamples = length(KnockDownVals)

    #Vector to save all IFN dynamics from simulation
    solKD =  Vector(undef,kdSamples)

    #Loop through the KD percents and solve ODEs
    for (i,percent) in enumerate(KnockDownVals)
        #Copy the parameters and change desired one
        θcurrent = deepcopy(prob.p)
        #Check if running ABM or ODE
        if θcurrent isa Array #Are the parameters from the ODE model?
            θcurrent[par2Change] = θcurrent[par2Change] * percent
            #Print to track progress
            println(θcurrent[[20,23]])
        else #Oh nope its the ABM model then
            θcurrent.par[par2Change] = θcurrent.par[par2Change] * percent
            #Print to track progress
            println(θcurrent.par[[20,23]])
        end

        #Redefine problem with the new parameters and solve
        probKD = remake(prob; p=θcurrent)
        sol = solve(probKD,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
        #Store the IFN concentrations
        solKD[i] = θcurrent isa Array ? vec(sol[7,:]) : vec(sol[:,:,7,:])
    end

    #Return the IFN dynamics
    return solKD
end

########################################################################
#Knockdown ABM simulations
########################################################################

#Define a problem
@everywhere prob = ModelSetup(:ISD,:notStochastic)
#Run the KD simulation
KDSols = pmap(x-> KnockDownSim(x,prob),[20,23])

#Save the simulation in a Dataframe to pass on to R (for plotting)
KDData = DataFrame()

#Seperate out the columns for each KD
KDData.IRF7KD = vcat(KDSols[1]...)
KDData.TREXKD = vcat(KDSols[2]...)

#Create columns that keep track of cell ID, time, and KD percent
timeLength = 481
KnockDownVals=1:-0.25:0
percentLabels = [string(i)*"%" for i=0:25:100]
kdSamples = length(KnockDownVals)

#Put everything into the DataFrame
KDData.Cell = repeat(1:nCells,timeLength*kdSamples)
KDData.Time = repeat(0:0.1:48,inner=nCells,outer=kdSamples)
KDData.Percent = repeat(KnockDownVals,inner=nCells*timeLength)

########################################################################
#Knockdown ODE simulations
########################################################################

#ODE IC
ODEu0 = zeros(species)
ODEu0[[1,2,3,5]] .= m2c([1e3,1e3,1e3, 1e4])

#ODE Pars
ODEpars = [2.6899, 4.8505, 0.0356, 7.487, 517.4056, 22328.3852, 11226.3682,0.9341,
         206.9446, 10305.461, 47639.70295,3.8474, 13.006, 78.2048, 0.0209,
         0.0059, 0.001, 0.0112, 0.001, 99.9466, 15.1436,0.0276, 237539.3249,
         61688.259, 0.96, 0.347, 12.2428736,1.2399, 1.5101, 0.347, 0.165, 6.9295,
         0.0178]
 #Append  cGAStot, Stingtot, IRF3tot to the parameters
 append!(ODEpars,ODEu0[[1,3,5]])

 #ODE Model
 function ODEmodel(du,u,p,t)
   #Species
   cGAS, DNA, Sting, cGAMP, IRF3, IFNβm, IFNβ, STAT, SOCSm, IRF7m, TREX1m, IRF7, TREX1 = u

   #Parameters
   k1f, k1r, k3f, k3r, k4f, kcat5, Km5, k5r, kcat6, Km6, kcat7, Km7, kcat8, Km8, k8f, k9f, k10f1, k10f2, k11f, k12f, k13f, k6f, kcat2, Km2, τ4, τ6, τ7, τ8, τ9, τ10, τ11, τ12, τ13, cGAStot, Stingtot, IRF3tot = p

   #Update derivatives for each species according to model
   du[1] = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS)
   du[2] = -k1f*cGAS*DNA + k1r*(cGAStot - cGAS) - kcat2*TREX1*DNA / (Km2 + DNA)
   du[3] = -k3f*cGAMP*Sting + k3r*(Stingtot - Sting)
   du[4] = k4f*(cGAStot - cGAS) - k3f*cGAMP*Sting + k3f*(Stingtot - Sting) - τ4*cGAMP
   du[5] = -kcat5*IRF3*(Stingtot - Sting) / (Km5 +IRF3) + k5r*(IRF3tot - IRF3)
   du[6] = kcat6*(IRF3tot - IRF3) / (Km6 + (IRF3tot - IRF3)) + k6f*IRF7 - τ6*IFNβm
   du[7] = kcat7*IFNβm / (Km7 + IFNβm) - τ7*IFNβ
   du[8] = kcat8*IFNβ / (Km8 + IFNβ) * 1.0/(1.0+k8f*SOCSm) - τ8*STAT
   du[9] = k9f*STAT - τ9*SOCSm
   du[10] = k10f1*STAT + k10f2*IRF7 - τ10*IRF7m
   du[11] = k11f*STAT - τ11*TREX1m
   du[12] = k12f*IRF7m - τ12*IRF7
   du[13] = k13f*TREX1m - τ13*TREX1
 end

 ODEprob = ODEProblem(ODEmodel,ODEu0,tspan,ODEpars)

 #Run the KD simulation (THis could also be parallelized, but its really fast)
 ODEsols = map(x-> KnockDownSim(x,ODEprob),[20,23])

#GEnerate the DataFrame for the ODE results
 KDDataODE = DataFrame()
 KDDataODE.IRF7KD = vcat(ODEsols[1]...)
 KDDataODE.TREXKD = vcat(ODEsols[2]...)
 KDDataODE.Time = repeat(0:0.1:48,outer=kdSamples)
 KDDataODE.Percent = repeat(KnockDownVals,inner=timeLength)

 ########################################################################
 #Plotting the ABM and ODE solutions
 ########################################################################

 @rput KDData
 @rput KDDataODE
 R"""
 library(ggplot2)
 library(ggpubr)

 stateAve = aggregate(KDData[,1:2], list(KDData$Percent,KDData$Time), mean)
 colnames(stateAve)[1] <- "Percent"
 colnames(stateAve)[2] <- "Time"


 stateSD = aggregate(KDData[,1:2], list(KDData$Percent,KDData$Time), sd)
 low = stateAve[,3:4] -  stateSD[,3:4]
 high = stateAve[,3:4] + stateSD[,3:4]

 commonFigureOptions <- list(scale_x_continuous(breaks=seq(0, 48, 12)),
   theme_pubr(border=TRUE),
   xlab("Time (hours)"),
   theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

 p1 <- ggplot(KDDataODE) + geom_line(aes(y=TREXKD, x=Time, group=factor(Percent), color = factor(Percent))) +
   ggtitle("ODE Model: \n TREX1 Knockdown") +
   ylab("IFN (nM)") +
   commonFigureOptions

 p2 <- ggplot(KDDataODE) + geom_line(aes(y=IRF7KD, x=Time, group=factor(Percent), color = factor(Percent))) +
   ggtitle("ODE Model: \n IRF7 Knockdown") +
   ylab("IFN (nM)") +
   commonFigureOptions

 p3 <- ggplot(stateAve) + geom_line(aes(y=TREXKD, x=Time, group=factor(Percent), color = factor(Percent))) +
   geom_ribbon(aes(ymin=low$TREXKD, ymax=high$TREXKD, x=Time,group=factor(Percent),fill = factor(Percent)), alpha = 0.2) +
   ggtitle("ABM: \n TREX1 Knockdown") +
   ylab("Average IFN (nM)") +
   commonFigureOptions

 p4 <- ggplot(stateAve) + geom_line(aes(y=IRF7KD, x=Time, group=factor(Percent), color = factor(Percent))) +
   geom_ribbon(aes(ymin=low$IRF7KD, ymax=high$IRF7KD, x=Time,group=factor(Percent),fill = factor(Percent)), alpha = 0.2) +
   ggtitle("ABM: \n IRF7 Knockdown") +
   ylab("Average IFN (nM)") +
   commonFigureOptions

   figure <- ggarrange(p1, p2, p3, p4,
                       labels = c("A", "B", "C", "D"),
                       common.legend = TRUE, legend = "right",
                       align = "hv",
                       ncol = 2, nrow = 2)
 
 ggsave("./Figures/Figure4.pdf")
 """
