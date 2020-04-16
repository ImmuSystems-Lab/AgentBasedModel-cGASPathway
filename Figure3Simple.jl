#Shows bulk interacellular concentration profiles
using RCall
include("ProblemGenerator.jl")
include("VirusCallBacks.jl")


#Create a talbe to hold all the simulation data


#Create a function to output a table
function TableMaker(prob,sol,probInfo)

  timeLength = length(sol.t)

  T = DataFrame()
  T.ModelType = fill(probInfo,nCells*timeLength)


  T.CellType = repeat([i > 0.0 ? "Primary" : "Secondary" for i in vec(prob.u0[:, :, 2])], timeLength)
  T.Time = repeat(sol.t, inner = nCells)

  wantedSpecies = ["DNA","IFNb","IRF7","TREX1"]
  wantedSpeciesIdx = [findfirst(isequal(str),statesNames) for str in wantedSpecies]

  for (i, st) in zip(wantedSpeciesIdx,wantedSpecies)
    T[!, Symbol(st)] = vec(sol[:, :, i, :])
  end

  return T
end

#--------------ISD (Det, homo)--------------
ISDprob = ModelSetup(:ISD, :notStochastic, :Homo)
ISDsol = @time solve(ISDprob, CVODE_BDF(linear_solver = :GMRES), saveat = 1.0)

probInfo = "ISD/Det/Homo"
fig3Table = TableMaker(ISDprob,ISDsol,probInfo)

#--------------Virus (Det, homo)--------------
Virprob = ModelSetup(:Virus, :notStochastic, :Homo)
Virsol = @time solve(Virprob, CVODE_BDF(linear_solver = :GMRES), saveat = 1.0,callback = cb)

probInfo = "Virus/Det/Homo"
fig3Table = vcat(fig3Table,TableMaker(Virprob,Virsol,probInfo))



#--------------ISD (Stoch, hetero)--------------
ISDprob2 = ModelSetup(:ISD, :Stochastic, :Hetero)
ISDsol2 = @time solve(ISDprob2, CVODE_BDF(linear_solver = :GMRES), saveat = 1.0)

probInfo = "ISD/Stoch/Hetero"
fig3Table = vcat(fig3Table,TableMaker(ISDprob2,ISDsol2,probInfo))

#--------------Virus (Stoch, hetero)--------------
Virprob2 = ModelSetup(:Virus, :Stochastic, :Hetero)
Virsol2 = @time solve(Virprob2, CVODE_BDF(linear_solver = :GMRES), saveat = 1.0,callback = cb)

probInfo = "Virus/Stoch/Hetero"
fig3Table = vcat(fig3Table,TableMaker(Virprob2,Virsol2,probInfo))

@rput fig3Table
R"""
library(tidyverse)
library(ggpubr)
library(scales)


df_tidy <- fig3Table %>%
  pivot_longer(c(DNA,IFNb,IRF7,TREX1),names_to = "SpeciesNames",values_to = "SpeciesValues") %>%
  group_by(ModelType,CellType,SpeciesNames,Time) %>%
  summarise(AveState = mean(SpeciesValues),
            q05 = quantile(SpeciesValues,0.05),
            q95 = quantile(SpeciesValues,0.95))

CustomPlotter <- function(dfIN,typeIn, nameIn){

  df <- dfIN %>%
    filter(ModelType == typeIn,SpeciesNames == nameIn)

  p <- ggplot(df, aes(x=Time, y=AveState)) +
      geom_line(aes(x=Time, y=AveState, color=CellType)) +
      geom_ribbon(aes(ymin=q05, ymax=q95, fill=CellType),alpha=0.3) +
      scale_x_continuous(breaks=seq(0, 48, 12)) +
      ggtitle(nameIn) +
      xlab("Time (hours)") +
      ylab("nM") +
      theme_pubr(border=TRUE) +
      theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1)


return(p)
}

# I know facet exists but it wasn't working how I wanted...

p1 <- CustomPlotter(df_tidy,"ISD/Det/Homo", "DNA")
p2 <- CustomPlotter(df_tidy,"ISD/Det/Homo", "IFNb")
p3 <- CustomPlotter(df_tidy,"ISD/Det/Homo", "IRF7")
p4 <- CustomPlotter(df_tidy,"ISD/Det/Homo", "TREX1")

p5 <- CustomPlotter(df_tidy,"Virus/Det/Homo", "DNA")
p6 <- CustomPlotter(df_tidy,"Virus/Det/Homo", "IFNb")
p7 <- CustomPlotter(df_tidy,"Virus/Det/Homo", "IRF7")
p8 <- CustomPlotter(df_tidy,"Virus/Det/Homo", "TREX1")

p9 <- CustomPlotter(df_tidy,"ISD/Stoch/Hetero", "DNA")
p10 <- CustomPlotter(df_tidy,"ISD/Stoch/Hetero", "IFNb")
p11 <- CustomPlotter(df_tidy,"ISD/Stoch/Hetero", "IRF7")
p12 <- CustomPlotter(df_tidy,"ISD/Stoch/Hetero", "TREX1")

p13 <- CustomPlotter(df_tidy,"Virus/Stoch/Hetero", "DNA")
p14 <- CustomPlotter(df_tidy,"Virus/Stoch/Hetero", "IFNb")
p15 <- CustomPlotter(df_tidy,"Virus/Stoch/Hetero", "IRF7")
p16 <- CustomPlotter(df_tidy,"Virus/Stoch/Hetero", "TREX1")

figure <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,
                    labels = "AUTO",
                    common.legend = TRUE, legend = "right",
                    align = "hv",
                    ncol = 4, nrow = 4)

ggsave("Figures/testF3.pdf",width = 12,height=9,units="in")
"""
