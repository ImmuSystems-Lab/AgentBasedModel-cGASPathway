#Shows what a typical ISD and Viral Infections looks like
using RCall
include("ProblemGenerator.jl")

########################################################################
#Heatmaps of IFN dynamics (Figure 2)
########################################################################

probISD = ModelSetup(:ISD,:notStochastic)
    solISD = @time solve(probISD,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
ISD10hours = DataFrame()
ISD10hours.x = repeat(1:N,N)
ISD10hours.y = repeat(1:N,inner=N)
ISD10hours.IFN = vec(solISD(10.0)[:,:,7])

probVirus = ModelSetup(:Virus,:Stochastic)
    solVirus = @time solve(probVirus,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
Virus10hours = DataFrame()
Virus10hours.x = repeat(1:N,N)
Virus10hours.y = repeat(1:N,inner=N)
Virus10hours.IFN = vec(solVirus(10.0)[:,:,7])

@rput ISD10hours
@rput Virus10hours
R"""
  library(ggplot2)
  library(ggpubr)

  #To draw an arc on the heatmap to show initial condition
  arc <- data.frame(
  x0 = seq(0, 50, length.out = 1000),
  y0 = sqrt(50^2 - (seq(0, 50, length.out = 1000)^2))
  )

  commonFigureOptions <- list(scale_fill_distiller(palette = "Spectral",guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")),
    scale_x_continuous(expand=c(0,0)),
    scale_y_continuous(expand=c(0,0)),
    theme_bw(base_size = 14),
    ylab("Cell"),
    xlab("Cell"),
    labs(fill="IFN (nM)"),
    theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

  p1 <- ggplot() +
    geom_raster(ISD10hours, mapping =aes(x, y, fill=IFN)) +
    geom_line(data=arc, aes(x=x0, y=y0),size=0.2,linetype="dashed") +
    labs(fill="IFN (nM)") +
    ggtitle("ISD Transfection: \n Time = 10 hours") +
    commonFigureOptions

  p2 <- ggplot(Virus10hours, aes(x, y, fill=IFN)) +
    geom_raster(aes(fill=IFN)) +
    ggtitle("Virus Infection: \n Time = 10 hours") +
    commonFigureOptions

  figure <- ggarrange(p1, p2,
                      labels = c("A", "B"),
                      common.legend = TRUE, legend = "right",
                      align = "h",
                      ncol = 2, nrow = 2)

  ggsave("./Figures/Figure2.pdf")
"""
