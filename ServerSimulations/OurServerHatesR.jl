using RCall, CSV, DataFrames

########################################################################
#Chronic Infected (Figure 6)
########################################################################

Det_Homo = CSV.read("noStochFigure6Data.csv")
Det_Hetero = CSV.read("VirusFigure6_Det_Hetero.csv")
Stoch_Homo = CSV.read("StochFigure6Data.csv")

@rput Det_Homo
@rput Det_Hetero
@rput Stoch_Homo
R"""
library(ggplot2)
library(ggpubr)

#Find min/max to put everything on common scale
IFNVals <- c(Det_Homo$IFN,Det_Hetero$IFN,Stoch_Homo$IFN)
lowIFN <- min(IFNVals)
highIFN <- max(IFNVals)

 commonFigureOptions <- list(
   theme_pubr(border=TRUE),
   xlab("Time (hours)"),
   labs(fill=expression(paste("IFN", beta, " (nM)"))),
   theme_pubr(border=TRUE,base_size = 14),
   scale_fill_distiller(palette = "Spectral",
    guide = guide_colorbar(frame.colour = "black",
    ticks.colour = "black"),
    limits = c(lowIFN,highIFN)),
   scale_x_continuous(expand=c(0,0)),
   scale_y_continuous(expand=c(0,0)),
   xlab(bquote("Interferon Degradation" ~ tau[7])),
   ylab(bquote("JAK/STAT Activity" ~ k[cat8])),
   theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1) )


p1 <- ggplot(Det_Homo, aes(kcat8, tau7, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  ggtitle("ISD: Deterministic + Homogeneous") +
  commonFigureOptions

p2 <- ggplot(Det_Hetero, aes(kcat8, tau7, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  ggtitle("ISD: Deterministic + Heterogeneous") +
  commonFigureOptions

p3 <- ggplot(Stoch_Homo, aes(kcat8, tau7, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  ggtitle("ISD: Stochastic + Homogeneous") +
  commonFigureOptions

figure <- ggarrange(p1, p2, p3,
                    labels = "AUTO",
                    common.legend = TRUE,legend = "right",
                    align = "hv",
                    ncol = 3, nrow = 1)

ggsave("../Figures/Figure6_All.pdf",width = 14,height=5)
"""


########################################################################
#Stochastic IFN (Figure 7)
########################################################################

#We need to pass 2 subsetted dataframes into R for plotting
timeDynamics = filter(row->row[:Percent]==0.2,VirusSim)

@rput VirusSim
R"""

timeDynamics

"""
