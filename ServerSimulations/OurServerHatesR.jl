using RCall, CSV, DataFrames

#Figure 6

chronicIFN = CSV.read("StochFigure6Data.csv")

@rput chronicIFN
R"""
library(ggplot2)
library(ggpubr)

p1 <- ggplot(chronicIFN, aes(kcat8, tau7, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  labs(fill="IFN (nM)") +
  ggtitle("Chronic Inflammation") +
  theme_pubr(border=TRUE) +
  scale_fill_distiller(palette = "Spectral",guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),limits=c(0.0,max(chronicIFN$IFN))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw(base_size = 14) +
  xlab(bquote("Interferon Degradation" ~ tau[7])) +
  ylab(bquote("JAK/STAT Activity" ~ k[cat8])) +
  labs(fill="IFN (nM)") +
  theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1)


ggsave("../Figures/Figure6Stoch.pdf")
"""
