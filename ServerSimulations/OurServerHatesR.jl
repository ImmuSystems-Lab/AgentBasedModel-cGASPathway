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

 commonFigureOptions <- list(
   theme_pubr(border=TRUE),
   xlab("Time (hours)"),
   labs(fill=expression(paste("IFN", beta, " (nM)"))),
   theme_pubr(border=TRUE,base_size = 14),
   scale_x_continuous(expand=c(0,0)),
   scale_y_continuous(expand=c(0,0)),
   xlab(bquote("Interferon Degradation" ~ tau[7])),
   ylab(bquote("JAK/STAT Activity" ~ k[cat8])),
   theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1) )


p1 <- ggplot(Det_Homo, aes(kcat8, tau7, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  ggtitle("ISD \n Deterministic \n Homogeneous") +
  scale_fill_distiller(palette = "Spectral",
   guide = guide_colorbar(frame.colour = "black",
   ticks.colour = "black")) +
  commonFigureOptions

p2 <- ggplot(Det_Hetero, aes(kcat8, tau7, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  ggtitle("ISD \n Deterministic \n Heterogeneous") +
  scale_fill_distiller(palette = "Spectral",
   guide = guide_colorbar(frame.colour = "black",
   ticks.colour = "black")) +
  commonFigureOptions

Stoch_Homo$IFN[Stoch_Homo$IFN<0] <- 0

p3 <- ggplot(Stoch_Homo, aes(kcat8, tau7, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  ggtitle("ISD \n Stochastic \n Homogeneous") +
  scale_fill_distiller(palette = "Spectral",
   guide = guide_colorbar(frame.colour = "black",
   ticks.colour = "black"),
   limits = c(0,max(Stoch_Homo$IFN))) +
  commonFigureOptions

figure <- ggarrange(p1, p2, p3,
                    labels = "AUTO",
                    legend = "right",
                    align = "hv",
                    ncol = 3, nrow = 1)

ggsave("../Figures/Figure6_All.pdf",width = 15,height=5)
"""


########################################################################
#Stochastic IFN (Figure 7)
########################################################################

R"""
library(tidyverse)
library(ggpubr)
library(scales)

df <- read.csv("VirusFigure7DataHomo50cells.csv")

#And that's when I discovered how to use tidyverse

df_tidy_median <- df %>%
  pivot_longer(c(Healthy,Infected,Dead),names_to = "CellState",values_to = "CellPercent") %>%
  filter(Percent==0.2) %>%
  group_by(Time,CellState) %>%
  summarise(medianPercentCells = median(CellPercent),
            q25 = quantile(CellPercent,0.25),
            q75 = quantile(CellPercent,0.75))

df_tidy_median$CellState = factor(df_tidy_median$CellState,ordered=TRUE,levels = c("Healthy","Infected","Dead"))

color_list <- c("#377EB8","#4DAF4A","#E41A1C")

p1 <- ggplot(df_tidy_median, aes(x=Time, y=medianPercentCells)) +
  geom_line(aes(x=Time, y=medianPercentCells, color=CellState)) +
  geom_ribbon(aes(ymin=q25, ymax=q75, fill=CellState),alpha=0.3) +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  scale_x_continuous(breaks=seq(0, 48, 12)) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Time (hours)") +
  ylab("Cell Percent") +
  theme_pubr(border=TRUE)


df_tidy_box <- df %>%
  filter(Time==48)

p2 <-ggplot(df_tidy_box, aes(x=as.factor(Percent), y=Dead)) +
  stat_boxplot(geom = "errorbar", width = 0.2,size=0.1) +
  geom_boxplot(fatten=1, fill = "#f6b9ba") +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_y_continuous(labels = scales::percent,limits=c(0,1)) +
  scale_x_discrete(labels = percent((0:10)/10)) +
  xlab(expression("IFN"*beta*" Producing Cells")) +
  ylab("Dead Cells") +
  theme_pubr(border=TRUE)


df <- read.csv("VirusFigure7DataHetero50cells.csv")

df_tidy_median <- df %>%
  pivot_longer(c(Healthy,Infected,Dead),names_to = "CellState",values_to = "CellPercent") %>%
  filter(Percent==0.2) %>%
  group_by(Time,CellState) %>%
  summarise(medianPercentCells = median(CellPercent),
            q25 = quantile(CellPercent,0.25),
            q75 = quantile(CellPercent,0.75))

df_tidy_median$CellState = factor(df_tidy_median$CellState,ordered=TRUE,levels = c("Healthy","Infected","Dead"))

color_list <- c("#377EB8","#4DAF4A","#E41A1C")

p3 <- ggplot(df_tidy_median, aes(x=Time, y=medianPercentCells)) +
  geom_line(aes(x=Time, y=medianPercentCells, color=CellState)) +
  geom_ribbon(aes(ymin=q25, ymax=q75, fill=CellState),alpha=0.3) +
  scale_fill_manual(values=color_list) +
  scale_color_manual(values=color_list) +
  scale_x_continuous(breaks=seq(0, 48, 12)) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Time (hours)") +
  ylab("Cell Percent") +
  theme_pubr(border=TRUE)


df_tidy_box <- df %>%
  filter(Time==48)

p4 <-ggplot(df_tidy_box, aes(x=as.factor(Percent), y=Dead)) +
  stat_boxplot(geom = "errorbar", width = 0.2,size=0.1) +
  geom_boxplot(fatten=1, fill = "#f6b9ba") +
  scale_y_continuous(labels = scales::percent,limits=c(0,1)) +
  scale_x_discrete(labels = percent((0:10)/10)) +
  xlab(expression("IFN"*beta*" Producing Cells")) +
  ylab("Dead Cells") +
  theme_pubr(border=TRUE)




figure <- ggarrange(p1, p2, p3, p4,
                    labels = "AUTO",
                    common.legend = TRUE, legend = "right",
                    widths = c(3,6),
                    align = "hv",
                    ncol = 2, nrow = 2)

ggsave("../Figures/Figure7.pdf",width = 9,height=6,units="in")
"""
