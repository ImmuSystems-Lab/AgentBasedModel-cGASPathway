using RCall
using Distributed
addprocs(2)
@everywhere include("ProblemGenerator.jl")


#Need a function that creates distributions of heterogenous states
@everywhere AddNoise2States(σ) =
  [TruncatedNormal(μ, σ * μ, 0, Inf) for μ in nonZeroSpeciesValues]

#Create a function that calculate the peak IFN distribution given:
 #The number of cells initially infected
 #The amount of noise desired for the initial protien concentrations
@everywhere function HeteroSimulation(primeCellPercent, varIC, prob)
    #Make a copy of the initial conditions
  u0 = copy(prob.u0)
    #Get number of primary cells wanted from percent
  primeCellNum = round(Int, nCells * primeCellPercent)
  secondaryCellNum = nCells - primeCellNum
    #Randomly place them in the ABM
  cellLocations = shuffle([ones(primeCellNum); zeros(secondaryCellNum)])
    #Update the number of initially infected cells (DNA conc.)
  u0[:, :, 2] = m2c.(1e3 .* reshape(cellLocations, N, N))

    @show varIC
    #Add noise to the nonzero initial conditions
  noiseDistributions = AddNoise2States(varIC)
  for (i, index) in enumerate(nonZeroSpeciesIdx)
    u0[:, :, index] = rand(noiseDistributions[i], N, N)
  end

    #Make a copy of the parameter container
  θ = deepcopy(prob.p)
    #Reset the number of cells infected initially
  θ.cellsInfected = fill(Inf, N, N)
  θ.cellsInfected[findall(u0[:, :, 2] .> 0.0), 1] .= 0.0
    #Reset the mass balances
  θ.mass = [u0[:, :, i] for i in nonZeroSpeciesIdx]

    #Define the new problem and solve
  probHetero = remake(prob; u0 = u0, p = θ)
  solHetero = solve(probHetero, CVODE_BDF(linear_solver = :GMRES))

    #Seperate by primary and secondary cells
  isPrimary = u0[:, :, 2] .> 0.0
  maxIFN = [maximum(solHetero[coord, 7, :]) for coord in cellIndicies]

  return [maxIFN[isPrimary], maxIFN[.!isPrimary]]
end


#Run the simulation as either Stochastic/Deterministic and Hetero
@everywhere probDet = ModelSetup(:ISD, :notStochastic, :Hetero)
@everywhere probStoch = ModelSetup(:ISD, :Stochastic, :Hetero)


varICRange = 0.0:0.1:1.0
varIClength = length(varICRange)
primeCellPercent = 0.63
primeCellNum = round(Int,nCells * primeCellPercent)

maxIFNDetHet = pmap(x -> HeteroSimulation(primeCellPercent, x, probDet), varICRange)
maxIFNStochHet = pmap(x -> HeteroSimulation(primeCellPercent, x, probStoch), varICRange)

function DataFrameMaker(sols)
  #Seperate out the primary and secondary cells
  primary = [sols[i][1] for i=1:varIClength]
  secondary = [sols[i][2] for i=1:varIClength]

  #Convert into a dataframe with 3 columns (percent,IFN, and celltype (primary or secondary))
  primary = DataFrame(percent=repeat(1:varIClength,inner=primeCellNum),
                      IFN = vcat(primary...),
                      CellState = repeat(["Primary"],length(vcat(primary...))) )
  secondary = DataFrame(percent=repeat(1:varIClength,inner= nCells-primeCellNum),
                        IFN = vcat(secondary...),
                        CellState = repeat(["Secondary"],length(vcat(secondary...))) )
  #Combine the two datadrames and return th result
  return vcat(primary,secondary)
end


ISDDetHet = DataFrameMaker(maxIFNDetHet)
ISDStochHet = DataFrameMaker(maxIFNStochHet)


@rput ISDDetHet
@rput ISDStochHet
R"""
library(ggplot2)
library(ggpubr)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

commonFigureOptions <- list(geom_split_violin(scale="width",trim=FALSE),
  theme_pubr(border=TRUE),
  xlab("Initial Condition Variance"),
  ylab("IFN (nM)"),
  scale_x_discrete(labels= unlist(lapply(seq(0.0,1.0,0.1),toString))),
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)),
  theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 0.5)
  )


p1 <- ggplot(ISDDetHet, aes(factor(percent), IFN, fill = CellState)) +
  ggtitle("ISD: Deterministic (63% of Cell Infected)") +
  commonFigureOptions

p2 <- ggplot(ISDStochHet, aes(factor(percent), IFN, fill = CellState)) +
  ggtitle("ISD: Stochastic (63% of Cell Infected)") +
  commonFigureOptions

figure <- ggarrange(p1,p2,
                    labels = "AUTO",
                    common.legend = TRUE, legend = "right",
                    align = "h",
                    ncol = 1, nrow = 2)

ggsave("./Figures/Figure5.pdf")
"""
