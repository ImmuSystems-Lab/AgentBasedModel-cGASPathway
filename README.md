#Agent-Based Modeling Reveals Benefits of Heterogeneous and Stochastic Cell Populations during cGAS-Mediated IFNβ Production

This repository implements a PDE/ABM model of interferon signaling across a population of 40,000 cells (200 by 200 grid). Each cell contains a set of ordinary differential equations that track intracellular concentrations overtime in response to interferon stimulatory DNA (ISD) or to viral infection. These stimuli cause cells to upregulate IFNβ which is excreted and allowed to diffuses across the cell population. As neighboring cells detect diffused IFNβ, they begin to express interferon stimulated genes (ISD) that feedback to regulate their own interferon response. The simulation terminates when all the cell die or all of the ISD/virus is eliminated. 

Here is a simple example of running a simulation:

```julia
include("ProblemGenerator.jl")
include("VirusCallBacks.jl")

prob = ModelSetup(:Virus,:notStochastic,:Homo)
    sol = @time solve(probISD_nS_Ho,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)
```
The `ModelSetup` function takes in 3 arguements: Virus/ISD, notStochastic/Stochastic, and Homo/Hetero. The simulation can be run with a viral infection or ISD transfection. The main difference being that a virus can replicate and ISD cannot. Stochasticity affects the number of cells capable of producing an interferon response (the default allowing ~20% of cells to respond with IFNβ) and heterogeneity adds variability to the initial concentrations of constitutively active species. 
