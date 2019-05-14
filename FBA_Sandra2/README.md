This repository holds the code for my final project for CHEME 7770.
To run the code, you will need to have Julia 1.0. (It may work on an older version, but I have not tried it myself.
Changes to TXTL parameters should be made in TXTLDictionary.jl unless the changes are for a Monte Carlo simulation, in which case the changes should be made by using include("TXTLDictionary.jl") in the beginning of Ensemble.jl.
Flux bounds, species bounds and the objective function(s) can be set in Bounds.jl.
To execute a single run, change the directory to wherever the clone of this repository is kept and type include("Solve.jl") in your Julia REPL or however you usually run Julia.
Solve.jl returns an array of the fluxes of interest.
To execute the Monte Carlo, make sure you are in the correct directory and then type include("Ensemble.jl") in the REPL.
Ensemble.jl saves a figure in .PNG format to disk and returns the maximum values of the fluxes of interest out of the entire ensemble.
