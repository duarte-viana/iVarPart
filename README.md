# iVarPart

"Partitioning environment and space in species-by-site matrices: a comparison of methods for community ecology and macroecology"

Authors: Duarte S. Viana, Petr Keil, Alienor Jeliazkov


The files contained in this repository provide the R code for simulating and analysing the data corresponding to the mentioned paper. Below we provide a brief description of each file.

“Methods_functions.R”: this script contains the code of the functions used to analyse the simulated data, i.e., the functions for the constrained ordination methods, distance-based methods, generalised linear models and tree-based methods.

“R2D2_functions.R”: this script contains the functions used to calculate the different R-squared metrics, as well as the variation partitioning function.

“Data_simulation.R”: this script contains the R code for the data simulations. Each simulation corresponds to a given combination of parameters provided in the table contained in the file “pars_GB.Rdata” (each row of the table corresponds to one combination and thus one simulation; note that each combination has 10 replicates).

"Data_analysis.R": this script contains the R code for analysing the simulated data (with code in “Data_simulation.R”) with the statistical methods coded in the files “Methods_functions.R” and “R2D2_functions.R”.
