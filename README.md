# iVarPart

"Partitioning environment and space in site-by-species matrices: a comparison of methods for community ecology and macroecology"
Duarte S. Viana, Petr Keil, Alienor Jeliazkov


The files contained in this repository provide the R code for simulating and analysing the data corresponding to the mentioned paper. Below we provide a brief description of each file.

“Methods_functions.R”: this script contains the code of the functions used to analyse the simulated data in "Exercise 2", i.e., the functions for the constrained ordination methods, distance-based methods, generalised linear models and tree-based methods.

“R2D2_functions.R”: this script contains the functions used to calculate the different R-squared metrics, as well as the variation partitioning function.

“Data_simulation.R”: this script contains the R code for the data simulations. Each simulation corresponds to a given combination of parameters provided in the table contained in the file “pars1.Rdata” (for Exercise 1) and “pars_VP.Rdata” (for Exercise 2) (each row of the table corresponds to one combination and thus one simulation).

"Exercise1.R": this script contains the R code for analysing the simulated data "sim1.Rdata". 

"Exercise2.R": this script contains the R code for analysing the simulated data "sim_data_vp.Rdata" with the statistical methods coded in the files “Methods_functions.R” and “R2D2_functions.R”. 
