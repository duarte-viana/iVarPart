# iVarPart

"Disentangling spatial and environmental effects: flexible methods for community ecology and macroecology"
Duarte S. Viana, Petr Keil, Alienor Jeliazkov


The files contained in this repository provide the R code for simulating and analysing the data corresponding to the mentioned paper. Below we provide a brief description of each file.

“Methods_functions.R”: this script contains the code of the functions used to fit the different models (Table 1 of the main paper).

“R2D2_functions.R”: this script contains the functions used to calculate the different R-squared metrics.

“Data_simulation.R”: this script contains the R code for the data simulations. Each simulation corresponds to a given combination of parameters provided in the table contained in the file “pars.Rdata” (each row of the table corresponds to one combination and thus one simulation).

"Model_fitting.R": this script contains the R code for fitting the different models in Table 1 (of the main paper), coded in “Methods_functions.R”, to the simulated data "sim_data.Rdata". 

"Analysis_results.R": this script contains the R code for analysing the performance of the different models to fit the data and partition explained variation. 
