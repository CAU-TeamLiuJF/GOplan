# INTRODUCTION
The design of breeding programs holds significant importance in achieving substantial economic gains. Simulation is the optimal approach because testing the effectiveness of a breeding program through trials is time-consuming and costly. GOplan is an R package tailored to help manage animal breeding programs of both single population and crossbreeding systems, as well as optimize breeding programs. GOplan has three wrapped functions, runCore() for evaluating breeding programs’ genetic progress and population diversity change, runWhole() for predicting crossbreed’s phenotype and evaluating the economic profit of crossbreeding systems, and runOpt() for optimizing breeding programs to get the greatest economic profit. Three demos were displayed to demonstrate its function. GOplan is a user-friendly and comprehensive R package with wrapped mainstream crossbreeding frameworks, which simplify the process of constructing crossbreeding systems. It also supports breeding program optimization by using the Bayesian optimization method.

# INSTALLATION
Compressed R package are stored in the branch "Rpackage".

For Windows user： install.packages("GOplan_0.1.0_windows.tar.gz",repos = NULL)

For Linux user： install.packages("GOplan_0.1.0_Linux.tar.gz",repos = NULL)

# demo
All the example files are stored in the folder "example_prm". Note: before you run the example, you should change the "out_path" in the example parameter files. The "out_path" must end with "/".

## Example 1: estimating the results of the nucleus breeding program with different productive lifetime of dams (the variable's symbol is "Yd").

setwd("./example_prm")

runCore(prm_path = "prm_Core.txt")

Note: To save time, we set QUICK as TRUE, which means we do not actually run breeding value estimation. In this case, we just concerned about one trait, you can consider multiple traits by modifying the parameter file.

## Example 2: estimating the results of three-way crossbreeding system breeding programs with different "SorP" (this variable means the male: female ratio of terminal cross).

runWhole(prm_path = "prm_Whole.txt")

Note: To save time, we set QUICK as TRUE. You could also compare the two-way or four-way crossbreeding systems.

## Example 3: running the breeding program optimization of a three-way crossbreeding system under a fixed number of final productions.

runOpt(prm_path = "prm_Opt.txt")

Note: To save time, we set QUICK as TRUE. This function uses the Bayesian Optimization method to find the optimized program. You can define the number of iterations and selected points in each iteration through the parameter file in the "Population Structure" section. The goal is to find the breeding program that yields more profit by optimizing the productive lifetime of each subpopulation and the female size of the nucleus.


If you successfully finish the running, you will see these:
Analyse finish!
Results generated in: xxxx 
Time difference of XXX mins




