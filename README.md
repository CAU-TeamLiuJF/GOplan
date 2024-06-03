# GOplan
R package for managing breeding programs。

# INSTALLATION
Compressed R package are stored in branch "Rpackage".

For Windows user： install.packages("GOplan_0.1.0_windows.tar.gz",repos = NULL)

For Linux user： install.packages("GOplan_0.1.0_Linux.tar.gz",repos = NULL)

# demo
All the example files are stored in folder "example_prm". Note: before you run the example, you should change the "out_path" in the example parameter files. The "out_path" must end with "/".

Example 1:

setwd("./example_pam")

runCore(prm_path = "prm_Core.txt")

This example shows how to running nucleus population breeding program with a varaible "Yd"(the productive lifetime of dams in population). To save time, we set QUICK as TRUE.

Example 2:

runWhole(prm_path = "prm_Whole.txt")

This example shows how to running the breeding program evaluation of a three-way crossbreeding system with a varaible "SorP"(the male:female ratio of terminal cross). To save time, we set QUICK as TRUE.

Example 3:

res = runOpt(prm_path = "prm_Opt.txt")

getOptRes(res = res, out_path = out_path) # transfer the results to xlsx file

This example shows how to running the breeding program optimization of a three-way crossbreeding system under a fixed number of final production. To save time, we set QUICK as TRUE.


If you successfully finish the running, you will see these:
Analyse finish!
Results generated in: xxxx 
Time difference of XXX mins





