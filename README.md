# GOplan
R package for managing breeding programs。

# INSTALLATION
install.packages("GOplan_0.1.0.tar.gz",repos = NULL)

# demo
All the example files are stored in folder "example_prm".

Example 1:

runCore(prm_path = "prm_Core.txt")

This example shows how to running nucleus population breeding program with a varaible "Yd"(the productive lifetime of dams in population). To save time, we set QUICK as TRUE.

Example 2:

runWhole(prm_path = "prm_Whole.txt")

This example shows how to running the breeding program evaluation of a three-way crossbreeding system with a varaible "SorP"(the male:female ratio of terminal cross). To save time, we set QUICK as TRUE.

Example 3:

runOpt(prm_path = "prm_Opt.txt")

This example shows how to running the breeding program optimization of a three-way crossbreeding system under a fixed number of final production. To save time, we set QUICK as TRUE.

If you successfully finish the running, you will see these:
Analyse finish!
Results generated in: xxxx 
Time difference of 30.48747 mins



