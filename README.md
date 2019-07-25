# code_Fauvernier_JRSSC_2019

Illustration code for the article 

## Multi-dimensional penalized hazard model with continuous covariates: applications for studying trends and social inequalities in cancer survival
Mathieu Fauvernier, Laurent Roche, Zoé Uhry, Laure Tron, Nadine Bossard, Laurent Remontet and the CENSUR Working Survival Group.

Corresponding author:
Mathieu FAUVERNIER
mathieu.fauvernier@chu-lyon.fr

The provided code is composed of two main code files :
- The file entitled "0000_Simulation_code.R" details all the necessary steps to fit the models used in the simulation study. 
In order to use it you need to change the « Dir » variable at the beginning, choose a scenario and a sample size and download the corresponding files starting by « ListDataSimulation » and « DataDesign » from : https://zenodo.org/record/3268769#.XR4ird9l_ct. We also added the code for the rstpm2 model (as a commentary) if you want to reproduce the results of the Supplementary Material C.2.

- The file "0000_Compar_packages.R"  contains code used in the Supplementary Material C.1 (only for sample size = 2,000) to compare survPen with existing approaches. For sample size = 20,000 you must use the following dataset :
« datCancer20000 <- do.call("rbind", ListDataSimulation.MI2 [1:10]) »
where ListDataSimulation.MI2 is the RData of simulated datasets for Cervix scenario and sample size 2000.


