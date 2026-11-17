
This repository contains all data, code and results for the article 
"Strong genotype-by-genotype and genotype-by-environment interactions exhibited by antimicrobial resistance mutations in the absence of drugs",
by Williams *et al.*

## Data

All data used in this study is contained in the folder "data". 
Four subfolders in this folder (with names of the form "species_author") contain the optical density (OD) data obtained across growth assays.
Each of these subfolders contains both the raw OD data obtained from BioTek platereaders (in subfolders "OD_xxx"),
and metadata linking each well in each growth assay to bacterial species, genotypes and temperatures (in subfolders "spec_xxx").
In adddition, subfolder "genetic_data" contains a table "muts.csv" that has genetic information for each genotype used in this study.
Finally, the file "global_key.csv" is a small table linking project IDs to species and corersponding *rpoB* reference sequences.

## Code

All code is written in R and contained in the folder "R". There are three files in this folder:

* master.R: This is the main file, and the only file that needs to be run to reproduce all results. The code is subdivided into several sections:
  * Step 0: setting up the project; includes sourcing the two other R script files, attaching all R packages etc.
  * Step 1: imports and wrangles all OD data.
  * Step 2: estimates growth parameters from the OD data.
  * Steps 3 through 6: fitting various Bayesian models to the growth parameters.
  * Sections for producing all plots and tables from the results.
  * A final section outputting the session info.
* plotting.R: contains functions for producing the figures. This script is run by master.R.
* util.R: contains various helper functions. This script is also run by master.R.

## Results

All results are found across three folders, "output", "plots" and "tables". 
If deleted, these three folders and all of these files therein are produced by the script main.R, using the data in the "data" folder.

The folder "output" contains the following files:

* all_OD_data.RData: list of tables containing the wrangled optical density data for all growth assays.
* growth_parameters_easylinear_h20.csv: fitted growth parameters (maximum growth rate and maximum OD) for all growth assays, using the easylinear method (main method reported in article).
* brms_fit_xxx.RData: five files containing all Bayesian model fits.
* model_comparison.csv: table containing diagnostics for model comparisons
* sessionInfo.txt: containing all details of the latest run of the master.R script. 

The folder "plots" contains the four main figures of the article as well as Fig. S1 and an alternative to Fig. 4 showing a model fit without 3-way interactions.

The folder "tables" contains the following three files, all including one or more tables in LaTex format:

* mutations.tex: Table 1 from the main text, showing all mutants used in this study.
* comparisons_table.tex: Table S1, showing diagnostics for model comparisons.
* model_fit_tables.tex: Tables S2 through S15 in the ms, showing parameters estimates for all model fits. 


