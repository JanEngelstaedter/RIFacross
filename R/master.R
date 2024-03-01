library(tidyverse)
library(grow96)
library(Biostrings)
library(ALJEbinf)

source("./R/plotting.R")

# wrangling and qc:
data_Maddie <- processODData(specPath = "./data/A_baylyi_Maddie/specs/",
                             dataPath = "./data/A_baylyi_Maddie/OD_raw") |>
  blankODs(groups = "RIFconc") |>
  mutate(RIFconc = as.double(RIFconc))
qcODData(data_Maddie, path = "./output/QC/A_baylyi_Maddie/", silent = FALSE)

# plotting:
plot_growthData(data_Maddie, "./plots/A_baylyi_Maddie_growthcurves.pdf")
