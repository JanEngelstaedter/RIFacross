library(tidyverse)
library(grow96)
library(Biostrings)
library(ALJEbinf)

source("./R/plotting.R")

projects <- c("E_coli_Alicia", "A_baylyi_Maddie", "A_baylyi_Ian", "L_plantarum_Edward")
dat <- vector(mode = "list", length = length(projects))
names(dat) <- projects

for(i in 1:length(projects)) {
  cat("--------------------------------------------------------\n")
  cat(paste0(">>> Working on project ", projects[i], "\n"))
  cat("--------------------------------------------------------\n\n")
  # wrangling:
  dat[[i]] <- processODData(specPath = paste0("./data/", projects[i], "/specs/"),
                            dataPath = paste0("./data/", projects[i], "/OD_raw"))
  if(projects[i] == "A_baylyi_Ian") {
    dat[[i]] <- dplyr::rename(dat[[i]], RIFconc = Concentration)  
  }
  dat[[i]] <- dat[[i]] |>
    blankODs(groups = "RIFconc") |>
    mutate(RIFconc = as.double(RIFconc)) 
  
  # quality control:
  qcODData(dat[[i]], path = paste0("./output/QC/", projects[i], "/"))
  
  # plotting:
  plot_growthData(dat[[i]], paste0("./plots/", projects[i], "_growthcurves.pdf"))
  cat("\n")
}
