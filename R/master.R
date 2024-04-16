library(tidyverse)
library(grow96)
library(Biostrings)
library(ALJEbinf)
library(googlesheets4)


source("./R/plotting.R")

projects <- read_csv("./data/global_key.csv")
dat <- vector(mode = "list", length = nrow(projects))
names(dat) <- projects$Project_ID
analysis <- dat

for(i in 1:nrow(projects)) {
  cat("--------------------------------------------------------\n")
  cat(paste0(">>> Working on project ", projects$Project_ID[i], "\n"))
  cat("--------------------------------------------------------\n\n")
  # wrangling:
  dat[[i]] <- processODData(specPath = paste0("./data/", projects$Project_ID[i], "/specs/"),
                            dataPath = paste0("./data/", projects$Project_ID[i], "/OD_raw"))
  if(projects$Project_ID[i] == "A_baylyi_Ian") {
    dat[[i]] <- dplyr::rename(dat[[i]], RIFconc = Concentration)  
  }
  dat[[i]] <- dat[[i]] |>
    blankODs(groups = "RIFconc") |>
    mutate(RIFconc = as.double(RIFconc)) 
  
  # quality control:
  qcODData(dat[[i]], path = paste0("./output/QC/", projects$Project_ID[i], "/"))
  
  # plotting:
  plot_growthData(dat[[i]], paste0("./plots/", projects$Project_ID[i], "_growthcurves.pdf"))
  analysis[[i]] <- analyseODData(dat[[i]])
  cat("\n")
}

#save(analysis, file = "./output/analysis.RData")
load("./output/analysis.RData")

parms <- NULL
# turn analysis list into single data frame:
for(i in 1:nrow(projects)) {
  analysis[[i]]$pars$Project_ID <- projects$Project_ID[i]
  analysis[[i]]$pars$Species <- projects$Species[i]
  parms <- rbind(parms, analysis[[i]]$pars)
}

parms <- parms |>
  select(Project_ID, Species, everything())

# read in file with identities of mutations:
muts <- read_csv("./data/Genetic_data/muts.csv")

# attach rpoB mutations to growth parameters:
parms <- parms |>
  left_join(muts, by = join_by(Species, Strain == Genotype)) |>
  select(Project_ID, Species, Strain,
         Nt_mut_name, Nt_mut_name_Ecoli, AA_mut_name, AA_mut_name_Ecoli, 
         RIFconc, Replicate,
         mumax, r2, maxOD) |>
  arrange(Species, Strain, RIFconc, Replicate)

parms$AA_mut_name_Ecoli <- ifelse(parms$Species == "Lactiplantibacillus plantarum" & parms$Strain != "Wild Type",
                                 parms$Strain,
                                 parms$AA_mut_name_Ecoli)

parms$AA_mut_name_Ecoli <- case_when(parms$Species == "Lactiplantibacillus plantarum" & parms$Strain != "Wild Type" 
                                    ~paste0("rpoB_", parms$AA_mut_name_Ecoli), TRUE~parms$AA_mut_name_Ecoli)

parms$AA_mut_name_Ecoli <- ifelse(parms$Strain %in% c("MG1655", "WT", "AB3", "Wild Type"), "WT", parms$AA_mut_name_Ecoli)

#max growth rate in absence of antibiotics
mean <- parms  %>% filter(RIFconc==0) %>% 
  group_by(Species, AA_mut_name_Ecoli) %>% 
  summarise(mean_mumax = mean(mumax, na.rm = TRUE), 
            mean_SD = sd(mumax, na.rm = TRUE))  %>% 
  ungroup()

