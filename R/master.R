library(tidyverse)
library(grow96)
library(Biostrings)
library(ALJEbinf)
library(googlesheets4)
library(UpSetR)

source("./R/plotting.R")
source("./R/util.R")

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
  
  # manual blanking of low_2 plates in Edward's project
  # (in these plates, no RIF was added to BLANK wells on the border so they can't be used)
  if (projects$Project_ID[i] == "L_plantarum_Edward") {
    meanBlankOD_low2 <- dat[[i]] |>
      filter(WellType == "BLANK") |>
      filter(Plate == "Mutant_L.plantarum_RIF_MIC_low_2") |>
      pull(OD) |>
      mean()
    
    dat[[i]] <- dat[[i]] |>
      mutate(blankedOD = if_else(Plate == "Mutant_L.plantarum_RIF_MIC_low_2",
                                 OD - meanBlankOD_low2,
                                 blankedOD))
  }
  
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

# read in files with identities of mutations:
muts <- read_csv("./data/Genetic_data/muts.csv")
muts_Edward <- read_csv("./data/Genetic_data/Edward_mutants.csv") |>
  filter(!is.na(Strain_Name)) |>
  mutate(New_name = str_replace(AA_mut_name_Ecoli, "rpoB_", ""))
muts_Edward_lookup <- muts_Edward$Strain_Name
names(muts_Edward_lookup) <- muts_Edward$New_name

# attach rpoB mutations to growth parameters:
parms <- parms |>
  # standardise names for suseptible wildtype:
  mutate(Strain = ifelse(Strain %in% c("AB3", "MG1655", "Wild Type"), 
                         "WT", Strain)) |>
  # fix Edward's strain names:
  mutate(Strain = ifelse((Project_ID == "L_plantarum_Edward") & (Strain != "WT"), 
                         muts_Edward_lookup[Strain],
                         Strain)) |>
  left_join(muts, by = join_by(Species, Strain == Genotype)) |>
  select(Project_ID, Species, Strain,
         Nt_mut_name, Nt_mut_name_Ecoli, AA_mut_name, AA_mut_name_Ecoli, 
         RIFconc, Replicate,
         mumax, r2, maxOD) |>
  arrange(Species, Strain, RIFconc, Replicate)

# Alicia, your code below shouldn't be necessary anymore, so I've commented it out:
# parms$AA_mut_name_Ecoli <- ifelse(parms$Species == "Lactiplantibacillus plantarum" & parms$Strain != "Wild Type",
#                                  parms$Strain,
#                                  parms$AA_mut_name_Ecoli)
# 
# parms$AA_mut_name_Ecoli <- case_when(parms$Species == "Lactiplantibacillus plantarum" & parms$Strain != "Wild Type" 
#                                     ~paste0("rpoB_", parms$AA_mut_name_Ecoli), TRUE~parms$AA_mut_name_Ecoli)
# 
# parms$AA_mut_name_Ecoli <- ifelse(parms$Strain %in% c("MG1655", "WT", "AB3", "Wild Type"), "WT", parms$AA_mut_name_Ecoli)

#max growth rate in absence of antibiotics
mean <- parms  %>% filter(RIFconc==0) %>% 
  group_by(Species,Strain, AA_mut_name_Ecoli) %>% 
  summarise(mean_mumax = mean(mumax, na.rm = TRUE), 
            mean_SD = sd(mumax, na.rm = TRUE), 
            observations = n())  %>% 
    ungroup()


MICs <- parms |>
  group_by(across(Species:AA_mut_name_Ecoli)) |>
  summarise(MIC = getMIC(RIFconc, maxOD))

write_csv(parms, "./output/growth_parameters.csv")
write_csv(MICs, "./output/MICs.csv")


MICs_overlap <- MICs |> group_by(AA_mut_name_Ecoli) |>
  filter(n() > 1) |>
  arrange(AA_mut_name_Ecoli)

growth_parameters |>
  select(Species, AA_mut_name_Ecoli) |>
  distinct() |>
  group_by(AA_mut_name_Ecoli) |>
  arrange(AA_mut_name_Ecoli) |>
  View()
