# Step 0: Loading all libraries and scripts             ----
# __________________________________________________________

library(tidyverse)
library(grow96)
library(Biostrings)
library(ALJEbinf)
library(googlesheets4)
library(UpSetR)
library(varhandle)
library(ggpubr)

source("./R/plotting.R")
source("./R/util.R")


# Step 1: Loading and processing all data               ----
# __________________________________________________________
#input= assay data and spec files
#output= growth parameters and growth curve plots


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
  analysis[[i]] <- analyseODData(dat[[i]], tmin = 1)
  cat("\n")
}
save(analysis, file = "./output/analysis.RData")

#empty working environment to keep everything clean:
rm.all()

# Step 2: Add genetic information to data               ----
# __________________________________________________________
#input= growth parameters and genetic data from Google Sheet
#output= growth parameters with corresponding mutation


projects <- read_csv("./data/global_key.csv")
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

parms$AA_pos <- as.numeric(gsub("[^0-9]","", parms$AA_mut_name_Ecoli))
parms <- parms %>% arrange(AA_pos) %>% 
  mutate(New_name = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) %>% 
  filter(!(Species == "Acinetobacter baylyi" & (Strain %in% c("RM08", "RM15")))) # remove two mutants that are the same as RM10 and RM12
parms$New_name[parms$Nt_mut_name =="rpoB_dup1611"] <- "dup1611"
parms$New_name[parms$Strain =="WT"] <- "WT"

write_csv(parms, "./output/growth_parameters.csv")

#empty working environment to keep everything clean:
rm.all()

# Step 3: Analysing data                                ----
# __________________________________________________________
#input= growth parameters 
#output= mumax, maxOD and MIC parameters and overlap between species


parms <- read_csv("./output/growth_parameters.csv")

#max growth rate and maxOD in absence of antibiotics
means <- parms  %>% filter(RIFconc==0) %>% 
  group_by(Species,Strain, AA_mut_name_Ecoli) %>% 
  summarise(mean_mumax = mean(mumax, na.rm = TRUE), 
            SD_mumax = sd(mumax, na.rm = TRUE),
            mean_maxOD = mean(maxOD, na.rm = TRUE),
            SD_maxOD = sd(maxOD, na.rm = TRUE),
            observations = sum(!is.na(mumax)))  %>% #to remove the NA parameters 
    ungroup()



write_csv(means, "./output/means.csv")


MICs <- parms |>
  group_by(across(Species:AA_mut_name_Ecoli)) |>
  summarise(MIC = getMIC(RIFconc, maxOD))

WT <- MICs %>% filter (Strain=="WT")
WT_MIC <- setNames(WT$MIC,WT$Species)
MICs <- MICs %>% mutate(xfold=MIC/WT_MIC[Species])

MICs <- MICs %>% 
  mutate(New_name = str_replace(AA_mut_name_Ecoli, "rpoB_", ""))

write_csv(MICs, "./output/MICs.csv")

# assessing overlaps between mutations:
overlaps <- muts |>
  arrange(AA_pos_Ecoli) |>
  select(Species, Genotype, AA_mut_name_Ecoli) |>
  dplyr::rename(Strain = Genotype) |>
  filter(!is.na(AA_mut_name_Ecoli)) |>
  filter(!(Species == "Acinetobacter baylyi" & (Strain %in% c("RM08", "RM15")))) |>  # remove two mutants that are the same as RM10 and RM12
  left_join(means, by = join_by(Species, Strain, AA_mut_name_Ecoli)) |>
  left_join(MICs, by = join_by(Species, Strain, AA_mut_name_Ecoli))|>
  mutate(Phenotyped = !is.na(mean_maxOD)) |>
  select(AA_mut_name_Ecoli, Species, Strain, Phenotyped) |>
  pivot_wider(names_from = Species, values_from = c(Strain, Phenotyped))

# more phenotyping needed for the following:
need_phenotyping <- overlaps |>
  mutate(n_mutants = 3 - (is.na(`Strain_Acinetobacter baylyi`) + is.na(`Strain_Escherichia coli`) + is.na(`Strain_Lactiplantibacillus plantarum`))) |>
  mutate(n_phenotyped = (!is.na(`Phenotyped_Acinetobacter baylyi`) & `Phenotyped_Acinetobacter baylyi`) +
                        (!is.na(`Phenotyped_Escherichia coli`) & `Phenotyped_Escherichia coli`) +
                        (!is.na(`Phenotyped_Lactiplantibacillus plantarum`) & `Phenotyped_Lactiplantibacillus plantarum`)) |>
  filter(n_mutants > 1) |>
  filter(n_mutants > n_phenotyped)

# MICs_overlap <- MICs |> group_by(AA_mut_name_Ecoli) |>
#   filter(n() > 1) |>
#   arrange(AA_mut_name_Ecoli)

growth_parameters |>
  select(Species, AA_mut_name_Ecoli) |>
  distinct() |>
  group_by(AA_mut_name_Ecoli) |>
  arrange(AA_mut_name_Ecoli) |>
  View()

#empty working environment to keep everything clean:
rm.all()

# Step 4: Graphing                               ----
# __________________________________________________________
#input= growth parameters 
#output= graph of each species parameters and comparison of parameters between species 



means <- read.csv("./output/means.csv")
MICs <- read.csv("./output/MICs.csv")
parms <- read_csv("./output/growth_parameters.csv")


#Make parameters relative to species WT 
WT <- means %>% filter (Strain=="WT")
WT_mumax <- setNames(WT$mean_mumax,WT$Species)
WT_maxOD <- setNames(WT$mean_maxOD, WT$Species)
means <- means %>% mutate(relative_mumax=mean_mumax/WT_mumax[Species]) %>% 
  mutate(relative_maxOD=mean_maxOD/WT_maxOD[Species])


#graph relative parameters 
means_filt <- means %>%  filter(Strain != "R21") %>%
  filter(Strain!= "RM15") #R21 is the A.baylyi Nt duplication and RM15 is same as RM12 
means_filt = subset(means_filt, select = c(Species, AA_mut_name_Ecoli, relative_mumax, relative_maxOD))

meanwide <- pivot_wider(means_filt, names_from = Species, values_from = c(relative_mumax, relative_maxOD))
names(meanwide) <- make.names(names(meanwide), unique = TRUE)


MICs_filt <- MICs %>%  filter(Strain != "R21") %>% 
  filter(Strain!= "RM15") #R21 is the A.baylyi Nt duplication and RM15 is same as RM12 
MICs_filt = subset(MICs_filt, select = c(Species, AA_mut_name_Ecoli, xfold))

MICwide <- pivot_wider(MICs_filt, names_from = Species, values_from = xfold)
names(MICwide) <- make.names(names(MICwide), unique = TRUE)


#E coli v A baylyi parameter plots

ecab_mu<- ggplot(meanwide, aes(x=relative_mumax_Escherichia.coli, 
                               y=relative_mumax_Acinetobacter.baylyi, 
                               label=AA_mut_name_Ecoli))+
geom_point() +
geom_text(hjust=0, vjust=0)+
labs(title = "E.coli v A.baylyi mumax")


ecab_OD <- ggplot(meanwide, aes(x=relative_maxOD_Escherichia.coli,
                                y=relative_maxOD_Acinetobacter.baylyi, 
                                label=AA_mut_name_Ecoli))+
geom_point()+
geom_text(hjust=0, vjust=0)+
labs(title = "E.coli v A.baylyi maxOD")

ecab_MIC <- ggplot(MICwide, aes(x=Escherichia.coli,
                                y=Acinetobacter.baylyi, 
                                label=AA_mut_name_Ecoli))+
  geom_point()+
  geom_text(hjust=0, vjust=0)+
  labs(title = "E.coli v A.baylyi MIC")+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')



#L.plantarum v A.baylyi
lpab_mu <- ggplot(meanwide, aes(x=relative_mumax_Lactiplantibacillus.plantarum, 
                                y=relative_mumax_Acinetobacter.baylyi, 
                                label=AA_mut_name_Ecoli))+
geom_point() +
geom_text(hjust=0, vjust=0)+
labs(title = "L.plantarum v A.baylyi mumax")

lpab_OD <- ggplot(meanwide, aes(x=relative_maxOD_Lactiplantibacillus.plantarum,
                                y=relative_maxOD_Acinetobacter.baylyi, 
                                label=AA_mut_name_Ecoli))+
geom_point() +
geom_text(hjust=0, vjust=0)+
labs(title = "L.plantarum v A.baylyi maxOD")

lpab_MIC <- ggplot(MICwide, aes(x=Lactiplantibacillus.plantarum,
                                y=Acinetobacter.baylyi, 
                                label=AA_mut_name_Ecoli))+
  geom_point() +
  geom_text(hjust=0, vjust=0)+
  labs(title = "L.plantarum v A.baylyi MIC")+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')


#E.coli v L.plantarum 
eclp_mu<- ggplot(meanwide, aes(x=relative_mumax_Escherichia.coli,
                               y=relative_mumax_Lactiplantibacillus.plantarum, 
                               label=AA_mut_name_Ecoli))+
geom_point() +
geom_text(hjust=0, vjust=0)+
labs(title = "E.coli v L.planturum mumax")


eclp_OD <- ggplot(meanwide, aes(x=relative_maxOD_Escherichia.coli, 
                                y=relative_maxOD_Lactiplantibacillus.plantarum,
                                label=AA_mut_name_Ecoli))+
geom_point()+
geom_text(hjust=0, vjust=0)+
labs(title = "E.coli v L.plantarum maxOD")

eclp_MIC <- ggplot(MICwide, aes(x=Escherichia.coli, 
                                y=Lactiplantibacillus.plantarum,
                                label=AA_mut_name_Ecoli))+
  geom_point()+
  geom_text(hjust=0, vjust=0)+
  labs(title = "E.coli v L.plantarum MIC")+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')


fig_comp <- ggarrange(ecab_mu,lpab_mu, eclp_mu, 
                      ecab_OD, lpab_OD, eclp_OD, 
                      ecab_MIC, lpab_MIC, eclp_MIC,
                    ncol=3, nrow=3)
fig_comp

ggexport(fig_comp, filename = "comparison_parameters.pdf")

##########################
#graph parameters of all mutations 


#mumax
abmu <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=New_name, y=mumax))+
  geom_boxplot(aes(alpha=0))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.5) +
  theme_classic()+ 
  labs(title = "Acinetobacter baylyi")+
  theme(legend.position= "none", 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("")+
  coord_cartesian(ylim=c(0,0.1))

ecmu <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=New_name, y=mumax))+
  geom_boxplot(aes(alpha=0))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.5) +
  theme_classic()+ 
  labs(title = "Escherichia coli")+
  theme(legend.position= "none", 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("") +
  ylab("")+
  coord_cartesian(ylim=c(0,0.05))

lpmu <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Lactiplantibacillus plantarum") %>% 
  ggplot(aes(x=New_name, y=mumax))+
  geom_boxplot(aes(alpha=0))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.5) +
  theme_classic()+ 
  labs(title = "Lactiplantibacillus plantarum")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("")+
  ylab("")+
  coord_cartesian(ylim=c(0,0.0225))


#maxOD
abmaxOD <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=New_name, y=maxOD))+
  geom_boxplot(aes(alpha=0))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.5) +
  theme_classic()+ 
  #labs(title = "Acinetobacter baylyi")+
  theme(legend.position= "none", 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("")+
  coord_cartesian(ylim=c(0.3,1.15))

ecmaxOD <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=New_name, y=maxOD))+
  geom_boxplot(aes(alpha=0))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.5) +
  theme_classic()+ 
  #labs(title = "Escherichia coli")+
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("")+
  ylab("")+
  coord_cartesian(ylim=c(0.3,1.15))

lpmaxOD <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Lactiplantibacillus plantarum") %>% 
  ggplot(aes(x=New_name, y=maxOD))+
  geom_boxplot(aes(alpha=0))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.5) +
  theme_classic()+ 
  #labs(title = "Lactiplantibacillus plantarum")+
  theme(legend.position= "none", 
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("")+
  ylab("")+
  coord_cartesian(ylim=c(0.3,1.15))

#MICs

abMIC <- MICs %>% 
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=New_name, y=MIC)) +
  geom_point()+
  theme_classic()+ 
  xlab("")+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  scale_y_continuous(trans='log10')

ecMIC <- MICs %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=New_name, y=MIC)) +
  geom_point()+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("")+
  ylab("")+
  scale_y_continuous(trans='log10')

lpMIC <- MICs %>% 
  filter(Species =="Lactiplantibacillus plantarum") %>% 
  ggplot(aes(x=New_name, y=MIC)) +
  geom_point()+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+
  xlab("")+
  ylab("")+
  scale_y_continuous(trans='log10')
  


fig_species <- ggarrange(abmu,ecmu, lpmu,
                         abmaxOD, ecmaxOD, lpmaxOD,
                         abMIC, ecMIC, lpMIC,
                         ncol=3, nrow=3)
fig_species

ggexport(fig_species, filename = "species_parameters.pdf")


#empty working environment to keep everything clean:
rm.all()


# Step 5: Statistical Analysis                               ----
# __________________________________________________________

parms <- read_csv("./output/growth_parameters.csv")

parms$Strain <- as.factor(parms$Strain)
res_lm <- lm(maxOD ~ relevel(Strain, ref = "WT"),
             data = parms %>% 
               filter(Species == "Escherichia coli" & RIFconc == 0 ))
summary(res_lm)