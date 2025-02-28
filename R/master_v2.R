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
library(brms)
library(broom)
library(tidybayes)
library(ggrepel)

source("./R/plotting.R")
source("./R/util.R")


# Step 1: Loading and processing MIC assay data         ----
# __________________________________________________________
#input= raw OD data and spec files 
#output= growth parameters (analysis_h20.RData) and growth curve plots


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
  analysis[[i]] <- analyseODData(dat[[i]], tmin = 1, h = 20)
  cat("\n")
}
save(analysis, file = "./output/analysis_h20.RData")

#empty working environment to keep everything clean:
rm.all()


# Step 2: Add genetic information to data               ----
# __________________________________________________________
#input= analysis_h20.RData and genetic data from Google Sheet
#output= growth parameters with corresponding mutation (growth_parameters_h20.csv)

projects <- read_csv("./data/global_key.csv")
load("./output/analysis_h20.RData")

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
         Nt_mut_name, Nt_mut_name_Ecoli, AA_mut_name, AA_mut_name_Ecoli, SetTemperature,
         RIFconc, Replicate,
         mumax, r2, maxOD) |>
  arrange(Species, Strain, RIFconc, Replicate)

parms$AA_pos <- as.numeric(gsub("[^0-9]","", parms$AA_mut_name_Ecoli))
parms <- parms %>% arrange(AA_pos) %>% 
  mutate(New_name = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) %>% 
  filter(!(Species == "Acinetobacter baylyi" & (Strain %in% c("RM08", "RM15")))) # remove two mutants that are the same as RM10 and RM12
parms$New_name[parms$Nt_mut_name =="rpoB_dup1611"] <- "dup1611"
parms$New_name[parms$Strain =="WT"] <- "WT"

write_csv(parms, "./output/growth_parameters_h20.csv")

#empty working environment to keep everything clean:
rm.all()



# Step 3: Loading and processing pleiotropy data        ----
# __________________________________________________________
#input= raw OD data, spec files and genetic information
#output= Pleiotropy assay growth parameters (growth_parameters_pleiotropy.csv)


projects <- read_csv("./data/global_key.csv")
projects <- projects %>%  filter(Project_ID== "A_baylyi_Ian" | Project_ID == "E_coli_Alicia")
dat <- vector(mode = "list", length = nrow(projects))
names(dat) <- projects$Project_ID
analysis <- dat

for(i in 1:nrow(projects)) {
  cat("--------------------------------------------------------\n")
  cat(paste0(">>> Working on project ", projects$Project_ID[i], "\n"))
  cat("--------------------------------------------------------\n\n")
  dat[[i]] <- processODData(specPath = paste0("./data/", projects$Project_ID[i], "/specs_temp/"),
                            dataPath = paste0("./data/", projects$Project_ID[i], "/OD_temp"))
  dat[[i]] <- dat[[i]] |>
    blankODs(groups = "Strain") 
  dat[[i]] <- dat[[i]] |>
    filter(!is.na(Media))
  
  analysis[[i]] <- analyseODData(dat[[i]], tmin = 1, h = 20)
  cat("\n")
}

save(analysis, file = "./output/analysis_pleiotropy.RData")

#empty working environment to keep everything clean:
rm.all()


#Add genetic mutations 

projects <- read_csv("./data/global_key.csv")
projects <- projects %>%  filter(Project_ID== "A_baylyi_Ian" | Project_ID == "E_coli_Alicia")
load("./output/analysis_pleiotropy.RData")

pleioParms <- NULL
# turn analysis list into single data frame:
for(i in 1:nrow(projects)) {
  analysis[[i]]$pars$Project_ID <- projects$Project_ID[i]
  analysis[[i]]$pars$Species <- projects$Species[i]
  pleioParms <- rbind(pleioParms, analysis[[i]]$pars)
}

muts <- read_csv("./data/Genetic_data/muts.csv")

pleioParms <- pleioParms |>
  left_join(muts, by = join_by(Species, Strain == Genotype)) |>
  select(Project_ID, Species, Strain,
         Nt_mut_name, Nt_mut_name_Ecoli, AA_mut_name, AA_mut_name_Ecoli, 
         SetTemperature,
         mumax, r2, maxOD) |>
  arrange(Species, Strain)

pleioParms$AA_pos <- as.numeric(gsub("[^0-9]","", pleioParms$AA_mut_name_Ecoli))
pleioParms <- pleioParms %>% arrange(AA_pos) 
pleioParms$AA_mut_name_Ecoli[pleioParms$Strain =="WT"] <- "WT"

write_csv(pleioParms, "./output/growth_parameters_pleiotropy.csv")


# Step 4: Bayesian modelling of parameters              ----
# __________________________________________________________
#input= growth_parameters_pleiotropy.csv and growth_parameters_h20.csv
#output= 


pleioParms <- read_csv("./output/growth_parameters_pleiotropy.csv")
parms <- read_csv("./output/growth_parameters_h20.csv")
conditions <- read_csv("./data/conditions_key.csv")

#combine all replicates in the absence of antibiotics 
parms_0 <- parms %>% filter(RIFconc == 0 ) %>% 
  filter(!is.na(Strain))

parms_all <- bind_rows(pleioParms, parms_0) %>% select(-c(RIFconc, Replicate))



#remove this when nas are editted upstream
parms_all$mumax[is.na(parms_all$mumax)] <- 0 


dat <- c()
ratio_maxOD <- c()
ratio_mumax <- c()
Strain_name <- c()
fit_maxOD <- c()
fit_mumax <- c()
maxOD <- c()
mumax <- c()


for (i in 1:nrow(conditions)) {
  dat <- parms_all %>% 
    filter(parms_all$Species== conditions$Species[[i]]) 
  
  dat <- dat %>% 
    filter(dat$SetTemperature == conditions$Temperature[[i]])%>% 
    mutate(New_name = fct_relevel(factor(New_name), "WT")) 
  
  Strains <- data.frame(New_name = unique(dat$New_name)) %>% filter(New_name != "WT")
  
  #calculate parameters
  cat(paste0(">>> Working on ", conditions$Species[i]," in ", conditions$Temperature[i]," degrees C " , "\n"))
  condition <- conditions$Condition[i]
  
  fit_maxOD[[condition]] <- brm(maxOD ~ New_name, data = dat)
  fit_mumax[[condition]] <- brm(mumax ~ New_name, data = dat)
  
  
  posterior_WT_maxOD <- as_draws_df(fit_maxOD[[condition]])%>% pull(b_Intercept)
  posterior_WT_mumax <- as_draws_df(fit_mumax[[condition]])%>% pull(b_Intercept)
  
  
  #if statement to get relative WT? need to make the name for b_intercept the same as the WT? or the WT the same as the b_intercept?
  
  
  for(j in 1:nrow(Strains)){
    Strain_name <- paste0("b_New_name", Strains$New_name[j])
    posterior_mut_maxOD <- as_draws_df(fit_maxOD[[condition]]) %>% pull(Strain_name)
    posterior_mut_mumax <- as_draws_df(fit_mumax[[condition]]) %>% pull(Strain_name) 
    
    ratio_maxOD[[condition]][[Strain_name]]<- (posterior_WT_maxOD + posterior_mut_maxOD) / posterior_WT_maxOD
    ratio_mumax[[condition]][[Strain_name]]<- (posterior_WT_mumax + posterior_mut_mumax) / posterior_WT_mumax
    
    maxOD[[condition]][[Strain_name]] <- mean_hdi(ratio_maxOD[[condition]][[Strain_name]])
    mumax[[condition]][[Strain_name]] <- mean_hdi(ratio_mumax[[condition]][[Strain_name]])
    
    
  }
  
}

t_OD <- enframe(maxOD, name = "Condition", value = "relative_maxOD") %>% 
  unnest_longer(relative_maxOD, indices_to = "mutant")%>% 
  unnest(relative_maxOD) %>% 
  rename_with(~ paste("relative_maxOD", ., sep = "_"), -c(Condition, mutant))

t_mu <- enframe(mumax, name = "Condition", value = "relative_mumax") %>% 
  unnest_longer(relative_mumax, indices_to = "mutant")%>% 
  unnest(relative_mumax) %>% 
  rename_with(~ paste("relative_mumax", ., sep = "_"), -c(Condition, mutant))


rel_parms <- left_join(t_OD, t_mu, by= c("mutant", "Condition"))
rel_parms <- left_join(rel_parms, conditions, by = "Condition") %>% 
  select(Species, Temperature, Condition, mutant, 
         relative_maxOD_y, relative_maxOD_ymin, relative_maxOD_ymax,
         relative_mumax_y, relative_mumax_ymin, relative_mumax_ymax
         ) %>% 
  mutate(mutant = str_replace(mutant, "b_New_name", ""))

write_csv(rel_parms,"./output/bayesian_relative_parameters.csv" )

#wide dataframe for graphing
wide <- rel_parms %>% select( Condition , mutant, 
                             relative_maxOD_y, relative_maxOD_ymin, relative_maxOD_ymax,
                             relative_mumax_y, relative_mumax_ymin, relative_mumax_ymax
)
wide <- pivot_wider(wide,  names_from = Condition, values_from = c(relative_maxOD_y, relative_maxOD_ymin, relative_maxOD_ymax,
                                                                relative_mumax_y, relative_mumax_ymin, relative_mumax_ymax))
write_csv(wide,"./output/bayesian_relative_parameters_wide.csv" )






######TO DO 

#get rid of 37deg baylyi strain replicates that didn't grow (probably should do this in step 3)
#get relative WT and include HDI






# Step 5: Calculating MICs                              ----
# __________________________________________________________
#input= MIC assay parameters (growth_parameters_h20.csv)
#output= MICs.csv

parms <- read_csv("./output/growth_parameters_h20.csv")


MICs <- parms |>
  group_by(across(Species:AA_mut_name_Ecoli)) |>
  summarise(MIC = getMIC(RIFconc, maxOD))

WT <- MICs %>% filter (Strain=="WT")
WT_MIC <- setNames(WT$MIC,WT$Species)
MICs <- MICs %>% mutate(xfold=MIC/WT_MIC[Species])

MICs <- MICs %>% 
  mutate(New_name = str_replace(AA_mut_name_Ecoli, "rpoB_", ""))

write_csv(MICs, "./output/MICs.csv")



# Step 6: Epistasis and pleiotropy model                ----
# __________________________________________________________
#input=  
#output= 

#read combined df from step 4
#find overlapping Ecoli and Abaylyi mutations 
#filter so Rifconc=0, not Lplantarum and the mutants are only the overlap


maxOD_mod <- brm(maxOD ~ AA_mut_name_Ecoli*Species + 
             SetTemperature*Species + 
             AA_mut_name_Ecoli*SetTemperature, data = maxOD_nonneg)
summary(maxOD_mod)
plot(maxOD_mod)

mumax_mod <- brm(mumax ~ AA_mut_name_Ecoli*Species + 
                   SetTemperature*Species + 
                   AA_mut_name_Ecoli*SetTemperature, data = maxOD_nonneg)
summary(mumax_mod)
plot(mumax_mod)


# Fig 1: pleiotropy scatterplots                        ----
# __________________________________________________________
#input=  bayesian_relative_parameters_wide.csv
#output= pleiotropy graphs 

wide <- read_csv("./output/bayesian_relative_parameters_wide.csv" )




library(ggrepel)
#pleiotropy E. coli for 2 parameters at 2 temperatures

#mumax
pleio_mumax_EC <- ggplot(wide, aes(x= relative_mumax_y_EC30 , y=relative_mumax_y_EC37))+
  geom_point()+
  theme_classic()+
  labs(title = "mumax E.coli", x= "30 degrees C" , y= "37 degrees C")+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_mumax_ymin_EC30,
                     xmax= relative_mumax_ymax_EC30)
  )+
  geom_errorbar(aes(ymin= relative_mumax_ymin_EC37,
                    ymax= relative_mumax_ymax_EC37)
  )+
  geom_text_repel(aes(label= mutant ))


#maxOD
pleio_maxOD_EC <- ggplot(wide, aes(x= relative_maxOD_y_EC30 , y=relative_maxOD_y_EC37))+
  geom_point()+
  theme_classic()+
  labs(title = "maxOD E.coli", x= "30 degrees C" , y= "37 degrees C")+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_maxOD_ymin_EC30,
                     xmax= relative_maxOD_ymax_EC30)
  )+
  geom_errorbar(aes(ymin= relative_maxOD_ymin_EC37,
                    ymax= relative_maxOD_ymax_EC37)
  )+
  geom_text_repel(aes(label= mutant ))


#pleiotropy A. baylyi for 2 parameters at 2 temperatures

#mumax
pleio_mumax_AB <- ggplot(wide, aes(x= relative_mumax_y_AB30 , y=relative_mumax_y_AB37))+
  geom_point()+
  theme_classic()+
  labs(title = "mumax A. baylyi", x= "30 degrees C" , y= "37 degrees C")+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_mumax_ymin_AB30,
                     xmax= relative_mumax_ymax_AB30)
  )+
  geom_errorbar(aes(ymin= relative_mumax_ymin_AB37,
                    ymax= relative_mumax_ymax_AB37)
  )+
  geom_text_repel(aes(label= mutant ))


#maxOD
pleio_maxOD_AB <- ggplot(wide, aes(x= relative_maxOD_y_AB30 , y=relative_maxOD_y_AB37))+
  geom_point()+
  theme_classic()+
  labs(title = "maxOD A. baylyi", x= "30 degrees C" , y= "37 degrees C")+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_maxOD_ymin_AB30,
                     xmax= relative_maxOD_ymax_AB30)
  )+
  geom_errorbar(aes(ymin= relative_maxOD_ymin_AB37,
                    ymax= relative_maxOD_ymax_AB37)
  )+
  geom_text_repel(aes(label= mutant ))




#pleio_temp <- ggarrange(pleio_mumax, pleio_maxOD, 
           #             ncol = 1, nrow = 2)


#ggexport(pleio_temp, filename = "./plots/pleiotropy_scatter.pdf", 
    #     width = 15, 
     #    height = 15)








# Fig 2: log transformed pleiotropy                     ----
# __________________________________________________________
#input=  
#output= 

#previous test code 

test <- pleiowide %>% mutate(Ec_mumax_transformed = 
                               log2(mean_mumax_37_Escherichia.coli/mean_mumax_30_Escherichia.coli))

test <- test %>% mutate(Ec_maxOD_transformed = 
                          log2(mean_maxOD_37_Escherichia.coli/mean_maxOD_30_Escherichia.coli))

transformed_Ec <- ggplot(test, aes(Ec_mumax_transformed, Ec_maxOD_transformed)) +
  geom_point()+
  geom_vline(xintercept = mu_wt, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = od_wt, linetype = "dashed", color = "blue") +
  labs(title = "E. coli", 
       x = "mumax",
       y = "maxOD") 


wt_values <- test %>% filter(AA_mut_name_Ecoli == "WT")
mu_wt <- wt_values$Ec_mumax_transformed
od_wt <- wt_values$Ec_maxOD_transformed


test <- test %>% mutate(Ab_mumax_transformed = 
                          log2(mean_mumax_37_Acinetobacter.baylyi/mean_mumax_30_Acinetobacter.baylyi))

test <- test %>% mutate(Ab_maxOD_transformed = 
                          log2(mean_maxOD_37_Acinetobacter.baylyi/mean_maxOD_30_Acinetobacter.baylyi))

wt_values <- test %>% filter(AA_mut_name_Ecoli == "WT")
mu_wt_Ab <- wt_values$Ab_mumax_transformed
od_wt_Ab <- wt_values$Ab_maxOD_transformed

transformed_Ab <- ggplot(test %>% filter(AA_mut_name_Ecoli != "rpoB_S512P"), aes(Ab_mumax_transformed, Ab_maxOD_transformed)) +
  geom_point()+
  geom_vline(xintercept = mu_wt_Ab, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = od_wt_Ab, linetype = "dashed", color = "blue") +
  labs(title = "A. baylyi", 
       x = "mumax",
       y = "maxOD") 



# Fig 3: epistasis parameter scatterplots               ----
# __________________________________________________________
#input=  
#output= 


#previous code
##TO DO 
#use brms relative parameters and edit code 
#make sure all are on the same axes and labels are readable
#add in vertical and horizonal lines to mark the quandrants?

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



fig_comp <- ggarrange(ecab_mu,lpab_mu, eclp_mu, 
                      ecab_OD, lpab_OD, eclp_OD, 
                      ecab_MIC, lpab_MIC, eclp_MIC,
                      ncol=3, nrow=3)
fig_comp

ggexport(fig_comp, filename = "./plots/comparison_parameters.pdf")






# SupFig 1: fitness boxplots                            ----
# __________________________________________________________
#input=  
#output= 

#previous code
#TO DO
#edit the dimensions of the PDf so that it increasing with increasing boxplots (each boxplot is the same width)

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

ggexport(fig_species, filename = "./plots/species_parameters.pdf")




# SupFig 2: MIC growth curves                           ----
# __________________________________________________________
#input=  
#output= 

