# Step 0: Loading all libraries and scripts             ----
# __________________________________________________________

library(tidyverse)
library(grow96)
#library(Biostrings)
library(ALJEbinf)
library(googlesheets4)
library(UpSetR)
library(varhandle)
library(ggpubr)
library(brms)
library(broom)
library(tidybayes)
library(ggrepel)
library(patchwork)

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
#input= analysis_h20.RData and genetic data
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
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", ""))
muts_Edward_lookup <- muts_Edward$Strain_Name
names(muts_Edward_lookup) <- muts_Edward$rpoB

# attach rpoB mutations to growth parameters:
parms <- parms |>
  # standardise names for susceptible wildtype:
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
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) %>% 
  filter(!(Species == "Acinetobacter baylyi" & (Strain %in% c("RM08", "RM15")))) # remove two mutants that are the same as RM10 and RM12
parms$rpoB[parms$Nt_mut_name =="rpoB_dup1611"] <- "dup1611"
parms$rpoB[parms$Strain =="WT"] <- "WT"

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


#add the last replicate of A. baylyi
extra <- processODData(specPath = "./data/A_baylyi_Ian/specs_temp", 
                       dataPath = "./data/A_baylyi_Ian/OD_temp", 
                       filePrefix = "extra_") %>% 
  blankODs(groups = "ONCtemp") %>% 
  filter(ONCtemp == "37") %>% 
  filter(!is.na(Strain)) %>% 
  filter(! Well %in% c("C6", "C7", "D7")) %>% #abnormal replicates
  filter(!(Well == "D6" & Time_h > 15)) %>%  #Seems like mutation at end
  select(-ONCtemp)


extraAn <- analyseODData(extra,  tmin = 1, h = 20)

analysis$A_baylyi_Ian$pars <- analysis$A_baylyi_Ian$pars %>% 
  bind_rows(extraAn$pars)


save(analysis, file = "./output/analysis_pleiotropy.RData")


#empty working environment to keep everything clean:
rm.all()




#Add genetic mutations 

projects <- read_csv("./data/global_key.csv")
projects <- projects %>%  filter(Project_ID== "A_baylyi_Ian" | Project_ID == "E_coli_Alicia")
load("./output/analysis_pleiotropy.RData")



#remove A. baylyi replicates that didn't grow at 37degrees
analysis$A_baylyi_Ian$pars <- analysis$A_baylyi_Ian$pars %>% 
  filter(!(Date == "20/09/2024" & Strain %in% c("RM14", "R24")))
 


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

pleioParms <- pleioParms %>%
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) 
pleioParms$rpoB[pleioParms$Strain =="WT"] <- "WT"

pleioParms$mumax[is.na(pleioParms$mumax)] <- 0 

write_csv(pleioParms, "./output/growth_parameters_pleiotropy.csv")


# Step 4: Bayesian modelling of parameters              ----
# __________________________________________________________
#input= growth_parameters_pleiotropy.csv and growth_parameters_h20.csv
#output= Bayesian relative parameters bayesian_relative_parameters.csv, also in wide format


pleioParms <- read_csv("./output/growth_parameters_pleiotropy.csv")
parms <- read_csv("./output/growth_parameters_h20.csv")
conditions <- read_csv("./data/conditions_key.csv")

#combine all replicates in the absence of antibiotics 
parms_0 <- parms %>% filter(RIFconc == 0 ) %>% 
  filter(!is.na(Strain))

parms_all <- bind_rows(pleioParms, parms_0) %>% select(-c(RIFconc, Replicate))
write_csv(parms_all,"./output/combined_parms_absence_antibiotics.csv" )


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
    mutate(rpoB = fct_relevel(factor(rpoB), "WT")) 
  
  Strains <- data.frame(rpoB = unique(dat$rpoB)) 
  
  #calculate parameters
  cat(paste0(">>> Working on ", conditions$Species[i]," in ", conditions$Temperature[i]," degrees C " , "\n"))
  condition <- conditions$Condition[i]
  
  fit_maxOD[[condition]] <- brm(maxOD ~ rpoB, data = dat)
  fit_mumax[[condition]] <- brm(mumax ~ rpoB, data = dat)
  
  
  posterior_WT_maxOD <- as_draws_df(fit_maxOD[[condition]])%>% pull(b_Intercept)
  posterior_WT_mumax <- as_draws_df(fit_mumax[[condition]])%>% pull(b_Intercept)
  
  

  
  for(j in 1:nrow(Strains)){
    if (Strains$rpoB[j] == "WT") {
      Strain_name <- paste0("WT")
      posterior_mut_maxOD <- as_draws_df(fit_maxOD[[condition]]) %>% pull(b_Intercept)
      posterior_mut_mumax <- as_draws_df(fit_mumax[[condition]]) %>% pull(b_Intercept) 
      
      ratio_maxOD[[condition]][[Strain_name]]<- posterior_mut_maxOD / posterior_WT_maxOD
      ratio_mumax[[condition]][[Strain_name]]<- posterior_mut_mumax / posterior_WT_mumax
      
    } else {
      Strain_name <- paste0("b_rpoB", Strains$rpoB[j])
      posterior_mut_maxOD <- as_draws_df(fit_maxOD[[condition]]) %>% pull(Strain_name)
      posterior_mut_mumax <- as_draws_df(fit_mumax[[condition]]) %>% pull(Strain_name) 
      
      ratio_maxOD[[condition]][[Strain_name]]<- (posterior_mut_maxOD + posterior_WT_maxOD)/ posterior_WT_maxOD #need to add the intercept
      ratio_mumax[[condition]][[Strain_name]]<- (posterior_mut_mumax + posterior_WT_mumax)/ posterior_WT_mumax
      
    }
    
    maxOD[[condition]][[Strain_name]] <- mean_hdi(ratio_maxOD[[condition]][[Strain_name]])
    mumax[[condition]][[Strain_name]] <- mean_hdi(ratio_mumax[[condition]][[Strain_name]])
    
    
  }
  
}

t_OD <- enframe(maxOD, name = "Condition", value = "relative_maxOD") %>% 
  unnest_longer(relative_maxOD, indices_to = "rpoB")%>% 
  unnest(relative_maxOD) %>% 
  rename_with(~ paste("relative_maxOD", ., sep = "_"), -c(Condition, rpoB))

t_mu <- enframe(mumax, name = "Condition", value = "relative_mumax") %>% 
  unnest_longer(relative_mumax, indices_to = "rpoB")%>% 
  unnest(relative_mumax) %>% 
  rename_with(~ paste("relative_mumax", ., sep = "_"), -c(Condition, rpoB))


rel_parms <- left_join(t_OD, t_mu, by= c("rpoB", "Condition"))
rel_parms <- left_join(rel_parms, conditions, by = "Condition") %>% 
  select(Species, Temperature, Condition, rpoB, 
         relative_maxOD_y, relative_maxOD_ymin, relative_maxOD_ymax,
         relative_mumax_y, relative_mumax_ymin, relative_mumax_ymax
         ) %>% 
  mutate(rpoB = str_replace(rpoB, "b_rpoB", ""))

write_csv(rel_parms,"./output/bayesian_relative_parameters.csv" )

#wide dataframe for graphing
wide <- rel_parms %>% select( Condition , rpoB, 
                             relative_maxOD_y, relative_maxOD_ymin, relative_maxOD_ymax,
                             relative_mumax_y, relative_mumax_ymin, relative_mumax_ymax
)
wide <- pivot_wider(wide,  names_from = Condition, values_from = c(relative_maxOD_y, relative_maxOD_ymin, relative_maxOD_ymax,
                                                                relative_mumax_y, relative_mumax_ymin, relative_mumax_ymax))
write_csv(wide,"./output/bayesian_relative_parameters_wide.csv" )







# Step 5: Calculating MICs                              ----
# __________________________________________________________
#input= MIC assay parameters (growth_parameters_h20.csv)
#output= MICs.csv

parms <- read_csv("./output/growth_parameters_h20.csv")


MICs <- parms |>
  group_by(across(Species:AA_mut_name_Ecoli)) |>
  summarise(MIC = getMIC(RIFconc, maxOD))



MICs <- MICs %>%
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) 
MICs$rpoB[MICs$Strain =="WT"] <- "WT"
MICs$rpoB[MICs$Nt_mut_name =="rpoB_dup1611"] <- "dup1611"

MICs <- MICs %>%
  mutate(MIC = ifelse(
    MIC == Inf & (Species == "Acinetobacter baylyi" | Species == "Escherichia coli"), 4096,
    ifelse(MIC == Inf & Species == "Lactiplantibacillus plantarum", 512, MIC)
  ))


WT <- MICs %>% filter (Strain=="WT")
WT_MIC <- setNames(WT$MIC,WT$Species)
MICs <- MICs %>% mutate(xfold=MIC/WT_MIC[Species])



write_csv(MICs, "./output/MICs.csv")



# Step 6: Epistasis and pleiotropy model                ----
# __________________________________________________________
#input=  parameters in the absence of antibiotics
#output= brms model of epistasis and pleiotropy


parms_all <-  read_csv("./output/combined_parms_absence_antibiotics.csv" )

#find overlapping mutations 
overlapping <- parms_all %>% 
  filter(Species != "Lactiplantibacillus plantarum") %>% 
  distinct(Species, rpoB, Strain) %>% 
  pivot_wider(names_from = Species, values_from = Strain) %>% 
  drop_na("Acinetobacter baylyi") %>% 
  drop_na("Escherichia coli")

overlap <- unique(overlapping$rpoB)

pleio_mod <- parms_all %>% 
  filter(rpoB %in% overlap) %>% 
  filter(Species != "Lactiplantibacillus plantarum")


#some A.  baylyi have NAs- this is be gone after the last temp assay 
pleio_mod <- pleio_mod %>%
  mutate(across(c(mumax, maxOD), ~replace_na(., 0)))

dat_maxOD <- pleio_mod %>% 
  select(c(Species, rpoB, SetTemperature, maxOD))

dat_mumax <- pleio_mod %>% 
  select(c(Species, rpoB, SetTemperature, mumax))


maxOD_mod <- brm(maxOD ~ rpoB*Species + 
             SetTemperature*Species + 
             rpoB*SetTemperature, 
             data = dat_maxOD)
summary(maxOD_mod)
plot(maxOD_mod)
#try with the 3 way interaction- ask Andrew 


mumax_mod <- brm(mumax ~ rpoB*Species + 
                   SetTemperature*Species + 
                   rpoB*SetTemperature, 
                 data = dat_mumax)
summary(mumax_mod)
plot(mumax_mod)



# Fig 1: pleiotropy scatterplots                        ----
# __________________________________________________________
#input=  Bayesian relative parameters
#output= pleiotropy graphs 

wide <- read_csv("./output/bayesian_relative_parameters_wide.csv" )

set.seed(23) #for geom_repel

#pleiotropy E. coli for 2 parameters at 2 temperatures

#mumax
pleio_mumax_EC <- ggplot(wide, aes(x= relative_mumax_y_EC30 , y=relative_mumax_y_EC37))+
  geom_point()+
  theme_classic()+
  labs(title = expression(italic("E.coli")), x= expression( mu[max] * "(30"* degree *"C)") , y= expression(mu[max] *"(37"* degree *"C)"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_mumax_ymin_EC30,
                     xmax= relative_mumax_ymax_EC30),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_mumax_ymin_EC37,
                    ymax= relative_mumax_ymax_EC37),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )



#maxOD
pleio_maxOD_EC <- ggplot(wide, aes(x= relative_maxOD_y_EC30 , y=relative_maxOD_y_EC37))+
  geom_point()+
  theme_classic()+
  labs(title = "",x= expression("maxOD (30"* degree *"C)") , y= expression("maxOD (37"* degree *"C)"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_maxOD_ymin_EC30,
                     xmax= relative_maxOD_ymax_EC30),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_maxOD_ymin_EC37,
                    ymax= relative_maxOD_ymax_EC37),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


#pleiotropy A. baylyi for 2 parameters at 2 temperatures

#mumax
pleio_mumax_AB <- ggplot(wide, aes(x= relative_mumax_y_AB30 , y=relative_mumax_y_AB37))+
  geom_point()+
  theme_classic()+
  labs(title = expression(italic("A. baylyi")) , x= expression( mu[max] * "(30"* degree *"C)") , y= expression(mu[max]* "(37"* degree *"C)"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_mumax_ymin_AB30,
                     xmax= relative_mumax_ymax_AB30),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_mumax_ymin_AB37,
                    ymax= relative_mumax_ymax_AB37),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


#maxOD
pleio_maxOD_AB <- ggplot(wide, aes(x= relative_maxOD_y_AB30 , y=relative_maxOD_y_AB37))+
  geom_point()+
  theme_classic()+
  labs(title = "",x= expression("maxOD (30"* degree *"C)") , y= expression("maxOD (37"* degree *"C)"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  geom_errorbarh(aes(xmin= relative_maxOD_ymin_AB30,
                     xmax= relative_maxOD_ymax_AB30),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_maxOD_ymin_AB37,
                    ymax= relative_maxOD_ymax_AB37),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )




Fig1 <- (pleio_mumax_EC + pleio_mumax_AB )/(pleio_maxOD_EC + pleio_maxOD_AB)
Fig1 <- Fig1 + plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14))

ggexport(Fig1, filename = "./plots/Fig1_pleiotropy.pdf",  width = 10, height = 10)

#




# Fig 2: log transformed pleiotropy                     ----
# __________________________________________________________
#input=  Bayesian relative parameters 
#output= 

wide <- read_csv("./output/bayesian_relative_parameters_wide.csv" )

transformed <- wide %>% mutate(Ec_mumax_transformed = 
                               log2(relative_mumax_y_EC37/relative_mumax_y_EC30)) %>% 
                        mutate(Ec_mumax_min_transformed = 
                                 log2(relative_mumax_ymin_EC37/relative_mumax_ymin_EC30)) %>% 
                        mutate(Ec_mumax_max_transformed = 
                                 log2(relative_mumax_ymax_EC37/relative_mumax_ymax_EC30)) %>%
                        mutate(Ec_maxOD_transformed = 
                               log2(relative_maxOD_y_EC37/relative_maxOD_y_EC30)) %>%
                        mutate(Ec_maxOD_min_transformed = 
                                 log2(relative_maxOD_ymin_EC37/relative_maxOD_ymin_EC30)) %>% 
                        mutate(Ec_maxOD_max_transformed = 
                                 log2(relative_maxOD_ymax_EC37/relative_maxOD_ymax_EC30)) %>%
                        mutate(Ab_mumax_transformed = 
                              log2(relative_mumax_y_AB37/relative_mumax_y_AB30)) %>% 
                        mutate(Ab_mumax_min_transformed = 
                                 log2(relative_mumax_ymin_AB37/relative_mumax_ymin_AB30)) %>% 
                        mutate(Ab_mumax_max_transformed = 
                                 log2(relative_mumax_ymax_AB37/relative_mumax_ymax_AB30)) %>%
                        mutate(Ab_maxOD_transformed = 
                               log2(relative_maxOD_y_AB37/relative_maxOD_y_AB30)) %>% 
                        mutate(Ab_maxOD_min_transformed = 
                                 log2(relative_maxOD_ymin_AB37/relative_maxOD_ymin_AB30)) %>% 
                         mutate(Ab_maxOD_max_transformed = 
                                  log2(relative_maxOD_ymax_AB37/relative_maxOD_ymax_AB30)) 


####get WT relative parameters?
mu_wt <- log2(1/1)
od_wt <- log2(1/1)


set.seed(23) #for geom_repel

transformed_Ec_graph <- ggplot(transformed, aes(Ec_mumax_transformed, Ec_maxOD_transformed)) +
  geom_point()+
  theme_classic()+
  labs(title = expression(italic("E. coli")), 
       x = expression("log"[2]* "(" *mu[max] * "37/30)") ,
       y = expression("log"[2]* "(maxOD 37/30)") )+
 geom_vline(xintercept = mu_wt, linetype = "dotted", color = "blue") +
 geom_hline(yintercept = od_wt, linetype = "dotted", color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA) +
  geom_errorbarh(aes(xmin= Ec_mumax_min_transformed,
                     xmax= Ec_mumax_max_transformed),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= Ec_maxOD_min_transformed,
                    ymax= Ec_maxOD_max_transformed),
                alpha = 0.3
  )+
  theme(
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
 

transformed_Ab_graph <- ggplot(transformed, aes(Ab_mumax_transformed, Ab_maxOD_transformed)) +
  geom_point()+
  theme_classic()+
  labs(title = expression(italic("A. baylyi")), 
       x = expression("log"[2]* "(" *mu[max] * "37/30)") ,
       y = expression("log"[2]* "(maxOD 37/30)") )+
  geom_vline(xintercept = mu_wt, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = od_wt, linetype = "dotted", color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA) +
  geom_errorbarh(aes(xmin= Ab_mumax_min_transformed,
                     xmax= Ab_mumax_max_transformed),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= Ab_maxOD_min_transformed,
                    ymax= Ab_maxOD_max_transformed),
                alpha = 0.3
  )+
  theme(
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10), 
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


Fig2 <- transformed_Ec_graph + transformed_Ab_graph
Fig2 <- Fig2 + plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14))
Fig2

ggexport(Fig2, filename = "./plots/Fig2_transformed_pleiotropy.pdf",  width = 10, height = 5)

#

# Fig 3: epistasis parameter scatterplots               ----
# __________________________________________________________
#input=  Bayesian relative parameters
#output= 

wide <- read_csv("./output/bayesian_relative_parameters_wide.csv" )
set.seed(23) #for geom_repel

##TO DO 
#make sure all are on the same axes and labels are readable
#add in vertical and horizonal lines to mark the quandrants?

ecab_mu <- ggplot(wide, aes(x=relative_mumax_y_EC37, 
                               y=relative_mumax_y_AB30, 
                               label=rpoB))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  coord_cartesian(xlim = c(0.425, 1.4), ylim =c(0.425, 1.4))+
  labs(x =expression(mu* "max " *italic("(E. coli)")) , y =expression(mu* "max " * italic( "(A. baylyi)")))+
  theme_classic() +
  geom_errorbarh(aes(xmin= relative_mumax_ymin_EC37,
                     xmax= relative_mumax_ymax_EC37),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_mumax_ymin_AB30,
                    ymax= relative_mumax_ymax_AB30),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )



ecab_OD <- ggplot(wide, aes(x=relative_maxOD_y_EC37, 
                            y=relative_maxOD_y_AB30, 
                            label=rpoB))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  coord_cartesian(xlim = c(0.8, 1.36), ylim =c(0.8, 1.36))+
  labs(x= expression("maxOD " * italic("(E. coli)")) , y=expression("maxOD " * italic( "(A. baylyi)")))+
  theme_classic() +
  geom_errorbarh(aes(xmin= relative_maxOD_ymin_EC37,
                     xmax= relative_maxOD_ymax_EC37),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_maxOD_ymin_AB30,
                    ymax= relative_maxOD_ymax_AB30),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )




#L.plantarum v A.baylyi
lpab_mu <- ggplot(wide, aes(x=relative_mumax_y_LP30, 
                            y=relative_mumax_y_AB30, 
                            label=rpoB))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  coord_cartesian(xlim = c(0.425, 1.4), ylim =c(0.425, 1.4))+
  labs(x= expression(mu * "max "* italic("(L. plantarum)")), y=expression(mu* "max "* italic("(A. baylyi)")))+
  theme_classic() +
  geom_errorbarh(aes(xmin= relative_mumax_ymin_LP30,
                     xmax= relative_mumax_ymax_LP30),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_mumax_ymin_AB30,
                    ymax= relative_mumax_ymax_AB30),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )



lpab_OD <- ggplot(wide, aes(x=relative_maxOD_y_LP30, 
                            y=relative_maxOD_y_AB30, 
                            label=rpoB))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  coord_cartesian(xlim = c(0.8, 1.36), ylim =c(0.8, 1.36))+
  labs(x= expression("maxOD " * italic("(L. plantarum)"))  , y= expression("maxOD "* italic("(A. baylyi)")))+
  theme_classic() +
  geom_errorbarh(aes(xmin= relative_maxOD_ymin_LP30,
                     xmax= relative_maxOD_ymax_LP30),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_maxOD_ymin_AB30,
                    ymax= relative_maxOD_ymax_AB30),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )





#E.coli v L.plantarum 
eclp_mu<- ggplot(wide, aes(x=relative_mumax_y_EC37, 
                           y=relative_mumax_y_LP30, 
                           label=rpoB))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  coord_cartesian(xlim = c(0.425, 1.4), ylim =c(0.425, 1.4))+
  labs(x= expression(mu* "max " *italic("(E. coli)")), y= expression(mu* "max " *italic("(L. planturum)")) )+
  theme_classic() +
  geom_errorbarh(aes(xmin= relative_mumax_ymin_EC37,
                     xmax= relative_mumax_ymax_EC37),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_mumax_ymin_LP30,
                    ymax= relative_mumax_ymax_LP30),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )





eclp_OD <- ggplot(wide, aes(x=relative_maxOD_y_EC37, 
                            y=relative_maxOD_y_LP30, 
                            label=rpoB))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+
  coord_cartesian(xlim = c(0.8, 1.36), ylim =c(0.8, 1.36))+
  labs(x= expression("maxOD " *italic("(E. coli)")), y= expression("maxOD " *italic("(L. planturum)")))+
  theme_classic() +
  geom_errorbarh(aes(xmin= relative_maxOD_ymin_EC37,
                     xmax= relative_maxOD_ymax_EC37),
                 alpha = 0.3
  )+
  geom_errorbar(aes(ymin= relative_maxOD_ymin_LP30,
                    ymax= relative_maxOD_ymax_LP30),
                alpha = 0.3
  )+
  geom_label_repel(aes(label = rpoB), 
                   size = 4,
                   fill = alpha("white", 0.3), 
                   label.size = NA,
                   segment.color = NA)+
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )




Fig3 <- (ecab_mu + eclp_mu + lpab_mu)/ (ecab_OD + eclp_OD + lpab_OD)
Fig3 <- Fig3 +plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 16))
Fig3

ggexport(Fig3, filename = "./plots/Fig3_epistasis.pdf",  width = 15, height = 10)

#

# SupFig 1: fitness boxplots                            ----
# __________________________________________________________
#input=  growth_parameters_h20.csv, parameters from MIC assay
#output= boxplots of species and mutants' fitness

parms_all <- read_csv("./output/combined_parms_absence_antibiotics.csv" )

parms <- read_csv("./output/growth_parameters_h20.csv")
parms <- parms %>% filter(!is.na(rpoB))
MICs <- read_csv("./output/MICs.csv")
MICs <- MICs %>% filter(!is.na(rpoB))

#mumax
abmu <- parms %>% 
  filter(RIFconc==0) %>%
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=rpoB, y=mumax))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
    labs(title = expression(italic("A. baylyi")), y= "" , x= "")+
    theme(legend.position= "none", 
        plot.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  size = 10))+
  coord_cartesian(ylim=c(0.002,0.035))

ecmu <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=rpoB, y=mumax))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  labs(title = expression(italic("E. coli")), y= expression( mu[max]) , x= "")+
  theme(legend.position= "none", 
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  size = 10))+
  coord_cartesian(ylim=c(0.002,0.035))

lpmu <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Lactiplantibacillus plantarum") %>% 
  ggplot(aes(x=rpoB, y=mumax))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  labs(title = expression(italic("L. plantarum")))+
  theme(legend.position= "none",
        plot.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  size = 10))+
  xlab("")+
  ylab("")+
  coord_cartesian(ylim=c(0.002,0.035))


#maxOD
abmaxOD <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=rpoB, y=maxOD))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  theme(legend.position= "none", 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size= 10))+
  xlab("")+
  ylab("")+
  coord_cartesian(ylim=c(0.3,1.15))

ecmaxOD <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=rpoB, y=maxOD))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10))+
  xlab("")+
    coord_cartesian(ylim=c(0.3,1.15))

lpmaxOD <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species =="Lactiplantibacillus plantarum") %>% 
  ggplot(aes(x=rpoB, y=maxOD))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  theme(legend.position= "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10))+
  xlab("")+
  ylab("")+
  coord_cartesian(ylim=c(0.3,1.15)) 




#MICS

abMIC <- MICs %>% 
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=rpoB, y=MIC)) +
  geom_point()+
  theme_classic()+ 
  xlab("")+
  ylab("")+
  theme(axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1, size = 10))+
  scale_y_continuous(trans='log10')

ecMIC <- MICs %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=rpoB, y=MIC)) +
  geom_point()+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1, size = 10))+
  xlab("")+
  ylab("MIC")+
  scale_y_continuous(trans='log10')

lpMIC <- MICs %>% 
  filter(Species =="Lactiplantibacillus plantarum") %>% 
  ggplot(aes(x=rpoB, y=MIC)) +
  geom_point()+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1, size = 10))+
  xlab("")+
  ylab("")+
  scale_y_continuous(trans='log10')





SupFig1 <- (ecmu + abmu + lpmu)/(ecmaxOD + abmaxOD + lpmaxOD) / (ecMIC + abMIC +lpMIC)
SupFig1 

ggexport(SupFig1, filename = "./plots/SupFig1_boxplots.pdf", width=15 ,height= 15 )




# SupFig 2: Temperature assay absolute fitness          ----
# __________________________________________________________
#input=  opposite temperature values
#output= absolute fitness boxplots

pleioParms <- read_csv("./output/growth_parameters_pleiotropy.csv")


abmutemp <- pleioParms %>% 
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=rpoB, y=mumax))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  labs(title = expression(italic("A. baylyi")), y= "" , x= "")+
  theme(legend.position= "none", 
        plot.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  size = 10))+
  coord_cartesian(ylim=c(0.00,0.0215))

ecmutemp <- pleioParms %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=rpoB, y=mumax))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  labs(title = expression(italic("E. coli")), y= expression( mu[max]) , x= "")+
  theme(legend.position= "none", 
        plot.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  size = 10))+
  coord_cartesian(ylim=c(0.00,0.0215))



abmaxODtemp <- pleioParms %>% 
  filter(Species =="Acinetobacter baylyi") %>% 
  ggplot(aes(x=rpoB, y=maxOD))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  theme(legend.position= "none", 
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size= 10))+
  xlab("")+
  ylab("")+
  coord_cartesian(ylim=c(0,1.16))

ecmaxODtemp <- pleioParms %>% 
  filter(Species =="Escherichia coli") %>% 
  ggplot(aes(x=rpoB, y=maxOD))+
  geom_boxplot(aes(alpha=0.5))+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL), alpha = 0.3) +
  theme_classic()+ 
  theme(legend.position= "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10))+
  xlab("")+
  coord_cartesian(ylim=c(0,1.16))



SupFig2 <- (ecmutemp + abmutemp) / (ecmaxODtemp + abmaxODtemp)
SupFig2

ggexport(SupFig2, filename = "./plots/SupFig2_boxplots_temp.pdf", width=10 ,height= 8 )








