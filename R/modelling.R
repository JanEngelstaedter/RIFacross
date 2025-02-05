library(tidyverse)

pleioParms <- read.csv("./output/growth_parameters_pleiotropy.csv")
parms <- read_csv("./output/growth_parameters_h20.csv")

overlap_muts <-  unique(pleioParms$AA_mut_name_Ecoli)

parms$AA_mut_name_Ecoli <- ifelse(parms$Strain == "WT", "WT" , parms$AA_mut_name_Ecoli)

modelling <- parms %>% 
  filter(RIFconc==0) %>% 
  filter(Species != "Lactiplantibacillus plantarum") %>% 
  filter(AA_mut_name_Ecoli %in% overlap_muts)

pleio_mod <- bind_rows(pleioParms, modelling)
pleio_mod$SetTemperature <- as.factor(pleio_mod$SetTemperature)
pleio_mod$AA_mut_name_Ecoli <- relevel(factor(pleio_mod$AA_mut_name_Ecoli), ref = "WT")



maxOD <- pleio_mod %>% 
  select(c(Species, AA_mut_name_Ecoli, SetTemperature, maxOD))
mumax <- pleio_mod %>% 
  select(c(Species, AA_mut_name_Ecoli, SetTemperature, mumax))

maxOD_nonneg <- maxOD %>% filter(maxOD > 0)
mod <- lm(maxOD ~ AA_mut_name_Ecoli*Species + 
            SetTemperature*Species + 
            AA_mut_name_Ecoli*SetTemperature, data = maxOD_nonneg)
summary(mod)

#############################################
# brms version
library(brms)
mod <- brm(maxOD ~ AA_mut_name_Ecoli*Species + 
            SetTemperature*Species + 
            AA_mut_name_Ecoli*SetTemperature, data = maxOD_nonneg)
summary(mod)
plot(mod)



#############################################

mod <- lm(maxOD ~ AA_mut_name_Ecoli*SetTemperature*Species, data = maxOD)
summary(mod)
plot(mod)

mod <- glm(maxOD ~ AA_mut_name_Ecoli*SetTemperature*Species, data = maxOD_nonneg, family = Gamma(link = "inverse"))
summary(mod)
plot(mod)


mod_mu <- glm(mumax ~ AA_mut_name_Ecoli*SetTemperature*Species, data = mumax, family = Gamma(link = "inverse"))
summary(mod_mu)

mod <- lm(maxOD ~ SetTemperature*AA_mut_name_Ecoli, data = maxOD)
summary(mod)


