library(tidyverse)
library(brms)
library(broom)

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

##########################################
#trying brms for propogated CI 

pleioParms <- read.csv("./output/growth_parameters_pleiotropy.csv")


Testmuts <- pleioParms %>% 
  filter(Species == "Escherichia coli") %>% 
  mutate(AA_mut_name_Ecoli = fct_relevel(factor(AA_mut_name_Ecoli), "WT")) 


# df_WT <- Testmuts %>% filter(AA_mut_name_Ecoli =="WT")
# df_mut <- Testmuts %>% filter(AA_mut_name_Ecoli =="rpoB_D516G")


fit_y <- brm(maxOD ~ 1, data = Testmuts %>% filter(AA_mut_name_Ecoli =="rpoB_D516G"))
summary(fit_y)
plot(fit_y)

fit_z <- brm(maxOD ~ AA_mut_name_Ecoli, data = Testmuts)
summary(fit_z)
tidy(fixef(fit_z))
plot(fit_z)
#try using tidy.draws instead of as.draws
posterior_WT <- as_draws_df(fit_z) %>% pull(b_Intercept)
posterior_mut <- as_draws_df(fit_z) %>% pull(b_AA_mut_name_EcolirpoB_I572L)

ratio_mut_WT_maxod <- (posterior_WT + posterior_mut) / posterior_WT
hist(ratio_mut_WT_maxod)

posterior_x <- posterior_y / posterior_z
posterior_x <- (intercept + posterior_mut) / intercept

# Compute posterior mean and 95% credible interval
x_mean <- mean(posterior_x)
x_CI <- quantile(posterior_x, probs = c(0.025, 0.975))

# Print results
cat("Posterior mean of x:", x_mean, "\n")
cat("95% Credible Interval:", x_CI[1], "-", x_CI[2], "\n")







#write a loop? have the WT and then loop through all the different mutants? or is there an easier way?
#will need loops inside loops for different temps and species? 


