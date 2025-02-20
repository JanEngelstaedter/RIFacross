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

#get one condition, set the level 
Testmuts <- pleioParms %>% 
  filter(Species == "Escherichia coli") +
  mutate(AA_mut_name_Ecoli = fct_relevel(factor(AA_mut_name_Ecoli), "WT")) 

#fit the model to all the replicates to get the mean and posteriors
fit <- brm(maxOD ~ AA_mut_name_Ecoli, data = Testmuts)
summary(fit)
tidy(fixef(fit))
plot(fit)



#will need to loop this to do for all mutants in the condition

#try using tidy.draws instead of as.draws
posterior_WT <- as_draws_df(fit) %>% pull(b_Intercept)
posterior_mut <- as_draws_df(fit) %>% pull(b_AA_mut_name_EcolirpoB_I572L)#in the loop this will be the "i", but remember that it has to have the b_ before the column and mutation name 

#will need to use these pulls to also get the absolute means 


#to get it relative 
ratio_mut_WT_maxod <- (posterior_WT + posterior_mut) / posterior_WT
hist(ratio_mut_WT_maxod)


#isn't the intercept the mean for the WT? why is it saying to mean the posteriors?
# Compute posterior mean and 95% credible interval
x_mean <- mean(posterior_x)
x_CI <- quantile(posterior_x, probs = c(0.025, 0.975))

# Print results
cat("Posterior mean of x:", x_mean, "\n")
cat("95% Credible Interval:", x_CI[1], "-", x_CI[2], "\n")







#this will be the way to calculate all means-
#start with the entire dataframe of parameters for all temps,species
#then use loops to filter through the species, temps
#write a function that does the brms for the mutants in each condition for the 2 parameters and the Credible intervals
#change to highest posterior density intervals


#when calling the function, be able to specify the species and temp, use this for the loops 
#edit the dataframe so that WT is the first level 
#then will need to extract the data 
#extract the absolute value (and the absolute CIs?)
#then calculate the relative value 
#also the relative CIs

#find a way to make this a nice table that the loop fills out with info 

pleioParms <- read.csv("./output/growth_parameters_pleiotropy.csv")
conditions <- read_csv("./data/conditions_key.csv")
dat <- c()


fit_maxOD <- c()
pull_maxOD <- c()
fit_mumax <- c()
pull_mumax <- c()

for (i in 1:nrow(conditions)) {
  cat(paste0(">>> Working on ", conditions$Species[i]," in ", conditions$Temperature[i]," degrees C maxOD" , "\n"))
  dat <- pleioParms %>% 
    filter(pleioParms$Species== conditions$Species[[i]]) 
  
  dat <- dat %>% 
    filter(dat$SetTemperature == conditions$Temperature[[i]])%>% 
    mutate(AA_mut_name_Ecoli = fct_relevel(factor(AA_mut_name_Ecoli), "WT")) 
  
  
  fit_maxOD[[conditions$Condition[i]]] <- brm(maxOD ~ AA_mut_name_Ecoli, data = dat)
  pull_maxOD[[conditions$Condition[i]]] <- fixef(fit_maxOD[[conditions$Condition[i]]])

  
  cat(paste0(">>> Working on ", conditions$Species[i]," in ", conditions$Temperature[i]," degrees C mumax" , "\n"))
  
  fit_mumax[[conditions$Condition[i]]] <- brm(mumax ~ AA_mut_name_Ecoli, data = dat)
  pull_mumax[[conditions$Condition[i]]] <- fixef(fit_mumax[[conditions$Condition[i]]])
 
  
  #above to get the absolute parameter values (need to format it to the intercept)
  

   
}


#getting absolute values 
bluh <- pull_maxOD[["EC30"]]

coef_adj <- c(bluh[1,1])
for (i in 2:nrow(bluh)){
  coef_adj[i] <- bluh[1,1] + bluh[i,1]
}

bluh[1,1] + bluh[i,1]




#use posteriors (fit) to get relative parameters and CIs


for (i in 1:nrow(conditions)){
  dat <- pleioParms %>% 
    filter(pleioParms$Species== conditions$Species[[i]]) 
  dat <- dat %>% 
    filter(dat$SetTemperature == conditions$Temperature[[i]])
  
  Strains <- unique(dat$AA_mut_name_Ecoli)
  
  for (i in 1:nrows(Strains))
  
  posterior_WT <- as_draws_df(fit_maxOD[[conditions$Condition[i]]]) %>% pull(b_Intercept)
  posterior_mut <- as_draws_df(fit_maxOD[[conditions$Condition[i]]]) %>% pull(#this needs to be the Strain in the loop fromt the condition)
  
  
  
  
}



posterior_WT <- as_draws_df(fit_maxOD$EC30) %>% pull(b_Intercept)
posterior_mut <quantile(ratio_mut_WT_maxod, probs = c(0.025, 0.975))- as_draws_df(fit_maxOD$EC30) %>% pull(b_AA_mut_name_EcolirpoB_I572L)#in the loop this will be the "i", but remember that it has to have the b_ before the column and mutation name 
ratio_mut_WT_maxod <- (posterior_WT + posterior_mut) / posterior_WT
hist(ratio_mut_WT_maxod)
x_mean <- mean(ratio_mut_WT_maxod)
x_CI <- 



summary(fit$EC30)

#SetTemp is an integer- maybe why not working cause in conditions is a number 




#try to get the fits out of being a sub-dataframe thing
#try to get the filtering for the temperature working too 

#when adding in all the data, will need to take the RIF concentrations into account!
#look at previous code, may have filtered out the RIf wells 

#then I can try to extract the values






#try to combine everything together 
pleioParms <- read.csv("./output/growth_parameters_pleiotropy.csv")
conditions <- read_csv("./data/conditions_key.csv")

pleioParms[is.na(pleioParms)] <- 0 #for the entire dataset, need to figure out how to deal with the NAs, if there are any



dat <- c()
ratio_maxOD <- c()
ratio_mumax <- c()
Strain_name <- c()
fit_maxOD <- c()
fit_mumax <- c()
maxOD <- c()
mumax <- c()

post_summaries <- pleioParms |>
  distinct(Species, Strain, SetTemperature) |>
  mutate(mean_maxOD = NA, mean_mumax = NA) |>
  arrange(Species, SetTemperature)

for (i in 1:nrow(conditions)) {
  dat <- pleioParms %>% 
    filter(pleioParms$Species== conditions$Species[i]) 
  
  dat <- dat %>% 
    filter(dat$SetTemperature == conditions$Temperature[[i]])%>% 
    mutate(AA_mut_name_Ecoli = fct_relevel(factor(AA_mut_name_Ecoli), "WT")) 
  
  Strains <- data.frame(AA_mut_name_Ecoli = unique(dat$AA_mut_name_Ecoli)) %>% filter(AA_mut_name_Ecoli != "WT")
  
  #calculate parameters
  cat(paste0(">>> Working on ", conditions$Species[i]," in ", conditions$Temperature[i]," degrees C " , "\n"))
  condition <- conditions$Condition[i]
  
  fit_maxOD[[condition]] <- brm(maxOD ~ AA_mut_name_Ecoli, data = dat)
  fit_mumax[[condition]] <- brm(mumax ~ AA_mut_name_Ecoli, data = dat)
  
  
  posterior_WT_maxOD[[condition]] <- as_draws_df(fit_maxOD[[condition]])%>% pull(b_Intercept)
  posterior_WT_mumax[[condition]] <- as_draws_df(fit_mumax[[condition]])%>% pull(b_Intercept)
    
  
  #if statement to get relative WT? need to make the name for b_intercept the same as the WT? or the WT the same as the b_intercept?
  
  
  for(j in 2:nrow(Strains)){
      Strain_name <- paste0("b_AA_mut_name_Ecoli", Strains$AA_mut_name_Ecoli[j])
      posterior_mut_maxOD <- as_draws_df(fit_maxOD[[condition]]) %>% pull(Strain_name) 
      posterior_mut_mumax <- as_draws_df(fit_mumax[[condition]]) %>% pull(Strain_name) 
        
      ratio_maxOD[[condition]][[Strain_name]]<- (posterior_WT_maxOD + posterior_mut_maxOD) / posterior_WT_maxOD
      ratio_mumax[[condition]][[Strain_name]]<- (posterior_WT_mumax + posterior_mut_mumax) / posterior_WT_mumax
  
      post_summaries$mean_maxOD[Strain == Strains[i] & 
                                Species == conditions$Species[i] & 
                                SetTemperature == conditions$Temperature[i]] <- mean(ratio_maxOD[[condition]][[Strain_name]])
      #maxOD[[condition]][[Strain_name]] <- mean(ratio_maxOD[[condition]][[Strain_name]])# add HDI here!
      #mumax[[condition]][[Strain_name]] <- mean(ratio_mumax[[condition]][[Strain_name]])
  }
  
  }



str(test)
test <- do.call(rbind, maxOD)

bluh <- test |> as_tibble() %>% mutate(condition = row.names(test)) %>% 
  pivot_longer(-condition, names_to = "mutation", values_to = "maxod") #add HDI columns!


#pivot wider table the rows are unique mutations? (check with format needed for each of the tables but I think this is needed for the comparison across species and temperature plots )
#bind maxOD and mumax tables together? (I think only needed for log transformed pleiotropy?)



#Error: Column `b_AA_mut_name_EcolirpoB_S509R` doesn't exist. for AB37- think this in because of mutants that don't grow and give NAs
#how to include the WT in making it relative 




x_mean <- mean(ratio_maxOD$EC30$b_AA_mut_name_EcolirpoB_S509R)
    
    
hist(posterior_mut_maxOD)
