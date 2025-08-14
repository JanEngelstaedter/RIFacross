# Step 0: Loading all libraries and scripts             ----
# __________________________________________________________

library(tidyverse)
library(grow96)
library(varhandle)
library(brms)
#options(brms.backend = "cmdstanr")
library(broom)
library(tidybayes)
library(ggrepel)
library(patchwork)
library(sjPlot)
library(xtable)
library(bridgesampling)
library(loo)
library(bayesplot)

source("./R/plotting.R")
source("./R/util.R")

globsets <- list(
  h = 20, # parameter for easylinear
  random_seed = 22,
  cred_int = 0.95)

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
  # adding additional information from project:
  analysis[[i]]$pars$Project_ID <- projects$Project_ID[i]
  analysis[[i]]$pars$Species <- projects$Species[i]
  cat("\n")
}
save(analysis, file = "./output/analysis_h20.RData")

# produce single table with parameter estimates for all species:
parms <- NULL
for(i in 1:nrow(projects)) {
  parms <- rbind(parms, analysis[[i]]$pars)
}

# add genetic information:
muts <- read_csv("./data/Genetic_data/muts.csv")
muts_Edward_lookup <- muts |>
  filter(Ref_code == "ALJElab (Edward)") |>
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) |>
  pull(Genotype, name = rpoB)

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
rm.all.but()



# Step 2: Loading and processing data for GxE analyses    ----
# __________________________________________________________
#input= raw OD data, spec files and genetic information
#output= GxE assay growth parameters (growth_parameters_GxE.csv)

projects <- read_csv("./data/global_key.csv") |>
  filter(Project_ID== "A_baylyi_Ian" | Project_ID == "E_coli_Alicia")
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

save(analysis, file = "./output/analysis_GxE.RData")

#Add genetic mutations 
load("./output/analysis_GxE.RData")

#remove A. baylyi replicates that didn't grow at 37degrees
analysis$A_baylyi_Ian$pars <- analysis$A_baylyi_Ian$pars %>% 
  filter(!(Date == "20/09/2024" & Strain %in% c("RM14", "R24")))

GxEParms <- NULL
# turn analysis list into single data frame:
for(i in 1:nrow(projects)) {
  analysis[[i]]$pars$Project_ID <- projects$Project_ID[i]
  analysis[[i]]$pars$Species <- projects$Species[i]
  GxEParms <- rbind(GxEParms, analysis[[i]]$pars)
}

muts <- read_csv("./data/Genetic_data/muts.csv")

GxEParms <- GxEParms |>
  left_join(muts, by = join_by(Species, Strain == Genotype)) |>
  select(Project_ID, Species, Strain,
         Nt_mut_name, Nt_mut_name_Ecoli, AA_mut_name, AA_mut_name_Ecoli, 
         SetTemperature,
         mumax, r2, maxOD) |>
  arrange(Species, Strain)

GxEParms$AA_pos <- as.numeric(gsub("[^0-9]","", GxEParms$AA_mut_name_Ecoli))
GxEParms <- GxEParms %>% arrange(AA_pos) 
GxEParms$AA_mut_name_Ecoli[GxEParms$Strain =="WT"] <- "WT"

GxEParms <- GxEParms %>%
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) 
GxEParms$rpoB[GxEParms$Strain =="WT"] <- "WT"

GxEParms$mumax[is.na(GxEParms$mumax)] <- 0 

write_csv(GxEParms, "./output/growth_parameters_GxE.csv")

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Step 3: Fitting separate models to E. coli and A. baylyi (for GxE analysis) ----
# ________________________________________________________________________________
#input= combined_parms_absence_antibiotics.csv
#output= brms_fit_GxE_models.RData containing four brms fits

parms_all <- read_csv("./output/combined_parms_absence_antibiotics.csv") |>
  filter(Species != "Lactiplantibacillus plantarum") |>
  mutate(Species = fct_relevel(factor(Species), "Escherichia coli")) |>
  mutate(SetTemperature = fct_relevel(factor(SetTemperature), "37")) |>
  mutate(rpoB = fct_relevel(factor(rpoB), "WT"))

# filter out non-overlapping mutations:
overlap <- intersect(parms_all |> filter(Species == "Escherichia coli") |> pull(rpoB),
                     parms_all |> filter(Species == "Acinetobacter baylyi") |> pull(rpoB))
parms_all <- parms_all |>
  filter(rpoB %in% overlap)

species <- levels(parms_all$Species)

all_fits <- c()

for (i in 1:length(species)) {    
  dat <- parms_all |>
    filter(parms_all$Species== species[i]) |>
    select(SetTemperature, rpoB, mumax, maxOD)
  all_fits[["mumax"]][[species[i]]] <- brm(mumax ~ rpoB * SetTemperature, data = dat)
  all_fits[["maxOD"]][[species[i]]] <- brm(maxOD ~ rpoB * SetTemperature, data = dat)
}

save(all_fits, file = "output/brms_fit_GxE_models.RData")

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Step 4: Fitting models to all species pairs (for GxG analysis) ----
# ___________________________________________________________________
#input=  combined_parms_absence_antibiotics.csv
#output= brms_fit_GxG_models.RData, containing six brms fits

parms_all <- read_csv("./output/combined_parms_absence_antibiotics.csv") |>
  mutate(rpoB = fct_relevel(factor(rpoB), "WT"))

species_pairs <- data.frame(Pair = c("EC_AB", "EC_LP", "AB_LP"),
                            Species1 = c(rep("Escherichia coli", 2), "Acinetobacter baylyi"),
                            Species2 = c("Acinetobacter baylyi", rep("Lactiplantibacillus plantarum", 2)),
                            SetTemperature1 = c("37", "37", "30"),
                            SetTemperature1 = c("30", "30", "30"))

all_fits <- c()

for (i in 1:nrow(species_pairs)) {
  # filter out non-overlapping mutations:
  overlap <- intersect(parms_all |> filter(Species == species_pairs$Species1[i]) |> pull(rpoB),
                       parms_all |> filter(Species == species_pairs$Species2[i]) |> pull(rpoB))
  dat <- parms_all |>
    filter(rpoB %in% overlap) |>
    filter(Species %in% species_pairs[i, c("Species1", "Species2")]) |>
    select(Species, rpoB, mumax, maxOD) |>
    mutate(Species = fct_relevel(factor(Species), species_pairs$Species1[i]))
  all_fits[["fits"]][["mumax"]][[species_pairs$Pair[i]]] <- brm(mumax ~ rpoB * Species, data = dat)
  all_fits[["fits"]][["maxOD"]][[species_pairs$Pair[i]]] <- brm(maxOD ~ rpoB * Species, data = dat)
}

all_fits[["legend"]] <- species_pairs

save(all_fits, file = "output/brms_fit_GxG_models.RData")


# Step 5: Fitting single model to E. coli & A. baylyi data (for GxGxE analysis) ----
# __________________________________________________________________________________
#input=  combined_parms_absence_antibiotics.csv
#output= brms_fit_GxGxE_model.RData containing two brms fits

parms_all <-  read_csv("./output/combined_parms_absence_antibiotics.csv" )

#find overlapping mutations 
overlapping <- parms_all %>% 
  filter(Species != "Lactiplantibacillus plantarum") %>% 
  distinct(Species, rpoB, Strain) %>% 
  pivot_wider(names_from = Species, values_from = Strain) %>% 
  drop_na("Acinetobacter baylyi") %>% 
  drop_na("Escherichia coli")

overlap <- unique(overlapping$rpoB)

GxE_mod <- parms_all %>% 
  filter(rpoB %in% overlap) %>% 
  filter(Species != "Lactiplantibacillus plantarum") |>
  mutate(rpoB = fct_relevel(factor(rpoB), "WT")) |>
  mutate(SetTemperature = fct_relevel(factor(SetTemperature), "37")) |>
  mutate(Species = fct_relevel(factor(Species), "Escherichia coli"))


# model fit without 3-way interactions:
maxOD_fit_no3way <- brm(maxOD ~ rpoB*Species + 
             SetTemperature*Species + 
             rpoB*SetTemperature, 
             data = GxE_mod)
summary(maxOD_fit_no3way)
#plot(maxOD_fit_no3way)

mumax_fit_no3way <- brm(mumax ~ rpoB*Species + 
                                  SetTemperature*Species + 
                                  rpoB*SetTemperature, 
                        data = GxE_mod)
summary(mumax_fit_no3way)
#plot(mumax_fit_no3way)

# model fit with 3-way interactions:
maxOD_fit <- brm(maxOD ~ rpoB*Species*SetTemperature, 
                 data = GxE_mod)
summary(maxOD_fit)
#plot(maxOD_fit)

mumax_fit <- brm(mumax ~ rpoB*Species*SetTemperature, 
                 data = GxE_mod)
summary(mumax_fit)
#plot(mumax_fit)

# comparing the simpler vs. the 3-way model:

pp_check(mumax_fit_no3way)
pp_check(mumax_fit)

bridge_no3way <- bridge_sampler(mumax_fit_no3way)
bridge <- bridge_sampler(mumax_fit)
bayes_factor(bridge, bridge_no3way)

loo_no3way <- loo(mumax_fit_no3way, moment_match = TRUE)
loo <- loo(mumax_fit, moment_match = TRUE)
loo_compare(loo_no3way, loo)

# saving fits
all_fits <- list(mumax = mumax_fit_no3way,
                 maxOD = maxOD_fit_no3way)
save(all_fits, file = "output/brms_fit_GxGxE_model_no3way.RData")

all_fits <- list(mumax = mumax_fit,
                 maxOD = maxOD_fit)
save(all_fits, file = "output/brms_fit_GxGxE_model.RData")

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Fig 1: Info plot explaining following figures ----
# __________________________________________________________
#input=  Bayesian relative parameters
#output= GxE graphs 

W1 <- 1
W2 <- 0.8
M1 <- 0.5
M2 <-0.6
arrow_offset <- 0.02
arrow_length <- unit(0.02, "npc")
label_a <- "italic(a)"
label_b2 <- "italic(b)[2]"
label_bM <- "italic(b)[M]"
label_bM2 <- "italic(b)[M2]"

p <- plot_epistasis_bg(W1, W2, 0, 1.2, 0, 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#F21A00") +
  scale_x_continuous(breaks = c(W1, M1), labels = c(expression(italic(W)[1]), expression(italic(M)[1])), 
                     limits = c(0, 1.2), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(W2, M2), labels = c(expression(italic(W)[2]), expression(italic(M)[2])),
                     limits = c(0, 1.2), expand = c(0, 0)) +
  labs(x = "Trait in species 1", y = "Trait in species 2") +
  annotate("point", W1, W2, size = 4) +
  annotate("text", x = W1 + 0.05, y = W2 - 0.03, label = "bold(W)", parse = TRUE, size = 5) +
  annotate("point", M1, M2, size = 4) +
  annotate("text", x = M1 - 0.04, y = M2 + 0.03, label = "bold(M)", parse = TRUE, size = 5) +
  annotate("segment", x = 0, y = W1, xend = W1, yend = W1,
           arrow = arrow(type = "closed", length = arrow_length)) +
  annotate("text", x = W1/2, y = W1 + 0.03, 
           label = label_a, parse = TRUE, size = 4) +
  annotate("segment", x = W1, y = W1, xend = W1, yend = W2 + arrow_offset,
           arrow = arrow(type = "closed", length = arrow_length)) +
  annotate("text", x = W1 + 0.04, y = W2 + (W1-W2)/2 + 0.01, 
           label = label_b2, parse = TRUE, size = 4) +
  annotate("segment", x = W1, y = W2 + 0.005, xend = M1, yend = W2 - (W1 - M1) + 0.005,
           arrow = arrow(type = "closed", length = arrow_length)) +
  annotate("text", x = M1 + (W1 - M1)/2 + 0.01, y = W2 - (W1 - M1) + (W1 - M1)/2 - 0.05, 
           label = label_bM, parse = TRUE, size = 4) +
  annotate("segment", x = M1, y = W2 - (W1 - M1), xend = M1, yend = M2 - arrow_offset,
           arrow = arrow(type = "closed", length = arrow_length)) +
  annotate("text", x = M1 + 0.04, y = W2 - (W1 - M1) + (M2 - W2 + (W1 - M1))/2, 
           label = label_bM2, parse = TRUE, size = 4) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 4 * 2.83),
        axis.text.y = element_text(size = 4 * 2.83), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

ggsave("plots/Fig1_explainer.png", plot = p, height = 5, width = 4.8)

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Fig 2: GxE scatterplots ----
# __________________________________________________________
#input=  Bayesian model fits
#output= GxE figure

set.seed(globsets$random_seed)

load("output/brms_fit_GxE_models.RData")
variables <- names(all_fits)
species <- names(all_fits[["mumax"]])
p <- c()

for(i in 1:length(variables)) {
  if (variables[i] == "mumax") {
    xlab <- expression(mu[max] * "(30"* degree *"C)")
    ylab <- expression(mu[max] *"(37"* degree *"C)")
  } else if (variables[i] == "maxOD") {
    xlab <- expression("maxOD (30"* degree *"C)")
    ylab <- expression("maxOD (37"* degree *"C)")
  }
  for(j in 1:length(species)) {
    fit <- all_fits[[variables[i]]][[species[j]]]
    newdata <- expand_grid(SetTemperature = levels(fit$data$SetTemperature),
                           rpoB = levels(fit$data$rpoB))
    post_preds <- fit |>
      epred_draws(newdata = newdata) |>
      group_by(across(all_of(colnames(newdata)))) |>
      mean_hdi(.epred, .width = globsets$cred_int)
    
    post_pred_1 <- post_preds |>
      filter(SetTemperature == "30") |>
      select(rpoB, .epred, .lower, .upper) |>
      rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)
    
    post_pred_2 <- post_preds |>
      filter(SetTemperature == "37") |>
      select(rpoB, .epred, .lower, .upper) |>
      rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)
    
    p[[variables[i]]][[species[j]]] <- plot_epistasis_abspars(post_pred_1, 
                                                              post_pred_2,
                                                              xlab, ylab, leglab = "GxE")
  }
}

Fig1 <- (p[["mumax"]][["Escherichia coli"]] + p[["mumax"]][["Acinetobacter baylyi"]])/
  (p[["maxOD"]][["Escherichia coli"]] + p[["maxOD"]][["Acinetobacter baylyi"]])
Fig1 <- Fig1 + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14),
        legend.position = "bottom")

ggsave(Fig1, filename = "./plots/Fig2_GxE_speciesfits.png",  width = 10, height = 10)

#empty working environment to keep everything clean:
rm.all.but("globsets")


# Fig 3: Epistasis scatterplots               ----
# __________________________________________________________
#input=  Bayesian relative parameters
#output= 

set.seed(globsets$random_seed)

load("output/brms_fit_GxG_models.RData")
variables <- names(all_fits$fits)
pairs <- names(all_fits$fits[["mumax"]])
species_pairs <- all_fits$legend
p <- c()

species_abbr <- c("Escherichia coli" = "E. coli",
                  "Acinetobacter baylyi" = "A. baylyi",
                  "Lactiplantibacillus plantarum" = "L. plantarum")

for(i in 1:length(variables)) {
  for(j in 1:length(all_fits$fits[[i]])) {
    if (variables[i] == "mumax") {
      xlab <- bquote(mu[max] * "(" * italic(.(species_abbr[species_pairs$Species1[j]])) * ")")
      ylab <- bquote(mu[max] * "(" * italic(.(species_abbr[species_pairs$Species2[j]])) * ")")
    } else if (variables[i] == "maxOD") {
      xlab <- bquote("maxOD(" * italic(.(species_abbr[species_pairs$Species1[j]])) * ")")
      ylab <- bquote("maxOD(" * italic(.(species_abbr[species_pairs$Species2[j]])) * ")")
    }
    
    fit <- all_fits$fits[[variables[i]]][[species_pairs$Pair[j]]]
    newdata <- expand_grid(Species = levels(fit$data$Species),
                           rpoB = levels(fit$data$rpoB))
    post_preds <- fit |>
      epred_draws(newdata = newdata) |>
      group_by(across(all_of(colnames(newdata)))) |>
      mean_hdi(.epred, .width = globsets$cred_int)
    
    post_pred_1 <- post_preds |>
      filter(Species == species_pairs$Species1[j]) |>
      select(rpoB, .epred, .lower, .upper) |>
      rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)
    
    post_pred_2 <- post_preds |>
      filter(Species == species_pairs$Species2[j]) |>
      select(rpoB, .epred, .lower, .upper) |>
      rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)
    
    p[[variables[i]]][[species_pairs$Pair[j]]] <- plot_epistasis_abspars(post_pred_1, 
                                                                         post_pred_2,
                                                                         xlab, ylab,
                                                                         leglab = "Epistasis")
  }
}

Fig3 <- (p[["mumax"]][["EC_AB"]] + p[["mumax"]][["EC_LP"]] + p[["mumax"]][["AB_LP"]])/
  (p[["maxOD"]][["EC_AB"]] + p[["maxOD"]][["EC_LP"]] + p[["maxOD"]][["AB_LP"]])
Fig3 <- Fig3 + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14),
        legend.position = "bottom")

ggsave(Fig3, filename = "./plots/Fig3_GxG_speciespairfits.png",  width = 15, height = 10)

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Fig 4: posterior distributions for single fit ----
# __________________________________________________________
#input=  Bayesian relative parameters
#output= GxE graphs 

# fits without 3-way interactions:
load("output/brms_fit_GxGxE_model_no3way.RData")

p_mumax <- plot_posteriors(all_fits$mumax)
p_maxOD <- plot_posteriors(all_fits$maxOD)

p <- (p_mumax | p_maxOD) + 
  plot_annotation(tag_levels = 'A')

ggsave("plots/GxGxE_no3way.pdf", p, width = 12, height = 8)


# fits with 3-way interactions:
load("output/brms_fit_GxGxE_model.RData")

p_mumax <- plot_posteriors(all_fits$mumax)
p_maxOD <- plot_posteriors(all_fits$maxOD)

p <- (p_mumax | p_maxOD) + 
  plot_annotation(tag_levels = 'A')

ggsave("plots/Fig4_GxGxE.pdf", p, width = 10, height = 8)

#empty working environment to keep everything clean:
rm.all.but("globsets")

# SupFig 1: fitness boxplots for all species and temperatures ----
# ________________________________________________________________
#input=  growth_parameters_h20.csv, parameters from MIC assay
#output= boxplots of species and mutants' fitness

set.seed(globsets$random_seed)

parms <- read_csv("./output/growth_parameters_h20.csv") |>
  filter(RIFconc < 1e-9) |>
  bind_rows(read_csv("./output/combined_parms_absence_antibiotics.csv" )) |>
  filter(!is.na(rpoB)) |>
  mutate(SetTemperature = fct_relevel(factor(SetTemperature), "30")) |>
  mutate(rpoB = fct_reorder(factor(rpoB), AA_pos, .na_rm = FALSE)) |>
  mutate(rpoB = fct_relevel(rpoB, "WT")) |>
  select(Species, SetTemperature, rpoB, mumax, maxOD)

cols <- c("30" = "#3B9AB2", "37" = "black")

p_mumax <- c()

for(species in unique(parms$Species)[c(2, 1, 3)]) {
  p_mumax[[species]] <- parms |>
    filter(Species == species) |>
    mutate(rpoB = droplevels(rpoB)) |>
    complete(SetTemperature, rpoB, fill = list(mumax = -11, maxOD = -11), explicit = FALSE) |>
    ggplot(aes(x=rpoB, y=mumax, col = SetTemperature, group = interaction(SetTemperature, rpoB))) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size = 0.5, alpha = 0.5) +
    theme_classic() +
    scale_color_manual(values = cols) +
    theme(legend.position= "none", 
          plot.title = element_text(size = 18, face = "italic"),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  size = 10)) +
    guides(col=guide_legend(title = "Temperature")) +
    coord_cartesian(ylim=c(-0.001,0.035))
  
  if (species == "Escherichia coli") {
    p_mumax[[species]] <- p_mumax[[species]] +
      labs(y= expression(mu[max]) , x= "", title = species)
  } else {
    p_mumax[[species]] <- p_mumax[[species]] +
      labs(y= "" , x= "", title = species)
  }
}

p_maxOD <- c()

for(species in unique(parms$Species)[c(2, 1, 3)]) {
  p_maxOD[[species]] <- parms |>
    filter(Species == species) |>
    mutate(rpoB = droplevels(rpoB)) |>
    complete(SetTemperature, rpoB, fill = list(mumax = -11, maxOD = -11), explicit = FALSE) |>
    ggplot(aes(x=rpoB, y=maxOD, col = SetTemperature, group = interaction(SetTemperature, rpoB))) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0), size = 0.5, alpha = 0.5) +
    theme_classic() +
    scale_color_manual(values = cols) +
    theme(legend.position= "bottom", 
          plot.title = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,  size = 10)) +
    guides(col=guide_legend(title = "Temperature")) +
    coord_cartesian(ylim=c(0.3, 1.15))
  
  if (species == "Escherichia coli") {
    p_maxOD[[species]] <- p_maxOD[[species]] +
      labs(y= "maxOD", x= "")
  } else {
    p_maxOD[[species]] <- p_maxOD[[species]] +
      labs(y= "" , x= "")
  }
}

SupFig1 <- p_mumax[["Escherichia coli"]] + p_mumax[["Acinetobacter baylyi"]] + p_mumax[["Lactiplantibacillus plantarum"]] +
  p_maxOD[["Escherichia coli"]] + p_maxOD[["Acinetobacter baylyi"]] + p_maxOD[["Lactiplantibacillus plantarum"]] +
  plot_layout(nrow = 2, widths = c(20, 20, 15.5), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(SupFig1, filename = "./plots/SupFig1_boxplots.pdf", width=15 , height= 8)

#empty working environment to keep everything clean:
rm.all.but("globsets")



# Tables for model fits                                 ----
# __________________________________________________________

# GxE models:
all_tables <- NULL
load("output/brms_fit_GxE_models.RData")
vars <- names(all_fits)
for(v in vars) {
  models <- names(all_fits[[v]])
  for(m in models) {
    all_tables <- brms::posterior_summary(all_fits[[v]][[m]], 
                          probs = c((1 - globsets$cred_int)/2, 0.5 + globsets$cred_int/2), 
                          central = FALSE) |>
    xtable(caption = posterior_table_caption("GxE", v, m, formula(all_fits[[v]][[m]])),
           digits = c(0, 4, 4, 4, 4),
           label = paste0("supptab:GxE_", v, "_", m)) |>
    print(type = "latex", file = "") |>
    capture.output() %>%
    c(all_tables, rep("", 3), .)
  }
}

# GxG models:
load("output/brms_fit_GxG_models.RData")
vars <- names(all_fits$fits)
for(v in vars) {
  models <- names(all_fits$fits[[v]])
  for(m in models) {
    all_tables <- brms::posterior_summary(all_fits$fits[[v]][[m]], 
                            probs = c((1 - globsets$cred_int)/2, 0.5 + globsets$cred_int/2), 
                            central = FALSE) |>
      xtable(caption = posterior_table_caption("GxG", v, m, formula(all_fits$fits[[v]][[m]])),
             digits = c(0, 4, 4, 4, 4),
             label = paste0("supptab:GxG_", v, "_", m)) |>
      print(type = "latex", file = "") |>
      capture.output() %>%
      c(all_tables, rep("", 3), .)
  }
}

# GxGxE models:
load("output/brms_fit_GxGxE_model.RData")
vars <- names(all_fits)
for(v in vars) {
  all_tables <- brms::posterior_summary(all_fits[[v]], 
                          probs = c((1 - globsets$cred_int)/2, 0.5 + globsets$cred_int/2), 
                          central = FALSE) |>
    xtable(caption = posterior_table_caption("GxGxE", v, formula = formula(all_fits[[v]])),
           digits = c(0, 4, 4, 4, 4),
           label = paste0("supptab:GxGxE_", v)) |>
    print(type = "latex", file = "", size = "\\footnotesize") |>
    capture.output() %>%
    c(all_tables, rep("", 3), .)
}

load("output/brms_fit_GxGxE_model_no3way.RData")
vars <- names(all_fits)
for(v in vars) {
  all_tables <- brms::posterior_summary(all_fits[[v]], 
                          probs = c((1 - globsets$cred_int)/2, 0.5 + globsets$cred_int/2), 
                          central = FALSE) |>
    xtable(caption = posterior_table_caption("GxGxE", v, formula = formula(all_fits[[v]])),
           digits = c(0, 4, 4, 4, 4),
           label = paste0("supptab:GxGxe_no3way_", v)) |>
    print(type = "latex", file = "", size = "\\footnotesize") |>
    capture.output() %>%
    c(all_tables, rep("", 3), .)
}

# polish and save as .tex file:
all_tables |>
  polish_posterior_tables() |>
  writeLines("tables/all_supp_tables.tex")



