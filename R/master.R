# Step 0: Loading all libraries and scripts             ----
# __________________________________________________________

here::i_am("R/master.R")
library(here)
library(tidyverse) 
library(grow96) # available from github via pak::pak("janengelstaedter/grow96")
library(varhandle)
library(brms)
options(brms.backend = "cmdstanr")
# (Install using install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos"))))
# or comment out to use rstan instead.)
library(broom)
library(tidybayes)
library(ggrepel)
library(patchwork)
library(sjPlot)
library(xtable)
library(bridgesampling)
library(loo)
library(bayesplot)

source(here("R", "plotting.R"))
source(here("R", "util.R"))

globsets <- list(
  h = 20, # parameter for easylinear
  random_seed = 22,
  cred_int = 0.95)

# create output folders if they don't exist yet:
for (folder in c("output", "plots", "tables")) {
  if (!dir.exists(here(folder))) dir.create(here(folder))
}

# Step 1: Importing and processing OD data             ----
# __________________________________________________________
# input= raw OD data and spec files 
# output= processed OD data

projects <- read_csv(here("data", "global_key.csv"))
dat <- vector(mode = "list", length = nrow(projects))
names(dat) <- projects$Project_ID

for(i in 1:nrow(projects)) {
  # wrangling:
  dat[[i]] <- processODData(specPath = here("data", projects$Project_ID[i], "specs"),
                            dataPath = here("data", projects$Project_ID[i], "OD_raw"))
  if(projects$Project_ID[i] == "A_baylyi_Ian") {
    dat[[i]] <- dplyr::rename(dat[[i]], RIFconc = Concentration)  
  }
  # manual correction of error where column in plate wasn't inoculated:
  if (projects$Project_ID[i] == "A_baylyi_Maddie") {
    dat[[i]] <- dat[[i]] |>
      mutate(WellType = if_else(Plate == "MIC_512_strains4" & Replicate == "1" & Row == "G", "EMPTY", WellType))
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
}

# additional data:
projects <- read_csv(here("data", "global_key.csv")) |>
  filter(Project_ID == "A_baylyi_Ian" | Project_ID == "E_coli_Alicia")
dat_add <- vector(mode = "list", length = nrow(projects))
names(dat_add) <- paste0(projects$Project_ID, "_add")

for(i in 1:nrow(projects)) {
  dat_add[[i]] <- processODData(specPath = here("data", projects$Project_ID[i], "specs_temp"),
                            dataPath = here("data", projects$Project_ID[i], "OD_temp"))
  dat_add[[i]] <- dat_add[[i]] |>
    blankODs(groups = "Strain") 
  dat_add[[i]] <- dat_add[[i]] |>
    filter(!is.na(Media))
}

#add the last replicate of A. baylyi
extra <- processODData(specPath = here("data", "A_baylyi_Ian", "specs_temp"), 
                       dataPath = here("data", "A_baylyi_Ian", "OD_temp")) %>% 
  filter(ONCtemp == "37") %>% 
  filter(!is.na(Strain)) %>% 
  blankODs() %>% 
  filter(! Well %in% c("C6", "C7", "D7")) %>% #abnormal replicates
  filter(!(Well == "D6" & Time_h > 15)) %>%  #Seems like mutation at end
  select(-ONCtemp)

# merge all data into single object and save it:
all_OD_data <- c(dat, dat_add, list(A_baylyi_Ian_extra = extra))

save(all_OD_data, file = here("output", "all_OD_data.RData"))

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Step 2: Analysing OD data                             ----
# __________________________________________________________
# input= project metadata, processed OD data, genetic data
# output= file with all estimated growth parameters

projects <- read_csv(here("data", "global_key.csv"))
load(file = here("output", "all_OD_data.RData"))
muts <- read_csv(here("data", "Genetic_data", "muts.csv"))

# merge all tables into single one and analyse:
OD_data_filtered <- all_OD_data |>
  bind_rows(.id = "Project_ID") |>
  filter(is.na(RIFconc) | RIFconc < 1e-9) |>
  filter(WellType == "DATA") |>
  mutate(Project_ID = str_replace(Project_ID, "_add", "")) |>
  mutate(Project_ID = str_replace(Project_ID, "_extra", "")) |>
  left_join(projects, by = join_by(Project_ID)) |>
  filter(!(Species == "Acinetobacter baylyi" & (Strain %in% c("RM08", "RM15")))) |> # remove two mutants that are the same as RM10 and RM12
  filter(!(Species == "Acinetobacter baylyi" & Date == "20/09/2024" & Strain %in% c("RM14", "R24"))) #remove A. baylyi replicates that didn't grow at 37degrees

growth_parms <- analyseODData(OD_data_filtered, method = "easylinear", tmin = 1, h = 20)$pars |> as_tibble()

# adding genetic information
muts_Edward_lookup <- muts |>
  filter(Ref_code == "ALJElab (Edward)") |>
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) |>
  pull(Genotype, name = rpoB)

growth_parms <- growth_parms |>
  # standardise names for susceptible wildtype:
  mutate(Strain = ifelse(Strain %in% c("AB3", "MG1655", "Wild Type"), 
                         "WT", Strain)) |>
  # fix Edward's strain names:
  mutate(Strain = ifelse((Project_ID == "L_plantarum_Edward") & (Strain != "WT"), 
                         muts_Edward_lookup[Strain],
                         Strain)) |>
  left_join(muts, by = join_by(Species, Strain == Genotype)) |>
  select(Project_ID, 
         Species, Strain, Nt_mut_name, Nt_mut_name_Ecoli, AA_mut_name, AA_mut_name_Ecoli, 
         SetTemperature, Replicate,
         mumax, r2, maxOD) |>
  mutate(AA_pos = as.numeric(gsub("[^0-9]","", AA_mut_name_Ecoli))) |>
  mutate(rpoB = str_replace(AA_mut_name_Ecoli, "rpoB_", "")) |>
  mutate(rpoB = if_else(Nt_mut_name == "rpoB_dup1611", "dup1611", rpoB)) |>
  mutate(rpoB = if_else(Strain == "WT", "WT", rpoB)) |>
  arrange(Species, Strain, Replicate)

write_csv(growth_parms, file = here("output", "growth_parameters_easylinear_h20.csv"))

#empty working environment to keep everything clean:
rm.all.but("globsets")


# Step 3: Fitting separate models to E. coli and A. baylyi (GxE analysis) ----
# ________________________________________________________________________________
# input= combined_parms_absence_antibiotics.csv
#o utput= brms_fit_GxE_models.RData containing four brms fits

growth_parms <- read_csv(here("output", "growth_parameters_easylinear_h20.csv")) |>
  filter(Species != "Lactiplantibacillus plantarum") |>
  mutate(Species = fct_relevel(factor(Species), "Escherichia coli")) |>
  mutate(SetTemperature = fct_relevel(factor(SetTemperature), "37")) |>
  mutate(rpoB = fct_relevel(factor(rpoB), "WT"))

# filter out non-overlapping mutations:
overlap <- intersect(growth_parms |> filter(Species == "Escherichia coli") |> pull(rpoB),
                     growth_parms |> filter(Species == "Acinetobacter baylyi") |> pull(rpoB))
growth_parms <- growth_parms |>
  filter(rpoB %in% overlap)

species <- levels(growth_parms$Species)

all_fits <- c()

for (i in 1:length(species)) {    
  dat <- growth_parms |>
    filter(growth_parms$Species== species[i]) |>
    select(SetTemperature, rpoB, mumax, maxOD)
  all_fits[["mumax"]][[species[i]]] <- brm(mumax ~ rpoB * SetTemperature, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
  all_fits[["maxOD"]][[species[i]]] <- brm(maxOD ~ rpoB * SetTemperature, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
}

save(all_fits, file = here("output", "brms_fit_GxE_models.RData"))

all_fits_noint <- c()

for (i in 1:length(species)) {    
  dat <- growth_parms |>
    filter(growth_parms$Species== species[i]) |>
    select(SetTemperature, rpoB, mumax, maxOD)
  all_fits_noint[["mumax"]][[species[i]]] <- brm(mumax ~ rpoB + SetTemperature, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
  all_fits_noint[["maxOD"]][[species[i]]] <- brm(maxOD ~ rpoB + SetTemperature, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
}

save(all_fits_noint, file = here("output", "brms_fit_GplusE_models.RData"))

# comparing the simpler vs. the interactions models:
load(file = here("output", "brms_fit_GplusE_models.RData"))
load(file = here("output", "brms_fit_GxE_models.RData"))

# E. coli
# mu_max 
loo_compare(all_fits_noint[["mumax"]][["Escherichia coli"]], all_fits[["mumax"]][["Escherichia coli"]])
bayes_R2(all_fits_noint[["mumax"]][["Escherichia coli"]])
bayes_R2(all_fits[["mumax"]][["Escherichia coli"]])

# max OD
loo_compare(all_fits_noint[["maxOD"]][["Escherichia coli"]], all_fits[["maxOD"]][["Escherichia coli"]])
bayes_R2(all_fits_noint[["maxOD"]][["Escherichia coli"]])
bayes_R2(all_fits[["maxOD"]][["Escherichia coli"]])

# A. baylyi
# mu_max 
loo_compare(all_fits_noint[["mumax"]][["Acinetobacter baylyi"]], all_fits[["mumax"]][["Acinetobacter baylyi"]])
bayes_R2(all_fits_noint[["mumax"]][["Acinetobacter baylyi"]])
bayes_R2(all_fits[["mumax"]][["Acinetobacter baylyi"]])

# max OD
loo_compare(all_fits_noint[["maxOD"]][["Acinetobacter baylyi"]], all_fits[["maxOD"]][["Acinetobacter baylyi"]])
bayes_R2(all_fits_noint[["maxOD"]][["Acinetobacter baylyi"]])
bayes_R2(all_fits[["maxOD"]][["Acinetobacter baylyi"]])

#empty working environment to keep everything clean:
rm.all.but("globsets")



# Step 4: Fitting models to all species pairs (GxG analysis) ----
# ___________________________________________________________________
# input=  combined_parms_absence_antibiotics.csv
# output= brms_fit_GxG_models.RData, containing six brms fits

growth_parms <- read_csv(here("output", "growth_parameters_easylinear_h20.csv")) |>
  mutate(rpoB = fct_relevel(factor(rpoB), "WT"))

species_pairs <- data.frame(Pair = c("EC_AB", "EC_LP", "AB_LP"),
                            Species1 = c(rep("Escherichia coli", 2), "Acinetobacter baylyi"),
                            Species2 = c("Acinetobacter baylyi", rep("Lactiplantibacillus plantarum", 2)),
                            SetTemperature1 = c("37", "37", "30"),
                            SetTemperature2 = c("30", "30", "30"))

all_fits <- c()

for (i in 1:nrow(species_pairs)) {
  # filter out non-overlapping mutations:
  overlap <- intersect(growth_parms |> filter(Species == species_pairs$Species1[i]) |> pull(rpoB),
                       growth_parms |> filter(Species == species_pairs$Species2[i]) |> pull(rpoB))
  dat <- growth_parms |>
    filter(rpoB %in% overlap) |>
    filter(((Species == species_pairs$Species1[i]) & (SetTemperature == species_pairs$SetTemperature1[i])) |
           ((Species == species_pairs$Species2[i]) & (SetTemperature == species_pairs$SetTemperature2[i]))) |>
    select(Species, rpoB, mumax, maxOD) |>
    mutate(Species = fct_relevel(factor(Species), species_pairs$Species1[i]))
  all_fits[["fits"]][["mumax"]][[species_pairs$Pair[i]]] <- brm(mumax ~ rpoB * Species, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
  all_fits[["fits"]][["maxOD"]][[species_pairs$Pair[i]]] <- brm(maxOD ~ rpoB * Species, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
  all_fits[["fits_noint"]][["mumax"]][[species_pairs$Pair[i]]] <- brm(mumax ~ rpoB + Species, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
  all_fits[["fits_noint"]][["maxOD"]][[species_pairs$Pair[i]]] <- brm(maxOD ~ rpoB + Species, data = dat) |>
    add_criterion(c("loo", "bayes_R2"), reloo = TRUE)
}

all_fits[["legend"]] <- species_pairs

save(all_fits, file = here("output", "brms_fit_GxG_models.RData"))

# comparing the simpler vs. the interactions models:
load(file = here("output", "brms_fit_GxG_models.RData"))
Pair = c("EC_AB", "EC_LP", "AB_LP")

R2_int_mumax <- list()
R2_int_maxOD <- list()
R2_noint_mumax <- list()
R2_noint_maxOD <- list()

compare_loo_maxOD <- list()
compare_loo_mumax <- list()

for (i in Pair) {
  R2_int_mumax[[i]] <- bayes_R2(all_fits[["fits"]][["mumax"]][[i]])
  R2_int_maxOD[[i]] <- bayes_R2(all_fits[["fits"]][["maxOD"]][[i]])
  R2_noint_mumax[[i]] <- bayes_R2(all_fits[["fits_noint"]][["mumax"]][[i]])
  R2_noint_maxOD[[i]] <- bayes_R2(all_fits[["fits_noint"]][["maxOD"]][[i]])

  compare_loo_mumax[[i]] <- loo_compare(all_fits[["fits"]][["mumax"]][[i]], all_fits[["fits_noint"]][["mumax"]][[i]])
  compare_loo_maxOD[[i]] <- loo_compare(all_fits[["fits"]][["maxOD"]][[i]], all_fits[["fits_noint"]][["maxOD"]][[i]])
}

# E. coli vs A. baylyi
# mu_max 
compare_loo_mumax[["EC_AB"]]
R2_int_mumax[["EC_AB"]]
R2_noint_mumax[["EC_AB"]]
# max_OD 
compare_loo_maxOD[["EC_AB"]]
R2_int_maxOD[["EC_AB"]]
R2_noint_maxOD[["EC_AB"]]

# E. coli vs L. plantarum
# mu_max 
compare_loo_mumax[["EC_LP"]]
R2_int_mumax[["EC_LP"]] 
R2_noint_mumax[["EC_LP"]]
# max_OD 
compare_loo_maxOD[["EC_LP"]]
R2_int_maxOD[["EC_LP"]]
R2_noint_maxOD[["EC_LP"]]

# A. baylyi vs L. plantarum
# mu_max 
compare_loo_mumax[["AB_LP"]]
R2_int_mumax[["AB_LP"]]
R2_noint_mumax[["AB_LP"]]
# max_OD 
compare_loo_maxOD[["AB_LP"]]
R2_int_maxOD[["AB_LP"]]
R2_noint_maxOD[["AB_LP"]]

#empty working environment to keep everything clean:
rm.all.but("globsets")




# Step 5: Fitting single model to E. coli & A. baylyi data (GxGxE analysis) ----
# __________________________________________________________________________________
# input=  combined_parms_absence_antibiotics.csv
# output= brms_fit_GxGxE_model.RData containing two brms fits

growth_parms <-  read_csv(here("output", "growth_parameters_easylinear_h20.csv"))

#find overlapping mutations 
overlapping <- growth_parms %>% 
  filter(Species != "Lactiplantibacillus plantarum") %>% 
  distinct(Species, rpoB, Strain) %>% 
  pivot_wider(names_from = Species, values_from = Strain) %>% 
  drop_na("Acinetobacter baylyi") %>% 
  drop_na("Escherichia coli")

overlap <- unique(overlapping$rpoB)

GxE_mod <- growth_parms %>% 
  filter(rpoB %in% overlap) %>% 
  filter(Species != "Lactiplantibacillus plantarum") |>
  mutate(rpoB = fct_relevel(factor(rpoB), "WT")) |>
  mutate(SetTemperature = fct_relevel(factor(SetTemperature), "37")) |>
  mutate(Species = fct_relevel(factor(Species), "Escherichia coli"))


# model fit without 3-way interactions:
maxOD_fit_no3way <- brm(maxOD ~ rpoB*Species + 
             SetTemperature*Species + 
             rpoB*SetTemperature, 
             data = GxE_mod) |>
  add_criterion("loo", reloo = TRUE) |>
  add_criterion("bayes_R2")
summary(maxOD_fit_no3way)
#plot(maxOD_fit_no3way)

mumax_fit_no3way <- brm(mumax ~ rpoB*Species + 
                                  SetTemperature*Species + 
                                  rpoB*SetTemperature, 
                        data = GxE_mod) |>
  add_criterion("loo", reloo = TRUE) |>
  add_criterion("bayes_R2")
summary(mumax_fit_no3way)
#plot(mumax_fit_no3way)

# model fit with 3-way interactions:
maxOD_fit <- brm(maxOD ~ rpoB*Species*SetTemperature, 
                 data = GxE_mod) |>
  add_criterion("loo", reloo = TRUE) |>
  add_criterion("bayes_R2")
summary(maxOD_fit)
#plot(maxOD_fit)

mumax_fit <- brm(mumax ~ rpoB*Species*SetTemperature, 
                 data = GxE_mod) |>
  add_criterion("loo", reloo = TRUE) |>
  add_criterion("bayes_R2")
summary(mumax_fit)
#plot(mumax_fit)


# saving fits
all_fits <- list(mumax = mumax_fit_no3way,
                 maxOD = maxOD_fit_no3way)
save(all_fits, file = here("output", "brms_fit_GxGxE_model_no3way.RData"))

all_fits <- list(mumax = mumax_fit,
                 maxOD = maxOD_fit)
save(all_fits, file = here("output", "brms_fit_GxGxE_model.RData"))

# comparing the simpler vs. the 3-way model:
all_fits_noint <- get(load(file = here("output", "brms_fit_GxGxE_model_no3way.RData")))
all_fits <- get(load(file = here("output", "brms_fit_GxGxE_model.RData")))

# mu_max

#pp_check(mumax_fit_no3way)
#pp_check(mumax_fit)

#bridge_no3way <- bridge_sampler(mumax_fit_no3way)
#bridge <- bridge_sampler(mumax_fit)
#bayes_factor(bridge, bridge_no3way)

loo_compare(all_fits_noint[["mumax"]], all_fits[["mumax"]])
bayes_R2(all_fits_noint[["mumax"]])
bayes_R2(all_fits[["mumax"]])

# max OD

loo_compare(all_fits_noint[["maxOD"]], all_fits[["maxOD"]])
bayes_R2(all_fits_noint[["maxOD"]])
bayes_R2(all_fits[["maxOD"]])


# Step 6: Model comparisons ----
# ___________________________________________________________________
# input=  combined_parms_absence_antibiotics.csv
# output= brms_fit_GxG_models.RData, containing six brms fits

lo <- (1 - globsets$cred_int) / 2
hi <- 1 - lo

# load all fit objects:
ge    <- get(load(here("output", "brms_fit_GxE_models.RData")))          # interaction
gen   <- get(load(here("output", "brms_fit_GplusE_models.RData")))       # additive
gg    <- get(load(here("output", "brms_fit_GxG_models.RData")))          # $fits / $fits_noint / $legend
gxe3  <- get(load(here("output", "brms_fit_GxGxE_model.RData")))         # $mumax / $maxOD  (3-way)
gxe3n <- get(load(here("output", "brms_fit_GxGxE_model_no3way.RData")))  # $mumax / $maxOD  (no 3-way)

# build the raw comparison table:
comparison_table <- bind_rows(
  
  # GxE: rpoB + Temp  vs  rpoB * Temp, per species
  map_dfr(names(ge$mumax), \(sp) bind_rows(
    extract_comparison(ge$mumax[[sp]], gen$mumax[[sp]], "GxE", "mumax", sp),
    extract_comparison(ge$maxOD[[sp]], gen$maxOD[[sp]], "GxE", "maxOD", sp)
  )),
  
  # GxG: rpoB + Species  vs  rpoB * Species, per species pair
  map_dfr(gg$legend$Pair, \(pr) bind_rows(
    extract_comparison(gg$fits$mumax[[pr]], gg$fits_noint$mumax[[pr]], "GxG", "mumax", pr),
    extract_comparison(gg$fits$maxOD[[pr]], gg$fits_noint$maxOD[[pr]], "GxG", "maxOD", pr)
  )),
  
  # GxGxE: all 2-way  vs  full 3-way
  extract_comparison(gxe3$mumax, gxe3n$mumax, "GxGxE", "mumax", "EC_AB"),
  extract_comparison(gxe3$maxOD, gxe3n$maxOD, "GxGxE", "maxOD", "EC_AB")
)

# tidy labels:
pretty_contrast <- c(
  "Escherichia coli"     = "E. coli",
  "Acinetobacter baylyi" = "A. baylyi",
  "EC_AB" = "E. coli / A. baylyi",
  "EC_LP" = "E. coli / L. plantarum",
  "AB_LP" = "A. baylyi / L. plantarum"
)

comparison_table <- comparison_table |>
  mutate(contrast = coalesce(pretty_contrast[contrast], contrast))

write_csv(comparison_table, here("output", "model_comparison.csv"))

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
M2_noepi <- W2 + M1 - W1

arrow_offset <- 0.02
arrow_length <- unit(0.02, "npc")
label_a <- "italic(a)"
label_b2 <- "italic(b)[2]"
label_bM <- "italic(b)[M]"
label_bM2 <- "italic(b)[M2]"

# plot A (traditional fitness landscape):
data_plot <- tibble(Genotype = c(expression(italic(W)[1]), 
                                 expression(italic(W)[2]), 
                                 expression(italic(M)[1]),
                                 expression(italic(M)[2]),
                                 expression(italic(M)[2])),
                    Steps = c(0, 1, 1, 2, 2),
                    Trait = c(W1, W2, M1, M2, M2_noepi))

plot_theme <-   theme_light() +
  theme(axis.text.x = element_text(size = 4 * 2.83),
        axis.text.y = element_text(size = 4 * 2.83), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


pA <- ggplot(data_plot) +
  geom_line(aes(Steps, Trait), data = data_plot[c(2, 5),], col = "grey") +
  geom_line(aes(Steps, Trait), data = data_plot[c(3, 5),], col = "grey") +
  geom_point(aes(Steps, Trait), size = 4, col = c(rep("black", 4), "grey")) +
  geom_line(aes(Steps, Trait), data = data_plot[c(-3, -5),]) +
  geom_line(aes(Steps, Trait), data = data_plot[c(-2, -5),]) +
  scale_x_continuous(labels = NULL, breaks = c(0, 1, 2), limits = c(0, 2)) +
  scale_y_continuous(labels = NULL, limits = c(0, 1.2), expand = c(0, 0)) +
  annotate("text", 
           x = data_plot$Steps, 
           y = data_plot$Trait + 0.06, 
           label = data_plot$Genotype, 
           parse = TRUE, 
           size = 5,
           col = c(rep("black", 4), "grey")) +
  labs(x = "Genotype", y = "Trait") +
  plot_theme

pA_legend <- ggplot(data_plot) +
  scale_x_continuous(labels = NULL, limits = c(0.9, 5)) +
  scale_y_continuous(labels = NULL, limits = c(0, 1.2), expand = c(0, 0)) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = W1),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_segment(aes(x = 2, y = W1, xend = 2, yend = W2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_segment(aes(x = 3, y = W1, xend = 3, yend = M1),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  geom_segment(aes(x = 4, y = W2 - W1 + M1, xend = 4, yend = M2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  annotate("text", x = 1, y = 1.06, label = label_a, parse = TRUE, size = 3) +
  annotate("text", x = 2, y = 1.06, label = label_b2, parse = TRUE, size = 3) +
  annotate("text", x = 3, y = 1.06, label = label_bM, parse = TRUE, size = 3) +
  annotate("text", x = 4, y = 1.06, label = label_bM2, parse = TRUE, size = 3) +
  theme_void()

# plot B (our own visualisation):
pB <- plot_epistasis_bg(W1, W2, 0, 1.2, 0, 1.2) +
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
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p <- (pA + labs(tag = "A")) + pA_legend + (pB + labs(tag = "B")) +
  plot_layout(widths = c(1, 0.21, 1))

ggsave(here("plots", "Fig1_explainer.png"), plot = p, height = 5, width = 9.6)

#empty working environment to keep everything clean:
rm.all.but("globsets")


# Fig 2: Epistasis scatterplots               ----
# __________________________________________________________
# input=  Bayesian relative parameters
# output= Fig2_GxG_speciespairfits.png

set.seed(globsets$random_seed)

load(here("output", "brms_fit_GxG_models.RData"))
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

fig <- (p[["mumax"]][["EC_AB"]] + p[["mumax"]][["EC_LP"]] + p[["mumax"]][["AB_LP"]])/
  (p[["maxOD"]][["EC_AB"]] + p[["maxOD"]][["EC_LP"]] + p[["maxOD"]][["AB_LP"]])
fig <- fig + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14),
        legend.position = "bottom")

ggsave(fig, filename = here("plots", "Fig2_GxG_speciespairfits.png"),  width = 15, height = 10)

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Fig 3: GxE scatterplots ----
# __________________________________________________________
# input=  Bayesian model fits
# output= GxE figure

set.seed(globsets$random_seed)

load(here("output", "brms_fit_GxE_models.RData"))
variables <- names(all_fits)
species <- names(all_fits[["mumax"]])
p <- c()

for(i in 1:length(variables)) {
  if (variables[i] == "mumax") {
    xlab <- expression(mu[max] * "(37"* degree *"C)")
    ylab <- expression(mu[max] *"(30"* degree *"C)")
  } else if (variables[i] == "maxOD") {
    xlab <- expression("maxOD (37"* degree *"C)")
    ylab <- expression("maxOD (30"* degree *"C)")
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
      filter(SetTemperature == "37") |>
      select(rpoB, .epred, .lower, .upper) |>
      rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)
    
    post_pred_2 <- post_preds |>
      filter(SetTemperature == "30") |>
      select(rpoB, .epred, .lower, .upper) |>
      rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)
    
    p[[variables[i]]][[species[j]]] <- plot_epistasis_abspars(post_pred_1, 
                                                              post_pred_2,
                                                              xlab, ylab, leglab = "GxE")
  }
}

fig <- (p[["mumax"]][["Escherichia coli"]] + p[["mumax"]][["Acinetobacter baylyi"]])/
  (p[["maxOD"]][["Escherichia coli"]] + p[["maxOD"]][["Acinetobacter baylyi"]])
fig <- fig + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14),
        legend.position = "bottom")

ggsave(fig, filename = here("plots", "Fig3_GxE_speciesfits.png"),  width = 10, height = 10)

#empty working environment to keep everything clean:
rm.all.but("globsets")


# Fig 4: posterior distributions for single fit ----
# __________________________________________________________
# input=  Bayesian relative parameters
# output= GxE graphs 

# fits without 3-way interactions:
load(here("output", "brms_fit_GxGxE_model_no3way.RData"))

p_mumax <- plot_posteriors(all_fits$mumax)
p_maxOD <- plot_posteriors(all_fits$maxOD)

p <- (p_mumax | p_maxOD) + 
  plot_annotation(tag_levels = 'A')

ggsave(here("plots", "GxGxE_no3way.pdf"), p, width = 12, height = 8)


# fits with 3-way interactions:
load(here("output", "brms_fit_GxGxE_model.RData"))

p_mumax <- plot_posteriors(all_fits$mumax)
p_maxOD <- plot_posteriors(all_fits$maxOD)

p <- (p_mumax | p_maxOD) + 
  plot_annotation(tag_levels = 'A')

ggsave(here("plots", "Fig4_GxGxE.pdf"), p, width = 10, height = 8)

#empty working environment to keep everything clean:
rm.all.but("globsets")

# SupFig 1: fitness boxplots for all species and temperatures ----
# ________________________________________________________________
# input=  growth_parameters_h20.csv, parameters from MIC assay
# output= boxplots of species and mutants' fitness

set.seed(globsets$random_seed)

parms <- read_csv(here("output", "growth_parameters_easylinear_h20.csv")) |>
  mutate(SetTemperature = fct_relevel(factor(SetTemperature), "30")) |>
  mutate(rpoB = fct_reorder(factor(rpoB), AA_pos, .na_rm = FALSE)) |>
  mutate(rpoB = fct_relevel(rpoB, "WT")) |>
  select(Species, SetTemperature, rpoB, mumax, maxOD)

# select only mutations found in >=2 species:
overlapping_muts <- parms |>
  distinct(Species, rpoB) |>
  count(rpoB) |>
  filter(n >= 2) |>
  pull(rpoB)

parms <- filter(parms, rpoB %in% overlapping_muts)

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
    coord_cartesian(ylim=c(0, 1.15))
  
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

ggsave(SupFig1, filename = here("plots", "SupFig1_boxplots.pdf"), width=15 , height= 8)

#empty working environment to keep everything clean:
rm.all.but("globsets")


# Table of all mutants                                  ----
# __________________________________________________________

# identify mutations found in >=2 species:
overlapping_muts <- read_csv(here("output", "growth_parameters_easylinear_h20.csv")) |>
  distinct(Species, rpoB) |>
  count(rpoB) |>
  filter(n >= 2) |>
  pull(rpoB)

read_csv(here("output", "growth_parameters_easylinear_h20.csv")) |>
  filter(rpoB %in% overlapping_muts) |>
  distinct(Species, rpoB, AA_mut_name, AA_mut_name_Ecoli)  |>
  mutate(Species = abbreviate_species_name(Species)) |>
  filter(rpoB != "WT") |>
  mutate(AA_mut_name = str_replace(AA_mut_name, "rpoB_", "")) |>
  mutate(AA_pos = str_extract(AA_mut_name, "\\d+")) |>
  mutate(Position = str_extract(AA_mut_name_Ecoli, "\\d+")) |>
  mutate(Mutation = str_replace(AA_mut_name, "^([A-Za-z])\\d+([A-Za-z])$", "\\1 → \\2")) |>
  select(-AA_mut_name) |>
  pivot_wider(names_from = Species, values_from = AA_pos) |>
  select(Position, Mutation, "E. coli", "A. baylyi", "L. plantarum") |>
  arrange(Position) |>
  xtable(caption = "Rifampicin resistance mutants used in this study.
         The first columns show the amino acid position in \\textit{E. coli} coordinates and the amino acid substitution. 
         The remaining columns show the mutations used in this study, with numbers indicating the amino acid positions in the three focal species. 
         Mutations shown in the same row are homologues.",
         label = "tab:mutations") |>
  print(type = "latex", file = here("tables", "mutations.tex"), include.rownames = FALSE, floating = TRUE,
        table.placement = "H",
        size = "\\small\\renewcommand{\\arraystretch}{0.8}")

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Tables for model fits                                 ----
# __________________________________________________________

all_tables <- NULL

# GxG models:
load(here("output", "brms_fit_GxG_models.RData"))
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

# GxE models:
load(here("output", "brms_fit_GxE_models.RData"))
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

# GxGxE models:
load(here("output", "brms_fit_GxGxE_model.RData"))
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

load(here("output", "brms_fit_GxGxE_model_no3way.RData"))
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
  writeLines(here("tables", "model_fit_tables.tex"))

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Table for model comparisons                          ----
# __________________________________________________________

tab <- read_csv(here("output", "model_comparison.csv")) 
tab |>
  mutate(Model = paste0(response, " TILDE ", analysis)) |>
  select(-analysis, -response) |>
  select(Model, everything()) |>
  rename(Contrast = "contrast", Favoured = "favoured") |>
  xtable(caption = "Comparisons between model fits including interaction terms (e.g., maxOD $\\sim$ GxE) and the corresponding models without interaction terms (e.g., maxOD $\\sim$ G+E). Columns indicate the model and response variable for each species or species pair, the difference in the Expected Log Pointwise Predictive Density (ELPD) between the standard errors (se), which model is favoured through leave-one-out cross validation, as well as the $R^2$ values for both models along with their 95\\% credible intervals.",
         digits = rep(3, 12),
         label = "supptab:model_comparisons") |>
  print(type = "latex", file = "", 
        size = "\\footnotesize",
        include.rownames = FALSE) |>
  capture.output() |>
  str_replace_all("elpd\\\\_diff", "$\\\\Delta$ELPD") |>
  str_replace_all("se\\\\_diff", "$\\\\Delta$se") |>
  str_replace_all("R2\\\\_additive\\\\_lo", "Q2.5") |>
  str_replace_all("R2\\\\_additive\\\\_hi", "Q97.5") |>
  str_replace_all("R2\\\\_interaction\\\\_lo", "Q2.5") |>
  str_replace_all("R2\\\\_interaction\\\\_hi", "Q97.5") |>
  str_replace_all("R2\\\\_additive", "$R\\^2_{\\\\text{additive}}$") |>
  str_replace_all("R2\\\\_interaction", "$R\\^2_{\\\\text{interaction}}$") |>
  str_replace_all("TILDE", "$\\\\sim$") |>
  writeLines(here("tables", "comparison_table.tex"))

#empty working environment to keep everything clean:
rm.all.but("globsets")

# Session info                                          ----
# __________________________________________________________
# record R and package versions used to produce all outputs:

writeLines(capture.output(sessionInfo()), here("output", "sessionInfo.txt"))
