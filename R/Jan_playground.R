
# Comparing global fit with individual fits to species x T conditions ----

# load data and use better names for fit objects:
load("output/brms_fit_single_model_mumax_threeway.RData")
load("output/brms_fit_separate_models_mumax.RData")

fit_mumax_ind <- fit_mumax # individual fits
fit_mumax_glob <- mumax_fit # global fit

rm(fit_mumax, mumax_fit)

# load estimated parameters:
parms_all <- read_csv("./output/combined_parms_absence_antibiotics.csv")

# example: mumax for E. coli mutation 526D at 37C

comparison <- tibble(WT_ind = as_draws_df(fit_mumax_ind$EC37) %>% pull(b_Intercept)) |>
  mutate(H526D_ind = WT_ind + as_draws_df(fit_mumax_ind$EC37) %>% pull(b_rpoBH526D)) |>
  mutate(WT_glob = as_draws_df(fit_mumax_glob) %>% pull(b_Intercept)) |>
  mutate(H526D_glob = WT_glob + as_draws_df(fit_mumax_glob) %>% pull(b_rpoBH526D))

mumax_estimates_H526D <- parms_all |>
  filter(Species == "Escherichia coli" & as.integer(SetTemperature) == 37L & rpoB == "H526D") |>
  pull(mumax)

mumax_estimates_WT <- parms_all |>
  filter(Species == "Escherichia coli" & as.integer(SetTemperature) == 37L & rpoB == "WT") |>
  pull(mumax)

ggplot(comparison) +
  geom_density(aes(WT_ind), col = "blue") +
  geom_density(aes(WT_glob), col = "blue", linetype = "dashed") +
  geom_vline(xintercept = mean(mumax_estimates_WT, na.rm = TRUE), col = "blue") +
  geom_density(aes(H526D_ind), col = "red") +
  geom_density(aes(H526D_glob), col = "red", linetype = "dashed") +
  geom_vline(xintercept = mean(mumax_estimates_H526D), col = "red") +
  theme_minimal()

GxE_mod |>
  filter(Species == "Escherichia coli" & SetTemperature == "30" & rpoB == "H526Y") |>
  pull(maxOD) |>
  mean(na.rm = TRUE)

GxE_mod |>
  filter(Species == "Escherichia coli" & SetTemperature == "37" & rpoB == "WT") |>
  pull(mumax) |>
  mean(na.rm = TRUE)

plot_posteriors(fit_mumax_glob)

# Fig 1_alt: GxE scatterplots but on a non-normalised scale ----
# __________________________________________________________
#input=  Bayesian relative parameters
#output= GxE graphs 

# preparing data:
set.seed(globsets$random_seed)

load(file = "output/brms_fit_separate_models_mumax.RData")
conditions <- read_csv("./data/conditions_key.csv")

as_draws_df(fit_maxOD[["EC37"]])$b_Intercept |> var()

abs_parms <- fit_mumax |>
  map(brms::posterior_summary, probs = c(0.025, 0.975), central = FALSE) |>
  map(as_tibble, rownames = "Coefficient") |>
  bind_rows(.id = "Condition") |>
  left_join(conditions, by = join_by(Condition)) |>
  filter(str_detect(Coefficient, "b_")) |>
  mutate(rpoB = if_else(Coefficient == "b_Intercept", "WT", str_replace(Coefficient, "b_rpoB", "")))

intercepts <- abs_parms |>
  filter(Coefficient == "b_intercept") |>
  select()



rel_parms <- read_csv("./output/bayesian_relative_parameters.csv" ) |>
  select(-Condition) |>
  pivot_wider(names_from = Temperature, 
              values_from = c(relative_maxOD_y, relative_maxOD_ymin, relative_maxOD_ymax,
                              relative_mumax_y, relative_mumax_ymin, relative_mumax_ymax))

# mumax plots:

xlab <- expression(mu[max] * "(30"* degree *"C)")
ylab <- expression(mu[max] *"(37"* degree *"C)")
p_data <- rel_parms |>
  mutate(x = relative_mumax_y_30, xmin = relative_mumax_ymin_30, xmax = relative_mumax_ymax_30,
         y = relative_mumax_y_37, ymin = relative_mumax_ymin_37, ymax = relative_mumax_ymax_37)

p_mumax_EC <- plot_epistasis_relpars(p_data |> filter(Species == "Escherichia coli"), 
                                     fit_mumax[["EC30"]], fit_mumax[["EC37"]],
                                     0.7, 1.15, 0.7, 1.15, xlab, ylab)
p_mumax_AB <- plot_epistasis_relpars(p_data |> filter(Species == "Acinetobacter baylyi"), 
                                     fit_mumax[["AB30"]], fit_mumax[["AB37"]],
                                     0.65, 1.1, -0.05, 1.1, xlab, ylab)

# maxOD plots:
load(file = "output/brms_fit_separate_models_maxOD.RData")
xlab <- expression("maxOD (30"* degree *"C)")
ylab <- expression("maxOD (37"* degree *"C)")
p_data <- rel_parms |>
  mutate(x = relative_maxOD_y_30, xmin = relative_maxOD_ymin_30, xmax = relative_maxOD_ymax_30,
         y = relative_maxOD_y_37, ymin = relative_maxOD_ymin_37, ymax = relative_maxOD_ymax_37)

p_maxOD_EC <- plot_epistasis_relpars(p_data |> filter(Species == "Escherichia coli"), 
                                     fit_maxOD[["EC30"]], fit_maxOD[["EC37"]],
                                     0.7, 1.4, 0.7, 1.05, xlab, ylab)
p_maxOD_AB <- plot_epistasis_relpars(p_data |> filter(Species == "Acinetobacter baylyi"), 
                                     fit_maxOD[["AB30"]], fit_maxOD[["AB37"]],
                                     0.7, 1.4, -0.05, 1.3, xlab, ylab)

Fig1 <- (p_mumax_EC + p_mumax_AB )/(p_maxOD_EC + p_maxOD_AB)
Fig1 <- Fig1 + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14),
        legend.position = "bottom")

ggsave(Fig1, filename = "./plots/Fig1_GxE.pdf",  width = 10, height = 10)

#empty working environment to keep everything clean:
rm.all.but("globsets")



# turning global fit into local fit coefficients

fit <- mumax_fit
newdata <- expand_grid(Species = levels(fit$data$Species),
                         SetTemperature = levels(fit$data$SetTemperature),
                         rpoB = levels(fit$data$rpoB))

post_preds <- fit |>
  epred_draws(newdata = newdata) |>
  group_by(across(all_of(colnames(newdata)))) |>
  mean_hdi(.epred, .width = 0.95)

post_pred_1 <- post_preds |>
  filter(Species == "Escherichia coli", SetTemperature == "30") |>
  select(rpoB, .epred, .lower, .upper) |>
  rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)

post_pred_2 <- post_preds |>
  filter(Species == "Escherichia coli", SetTemperature == "37") |>
  select(rpoB, .epred, .lower, .upper) |>
  rename(mean = .epred, hdi_lower = .lower, hdi_upper = .upper)

xlab <- expression(mu[max] * "(30"* degree *"C)")
ylab <- expression(mu[max] *"(37"* degree *"C)")
p_mumax_EC <- plot_epistasis_abspars(post_pred_1, 
                                     post_pred_2,
                                     0.0132, 0.0192, 0.0212, 0.0305, xlab, ylab, leglab = "GxE")
# --> different fit, not surprising?

########################################################################
# Fit model to individual species (--> GxE plots)
########################################################################

# fitting:
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

# plotting:
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
      mean_hdi(.epred, .width = 0.95)
    
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

ggsave(Fig1, filename = "./plots/Fig1_GxE_speciesfits.png",  width = 10, height = 10)


########################################################################
# Fit model to species pairs (--> Epistasis plots)
########################################################################

# fitting:
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

# plotting:
load("output/brms_fit_GxG_models.RData")
variables <- names(all_fits$fits)
pairs <- names(all_fits$fits[["mumax"]])
species_pairs <- all_fits$legend
p <- c()

for(i in 1:length(variables)) {
  for(j in 1:length(all_fits$fits[[i]])) {
    if (variables[i] == "mumax") {
      xlab <- expression(mu[max] * eval(paste0("(", species_pairs$Species1[j], ")")))
      ylab <- expression(mu[max] * eval(paste0("(", species_pairs$Species2[j], ")")))
    } else if (variables[i] == "maxOD") {
      xlab <- paste0("maxOD (", species_pairs$Species1[j], ")")
      ylab <- paste0("maxOD (", species_pairs$Species2[j], ")")
    }
    
    fit <- all_fits$fits[[variables[i]]][[species_pairs$Pair[j]]]
    newdata <- expand_grid(Species = levels(fit$data$Species),
                           rpoB = levels(fit$data$rpoB))
    post_preds <- fit |>
      epred_draws(newdata = newdata) |>
      group_by(across(all_of(colnames(newdata)))) |>
      mean_hdi(.epred, .width = 0.95)
    
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
                                                                         xlab, ylab, leglab = "Epistasis")
  }
}

Fig3 <- (p[["mumax"]][["EC_AB"]] + p[["mumax"]][["EC_LP"]] + p[["mumax"]][["AB_LP"]])/
  (p[["maxOD"]][["EC_AB"]] + p[["maxOD"]][["EC_LP"]] + p[["maxOD"]][["AB_LP"]])
Fig3 <- Fig3 + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14),
        legend.position = "bottom")

ggsave(Fig3, filename = "./plots/Fig3_GxG_speciespairfits.png",  width = 15, height = 10)
