
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
