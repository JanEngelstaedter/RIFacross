# Fig 1_alt: GxE scatterplots but on a non-normalised scale ----
# __________________________________________________________
#input=  Bayesian relative parameters
#output= GxE graphs 

# preparing data:
set.seed(globsets$random_seed)

load(file = "output/brms_fit_separate_models_mumax.RData")
conditions <- read_csv("./data/conditions_key.csv")

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
