# facet plot of growth data for all strains and RIF concentrations in the data set:
plot_growthData <- function(dat, fileName = NULL) {
  dat_plot <- dat |>
    mutate(RIFconc = factor(RIFconc |> 
                              round(6) |> 
                              as.character(), 
                            levels = RIFconc |> 
                              round(6) |>
                              sort() |>
                              as.character() |>
                              unique())) |>
    filter(WellType == "DATA", !is.na(Strain), !is.na(RIFconc))
  
  dat_missing <- dat_plot |>
    select(Strain, RIFconc) |>
    distinct() |>
    tidyr::expand(Strain, RIFconc) |>
    anti_join(dat_plot, by = join_by(Strain, RIFconc))
  
  p <- ggplot(data = dat_plot) +
    geom_line(mapping = aes(x = Time_h, y = blankedOD, group = interaction(Plate, Replicate))) +
    geom_rect(data = dat_missing, 
               fill = "grey", colour = NA, xmin = -Inf, xmax = Inf,
               ymin = -Inf,ymax = Inf) +
    scale_x_continuous(breaks = c(0, 10, 20)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), 
                       limits = c(-0.01, round(max(dat_plot$blankedOD) + 0.05, 1))) +
    facet_grid(rows = vars(Strain), cols = vars(RIFconc)) +
    labs(x = "Time (hours)", y = "Optical Densiy (A600)") +
    theme_bw()
  
  if (!is.null(fileName)) {
    w <- 2 * (dat_plot |> pull(RIFconc) |> unique() |> length() + 1) / 3
    h <- 2 * (dat_plot |> pull(Strain) |> unique() |> length() + 1) / 3
    ggsave(fileName, p, width = w, height = h)
  }
  return(p)
}


# plots the posterior distribution of all parameters from a brms fit:
plot_posteriors <- function(fit) {
  fit_for_plot <- brms::posterior_summary(fit, probs = c(0.025, 0.975)) |>
    as.data.frame() |>
    rownames_to_column("Parameter") |>
    filter(grepl("^b_", Parameter)) |>
    mutate(Parameter = gsub("^b_", "", Parameter)) |>
    mutate(Parameter = gsub("SpeciesAcinetobacterbaylyi", "AB", Parameter)) |>
    mutate(Parameter = gsub("SetTemperature30", "T30", Parameter)) |>
    mutate(Parameter = gsub("rpoB.", "", Parameter)) |>
    filter(Parameter != "Intercept") |>
    mutate(Parameter_group = case_when(
      grepl("^[0-9]", Parameter) & grepl("AB$", Parameter) ~ "Mutation:Species",
      grepl("^[0-9]", Parameter) & grepl("T30$", Parameter) ~ "Mutation:Temperature",
      grepl("^[0-9]", Parameter) ~ "Mutation",
      .default = "Species & Temperature")) |>
    mutate(relation_to_zero = case_when(
      Q2.5 > 0 ~ "above",
      Q97.5 < 0 ~ "below",
      TRUE ~ "includes")) |>
    group_by(Parameter_group) |>
    arrange(Parameter) |>
    ungroup() |>
    arrange(Parameter_group) |>
    mutate(Parameter = as_factor(Parameter)) |>
    mutate(Parameter = fct_relevel(Parameter, "AB", "T30", "AB:T30")) |>
    mutate(Parameter_group = as_factor(Parameter_group)) |>
    mutate(Parameter_group = fct_relevel(Parameter_group, "Species & Temperature"))
  
  p <- ggplot(fit_for_plot, aes(x = Parameter, y = Estimate, color = relation_to_zero)) +
    geom_point() +
    geom_errorbar(aes(ymin = `Q2.5`, ymax = `Q97.5`), width = 0.2) +
    coord_flip() +
    scale_x_discrete(limits = rev)  +
    scale_color_manual(
      values = c("above" = "blue", "below" = "red", "includes" = "gray40"),
      labels = NULL
    ) +
    ggforce::facet_col(~ Parameter_group, scales = "free_y", space = "free") +
    #facet_wrap(~ Parameter_group, scales = "free_y", ncol = 1) +
    labs(y = "Estimate", x = "Parameter") +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(p)
}

