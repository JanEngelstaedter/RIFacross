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
  fit_for_plot <- brms::posterior_summary(fit, probs = c(0.025, 0.975), central = FALSE) |>
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


#' An area plot showing different types of epistasis
#'
#' @param W1 Trait of Wt in environment or background 1
#' @param W2 Trait of Wt in environment or background 2
#' @param min1 Lower range on x axis
#' @param max1 Upper range on x axis
#' @param min2 Lower range on y axis
#' @param max2 Upper range on y axis
#' @param scale Scaling of axes. Options are "ID" (default, no scaling),
#' "normalised" (everything scaled multiplicatively to WT values), and
#' "lognormalised" (scaled multiplicatively to WT values, and then log2 applied).
#' @param n Resolution (number of tiles) on each axis
#'
#' @returns A ggplot object
#' 
plot_epistasis_bg <- function(W1, W2, min1, max1, min2, max2, scale = "ID", n = 201, leglab = "Epistasis") {
  epi_types <- c("magnitude", "sign", "reciprocal")
  if (scale == "ID") {
    M1 <- seq(min1, max1, length.out = n)
    M2 <- seq(min2, max2, length.out = n)
    
    epi <- expand.grid(M1 = M1, M2 = M2) |>
      rowwise() |>
      mutate(E = epistasis(W1, W2, M1, M2)) |>
      mutate(type = epistasis_type(W1, W2, M1, M2)) |>
      mutate(sign = if_else(E<0, "-", if_else(E>0, "+", ""))) |>
      mutate(type = factor(type, levels = epi_types)) |>
      add_dummy_levels("type", epi_types)
    
    zeroE_line <- data.frame(M1 = M1) |>
      rowwise() |>
      mutate(M2 = zeroE_AB(W1, W2, M1))
    
    p <- ggplot(epi) +
      geom_tile(aes(x = M1, y = M2, fill = type), alpha = 0.3) +
      geom_line(data = zeroE_line, mapping = aes(x = M1, y = M2), col = "#3B9AB2") 
  } else if (scale == "normalised") {
    M1_rel <- seq(min1, max1, length.out = n)
    M2_rel <- seq(min2, max2, length.out = n)
    
    epi <- expand.grid(M1_rel = M1_rel, M2_rel = M2_rel) |>
      rowwise() |>
      mutate(M1 = M1_rel * W1, M2 = M2_rel * W2) |>
      mutate(E = epistasis(W1, W2, M1, M2)) |>
      mutate(type = epistasis_type(W1, W2, M1, M2)) |>
      mutate(sign = if_else(E<0, "-", if_else(E>0, "+", ""))) |>
      mutate(type = factor(type, levels = epi_types)) |>
      add_dummy_levels("type", epi_types)
    
    zeroE_line <- data.frame(M1_rel = M1_rel) |>
      rowwise() |>
      mutate(M1 = M1_rel * W1) |>
      mutate(M2 = zeroE_AB(W1, W2, M1)) |>
      mutate(M2_rel = M2/W2)
    
    p <- ggplot(epi) +
      geom_tile(aes(x = M1_rel, y = M2_rel, fill = type), alpha = 0.3) +
      geom_line(data = zeroE_line, mapping = aes(x = M1_rel, y = M2_rel), col = "#3B9AB2")
  } else if (scale == "lognormalised") {
    log2_M1_rel <- seq(min1, max1, length.out = n)
    log2_M2_rel <- seq(min2, max2, length.out = n)
    
    epi <- expand.grid(log2_M1_rel = log2_M1_rel, log2_M2_rel = log2_M2_rel) |>
      rowwise() |>
      mutate(M1 = 2^log2_M1_rel * W1, M2 = 2^log2_M2_rel * W2) |>
      mutate(E = epistasis(W1, W2, M1, M2)) |>
      mutate(type = epistasis_type(W1, W2, M1, M2)) |>
      mutate(sign = if_else(E<0, "-", if_else(E>0, "+", ""))) |>
      mutate(type = factor(type, levels = epi_types)) |>
      add_dummy_levels("type", epi_types)
    
    zeroE_line <- data.frame(log2_M1_rel = log2_M1_rel) |>
      rowwise() |>
      mutate(M1 = 2^log2_M1_rel * W1) |>
      mutate(M2 = zeroE_AB(W1, W2, M1)) |>
      mutate(log2_M2_rel = log2(M2/W2))
    
    p <- ggplot() +
      geom_tile(data = epi, mapping = aes(x = log2_M1_rel, y = log2_M2_rel, fill = type), alpha = 0.3) +
      geom_line(data = zeroE_line, mapping = aes(x = log2_M1_rel, y = log2_M2_rel), col = "#3B9AB2")
  } else {
    stop("Unsupported value of scale argument.")
  }
  
  cols <- c("magnitude" = "#3B9AB2", "sign" = "#EBCC2A", "reciprocal" = "#F21A00")
  p <- p + 
    scale_x_continuous(limits = c(min1, max1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min2, max2), expand = c(0, 0)) +
    scale_fill_manual(values = cols, drop = FALSE) +
    theme_light() +
    guides(fill=guide_legend(title = leglab))
  return(p)
}

plot_epistasis_relpars <- function(rel_parms, fit1, fit2, min1, max1, min2, max2, xlab, ylab, leglab = "Epistasis") {
  W1 <- posterior_summary(fit1, probs = c(0.025, 0.975), central = FALSE)["b_Intercept", "Estimate"]
  W2 <- posterior_summary(fit2, probs = c(0.025, 0.975), central = FALSE)["b_Intercept", "Estimate"]
  
  p <- plot_epistasis_bg(W1, W2, min1, max1, min2, max2, scale = "normalised", leglab = leglab)
  p <- p +
    geom_point(aes(x = x, y = y), data = rel_parms) +
    labs(#title = expression(italic("E.coli")), 
      x = xlab, y = ylab) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#F21A00") +
    geom_errorbarh(aes(y = y, xmin = xmin, xmax = xmax),
                   data = rel_parms,
                   alpha = 0.3) +
    geom_errorbar(aes(x = x, ymin = ymin, ymax = ymax),
                  data = rel_parms,
                  alpha = 0.3) +
    geom_label_repel(aes(x = x, y = y, label = rpoB),
                     data = rel_parms,
                     size = 4,
                     fill = alpha("white", 0), 
                     label.size = NA,
                     segment.color = NA) +
    theme(
      plot.title = element_text(size = 18),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10), 
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  return(p)
}
