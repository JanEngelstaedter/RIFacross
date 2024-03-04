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
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.01, 1)) +
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
