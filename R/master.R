library(tidyverse)
library(grow96)
library(Biostrings)
library(ALJEbinf)

# wrangling and qc:
data_Maddie <- processODData(specPath = "./data/A_baylyi_Maddie/specs/",
                             dataPath = "./data/A_baylyi_Maddie/OD_raw") |>
  blankODs(groups = "RIFconc")
qcODData(data_Maddie, path = "./output/QC/A_baylyi_Maddie/", silent = FALSE)

# plotting:
data_Maddie_plot <- data_Maddie |>
  mutate(RIFconc = factor(RIFconc, 
                          levels = RIFconc |> 
                            as.double() |> 
                            round(6) |>
                            sort() |>
                            as.character() |>
                            unique())) |>
  filter(WellType == "DATA", !is.na(Strain), !is.na(RIFconc))

p <- ggplot(data = data_Maddie_plot) +
  geom_line(mapping = aes(x = Time_h, y = blankedOD, group = interaction(Plate, Replicate))) +
  scale_x_continuous(breaks = c(0, 10, 20)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(-0.01, 1)) +
  facet_grid(rows = vars(Strain), cols = vars(RIFconc)) +
  labs(x = "Time (hours)", y = "Optical Densiy (A600)") +
  theme_bw()

ggsave("./plots/A_baylyi_Maddie_growthcurves.pdf", p, width = 18, height = 12)

data_Maddie |>
  filter(is.na(RIFconc))

data_Maddie_plot |>
  filter(is.na(RIFconc))

data_Maddie |>
  pull(RIFconc) |>
  unique() |>
  as.double() |> 
  round(6) |>
  sort() |>
  as.character()
