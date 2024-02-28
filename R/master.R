library(tidyverse)
library(grow96)
library(Biostrings)
library(ALJEbinf)

data_Maddie <- processODData(specPath = "./data/A_baylyi_Maddie/specs/",
                             dataPath = "./data/A_baylyi_Maddie/OD_raw") |>
  blankODs(groups = "RIFconc")

qcODData(data_Maddie, path = "./output/QC/A_baylyi_Maddie/", silent = FALSE)

ggplot(data = filter(data_Maddie) +
  geom_point(mapping = aes(x = Time_h, y = blankedOD), size = 1) +
  scale_x_continuous(breaks = c(0, 10, 20)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_grid(rows = vars(Strain), cols = vars(round(RIFconc, 3))) +
  labs(x = "Time (hours)", y = "Optical Densiy (A600)") +
  theme_bw()
