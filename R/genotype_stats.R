#!/usr/bin/env Rscript
#from Christopher 17/5/24

library(SuperExactTest)

n <- 385 #Thanks to Negin.

#total-mutations-ALJE.csv is muts from master 

dat <- read.csv("total-mutations-ALJE.csv")

dat <- dat |>
	dplyr::filter(Gene == "rpoB") |>
	dplyr::mutate(AA_mut_name_Ecoli = ifelse(is.na(AA_mut_name_Ecoli), Nt_mut_name, AA_mut_name_Ecoli))

mutation_lists <- lapply(unique(dat$Species), function(sp) {
	dplyr::filter(dat, Species == sp) |> dplyr::pull(AA_mut_name_Ecoli)
})
names(mutation_lists) <- unique(dat$Species)

intersections <- SuperExactTest::supertest(mutation_lists, n = n, degree = 2:length(mutation_lists))

print(summary(intersections))
