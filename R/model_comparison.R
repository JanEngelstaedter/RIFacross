# Model comparison table: additive vs. interaction models
# ______________________________________________________________________
# Reads the saved brms fits (with cached `loo` and `bayes_R2` criteria)
# and builds one manuscript-ready table of LOO comparisons + Bayes R2.
#
# Assumes working directory = project root, and that Steps 3-5 of
# master.R have been run so the following files exist in output/:
#   brms_fit_GxE_models.RData          (all_fits       - GxE, interaction)
#   brms_fit_GplusE_models.RData       (all_fits_noint - GxE, additive)
#   brms_fit_GxG_models.RData          (all_fits       - $fits / $fits_noint / $legend)
#   brms_fit_GxGxE_model.RData         (all_fits       - $mumax / $maxOD, 3-way)
#   brms_fit_GxGxE_model_no3way.RData  (all_fits       - $mumax / $maxOD, no 3-way)

library(tidyverse)
library(brms)
library(loo)
library(xtable)

cred_int <- 0.95
lo <- (1 - cred_int) / 2
hi <- 1 - lo

# If a fit has no usable cached loo, recompute with reloo (slow) instead
# of a plain loo() that may return NA pointwise contributions.
RELOO_FALLBACK <- TRUE

# ---- helpers --------------------------------------------------------

get_R2_draws <- function(fit) {
  r2 <- fit$criteria$bayes_R2
  if (is.null(r2)) r2 <- bayes_R2(fit, summary = FALSE)
  as.numeric(r2)
}

get_loo <- function(fit) {
  l <- fit$criteria$loo
  if (is.null(l) || !all(is.finite(l$estimates[, "Estimate"]))) {
    message("  recomputing loo with reloo = TRUE ...")
    l <- loo(fit, reloo = TRUE)
  }
  l
}

extract_comparison <- function(fit_full, fit_reduced, analysis, response, contrast) {
  L_full <- get_loo(fit_full)
  L_red  <- get_loo(fit_reduced)

  e_full <- L_full$pointwise[, "elpd_loo"]
  e_red  <- L_red$pointwise[, "elpd_loo"]
  stopifnot(length(e_full) == length(e_red))

  d         <- e_full - e_red                # pointwise: interaction - additive
  elpd_diff <- sum(d)                        # = elpd_loo(full) - elpd_loo(reduced)
  se_diff   <- sqrt(length(d) * var(d))      # identical to loo_compare's se_diff

  if (!is.finite(elpd_diff))
    warning("non-finite elpd_diff for ", analysis, " / ", response, " / ", contrast,
            call. = FALSE)

  r2f <- get_R2_draws(fit_full)
  r2r <- get_R2_draws(fit_reduced)

  tibble(
    analysis, response, contrast,
    elpd_diff, se_diff,
    favoured = case_when(!is.finite(elpd_diff) ~ NA_character_,
                         elpd_diff > 0         ~ "interaction",
                         TRUE                  ~ "additive"),
    R2_additive       = mean(r2r),
    R2_additive_lo    = unname(quantile(r2r, lo)),
    R2_additive_hi    = unname(quantile(r2r, hi)),
    R2_interaction    = mean(r2f),
    R2_interaction_lo = unname(quantile(r2f, lo)),
    R2_interaction_hi = unname(quantile(r2f, hi))
  )
}


# ---- load all fit objects ------------------------------------------

ge    <- get(load("output/brms_fit_GxE_models.RData"))          # interaction
gen   <- get(load("output/brms_fit_GplusE_models.RData"))       # additive
gg    <- get(load("output/brms_fit_GxG_models.RData"))          # $fits / $fits_noint / $legend
gxe3  <- get(load("output/brms_fit_GxGxE_model.RData"))         # $mumax / $maxOD  (3-way)
gxe3n <- get(load("output/brms_fit_GxGxE_model_no3way.RData"))  # $mumax / $maxOD  (no 3-way)

# ---- build the raw comparison table -------------------------------

comparison_table <- bind_rows(

  # Step 3 - GxE: rpoB + Temp  vs  rpoB * Temp, per species
  map_dfr(names(ge$mumax), \(sp) bind_rows(
    extract_comparison(ge$mumax[[sp]], gen$mumax[[sp]], "GxE", "mumax", sp),
    extract_comparison(ge$maxOD[[sp]], gen$maxOD[[sp]], "GxE", "maxOD", sp)
  )),

  # Step 4 - GxG: rpoB + Species  vs  rpoB * Species, per species pair
  map_dfr(gg$legend$Pair, \(pr) bind_rows(
    extract_comparison(gg$fits$mumax[[pr]], gg$fits_noint$mumax[[pr]], "GxG", "mumax", pr),
    extract_comparison(gg$fits$maxOD[[pr]], gg$fits_noint$maxOD[[pr]], "GxG", "maxOD", pr)
  )),

  # Step 5 - GxGxE: all 2-way  vs  full 3-way
  extract_comparison(gxe3$mumax, gxe3n$mumax, "GxGxE", "mumax", "EC_AB"),
  extract_comparison(gxe3$maxOD, gxe3n$maxOD, "GxGxE", "maxOD", "EC_AB")
)

# ---- tidy labels -------------------------------------------------

pretty_contrast <- c(
  "Escherichia coli"     = "E. coli",
  "Acinetobacter baylyi" = "A. baylyi",
  "EC_AB" = "E. coli / A. baylyi",
  "EC_LP" = "E. coli / L. plantarum",
  "AB_LP" = "A. baylyi / L. plantarum"
)

comparison_table <- comparison_table |>
  mutate(contrast = coalesce(pretty_contrast[contrast], contrast))

comparison_table   # inspect in console

# raw values (optional):
write_csv(comparison_table, "output/model_comparison.csv")