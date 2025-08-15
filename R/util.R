
#' Estimation of minimum inhibitory concentrations (MICs)
#'
#' @param conc Vector of antibiotic concentrations
#' @param maxOD Corresponding vector of measured maximum optical densities (OD)
#' @param ICx Percentage by which maximum OD must be reduced to assign a 
#' population as non-growing for the purpose of MIC estimation
#'
#' @return A single number giving the MIC.
#' @export
#'
getMIC <- function(conc, maxOD, ICx = 90) {
  # calculate mean maxODs:
  result <- data.frame(conc = conc, maxOD = maxOD) |>
    group_by(conc) |>
    summarise(maxOD = mean(maxOD, na.rm = TRUE))
  
  # max OD at drug concentration zero:
  ODmax0 <- result |>
    filter(near(conc, 0)) |>
    pull(maxOD)
  
  # calculate MIC:
  result <- result |>  
    mutate(relativeOD = maxOD/ODmax0) |>
    filter(relativeOD < (1 - ICx/100)) |>
    pull(conc) |> 
    min(na.rm = TRUE)
  
  return(result)
}

# calculates additive epistasis:
epistasis <- function(ab, aB, Ab, AB) {
  return(ab + AB - Ab  - aB)
}

# Determines the type of epistasis (magnitude, sign or reciprocal sign).
# By default, no additivity (absence of epistasis) is reported as this is a special case.
epistasis_type <- function(ab, aB, Ab, AB,
                           report_additive = FALSE,
                           tol = 1e-2 * max(ab, aB, Ab, AB)) {
  if (report_additive && (abs(ab + AB - Ab - aB) < tol)) {
    E <- "additive"
  } else {
    # determine normalised ordering of values where,
    # after relabelling of the alleles, ab becomes the smallest:
    o <- rank(c(ab, aB, Ab, AB), ties.method = "first")
    if (o[1] == 1L) { # correct order
      n_o <- o # --> ab, aB, Ab, AB
    } else if (o[2] == 1L) {
      n_o <- c(o[2], o[1], o[4], o[3]) # --> aB, ab, AB, Ab
    } else if (o[3] == 1L) {
      n_o <- c(o[3], o[4], o[1], o[2]) # --> Ab, AB, ab, aB
    } else {
      n_o <- c(o[4], o[3], o[2], o[1]) # --> AB, aB, Ab, ab
    }
    if (n_o[4] == 2) {
      E <- "reciprocal"
    } else if (n_o[4] == 3) {
      E <- "sign"
    } else {
      E <- "magnitude"
    }
  }
  return(E)
}

# function to calculate AB trait that yields zero epistasis:
zeroE_AB <- function(ab, aB, Ab) {
  return(aB + Ab - ab)
}


# helper function (written by chatGPT) to add empty rows to a data frame so
# as to always plot all factor levels with the correct colour:
add_dummy_levels <- function(df, var, all_levels) {
  missing_levels <- setdiff(all_levels, unique(df[[var]]))
  if (length(missing_levels) == 0) return(df)
  
  dummy <- data.frame(matrix(NA, nrow = length(missing_levels), ncol = ncol(df)))
  colnames(dummy) <- colnames(df)
  dummy[[var]] <- factor(missing_levels, levels = all_levels)
  
  df[[var]] <- factor(df[[var]], levels = all_levels)
  rbind(df, dummy)
}


# helper function to (safely) shuffle the order in a vector:
shuffle <- function(x) {
  if (length(x) == 1L) {
    return(x) 
  } else {
    return(sample(x, length(x)))
  }
}


# function to prepare posterior predictions from model fit:
get_post_predictions <- function(fit) {
  post_pred <- fitted(fit, newdata = data.frame(rpoB = levels(fit$data$rpoB)), summary = FALSE) |>
    as_tibble(.name_repair = "minimal") %>%
    setNames(levels(fit$data$rpoB)) %>%
    pivot_longer(cols = everything(), names_to = "rpoB", values_to = "y_pred") |>
    group_by(rpoB) |>
    summarise(
      mean = mean(y_pred),
      hdi_lower = hdi(y_pred, prob = 0.95)[1],
      hdi_upper = hdi(y_pred, prob = 0.95)[2],
      .groups = "drop"
    )
  return(post_pred)
}


# function to determine plotting range from data:
range_with_padding <- function(x, mult = 0.05) {
  r <- range(x, na.rm = TRUE)
  padding <- diff(r) * mult
  r + c(-padding, padding)
}


posterior_table_caption <- function(model_type, variable, condition = NULL, formula = NULL) {
  cap <- paste0("Summary of the posterior distribution for the ",
                model_type,
                " model for ",
                variable,
                " as the reponse variable, and mutation")
  if (model_type == "GxE") {
    cap <- paste0(cap, " and temperature")
  } else if (model_type == "GxG") {
    cap <- paste0(cap, " and species") 
  } else if (model_type == "GxGxE") {
    cap <- paste0(cap, ", temperature and species")
  }
  cap <- paste0(cap, " as independent variables.")
  cap <- paste0(cap, " The columns show the mean, standard deviation, and credibile interval of the posterior distribution for each estimated parameter.")
  if (model_type == "GxE") {
    cap <- paste0(cap, " This model fit is for species ", condition, ". ")
  } else if (model_type == "GxG") {
    cap <- paste0(cap, " This model fit is for species pair ", condition, ". ")
  } else if (model_type == "GxGxE") {
    cap <- paste0(cap, " This model fit is for species pair \\textit{E.~coli} / \\textit{A.~baylyi}.")
  }
  cap <- paste0(cap, " Formula used for fit: ", as.character(formula[1]), ".")
  return(cap)
}

polish_posterior_tables <- function(tex) {
  tex_polished <- tex |>
    str_replace_all(" ~ ", " $\\\\sim$ ") |>
    str_replace_all("SetTemperature", "Temperature") |>
    str_replace_all("rpoB", "Mutation") |>
    str_replace_all("SpeciesAcinetobacterbaylyi", "SpeciesAbaylyi") |>
    str_replace_all("SpeciesLactiplantibacillusplantarum", "SpeciesLplantarum") |>
    str_replace_all(" EC_AB", " \\\\textit{E.~coli} / \\\\textit{A.~baylyi}") |>
    str_replace_all(" EC_LP", " \\\\textit{E.~coli} / \\\\textit{L.~plantarum}") |>
    str_replace_all(" AB_LP", " \\\\textit{A.~baylyi} / \\\\textit{L.~plantarum}") |>
    str_replace_all(" Acinetobacter baylyi", " \\\\textit{A.~baylyi}") |>
    str_replace_all(" Escherichia coli", " \\\\textit{E.~coli}") |>
    str_replace_all("_Acinetobacter baylyi", "_Abaylyi") |>
    str_replace_all("_Escherichia coli", "_Ecoli") |>
    str_replace_all("mumax", "$\\\\mu_{\\\\text{max}}$")
  return(tex_polished)
}

abbreviate_species_name <- function(s) {
  return(str_replace(s, "^(\\p{L})\\p{L}*", "\\1."))
}
