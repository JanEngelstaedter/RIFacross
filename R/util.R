
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
