
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
