#' Example dataset: dat_JMP
#'
#' A dataset containing the longitudinal measurements of quality of life and survival time
#'
#' @format A data frame with 118 rows and 36 variables:
#' \describe{
#'   \item{id}{subject ID}
#'   \item{trt}{treatment status (1 for treated, 0 for control)}
#'   \item{sex}{sex}
#'   \item{death}{censoring status (1 for death, 0 for censored)}
#'   \item{survival_time}{survival time in months}
#'   \item{time_}{time for each longitudinal measurement of QOL}
#'   \item{qol_}{quality of life at different time point}
#'   ...
#' }

"dat_JMP"
