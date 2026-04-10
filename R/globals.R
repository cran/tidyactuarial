# R/globals.R

#' @importFrom rlang .data :=
#' @importFrom utils head tail
NULL

# Evitar notas de R CMD check por variables usadas en NSE (dplyr/ggplot)
utils::globalVariables(c(
  "scenario_time",
  "scenario_id",
  "cashflow_id"
  # Si te salen mas en el futuro, agregalas a esta lista
))
