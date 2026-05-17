#' @importFrom utils head tail
#' @importFrom rlang .data
NULL

# Global variables used in tidy evaluation and dplyr/data-masking contexts.
# This avoids R CMD check notes for column names and special tidy-eval syntax.

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".data",
    ":=",
    "scenario_time",
    "scenario_id",
    "cashflow_id",
    "q_total",
    "p_total",
    "d_total"
  ))
}
