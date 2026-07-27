# Global variables used with tidy evaluation
#
# These names correspond to columns created and manipulated internally by
# tidyactuarial functions. Declaring them here prevents false-positive
# "no visible binding for global variable" notes during R CMD check.
#
# This file does not create objects in the package namespace and does not
# modify the behavior of any function.

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "cash_flow",
      "cashflow_id",
      "d_total",
      "duration_contribution",
      "p_total",
      "q_total",
      "remaining_time",
      "residual",
      "solution",
      "time",
      "value"
    )
  )
}
