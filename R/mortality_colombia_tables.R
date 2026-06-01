#' Colombian mortality tables
#'
#' @description
#' A tidy collection of Colombian mortality tables for life-contingency
#' examples. The dataset includes regulatory and pedagogical mortality tables
#' used in Colombian actuarial applications.
#'
#' @format A tibble with 12 variables:
#' \describe{
#'   \item{table_id}{Mortality table identifier.}
#'   \item{sex}{Sex category, typically \code{"male"} or \code{"female"}.}
#'   \item{x}{Integer actuarial age.}
#'   \item{lx}{Number of survivors at exact age \code{x}.}
#'   \item{dx}{Expected number of deaths between ages \code{x} and
#'   \code{x + 1}.}
#'   \item{qx}{One-year death probability between ages \code{x} and
#'   \code{x + 1}.}
#'   \item{px}{One-year survival probability between ages \code{x} and
#'   \code{x + 1}.}
#'   \item{mu_x}{Force of mortality at age \code{x}, when available.}
#'   \item{ex}{Complete life expectancy at age \code{x}, when available.}
#'   \item{source}{Source identifier or reference.}
#'   \item{qx_calc}{Death probability recalculated from \code{lx} and
#'   \code{dx}, when available.}
#'   \item{qx_diff}{Difference between reported \code{qx} and recalculated
#'   \code{qx}, when available.}
#' }
#'
#' @details
#' The dataset is intended for actuarial examples involving Colombian mortality
#' tables, survival probabilities, life annuities, life insurance present values,
#' and validation of life-table calculations.
#'
#' The variable names follow the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{x} denotes age, \code{lx} the number of lives,
#' \code{dx} the number of deaths, \code{qx} the one-year death probability,
#' \code{px} the one-year survival probability, \code{mu_x} the force of
#' mortality, and \code{ex} the life expectancy at age \code{x}.
#'
#' The variables \code{qx_calc} and \code{qx_diff} are included as validation
#' aids. They allow users to compare reported death probabilities against
#' probabilities reconstructed from \code{lx} and \code{dx}.
#'
#' Some tables may start at different initial ages, use different terminal ages,
#' or use different radix values. Users should filter the desired table and sex
#' before passing the data to life-contingency functions.
#'
#' @source Colombian mortality tables cleaned for tidyactuarial examples. Source
#' identifiers are provided in the \code{source} column.
#'
#' @examples
#' data(mortality_colombia_tables)
#'
#' head(mortality_colombia_tables)
#'
#' mortality_colombia_tables |>
#'   dplyr::count(table_id, sex)
#'
#' rv08_male <- mortality_colombia_tables |>
#'   dplyr::filter(table_id == "RV08_Rentistas_2005_2008", sex == "male") |>
#'   dplyr::select(x, lx, dx, qx, px, mu_x, ex)
#'
#' head(rv08_male)
#'
#' @docType data
#' @keywords datasets
"mortality_colombia_tables"
