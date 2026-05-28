#' SOA Illustrative Life Table
#'
#' A tidy version of the Society of Actuaries illustrative life table commonly
#' used in life-contingencies examples and benchmark calculations.
#'
#' @description
#' This dataset is intended for reproducible examples, internal validation, and
#' benchmark tests for life-contingency functions such as [annuity_x()],
#' [insurance_x()], [premium_x()], [reserve_x()], [annuity_xy()],
#' [insurance_xy()], [premium_xy()], and [reserve_xy()].
#'
#' @format A data frame with one row per integer age and the following columns:
#' \describe{
#'   \item{x}{Integer age.}
#'   \item{lx}{Number of lives surviving to exact age `x`.}
#'   \item{dx}{Number of deaths between exact ages `x` and `x + 1`.}
#'   \item{qx}{One-year probability of death between exact ages `x` and `x + 1`.}
#'   \item{px}{One-year probability of survival from exact age `x` to `x + 1`.}
#' }
#'
#' @details
#' The table is included as a convenient benchmark table for actuarial
#' calculations. The build script in `data-raw/soa08lt.R` constructs this tidy
#' dataset from the SOA illustrative actuarial table object distributed in the
#' `lifecontingencies` package.
#'
#' @source Society of Actuaries illustrative life table, commonly referenced in
#' Bowers et al. (1997), *Actuarial Mathematics*, Appendix 2A.
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and Nesbitt,
#' C. J. (1997). *Actuarial Mathematics*. Second edition. Society of Actuaries.
#'
#' Spedicato, G. A. (2013). The `lifecontingencies` package: performing
#' financial and actuarial mathematics calculations in R.
#'
#' @examples
#' data(soa08lt)
#' head(soa08lt)
#'
#' annuity_x(
#'   mortality_table = soa08lt,
#'   age = 65,
#'   rate = 0.06,
#'   timing = "due"
#' )
#'
#' @keywords datasets
"soa08lt"
