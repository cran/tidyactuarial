#' Discount factor under conventional or time-varying interest
#'
#' Computes the discount factor from time \code{t} back to time \code{s}
#' using exactly one of three interest specifications:
#' \itemize{
#'   \item a conventional interest rate \code{i};
#'   \item an accumulation function \code{a(t)};
#'   \item a time-varying force of interest \code{delta(t)}.
#' }
#'
#' The discount factor is the reciprocal of the corresponding accumulation
#' factor:
#' \deqn{v(s,t) = \frac{1}{A(s,t)}.}
#'
#' Therefore,
#' \deqn{v(s,t) = (1+i)^{-(t-s)}}
#' for a constant annual effective rate,
#' \deqn{v(s,t) = \frac{a(s)}{a(t)}}
#' for an accumulation function, and
#' \deqn{v(s,t) =
#' \exp\left(-\int_s^t \delta(u)\,du\right)}
#' for a time-varying force of interest.
#'
#' @param s Numeric vector of valuation times. Defaults to 0.
#' @param t Numeric vector of payment or accumulation times.
#' @param i Optional numeric vector of conventional interest-rate values.
#' @param i_type Character vector indicating the conventional interest-rate
#'   type. Allowed values are \code{"effective"},
#'   \code{"nominal_interest"}, \code{"nominal_discount"}, and
#'   \code{"force"}. Used only when \code{i} is supplied.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Used only when \code{i} is supplied.
#' @param a Optional accumulation function of one numeric time argument.
#' @param delta Optional force-of-interest function of one numeric time
#'   argument.
#' @param subdivisions Positive integer giving the maximum number of
#'   subintervals used by \code{stats::integrate()} when \code{delta} is
#'   supplied.
#' @param rel.tol Positive relative tolerance passed to
#'   \code{stats::integrate()} when \code{delta} is supplied.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric vector. If
#'   \code{TRUE}, returns a tibble with the inputs, accumulation factor, and
#'   discount factor.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of discount factors.
#'
#' If \code{tidy = TRUE}, a tibble containing the interval, interest model,
#' relevant intermediate quantities, accumulation factor, and discount factor.
#'
#' @details
#' Exactly one of \code{i}, \code{a}, or \code{delta} must be supplied.
#'
#' Numeric arguments may have length 1 or a common length. Scalars are
#' recycled over the remaining scenarios, making the function suitable for
#' use inside \code{dplyr::mutate()} pipelines.
#'
#' This function handles general discounting over an interval \eqn{[s,t]}.
#' For discount factors obtained specifically from a spot-rate curve, see
#' \code{\link{discount_factor_spot}}.
#'
#' @seealso \code{\link{accumulation_factor}},
#'   \code{\link{discount_factor_spot}},
#'   \code{\link{present_value}}, \code{\link{pv_flow}}
#'
#' @family interest
#' @family time-value
#'
#' @examples
#' # Constant annual effective interest
#' discount_factor(
#'   t = 5,
#'   i = 0.07
#' )
#'
#' # One rate recycled over several times
#' discount_factor(
#'   t = c(1, 2, 5),
#'   i = 0.07
#' )
#'
#' # Accumulation function
#' a_fun <- function(t) {
#'   exp(0.03 * t + 0.002 * t^2)
#' }
#'
#' discount_factor(
#'   s = 2,
#'   t = 5,
#'   a = a_fun
#' )
#'
#' # Time-varying force of interest
#' delta_fun <- function(t) {
#'   0.03 + 0.004 * t
#' }
#'
#' discount_factor(
#'   s = 2,
#'   t = 5,
#'   delta = delta_fun
#' )
#'
#' # Pipe-friendly use
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   scenarios <- tibble::tibble(
#'     i = 0.07,
#'     t = c(1, 2, 3, 5)
#'   )
#'
#'   scenarios |>
#'     dplyr::mutate(
#'       discount = discount_factor(
#'         t = t,
#'         i = i
#'       )
#'     )
#' }
#'
#' @export
discount_factor <- function(
    s = 0,
    t,
    i = NULL,
    i_type = "effective",
    m = 1,
    a = NULL,
    delta = NULL,
    subdivisions = 100L,
    rel.tol = 1e-8,
    tidy = FALSE
) {
  accumulation_out <- accumulation_factor(
    s = s,
    t = t,
    i = i,
    i_type = i_type,
    m = m,
    a = a,
    delta = delta,
    subdivisions = subdivisions,
    rel.tol = rel.tol,
    tidy = tidy
  )

  if (!tidy) {
    return(1 / accumulation_out)
  }

  accumulation_out |>
    dplyr::mutate(
      discount_factor = 1 / .data[["accumulation_factor"]]
    )
}
