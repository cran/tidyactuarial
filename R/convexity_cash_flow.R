#' Convexity of a general cash-flow stream
#'
#' Computes the yield convexity of a general cash-flow stream at a
#' specified valuation time.
#'
#' @param cf Numeric vector of cash-flow amounts.
#' @param t Numeric vector of cash-flow times. Times must be expressed in
#'   the same units used by `rate`.
#' @param rate Effective yield per unit of time. It must be a finite numeric
#'   scalar greater than `-1`.
#' @param valuation_time Finite numeric scalar indicating the valuation time.
#'   Defaults to `0`.
#' @param output Character scalar indicating the returned object:
#'   `"value"` returns convexity as a numeric scalar;
#'   `"audit"` returns a detailed tibble.
#'
#' @return
#' If `output = "value"`, a numeric scalar.
#'
#' If `output = "audit"`, a tibble with six columns:
#' `time`, `remaining_time`, `cash_flow`, `discount_factor`,
#' `present_value`, and `convexity_contribution`.
#'
#' @details
#' Let \eqn{\tau} be the valuation time. Under an effective yield \eqn{i},
#' the value of the remaining cash-flow stream is
#'
#' \deqn{
#' P_{\tau}(i)
#' =
#' \sum_{j:t_j>\tau}
#' C_j(1+i)^{-(t_j-\tau)}.
#' }
#'
#' The convexity returned by this function is the relative second
#' derivative of value with respect to the effective yield:
#'
#' \deqn{
#' \mathcal{C}_{\tau}
#' =
#' \frac{1}{P_{\tau}(i)}
#' \frac{d^2P_{\tau}(i)}{di^2}.
#' }
#'
#' Equivalently,
#'
#' \deqn{
#' \mathcal{C}_{\tau}
#' =
#' \frac{
#' \sum_{j:t_j>\tau}
#' (t_j-\tau)(t_j-\tau+1)
#' C_j(1+i)^{-(t_j-\tau+2)}
#' }{
#' P_{\tau}(i)
#' }.
#' }
#'
#' Together with modified duration, this convexity gives the second-order
#' approximation
#'
#' \deqn{
#' \frac{\Delta P_{\tau}}{P_{\tau}}
#' \approx
#' -D_{\mathrm{mod},\tau}\Delta i
#' +
#' \frac{1}{2}\mathcal{C}_{\tau}(\Delta i)^2.
#' }
#'
#' The function follows an ex-cash-flow convention: cash flows occurring
#' exactly at `valuation_time` are excluded. Therefore, the result represents
#' convexity immediately after any cash flow paid at that time.
#'
#' For streams containing both positive and negative cash flows, convexity
#' remains a yield-sensitivity measure, but it may be negative and may not
#' have the usual interpretation associated with conventional bonds.
#'
#' @examples
#' convexity_cash_flow(
#'   cf = c(100, 100, 1100),
#'   t = c(1, 2, 3),
#'   rate = 0.05
#' )
#'
#' convexity_cash_flow(
#'   cf = c(100, 100, 1100),
#'   t = c(1, 2, 3),
#'   rate = 0.05,
#'   valuation_time = 1
#' )
#'
#' convexity_cash_flow(
#'   cf = c(100, 100, 1100),
#'   t = c(1, 2, 3),
#'   rate = 0.05,
#'   output = "audit"
#' )
#'
#' @seealso [duration_cash_flow()]
#'
#' @export
convexity_cash_flow <- function(
  cf,
  t,
  rate,
  valuation_time = 0,
  output = c("value", "audit")
) {

  output <- match.arg(output)

  is_numeric_vector <- function(x) {
    is.numeric(x) && is.null(dim(x))
  }

  if (!is_numeric_vector(cf) || !is_numeric_vector(t)) {
    stop(
      "`cf` and `t` must be numeric vectors.",
      call. = FALSE
    )
  }

  if (length(cf) != length(t)) {
    stop(
      "`cf` and `t` must have the same length.",
      call. = FALSE
    )
  }

  if (length(cf) == 0L) {
    stop(
      "`cf` and `t` cannot be empty.",
      call. = FALSE
    )
  }

  if (any(!is.finite(cf)) || any(!is.finite(t))) {
    stop(
      "`cf` and `t` must contain only finite values.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(rate) ||
    length(rate) != 1L ||
    !is.finite(rate)
  ) {
    stop(
      "`rate` must be a finite numeric scalar.",
      call. = FALSE
    )
  }

  if (rate <= -1) {
    stop(
      "`rate` must be greater than -1.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(valuation_time) ||
    length(valuation_time) != 1L ||
    !is.finite(valuation_time)
  ) {
    stop(
      "`valuation_time` must be a finite numeric scalar.",
      call. = FALSE
    )
  }

  future_index <- t > valuation_time

  if (!any(future_index)) {
    stop(
      paste0(
        "No cash flows occur after `valuation_time` = ",
        valuation_time,
        "."
      ),
      call. = FALSE
    )
  }

  future_time <- t[future_index]
  future_cash_flow <- cf[future_index]

  order_index <- order(
    future_time,
    seq_along(future_time)
  )

  future_time <- future_time[order_index]
  future_cash_flow <- future_cash_flow[order_index]

  remaining_time <- future_time - valuation_time

  discount_factor <- exp(
    -remaining_time * log1p(rate)
  )

  present_value <- future_cash_flow *
    discount_factor

  if (
    any(!is.finite(discount_factor)) ||
    any(!is.finite(present_value))
  ) {
    stop(
      paste0(
        "Discounting produced non-finite values. Check `rate`, `t`, ",
        "and `valuation_time`."
      ),
      call. = FALSE
    )
  }

  total_present_value <- sum(present_value)
  present_value_scale <- sum(abs(present_value))

  if (
    !is.finite(total_present_value) ||
    !is.finite(present_value_scale)
  ) {
    stop(
      paste0(
        "Aggregating the present values produced a non-finite result. ",
        "Check the magnitude of `cf`, `rate`, and `t`."
      ),
      call. = FALSE
    )
  }

  near_zero_present_value <-
    present_value_scale == 0 ||
    abs(total_present_value) <=
      sqrt(.Machine$double.eps) *
        present_value_scale

  if (near_zero_present_value) {
    stop(
      paste0(
        "Convexity is undefined because the present value of the ",
        "remaining cash flows is zero or numerically close to zero."
      ),
      call. = FALSE
    )
  }

  convexity_weight <- (
    remaining_time *
      (remaining_time + 1)
  ) / (1 + rate)^2

  if (any(!is.finite(convexity_weight))) {
    stop(
      paste0(
        "The convexity weights are non-finite. Check the magnitude ",
        "of `rate`, `t`, and `valuation_time`."
      ),
      call. = FALSE
    )
  }

  convexity_numerator_component <- present_value *
    convexity_weight

  if (any(!is.finite(convexity_numerator_component))) {
    stop(
      paste0(
        "The convexity numerator is non-finite. Check the magnitude ",
        "of `cf`, `rate`, `t`, and `valuation_time`."
      ),
      call. = FALSE
    )
  }

  convexity_numerator <- sum(
    convexity_numerator_component
  )

  if (!is.finite(convexity_numerator)) {
    stop(
      paste0(
        "Aggregating the convexity contributions produced a ",
        "non-finite result. Check the magnitude of `cf`, `rate`, ",
        "`t`, and `valuation_time`."
      ),
      call. = FALSE
    )
  }

  convexity_value <- convexity_numerator /
    total_present_value

  if (!is.finite(convexity_value)) {
    stop(
      paste0(
        "Convexity is non-finite. Check the supplied cash flows ",
        "and valuation assumptions."
      ),
      call. = FALSE
    )
  }

  if (output == "value") {
    return(as.numeric(convexity_value))
  }

  convexity_contribution <-
    convexity_numerator_component /
      total_present_value

  tibble::tibble(
    time = future_time,
    remaining_time = remaining_time,
    cash_flow = future_cash_flow,
    discount_factor = discount_factor,
    present_value = present_value,
    convexity_contribution = convexity_contribution
  )
}
