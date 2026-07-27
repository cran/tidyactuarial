#' Duration of a general cash-flow stream
#'
#' Computes the Macaulay or modified duration of a general cash-flow
#' stream at a specified valuation time.
#'
#' @param cf Numeric vector of cash-flow amounts.
#' @param t Numeric vector of cash-flow times. Times must be expressed in
#'   the same units used by `rate`.
#' @param rate Effective yield per unit of time. It must be a finite numeric
#'   scalar greater than `-1`.
#' @param valuation_time Finite numeric scalar indicating the valuation time.
#'   Defaults to `0`.
#' @param duration Character scalar indicating the duration measure:
#'   `"macaulay"` or `"modified"`.
#' @param output Character scalar indicating the returned object:
#'   `"value"` returns the selected duration as a numeric scalar;
#'   `"audit"` returns a detailed tibble.
#'
#' @return
#' If `output = "value"`, a numeric scalar.
#'
#' If `output = "audit"`, a tibble with six columns:
#' `time`, `remaining_time`, `cash_flow`, `discount_factor`,
#' `present_value`, and `duration_contribution`.
#'
#' @details
#' Let \eqn{\tau} be the valuation time. Under an effective yield \eqn{i},
#' the value of the remaining cash-flow stream is
#'
#' \deqn{
#' P_{\tau}
#' =
#' \sum_{j:t_j>\tau}
#' C_j(1+i)^{-(t_j-\tau)}.
#' }
#'
#' The Macaulay duration measured from \eqn{\tau} is
#'
#' \deqn{
#' D_{\mathrm{Mac},\tau}
#' =
#' \frac{
#' \sum_{j:t_j>\tau}
#' (t_j-\tau)C_j(1+i)^{-(t_j-\tau)}
#' }{
#' P_{\tau}
#' }.
#' }
#'
#' The modified duration with respect to the effective yield is
#'
#' \deqn{
#' D_{\mathrm{mod},\tau}
#' =
#' \frac{D_{\mathrm{Mac},\tau}}{1+i}.
#' }
#'
#' The function follows an ex-cash-flow convention: cash flows occurring
#' exactly at `valuation_time` are excluded. Therefore, the result represents
#' duration immediately after any cash flow paid at that time.
#'
#' For streams containing both positive and negative cash flows, duration
#' remains a yield-sensitivity measure, but it may not have the usual
#' interpretation as a weighted-average payment time.
#'
#' @examples
#' duration_cash_flow(
#'   cf = c(100, 100, 1100),
#'   t = c(1, 2, 3),
#'   rate = 0.05
#' )
#'
#' duration_cash_flow(
#'   cf = c(100, 100, 1100),
#'   t = c(1, 2, 3),
#'   rate = 0.05,
#'   duration = "modified"
#' )
#'
#' duration_cash_flow(
#'   cf = c(100, 100, 1100),
#'   t = c(1, 2, 3),
#'   rate = 0.05,
#'   valuation_time = 1,
#'   output = "audit"
#' )
#'
#' @export
duration_cash_flow <- function(
  cf,
  t,
  rate,
  valuation_time = 0,
  duration = c("macaulay", "modified"),
  output = c("value", "audit")
) {

  duration <- match.arg(duration)
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
        "Duration is undefined because the present value of the ",
        "remaining cash flows is zero or numerically close to zero."
      ),
      call. = FALSE
    )
  }

  weighted_present_value <- remaining_time *
    present_value

  if (any(!is.finite(weighted_present_value))) {
    stop(
      paste0(
        "The duration numerator is non-finite. Check the magnitude ",
        "of `cf`, `t`, and `valuation_time`."
      ),
      call. = FALSE
    )
  }

  duration_numerator <- sum(weighted_present_value)

  if (!is.finite(duration_numerator)) {
    stop(
      paste0(
        "Aggregating the duration contributions produced a non-finite ",
        "result. Check the magnitude of `cf`, `t`, and `valuation_time`."
      ),
      call. = FALSE
    )
  }

  macaulay_duration <- duration_numerator /
    total_present_value

  modified_duration <- macaulay_duration /
    (1 + rate)

  duration_value <- switch(
    duration,
    macaulay = macaulay_duration,
    modified = modified_duration
  )

  if (!is.finite(duration_value)) {
    stop(
      paste0(
        "The selected duration is non-finite. Check the supplied ",
        "cash flows and valuation assumptions."
      ),
      call. = FALSE
    )
  }

  if (output == "value") {
    return(as.numeric(duration_value))
  }

  duration_contribution <- switch(
    duration,
    macaulay = weighted_present_value /
      total_present_value,
    modified = weighted_present_value /
      ((1 + rate) * total_present_value)
  )

  tibble::tibble(
    time = future_time,
    remaining_time = remaining_time,
    cash_flow = future_cash_flow,
    discount_factor = discount_factor,
    present_value = present_value,
    duration_contribution = duration_contribution
  )
}
