#' Tail Value at Risk (Expected Shortfall)
#'
#' Calculates TVaR_p(X) = [E(X 1(X > v)) + v(F(v) - p)] / (1-p), where
#' v = VaR_p(X). This correctly allocates only the necessary fraction of
#' a probability atom at v to the upper tail.
#'
#' @param x Finite nonnegative losses, a tibble of simulations, a supported
#'   severity specification, or an aggregate-distribution result.
#' @param level Probabilities strictly between zero and one.
#' @param value Optional loss column name for data frames.
#' @return A tibble with `level`, `tail_value_at_risk`, and `method`.
#' @details An aggregate distribution with unresolved probability mass cannot
#'   yield exact TVaR, even when its VaR is known: the missing loss
#'   expectation is not assumed to be zero. For continuous severity families
#'   with implemented stop-loss formulas, TVaR = VaR + E[(X-VaR)+]/(1-p).
#' @examples
#' tail_value_at_risk(c(0, 100, 100, 300), level = c(0.5, 0.9))
#' tail_value_at_risk(severity_spec("exponential", rate = 0.01), level = 0.95)
#' @export
tail_value_at_risk <- function(x, level = 0.95, value = NULL) {
  level <- .risk_validate_levels(level)
  if (inherits(x, "tidyact_severity_spec") &&
      !x$dist %in% c("discrete", "empirical", "degenerate")) {
    if (!is.null(value)) stop("`value` only selects a data-frame column.",
                              call. = FALSE)
    v <- .risk_exact_quantile(x, level)
    stop_loss <- stop_loss_premium(x, retention = v)
    tvar <- v + stop_loss$premium / (1 - level)
    method <- stop_loss$method[1L]
  } else {
    masses <- .risk_masses(x, value)
    if (masses$tail > 0) {
      stop("TVaR is unavailable: the aggregate distribution has an unresolved upper tail.",
           call. = FALSE)
    }
    v <- .risk_discrete_var(masses, level)
    tvar <- vapply(seq_along(level), function(k) {
      var <- v[k]
      p <- level[k]
      above <- masses$loss > var
      first_moment_above <- sum(masses$loss[above] *
                                  masses$probability[above])
      f_var <- sum(masses$probability[masses$loss <= var])
      (first_moment_above + var * max(0, f_var - p)) / (1 - p)
    }, numeric(1))
    method <- masses$method
  }
  tibble::tibble(level = level, tail_value_at_risk = tvar,
                 method = rep(method, length(level)))
}
