#' Simulate claim counts from a frequency model
#'
#' Draws independent realizations of the claim frequency random variable N
#' using the RNG registered for its distribution. A draw represents the
#' number of claims on one exposure or observation unit, not the portfolio
#' total unless the model itself is specified at portfolio level.
#'
#' @param x A [frequency_spec()] or [empirical_frequency()] specification.
#' @param n_sim Positive integer number of simulations.
#' @param seed Optional nonnegative integer seed. When supplied, the caller's
#'   random-number state is restored after the call.
#' @return A tibble with `simulation` and `claim_count` columns.
#' @details A missing RNG provider is an error, not a reason to silently
#'   approximate the model. For an empirical frequency, the draws follow the
#'   observed mass function exactly.
#' @examples
#' simulate_frequency(frequency_spec("poisson", lambda = 3), n_sim = 5, seed = 10)
#' @export
simulate_frequency <- function(x, n_sim = 1000L, seed = NULL) {
  .collective_validate_frequency(x)
  n_sim <- .risk_check_n_sim(n_sim)
  .risk_with_seed(seed, {
    tibble::tibble(
      simulation = seq_len(n_sim),
      claim_count = .risk_draw(x, n_sim, "frequency")
    )
  })
}
