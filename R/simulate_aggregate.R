#' Simulate a collective risk model
#'
#' Samples independent realizations of S = X[1] + ... + X[N], drawing
#' frequencies first and then independent severities conditional on those
#' counts. A zero-frequency realization has aggregate loss zero.
#' The model assumes iid severities independent of frequency. Whenever
#' the moments exist, E[S] = E[N] E[X], and
#' Var(S) = E[N] Var(X) + Var(N) E[X]^2.
#'
#' @param x A [collective_risk_spec()] object. Zero-payment specifications
#'   produced by [to_payment_basis()] are also accepted.
#' @param n_sim Positive integer number of independent realizations.
#' @param seed Optional nonnegative integer random seed, restored after use.
#' @param max_claims Maximum total number of individual claim draws in one
#'   invocation. Exceeding this protective allocation limit is an error, not
#'   a truncation of the frequency distribution.
#' @return A tibble with `simulation`, `claim_count`, and `aggregate_loss`.
#' @details The resulting aggregate distribution follows the frequency and
#'   severity specifications on their explicitly declared count basis.
#'   All individual severities are iid and independent of frequency.
#' @examples
#' model <- collective_risk_spec(frequency_spec("poisson", lambda = 2),
#'                               severity_spec("degenerate", value = 100))
#' simulate_aggregate(model, n_sim = 5, seed = 30)
#' @export
simulate_aggregate <- function(x, n_sim = 1000L, seed = NULL,
                               max_claims = 1000000L) {
  n_sim <- .risk_check_n_sim(n_sim)
  .risk_check_seed(seed)
  max_claims <- .risk_check_max_claims(max_claims)
  if (inherits(x, "tidyact_zero_payment_basis_spec")) {
    return(tibble::tibble(simulation = seq_len(n_sim),
                          claim_count = rep(0, n_sim),
                          aggregate_loss = rep(0, n_sim)))
  }
  if (!inherits(x, "tidyact_collective_risk_spec")) {
    stop("`x` must be a collective risk specification.", call. = FALSE)
  }
  .collective_validate_frequency(x$frequency)
  .collective_validate_severity(x$severity)
  .risk_with_seed(seed, {
    counts <- .risk_draw(x$frequency, n_sim, "frequency")
    if (any(counts > max_claims)) {
      stop("Simulated claim count exceeds `max_claims`.", call. = FALSE)
    }
    total <- sum(counts)
    if (!is.finite(total) || total > max_claims ||
        total > .Machine$integer.max) {
      stop("Total simulated claims exceed `max_claims`; no counts were truncated.",
           call. = FALSE)
    }
    aggregate <- numeric(n_sim)
    if (total > 0) {
      positive <- which(counts > 0)
      # rep.int's times must be integer; already protected by max_claims.
      groups <- rep(positive, times = as.integer(counts[positive]))
      amounts <- .risk_draw(x$severity, as.integer(total), "severity")
      by_simulation <- rowsum(matrix(amounts, ncol = 1L), group = groups, reorder = FALSE)
      aggregate[positive] <- as.numeric(by_simulation[, 1L])
      if (any(!is.finite(aggregate))) {
        stop("Aggregate claim amounts overflowed during simulation.",
             call. = FALSE)
      }
    }
    tibble::tibble(simulation = seq_len(n_sim), claim_count = counts,
                   aggregate_loss = aggregate)
  })
}
