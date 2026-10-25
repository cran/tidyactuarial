#' Specify a change from per-loss to per-payment basis
#' @md
#'
#' @description
#' Records the exact actuarial change of basis for a collective risk model
#' whose severity is a per-loss payment specification from [apply_coverage()].
#' The conversion uses binomial thinning of the loss count and conditioning
#' the payment amount on a strictly positive insurer payment.
#'
#' @param x A list-backed collective risk specification containing `frequency`
#'   (a `tidyact_frequency_spec`) and `severity` (a
#'   `tidyact_payment_spec` with `basis = "per_loss"`). The collective class
#'   `tidyact_collective_risk_spec`, as provided by [collective_risk_spec()].
#'   The input must have `count_basis = "loss"`.
#' @param label Optional nonempty character scalar identifying the transformed
#'   model. If `NULL`, the collective model label is inherited when present.
#'
#' @details
#' Write `N_L` for the number of losses and `Y_L` for the insurer payment
#' **per loss**, including zero payments. Put
#' \deqn{q = P(Y_L > 0) = 1-F_{Y_L}(0).}
#' Assuming independent, identically distributed loss severities independent
#' of `N_L`, the number of payments has the binomial-thinning law
#' \deqn{N_P \mid N_L=n \sim Binomial(n,q).}
#' For `q > 0`, the **per-payment** severity is
#' \deqn{Y_P \stackrel{d}{=} (Y_L \mid Y_L>0)}
#' and its CDF is `0` for `y <= 0` and
#' \deqn{F_{Y_P}(y)=\frac{F_{Y_L}(y)-(1-q)}{q}, \quad y>0.}
#' Consequently, the probability generating functions satisfy
#' \deqn{G_{N_P}(z)=G_{N_L}(1-q+qz),}
#' and the aggregate payments have the same distribution on either basis.
#'
#' Poisson counts remain Poisson with mean `q * lambda`; binomial counts
#' remain binomial with unchanged number of trials and probability multiplied
#' by `q`. A negative binomial frequency with `size = r` and `mu` transforms
#' to the same `size` and mean `q * mu`. For mixed Poisson models,
#' the conditional intensity becomes `q * Lambda` **without** discarding
#' the mixing variable. Other frequency families must be thinned via the
#' general PGF rule rather than silently assumed to stay in the same family.
#'
#' When `q = 1`, both basis representations coincide. When `q = 0`,
#' payment frequency and aggregate payment are identically zero, while
#' `Y_L | Y_L > 0` is **undefined**. A fabricated point mass at zero must
#' not be used as a per-payment severity.
#'
#' The conversion is executable for registered Poisson, binomial, negative
#' binomial, geometric and finite empirical frequencies. No unknown/custom
#' frequency family is silently approximated. The probability of a positive
#' payment is resolved from the exact source CDF or finite discrete masses.
#' Continuous per-payment severities inherit only the registered exact
#' operations of their source; unavailable moment orders remain unavailable.
#' For negative binomial and geometric frequencies, thinning uses
#' \eqn{p_P=p/(p+q(1-p))}. If a positive payment frequency would round to
#' the degenerate parameter `prob = 1`, the conversion stops and the
#' per-loss model should be retained; positive risk is never silently removed.
#'
#' @return
#' For `q > 0`, a `tidyact_payment_basis_spec` that also inherits from
#' `tidyact_collective_risk_spec`, with thinned frequency, conditional positive
#' severity, `count_basis = "payment"`, the unchanged source, and metadata.
#' For `q = 0`, an explicit `tidyact_zero_payment_basis_spec` with no
#' per-payment severity, since conditioning on zero-probability events is
#' undefined. This degenerate result is not a collective risk specification.
#'
#' @references
#' Klugman, S. A., Panjer, H. H., and Willmot, G. E. (2012).
#' *Loss Models: From Data to Decisions*, 4th ed., Appendix C.
#'
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and
#' Nesbitt, C. J. (1997). *Actuarial Mathematics*, 2nd ed., Chapter 12.
#'
#' @seealso [coverage_spec()], [apply_coverage()]
#'
#' @examples
#' # Fully executable Poisson collective model with an ordinary deductible.
#' f <- frequency_spec("poisson", lambda = 4)
#' x <- severity_spec("discrete", values = c(50, 150, 250),
#'                    probs = c(0.2, 0.5, 0.3))
#' paid <- x |> apply_coverage(coverage_spec(deductible = 100))
#' model <- collective_risk_spec(f, paid)
#' by_payment <- model |> to_payment_basis()
#' by_payment$count_basis
#' by_payment$thinning$probability
#' aggregate_moments(by_payment, order = 1:2)
#'
#' @export
to_payment_basis <- function(x, label = NULL) {
  if (!inherits(x, "tidyact_collective_risk_spec") || !is.list(x)) {
    stop("`x` must be a collective risk specification.", call. = FALSE)
  }
  if (!identical(x$count_basis, "loss")) {
    stop("`x$count_basis` must be 'loss' before conversion.", call. = FALSE)
  }
  if (!inherits(x$frequency, "tidyact_frequency_spec")) {
    stop("`x$frequency` must be a frequency specification.", call. = FALSE)
  }
  payment <- x$severity
  if (!inherits(payment, "tidyact_payment_spec") ||
      !identical(payment$basis, "per_loss") ||
      !inherits(payment$coverage, "tidyact_coverage_spec")) {
    stop("`x$severity` must be a per-loss payment from apply_coverage().",
         call. = FALSE)
  }
  .collective_validate_label(label)
  if (is.null(label)) label <- x$label

  q <- .payment_basis_exact_q(payment)
  if (q == 0) {
    return(structure(
      list(schema_version = 1L, source = x, count_basis = "payment",
           basis = "per_payment", thinning = list(method = "binomial",
           probability = 0, source = x$frequency),
           positive_severity = NULL, resolution_status = "zero_payment",
           label = label,
           metadata = list(aggregate_loss = 0,
                           reason = "Conditioning on a zero-probability event is undefined")),
      class = c("tidyact_zero_payment_basis_spec", "tidyact_payment_basis_spec")
    ))
  }

  frequency <- .payment_basis_thin_frequency(x$frequency, q)
  positive <- .payment_basis_positive_severity(payment, q)
  result <- collective_risk_spec(
    frequency = frequency,
    severity = positive,
    count_basis = "payment",
    label = label
  )
  result$source <- x
  result$basis <- "per_payment"
  result$thinning <- list(method = "binomial", probability = q,
                          source = x$frequency)
  result$positive_severity <- positive
  result$resolution_status <- "resolved"
  result$metadata$aggregate_invariance <- "in_distribution"
  class(result) <- c("tidyact_payment_basis_spec", class(result))
  result
}

.payment_basis_exact_q <- function(payment) {
  if (identical(payment$dist, "discrete")) {
    q <- sum(payment$params$probs[payment$params$values > 0])
  } else {
    coverage <- payment$coverage
    terms <- .coverage_effective_terms(coverage)
    if (terms$coinsurance == 0 || terms$limit == 0) return(0)
    cutoff <- terms$deductible / terms$inflation
    # Compute the upper tail directly. 1 - F(cutoff) loses small q.
    q <- .payment_source_probability(payment$source, cutoff,
                                     lower.tail = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) ||
      !is.finite(q) || q < 0 || q > 1) {
    stop("The source severity did not provide a valid exact payment probability.",
         call. = FALSE)
  }
  as.numeric(q)
}

.payment_basis_thin_frequency <- function(x, q) {
  if (q == 1) return(x)
  p <- x$params
  label <- x$label
  if (x$dist %in% c("negative_binomial", "geometric")) {
    # p/(p+q(1-p)) avoids overflow in q(1-p)/p for tiny original p.
    thinned_prob <- p$prob / (p$prob + q * (1 - p$prob))
    if (thinned_prob == 1 && p$prob < 1 && q > 0) {
      stop("The thinned frequency cannot be represented by a `prob` below 1 ",
           "at machine precision. Retain the per-loss model instead.",
           call. = FALSE)
    }
  }
  switch(
    x$dist,
    poisson = frequency_spec("poisson", lambda = p$lambda * q, label = label),
    binomial = frequency_spec("binomial", size = p$size,
                              prob = p$prob * q, label = label),
    negative_binomial = frequency_spec("negative_binomial", size = p$size,
                         prob = thinned_prob,
                         label = label),
    geometric = frequency_spec("geometric",
                               prob = thinned_prob,
                               label = label),
    empirical = {
      tab <- .payment_basis_thin_pmf(p$values, p$probs, q)
      tab <- tab[tab$probability > 0, , drop = FALSE]
      empirical_frequency(tab[["count"]], weight = tab[["probability"]],
                          label = label)
    },
    stop("Exact thinning is not yet registered for frequency family '",
         x$dist, "'. No mixing model is approximated or discarded.",
         call. = FALSE)
  )
}

.payment_basis_positive_severity <- function(payment, q) {
  label <- payment$label
  t <- .coverage_effective_terms(payment$coverage)
  # Under a franchise whose cap is below the deductible, every positive
  # payment equals the cap, regardless of the ground-up severity law.
  constant_franchise <- t$deductible_type == "franchise" &&
    is.finite(t$limit) && t$limit <= t$deductible
  if (constant_franchise) {
    severity <- severity_spec("degenerate", value = t$coinsurance * t$limit,
                              label = label)
  } else if (identical(payment$dist, "discrete")) {
    p <- payment$params
    positive <- p$values > 0
    if (!any(positive)) stop("No strictly positive payment values.", call. = FALSE)
    probs <- p$probs[positive] / q
    probs <- probs / sum(probs)
    severity <- severity_spec("discrete", values = p$values[positive],
                              probs = probs, label = label)
  } else {
    raw_type <- payment$support$type
    # Removing the zero atom from a capped-free continuous source leaves a
    # continuous positive-payment law. The cap creates an upper atom otherwise.
    if (identical(payment$source$support$type, "continuous") &&
        is.infinite(t$limit)) {
      raw_type <- "continuous"
    }
    custom <- list(
      cdf = function(y, params, lower.tail = TRUE, log.p = FALSE) {
        if (!is.numeric(y)) stop("`y` must be numeric.", call. = FALSE)
        tail <- .payment_source_call(params$payment, "cdf", y,
                     lower.tail = FALSE, log.p = FALSE)
        tail <- tail / params$q
        tail[!is.na(y) & y <= 0] <- 1
        tail[!is.na(y) & y == Inf] <- 0
        if (any(!is.na(tail) & (tail < -1e-10 | tail > 1 + 1e-10))) {
          stop("Payment survival is inconsistent with the conditioning probability.",
               call. = FALSE)
        }
        tail <- pmax(0, pmin(1, tail))
        out <- if (lower.tail) 1 - tail else tail
        if (log.p) out <- log(out)
        out
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        p <- .severity_valid_prob(prob, log.p = log.p)
        if (log.p) p <- exp(p)
        if (!lower.tail) p <- 1 - p
        .payment_basis_positive_quantile(p, params$payment, params$q)
      },
      rng = function(n, params) {
        .payment_basis_positive_quantile(stats::runif(n),
                                         params$payment, params$q)
      },
      pdf = function(y, params, log = FALSE) {
        ans <- .payment_source_call(params$payment, "pdf", y,
                                    log = FALSE) / params$q
        if (log) log(ans) else ans
      },
      pmf = function(y, params, log = FALSE) {
        ans <- .payment_source_call(params$payment, "pmf", y,
                                    log = FALSE) / params$q
        ans[!is.na(y) & y <= 0] <- 0
        if (log) log(ans) else ans
      },
      raw_moment = function(order, params) {
        raw <- .payment_source_call(params$payment, "raw_moment", order)
        if (order == 0L) 1 else raw / params$q
      }
    )
    if (identical(raw_type, "continuous")) {
      custom$pmf <- NULL
    } else if (identical(raw_type, "discrete")) {
      custom$pdf <- NULL
    }
    severity <- structure(
      list(schema_version = 1L, dist = "custom",
           params = list(payment = payment, q = q),
           support = list(type = raw_type, lower = 0,
                          upper = payment$support$upper),
           data = NULL, custom = custom,
           label = label, metadata = list(count_basis = "payment")),
      class = c("tidyact_severity_spec", "tidyact_distribution_spec")
    )
  }
  severity$metadata$count_basis <- "payment"
  severity$source <- payment
  severity$basis <- "per_payment"
  class(severity) <- c("tidyact_positive_payment_spec", class(severity))
  severity
}

.payment_basis_positive_quantile <- function(p, payment, q) {
  ans <- rep(NA_real_, length(p))
  valid <- !is.na(p)
  if (!any(valid)) return(ans)
  if (any(p[valid] < 0 | p[valid] > 1)) return(rep(NaN, length(p)))
  src <- payment$source
  cov <- payment$coverage
  # Survival-scale quantiles avoid rounding 1-q+qp to one when q is small.
  inside <- which(valid & p > 0)
  if (length(inside)) {
    original <- .payment_source_call(src, "quantile",
                   q * (1 - p[inside]), lower.tail = FALSE,
                   log.p = FALSE)
    ans[inside] <- .payment_extended_amount(original, cov)
  }
  at_zero <- which(valid & p == 0)
  if (length(at_zero)) {
    terms <- .coverage_effective_terms(cov)
    if (identical(src$dist, "custom")) {
      ans[at_zero] <- NA_real_  # Exact positive-support infimum may be unavailable.
    } else if (terms$deductible_type == "franchise") {
      ans[at_zero] <- terms$coinsurance * pmin(terms$deductible, terms$limit)
    } else {
      ans[at_zero] <- 0
    }
  }
  ans
}

# Exact, pure mathematical adapter. The source CDF must not be approximated.
.payment_basis_positive_probability <- function(payment, loss_cdf = NULL) {
  if (!inherits(payment, "tidyact_payment_spec") ||
      !identical(payment$basis, "per_loss") ||
      !inherits(payment$coverage, "tidyact_coverage_spec")) {
    stop("`payment` must be a per-loss payment specification.",
         call. = FALSE)
  }
  coverage <- payment$coverage
  if (coverage$coinsurance == 0 || coverage$limit == 0) {
    return(0)
  }
  if (!is.function(loss_cdf)) {
    return(.payment_basis_exact_q(payment))
  }
  terms <- .coverage_effective_terms(coverage)
  cutoff <- terms$deductible / terms$inflation
  p0 <- loss_cdf(cutoff)
  if (!is.numeric(p0) || length(p0) != 1L || is.na(p0) ||
      !is.finite(p0) || p0 < 0 || p0 > 1) {
    stop("`loss_cdf` must return one exact probability in [0, 1].",
         call. = FALSE)
  }
  1 - p0
}

# Exact conditional CDF given a per-loss payment CDF and q > 0.
.payment_basis_conditional_cdf <- function(y, payment_cdf, q) {
  if (!is.numeric(q) || length(q) != 1L || is.na(q) ||
      !is.finite(q) || q <= 0 || q > 1) {
    stop("`q` must be strictly positive and at most one.", call. = FALSE)
  }
  if (!is.numeric(y) || !is.function(payment_cdf)) {
    stop("`y` must be numeric and `payment_cdf` must be a function.",
         call. = FALSE)
  }
  out <- rep(NA_real_, length(y))
  valid <- !is.na(y)
  out[valid & y <= 0] <- 0
  out[valid & y == Inf] <- 1
  idx <- which(valid & y > 0 & y < Inf)
  if (length(idx)) {
    p <- payment_cdf(y[idx])
    if (!is.numeric(p) || length(p) != length(idx) || anyNA(p) ||
        any(!is.finite(p)) || any(p < 0 | p > 1)) {
      stop("`payment_cdf` must return exact CDF probabilities in [0, 1].",
           call. = FALSE)
    }
    z <- (p - (1 - q)) / q
    # Bound only roundoff, not an invalid distribution or conditioning event.
    if (any(z < -1e-10 | z > 1 + 1e-10)) {
      stop("The payment CDF is inconsistent with `q`.", call. = FALSE)
    }
    out[idx] <- pmax(0, pmin(1, z))
  }
  out
}

# Finite-support thinning identity for rigorous discrete reference tests.
# This is not the algorithm for an unbounded frequency distribution.
.payment_basis_thin_pmf <- function(n, probability, q, tol = 1e-10) {
  if (!is.numeric(n) || length(n) == 0L || anyNA(n) ||
      any(!is.finite(n)) || any(n < 0) || any(n != floor(n)) ||
      !is.numeric(probability) || length(probability) != length(n) ||
      anyNA(probability) || any(!is.finite(probability)) ||
      any(probability < 0) || abs(sum(probability) - 1) > tol) {
    stop("Provide finite nonnegative integer counts and probabilities summing to one.",
         call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) ||
      !is.finite(q) || q < 0 || q > 1) {
    stop("`q` must be a probability in [0, 1].", call. = FALSE)
  }
  max_count <- max(n)
  if (max_count > 1e5) {
    stop("Finite-support PMF expansion is limited to 100000 counts.",
         call. = FALSE)
  }
  k <- seq.int(0, max_count)
  p <- vapply(k, function(ki) {
    sum(probability * stats::dbinom(ki, size = n, prob = q))
  }, numeric(1))
  tibble::tibble(count = k, probability = p)
}

# General exact frequency PGF identity; preserves hierarchical mixing.
.payment_basis_thin_pgf <- function(z, source_pgf, q) {
  if (!is.numeric(z) || anyNA(z) || any(!is.finite(z)) ||
      any(z < 0 | z > 1) || !is.function(source_pgf)) {
    stop("`z` must be in [0, 1] and `source_pgf` a function.",
         call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) ||
      !is.finite(q) || q < 0 || q > 1) {
    stop("`q` must be a probability in [0, 1].", call. = FALSE)
  }
  answer <- source_pgf(1 - q + q * z)
  if (!is.numeric(answer) || length(answer) != length(z) || anyNA(answer) ||
      any(!is.finite(answer)) || any(answer < -1e-10 | answer > 1 + 1e-10)) {
    stop("`source_pgf` must return valid PGF probabilities.",
         call. = FALSE)
  }
  pmax(0, pmin(1, answer))
}
