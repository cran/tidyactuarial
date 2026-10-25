#' Apply an insurance contract to a loss severity specification
#' @md
#'
#' @description
#' Creates a per-loss payment severity specification by combining an existing
#' ground-up severity model with a contract created by [coverage_spec()]. The
#' original model is retained and not discretized, fitted, or simulated.
#'
#' @param x A severity specification inheriting from `tidyact_severity_spec`.
#'   It must represent nonnegative ground-up losses, not an existing payment
#'   specification.
#' @param coverage A `tidyact_coverage_spec` object from [coverage_spec()].
#' @param label Optional nonempty character scalar. If `NULL`, the source
#'   severity label is inherited when present; otherwise the coverage label
#'   is used when present.
#'
#' @details
#' The payment variable is defined on a **per-loss** basis, retaining every
#' loss event even when the associated insurer payment is zero. This is not
#' the conditional payment-per-payment severity distribution.
#'
#' Let `X >= 0` be the ground-up severity, `f > 0` the inflation factor,
#' `d_f` the indexed (or fixed) deductible, `L_f` the indexed (or fixed)
#' layer limit, and `c` the insurer coinsurance fraction. The exact
#' transformation is
#' \deqn{Y_L = c \min\{(fX-d_f)_+, L_f\}}
#' for an ordinary deductible and
#' \deqn{Y_L = c \min\{fX,L_f\} 1(fX>d_f)}
#' for a franchise deductible. The `limit` caps the covered layer **before**
#' coinsurance; it is not a final insurer-payment limit.
#'
#' A continuous ground-up distribution can become mixed: probability mass
#' arises at zero and, when the limit is finite and binding, at the maximum
#' payment `c * L_f`. Pre-existing atoms are preserved under the pushforward
#' and may merge. When `coinsurance = 0` or `limit = 0`, the payment is
#' identically zero. No smoothing or implicit discretization is performed.
#'
#' The result is an executable per-loss payment distribution specification;
#' its `source` and `coverage` fields contain everything necessary to
#' recover the payment law from the underlying severity provider. The
#' internal pure functions `.coverage_payment_amount()`,
#' `.payment_cdf_from_loss()`, `.payment_boundary_atoms()`, and
#' `.payment_discrete_pushforward()` implement the transformation exactly
#' given loss values or exact source probabilities. Finite discrete sources
#' are transformed exactly into a finite discrete distribution. Continuous or mixed sources use the existing custom-provider
#' contract for CDF, quantiles, random generation and available exact moments.
#' Moment orders without an exact registered formula are reported unavailable;
#' no numerical integration or discretization is introduced silently.
#'
#' @return
#' An S3 object with classes `tidyact_payment_spec`,
#' `tidyact_severity_spec`, and `tidyact_distribution_spec`. The list stores
#' `source`, `coverage`, `basis = "per_loss"`, `label`, and a normalized
#' severity specification (`dist = "discrete"` for finite discrete losses or
#' `dist = "custom"` for continuous and mixed losses). The support bound
#' is not asserted to be the smallest possible bound under the original
#' distribution.
#'
#' @references
#' Klugman, S. A., Panjer, H. H., and Willmot, G. E. (2012).
#' *Loss Models: From Data to Decisions*, 4th ed., Chapter 8.
#'
#' Frees, E. W., Lee, G., and Yang, L., *Loss Data Analytics*,
#' chapters on claim severity and coverage modifications.
#'
#' @seealso [coverage_spec()]
#'
#' @examples
#' # Full, registered severity model.
#' ground_up <- severity_spec("gamma", shape = 2, rate = 0.001,
#'                            label = "Ground-up severity")
#' contract <- coverage_spec(deductible = 500, limit = 10000,
#'                           coinsurance = 0.8)
#' payments <- ground_up |> apply_coverage(contract)
#' inherits(payments, "tidyact_payment_spec")
#' payments$basis
#' payments$coverage$deductible
#'
#' @export
apply_coverage <- function(x, coverage, label = NULL) {
  if (!inherits(x, "tidyact_severity_spec") || !is.list(x) ||
      !is.list(x$support)) {
    stop("`x` must be a valid severity specification.", call. = FALSE)
  }
  if (inherits(x, "tidyact_payment_spec") ||
      inherits(x, "tidyact_positive_payment_spec")) {
    stop("`x` is already a payment specification; nesting is not supported.",
         call. = FALSE)
  }
  .severity_validate_support(x$support)
  if (x$support$lower < 0) {
    stop("`x` must represent nonnegative losses.", call. = FALSE)
  }
  if (!inherits(coverage, "tidyact_coverage_spec")) {
    stop("`coverage` must be created by coverage_spec().", call. = FALSE)
  }
  .severity_validate_label(label)
  if (is.null(label)) {
    label <- if (!is.null(x$label)) x$label else coverage$label
  }
  terms <- .coverage_effective_terms(coverage)
  certainly_zero <- terms$coinsurance == 0 || terms$limit == 0
  finite_discrete <- x$dist %in% c("discrete", "empirical", "degenerate")
  if (finite_discrete) {
    values <- if (identical(x$dist, "degenerate")) x$params$value else x$params$values
    probs <- if (identical(x$dist, "degenerate")) 1 else x$params$probs
  }
  if (certainly_zero) {
    out_dist <- "discrete"
    params <- list(values = 0, probs = 1)
    support <- list(type = "discrete", lower = 0, upper = 0)
    custom <- NULL
  } else if (finite_discrete) {
    table <- .payment_discrete_pushforward(values, probs, coverage)
    out_dist <- "discrete"
    params <- list(values = table$payment, probs = table$probability)
    support <- list(type = "discrete", lower = min(params$values),
                    upper = max(params$values))
    custom <- NULL
  } else {
    cdf <- .payment_source_operation(x, "cdf")
    quantile <- .payment_source_operation(x, "quantile")
    rng <- .payment_source_operation(x, "rng")
    if (!is.function(cdf) || !is.function(quantile) || !is.function(rng)) {
      stop("The source severity requires exact CDF, quantile and RNG providers.",
           call. = FALSE)
    }
    out_dist <- "custom"
    params <- list(source = x, coverage = coverage)
    custom <- .payment_custom_callbacks()
    upper <- if (is.finite(terms$limit)) {
      terms$coinsurance * terms$limit
    } else if (is.finite(x$support$upper)) {
      .coverage_payment_amount(x$support$upper, coverage)
    } else {
      Inf
    }
    if (is.na(upper) || (is.infinite(upper) && is.finite(terms$limit))) {
      stop("The insured payment bound cannot be represented numerically.",
           call. = FALSE)
    }
    kind <- if (identical(x$support$type, "continuous") &&
                terms$deductible == 0 && is.infinite(terms$limit)) {
      "continuous"
    } else if (x$support$type %in% c("discrete", "lattice")) {
      "discrete"
    } else {
      "mixed"
    }
    support <- list(type = kind, lower = 0, upper = upper)
    if (identical(kind, "continuous")) custom$pmf <- NULL
    if (identical(kind, "discrete")) custom$pdf <- NULL
  }
  support <- .severity_validate_support(support)
  structure(
    list(schema_version = 1L, dist = out_dist, params = params,
         support = support, data = NULL, custom = custom,
         source = x, coverage = coverage, basis = "per_loss",
         label = label,
         metadata = list(count_basis = "loss", transformation = "coverage",
                         source_distribution = x$dist)),
    class = c("tidyact_payment_spec", "tidyact_severity_spec",
              "tidyact_distribution_spec")
  )
}

#' Print a payment severity specification
#' @param x A `tidyact_payment_spec` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.tidyact_payment_spec <- function(x, ...) {
  cat("<tidyact_payment_spec>\n")
  cat("Basis: per loss\n")
  cat("Source: ", x$source$dist, "\n", sep = "")
  cat("Law: ", x$support$type, " [0, ", x$support$upper,"]\n", sep = "")
  if (!is.null(x$label)) cat("Label: ", x$label, "\n", sep = "")
  invisible(x)
}

# Internal mathematical helper; not exported.
.coverage_effective_terms <- function(coverage) {
  if (!inherits(coverage, "tidyact_coverage_spec")) {
    stop("`coverage` must be created by coverage_spec().", call. = FALSE)
  }
  f <- coverage$inflation_factor
  d <- coverage$deductible * if (coverage$index_deductible) f else 1
  L <- coverage$limit * if (coverage$index_limit) f else 1
  if (is.nan(L) || is.nan(d) || !is.finite(d) ||
      (!is.finite(L) && L != Inf)) {
    stop("Indexed contract amounts cannot be represented numerically.",
         call. = FALSE)
  }
  list(inflation = f, deductible = d, limit = L,
       coinsurance = coverage$coinsurance,
       deductible_type = coverage$deductible_type)
}

# Internal mathematical helper; not exported.
.coverage_payment_amount <- function(loss, coverage) {
  if (!is.numeric(loss) || anyNA(loss) || any(!is.finite(loss)) ||
      any(loss < 0)) {
    stop("`loss` must be a finite nonnegative numeric vector.", call. = FALSE)
  }
  t <- .coverage_effective_terms(coverage)
  if (t$coinsurance == 0 || t$limit == 0) {
    return(rep(0, length(loss)))
  }
  z <- t$inflation * loss
  if (any(!is.finite(z))) {
    stop("Inflated losses cannot be represented numerically.", call. = FALSE)
  }
  if (t$deductible_type == "ordinary") {
    layer <- pmin(pmax(z - t$deductible, 0), t$limit)
  } else {
    layer <- ifelse(z > t$deductible, pmin(z, t$limit), 0)
  }
  t$coinsurance * layer
}

# Internal mathematical helper; not exported.
.payment_cdf_from_loss <- function(payment, coverage, loss_cdf) {
  if (!is.function(loss_cdf)) {
    stop("`loss_cdf` must be an exact CDF function.", call. = FALSE)
  }
  if (!is.numeric(payment)) {
    stop("`payment` must be numeric.", call. = FALSE)
  }
  t <- .coverage_effective_terms(coverage)
  out <- rep(NA_real_, length(payment))
  valid <- !is.na(payment)
  out[valid & payment < 0] <- 0
  out[valid & payment == Inf] <- 1

  if (t$coinsurance == 0 || t$limit == 0) {
    out[valid & payment >= 0] <- 1
    return(out)
  }

  upper <- t$coinsurance * t$limit
  out[valid & payment >= upper] <- 1
  idx <- which(valid & payment >= 0 & payment < upper)
  if (length(idx)) {
    v <- payment[idx]
    if (t$deductible_type == "ordinary") {
      threshold <- (t$deductible + v / t$coinsurance) / t$inflation
    } else {
      threshold <- pmax(t$deductible, v / t$coinsurance) / t$inflation
    }
    p <- loss_cdf(threshold)
    if (!is.numeric(p) || length(p) != length(idx) || anyNA(p) ||
        any(!is.finite(p)) || any(p < 0 | p > 1)) {
      stop("`loss_cdf` must return exact probabilities in [0, 1].",
           call. = FALSE)
    }
    out[idx] <- p
  }
  out
}

# Internal mathematical helper; not exported.
.payment_boundary_atoms <- function(coverage, loss_cdf, loss_left_cdf = NULL) {
  if (!is.function(loss_cdf)) {
    stop("`loss_cdf` must be an exact CDF function.", call. = FALSE)
  }
  if (!is.null(loss_left_cdf) && !is.function(loss_left_cdf)) {
    stop("`loss_left_cdf` must be an exact left-limit CDF function.",
         call. = FALSE)
  }
  t <- .coverage_effective_terms(coverage)
  if (t$coinsurance == 0 || t$limit == 0) {
    return(tibble::tibble(payment = 0, probability = 1))
  }
  probability <- function(fun, value) {
    ans <- fun(value)
    if (!is.numeric(ans) || length(ans) != 1L || is.na(ans) ||
        !is.finite(ans) || ans < 0 || ans > 1) {
      stop("CDF callbacks must return scalar probabilities in [0, 1].",
           call. = FALSE)
    }
    ans
  }

  d <- t$deductible / t$inflation
  p0 <- probability(loss_cdf, d)
  cap <- t$coinsurance * t$limit
  if (is.infinite(cap)) {
    return(tibble::tibble(payment = 0, probability = p0))
  }
  if (t$deductible_type == "franchise" && t$limit <= t$deductible) {
    pcap <- 1 - p0
  } else {
    if (is.null(loss_left_cdf)) {
      stop("An exact `loss_left_cdf` is required for the cap atom.",
           call. = FALSE)
    }
    threshold <- if (t$deductible_type == "ordinary") {
      (t$deductible + t$limit) / t$inflation
    } else {
      t$limit / t$inflation
    }
    pcap <- 1 - probability(loss_left_cdf, threshold)
  }
  tibble::tibble(payment = c(0, cap), probability = c(p0, pcap))
}

# Internal mathematical helper; not exported.
.payment_discrete_pushforward <- function(loss, probability, coverage,
                                          tolerance = 1e-10) {
  if (!is.numeric(loss) || anyNA(loss) || any(!is.finite(loss)) ||
      any(loss < 0) || !is.numeric(probability) ||
      length(loss) != length(probability) || length(loss) == 0L ||
      anyNA(probability) || any(!is.finite(probability)) ||
      any(probability < 0) ||
      abs(sum(probability) - 1) > tolerance) {
    stop("Provide finite nonnegative losses and probabilities summing to one.",
         call. = FALSE)
  }
  payment <- .coverage_payment_amount(loss, coverage)
  distinct <- sort(unique(payment))
  weights <- vapply(distinct, function(v) sum(probability[payment == v]),
                    numeric(1))
  tibble::tibble(payment = distinct, probability = weights)
}


# Source provider resolution uses the established severity-spec registry.
.payment_source_operation <- function(x, operation) {
  if (identical(x$dist, "custom")) return(x$custom[[operation]])
  entry <- .severity_registry()[[x$dist]]
  if (is.null(entry)) return(NULL)
  entry[[operation]]
}

.payment_source_call <- function(x, operation, ..., error = TRUE) {
  fun <- .payment_source_operation(x, operation)
  if (!is.function(fun)) {
    if (error) stop("No exact source severity `", operation,
                    "` provider is available.", call. = FALSE)
    return(NULL)
  }
  fun(..., x$params)
}

.payment_source_probability <- function(x, threshold, lower.tail = TRUE) {
  answer <- .payment_source_call(x, "cdf", threshold,
                                 lower.tail = lower.tail, log.p = FALSE)
  if (!is.numeric(answer) || length(answer) != length(threshold) ||
      anyNA(answer) || any(!is.finite(answer)) ||
      any(answer < 0 | answer > 1)) {
    stop("The source CDF returned invalid probabilities.", call. = FALSE)
  }
  answer
}

.payment_source_atom <- function(x, point) {
  pmf <- .payment_source_operation(x, "pmf")
  if (!is.function(pmf)) return(rep(0, length(point)))
  result <- pmf(point, x$params, log = FALSE)
  if (!is.numeric(result) || length(result) != length(point) ||
      anyNA(result) || any(!is.finite(result)) ||
      any(result < 0 | result > 1)) {
    stop("The source PMF returned invalid probabilities.", call. = FALSE)
  }
  result
}

.payment_custom_callbacks <- function() {
  list(
    cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
      x <- params$source
      F <- function(t) .payment_source_probability(x, t)
      result <- if (lower.tail) {
        .payment_cdf_from_loss(q, params$coverage, F)
      } else {
        .payment_survival_from_loss(q, params$coverage,
          function(t) .payment_source_probability(x, t, lower.tail = FALSE))
      }
      if (log.p) result <- log(result)
      result
    },
    quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
      p <- .severity_valid_prob(prob, log.p = log.p)
      if (log.p) p <- exp(p)
      if (!lower.tail) p <- 1 - p
      source <- params$source
      original <- .payment_source_call(source, "quantile", p,
                                       lower.tail = TRUE, log.p = FALSE)
      .payment_extended_amount(original, params$coverage)
    },
    rng = function(n, params) {
      x <- .payment_source_call(params$source, "rng", n)
      .coverage_payment_amount(x, params$coverage)
    },
    pdf = function(x, params, log = FALSE) {
      ans <- .payment_continuous_density(x, params$source, params$coverage)
      if (log) log(ans) else ans
    },
    pmf = function(x, params, log = FALSE) {
      ans <- .payment_point_masses(x, params$source, params$coverage)
      if (log) log(ans) else ans
    },
    raw_moment = function(order, params) {
      .payment_exact_raw_moment(order, params$source, params$coverage)
    }
  )
}

.payment_extended_amount <- function(loss, coverage) {
  if (!is.numeric(loss)) stop("Source quantile must be numeric.", call. = FALSE)
  out <- rep(NA_real_, length(loss))
  finite <- !is.na(loss) & is.finite(loss)
  if (any(finite)) out[finite] <- .coverage_payment_amount(loss[finite], coverage)
  pos_inf <- !is.na(loss) & loss == Inf
  terms <- .coverage_effective_terms(coverage)
  if (any(pos_inf)) {
    out[pos_inf] <- if (terms$coinsurance == 0 || terms$limit == 0) 0 else
      terms$coinsurance * terms$limit
  }
  out
}

.payment_point_masses <- function(y, source, coverage) {
  if (!is.numeric(y)) stop("Payment points must be numeric.", call. = FALSE)
  ans <- numeric(length(y))
  ans[is.na(y)] <- NA_real_
  terms <- .coverage_effective_terms(coverage)
  if (terms$coinsurance == 0 || terms$limit == 0) {
    ans[!is.na(y) & y == 0] <- 1
    return(ans)
  }
  zero <- !is.na(y) & y == 0
  cutoff <- terms$deductible / terms$inflation
  if (any(zero)) ans[zero] <- .payment_source_probability(source, cutoff)
  cap <- terms$coinsurance * terms$limit
  capped <- !is.na(y) & y == cap & is.finite(cap)
  if (any(capped)) {
    if (terms$deductible_type == "franchise" &&
        terms$limit <= terms$deductible) {
      atom <- .payment_source_probability(source, cutoff, lower.tail = FALSE)
    } else {
      threshold <- if (terms$deductible_type == "ordinary") {
        (terms$deductible + terms$limit) / terms$inflation
      } else {
        terms$limit / terms$inflation
      }
      atom <- .payment_source_probability(source, threshold,
                                          lower.tail = FALSE) +
        (if (is.finite(threshold)) .payment_source_atom(source, threshold) else 0)
    }
    ans[capped] <- pmax(0, pmin(1, atom))
  }
  middle <- which(!is.na(y) & y > 0 & y < cap)
  if (length(middle)) {
    values <- y[middle]
    threshold <- if (terms$deductible_type == "ordinary") {
      (terms$deductible + values / terms$coinsurance) / terms$inflation
    } else {
      values / (terms$coinsurance * terms$inflation)
    }
    active <- if (terms$deductible_type == "franchise") {
      values > terms$coinsurance * terms$deductible
    } else rep(TRUE, length(values))
    prob <- numeric(length(values))
    if (any(active)) prob[active] <- .payment_source_atom(source, threshold[active])
    ans[middle] <- prob
  }
  ans
}

.payment_continuous_density <- function(y, source, coverage) {
  if (!is.numeric(y)) stop("Payment points must be numeric.", call. = FALSE)
  ans <- numeric(length(y))
  ans[is.na(y)] <- NA_real_
  density <- .payment_source_operation(source, "pdf")
  if (!is.function(density)) return(ans)
  terms <- .coverage_effective_terms(coverage)
  if (terms$coinsurance == 0 || terms$limit == 0) return(ans)
  cap <- terms$coinsurance * terms$limit
  idx <- which(!is.na(y) & y > 0 & y < cap)
  if (!length(idx)) return(ans)
  v <- y[idx]
  threshold <- if (terms$deductible_type == "ordinary") {
    (terms$deductible + v / terms$coinsurance) / terms$inflation
  } else {
    v / (terms$coinsurance * terms$inflation)
  }
  active <- if (terms$deductible_type == "franchise") {
    v > terms$coinsurance * terms$deductible
  } else rep(TRUE, length(v))
  z <- numeric(length(v))
  if (any(active)) {
    z[active] <- density(threshold[active], source$params, log = FALSE) /
      (terms$coinsurance * terms$inflation)
  }
  ans[idx] <- z
  ans
}

# Only mathematically exact formulas are offered. `NA_real_` means the
# installed severity registry does not supply the needed exact moment.
.payment_exact_raw_moment <- function(order, source, coverage) {
  if (order == 0L) return(1)
  terms <- .coverage_effective_terms(coverage)
  if (terms$coinsurance == 0 || terms$limit == 0) return(0)
  if (identical(source$dist, "exponential")) {
    rate <- source$params$rate / terms$inflation
    d <- terms$deductible
    L <- terms$limit
    if (identical(terms$deductible_type, "ordinary")) {
      prob <- if (is.infinite(L)) 1 else stats::pgamma(L, shape = order, rate = rate)
      return(terms$coinsurance^order * exp(-rate * d) *
               exp(lgamma(order + 1) - order * log(rate)) * prob)
    }
    if (L <= d) return((terms$coinsurance * L)^order * exp(-rate * d))
    upper_diff <- if (is.infinite(L)) {
      stats::pgamma(d, shape = order + 1L, rate = rate, lower.tail = FALSE)
    } else {
      stats::pgamma(d, shape = order + 1L, rate = rate, lower.tail = FALSE) -
        stats::pgamma(L, shape = order + 1L, rate = rate, lower.tail = FALSE)
    }
    truncated <- exp(lgamma(order + 1) - order * log(rate)) *
      pmax(0, upper_diff)
    capped <- if (is.infinite(L)) 0 else L^order * exp(-rate * L)
    return(terms$coinsurance^order * (truncated + capped))
  }
  if (terms$deductible == 0 && is.infinite(terms$limit)) {
    raw_fn <- .payment_source_operation(source, "raw_moment")
    if (is.function(raw_fn)) {
      moment <- raw_fn(order, source$params)
      return((terms$coinsurance * terms$inflation)^order * moment)
    }
  }
  NA_real_
}

# Strict survival P(Y_L > y), evaluated directly from source upper tails.
# This is more stable than subtracting the payment CDF from one.
.payment_survival_from_loss <- function(payment, coverage, loss_survival) {
  if (!is.numeric(payment) || !is.function(loss_survival)) {
    stop("A numeric payment vector and exact loss-survival callback are required.",
         call. = FALSE)
  }
  terms <- .coverage_effective_terms(coverage)
  ans <- rep(NA_real_, length(payment))
  valid <- !is.na(payment)
  ans[valid & payment < 0] <- 1
  if (terms$coinsurance == 0 || terms$limit == 0) {
    ans[valid & payment >= 0] <- 0
    return(ans)
  }
  cap <- terms$coinsurance * terms$limit
  ans[valid & payment >= cap] <- 0
  idx <- which(valid & payment >= 0 & payment < cap)
  if (length(idx)) {
    v <- payment[idx]
    threshold <- if (terms$deductible_type == "ordinary") {
      (terms$deductible + v / terms$coinsurance) / terms$inflation
    } else {
      pmax(terms$deductible, v / terms$coinsurance) / terms$inflation
    }
    z <- loss_survival(threshold)
    if (!is.numeric(z) || length(z) != length(idx) || anyNA(z) ||
        any(!is.finite(z)) || any(z < 0 | z > 1)) {
      stop("The source survival callback must return probabilities in [0, 1].",
           call. = FALSE)
    }
    ans[idx] <- z
  }
  ans
}
