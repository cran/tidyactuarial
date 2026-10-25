# Shared validation and exact finite-support transformations for risk measures.

.risk_validate_levels <- function(level) {
  if (!is.numeric(level) || !length(level) || anyNA(level) ||
      any(!is.finite(level)) || any(level <= 0 | level >= 1)) {
    stop("`level` must contain probabilities strictly between zero and one.",
         call. = FALSE)
  }
  as.numeric(level)
}

.risk_validate_retention <- function(retention) {
  if (!is.numeric(retention) || !length(retention) || anyNA(retention) ||
      any(!is.finite(retention)) || any(retention < 0)) {
    stop("`retention` must contain finite nonnegative amounts.",
         call. = FALSE)
  }
  as.numeric(retention)
}

.risk_values <- function(x, value = NULL) {
  if (is.numeric(x) && is.null(dim(x))) {
    if (!is.null(value)) {
      stop("`value` is only used to select a data-frame column.",
           call. = FALSE)
    }
    out <- as.numeric(x)
  } else if (is.data.frame(x)) {
    if (is.null(value)) {
      candidates <- c("aggregate_loss", "claim_amount", "claim_count",
                      "loss", "payment")
      matches <- intersect(candidates, names(x))
      if (!length(matches)) {
        stop("Supply `value` as the name of a numeric loss column.",
             call. = FALSE)
      }
      value <- matches[1L]
    }
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !value %in% names(x)) {
      stop("`value` must be the name of a column in `x`.", call. = FALSE)
    }
    out <- x[[value]]
  } else {
    stop("`x` must contain observed losses or a supported risk distribution.",
         call. = FALSE)
  }
  if (!is.numeric(out) || !length(out) || anyNA(out) ||
      any(!is.finite(out)) || any(out < 0)) {
    stop("Loss observations must be finite nonnegative numbers without NA.",
         call. = FALSE)
  }
  as.numeric(out)
}

.risk_masses <- function(x, value = NULL) {
  count_tail <- 0
  count_tail_lower <- Inf
  if (inherits(x, "tidyact_zero_payment_basis_spec")) {
    return(list(loss = 0, probability = 1, tail = 0,
                method = "degenerate"))
  }
  if (inherits(x, "tidyact_aggregate_distribution")) {
    if (!is.null(value)) {
      stop("`value` only selects a data-frame column.", call. = FALSE)
    }
    loss <- x$grid$loss
    probability <- x$grid$probability
    tail <- x$tail_probability
    method <- if (isTRUE(x$metadata$exact)) "aggregate_grid_exact" else
      "aggregate_grid_approximate_or_unverified"
    if (identical(x$method, "fft") &&
        isTRUE(x$metadata$count_tail_probability > 0)) {
      count_tail <- x$metadata$count_tail_probability
      count_tail_lower <- (x$metadata$count_truncation + 1) *
        x$model$severity$support$lower
    }
  } else if (inherits(x, "tidyact_severity_spec") &&
             x$dist %in% c("discrete", "empirical", "degenerate")) {
    if (!is.null(value)) stop("`value` is only for data frames.", call. = FALSE)
    if (identical(x$dist, "degenerate")) {
      loss <- x$params$value
      probability <- 1
    } else {
      loss <- x$params$values
      probability <- x$params$probs
    }
    tail <- 0
    method <- "finite_distribution"
  } else if (is.numeric(x) || is.data.frame(x)) {
    values <- .risk_values(x, value)
    loss <- sort(unique(values))
    probability <- tabulate(match(values, loss), nbins = length(loss)) /
      length(values)
    tail <- 0
    method <- "empirical_simulation"
  } else {
    stop("An exact finite mass representation is not available for `x`.",
         call. = FALSE)
  }
  if (!is.numeric(loss) || !is.numeric(probability) ||
      !length(loss) || length(loss) != length(probability) ||
      anyNA(loss) || anyNA(probability) ||
      any(!is.finite(loss)) || any(loss < 0) ||
      any(!is.finite(probability)) || any(probability < 0) ||
      !is.numeric(tail) || length(tail) != 1L || is.na(tail) ||
      !is.finite(tail) || tail < 0 || tail > 1) {
    stop("Invalid retained probability mass or unobserved tail in `x`.",
         call. = FALSE)
  }
  ordered <- order(loss)
  loss <- loss[ordered]
  probability <- probability[ordered]
  if (anyDuplicated(loss)) {
    unique_loss <- unique(loss)
    probability <- as.numeric(rowsum(matrix(probability, ncol = 1L), match(loss, unique_loss),
                                     reorder = TRUE)[, 1L])
    loss <- unique_loss
  }
  total <- sum(probability) + tail
  if (!is.finite(total) || abs(total - 1) > 1e-7) {
    stop("Retained probabilities and omitted tail must sum to one.",
         call. = FALSE)
  }
  list(loss = loss, probability = probability, tail = tail, method = method,
       count_tail = count_tail, count_tail_lower = count_tail_lower)
}

.risk_discrete_var <- function(masses, level) {
  cumulative <- cumsum(masses$probability)
  max_captured <- cumulative[length(cumulative)]
  vapply(level, function(p) {
    if (p > max_captured) {
      stop("Requested VaR quantile lies in the uncomputed aggregate tail.",
           call. = FALSE)
    }
    k <- which(cumulative >= p)[1L]
    if (is.na(k)) {
      stop("Requested VaR quantile lies in the uncomputed aggregate tail.",
           call. = FALSE)
    }
    previous <- if (k == 1L) 0 else cumulative[k - 1L]
    if (isTRUE(masses$count_tail > 0) &&
        masses$loss[k] > masses$count_tail_lower &&
        previous + masses$count_tail >= p) {
      stop("Requested VaR is unresolved because FFT count truncation can ",
           "change the quantile. Reduce `tol` in aggregate_distribution().",
           call. = FALSE)
    }
    masses$loss[k]
  }, numeric(1))
}

.risk_exact_quantile <- function(x, level) {
  if (!inherits(x, "tidyact_severity_spec")) {
    stop("`x` has no registered continuous quantile provider.", call. = FALSE)
  }
  fn <- if (identical(x$dist, "custom")) x$custom$quantile else {
    registry <- .severity_registry()
    entry <- registry[[x$dist]]
    if (is.null(entry)) NULL else entry$quantile
  }
  if (!is.function(fn)) {
    stop("This severity has no exact quantile provider.", call. = FALSE)
  }
  quantile <- fn(level, x$params, lower.tail = TRUE, log.p = FALSE)
  if (!is.numeric(quantile) || length(quantile) != length(level) ||
      anyNA(quantile) || any(quantile < 0)) {
    stop("The severity quantile provider returned invalid results.",
         call. = FALSE)
  }
  as.numeric(quantile)
}

.risk_stop_loss_exact <- function(x, d) {
  if (.risk_exponential_payment(x)) {
    t <- .coverage_effective_terms(x$coverage)
    if (t$coinsurance == 0 || t$limit == 0) return(0)
    a <- d / t$coinsurance
    if (a >= t$limit) return(0)
    rate <- x$source$params$rate / t$inflation
    if (t$deductible_type == "ordinary") {
      threshold <- t$deductible + a
      width <- t$limit - a
      base <- 0
    } else {
      threshold <- max(t$deductible, a)
      width <- max(0, t$limit - threshold)
      base <- max(0, min(threshold, t$limit) - a)
    }
    # Integrate the exponential survival analytically. expm1 preserves
    # thin layers; no quadrature over a possibly enormous contract limit.
    return(t$coinsurance * exp(-rate * threshold) *
             (base + (-expm1(-rate * width)) / rate))
  }
  p <- x$params
  if (identical(x$dist, "exponential")) {
    return(exp(-p$rate * d) / p$rate)
  }
  if (identical(x$dist, "gamma")) {
    mean <- p$shape / p$rate
    return(pmax(0, mean * stats::pgamma(d, shape = p$shape + 1,
               rate = p$rate, lower.tail = FALSE) -
               d * stats::pgamma(d, shape = p$shape,
               rate = p$rate, lower.tail = FALSE)))
  }
  if (identical(x$dist, "lognormal")) {
    if (d == 0) return(exp(p$meanlog + p$sdlog^2 / 2))
    upper_first_moment <- exp(p$meanlog + p$sdlog^2 / 2) *
      stats::pnorm(log(d), mean = p$meanlog + p$sdlog^2,
                   sd = p$sdlog, lower.tail = FALSE)
    tail <- stats::plnorm(d, meanlog = p$meanlog, sdlog = p$sdlog,
                          lower.tail = FALSE)
    return(pmax(0, upper_first_moment - d * tail))
  }
  if (identical(x$dist, "pareto")) {
    if (p$shape <= 1) return(Inf)
    return(p$scale / (p$shape - 1) *
             exp((1 - p$shape) * log1p(d / p$scale)))
  }
  if (identical(x$dist, "weibull")) {
    shape <- p$shape
    scale <- p$scale
    z <- (d / scale)^shape
    upper_first_moment <- scale * gamma(1 + 1 / shape) *
      stats::pgamma(z, shape = 1 + 1 / shape, lower.tail = FALSE)
    tail <- stats::pweibull(d, shape = shape, scale = scale,
                           lower.tail = FALSE)
    return(pmax(0, upper_first_moment - d * tail))
  }
  if (x$dist %in% c("burr", "loglogistic")) {
    tail_index <- if (identical(x$dist, "burr")) p$shape1 * p$shape2 else
      p$shape
    if (tail_index <= 1) return(Inf)
    return(.risk_integrated_stop_loss(x, d))
  }
  if (identical(x$dist, "custom") && is.finite(x$support$upper) &&
      is.function(x$custom$cdf)) {
    return(.risk_integrated_stop_loss(x, d))
  }
  stop("No supported stop-loss formula or bounded-tail integration for '",
       x$dist, "'. Use simulations or a finite distribution.",
       call. = FALSE)
}


.risk_integrated_stop_loss <- function(x, retention) {
  upper <- x$support$upper
  if (is.finite(upper) && retention >= upper) return(0)
  callback <- if (identical(x$dist, "custom")) x$custom$cdf else
    .severity_registry()[[x$dist]]$cdf
  if (!is.function(callback)) {
    stop("Numerical stop-loss integration requires a survival provider.",
         call. = FALSE)
  }
  survival <- function(y) {
    z <- callback(y, x$params, lower.tail = FALSE, log.p = FALSE)
    if (!is.numeric(z) || length(z) != length(y) || anyNA(z) ||
        any(!is.finite(z)) || any(z < -1e-10 | z > 1 + 1e-10)) {
      stop("Invalid survival values during stop-loss integration.",
           call. = FALSE)
    }
    pmax(0, pmin(1, z))
  }
  answer <- tryCatch(
    stats::integrate(survival, lower = retention, upper = upper,
                     rel.tol = 1e-8, subdivisions = 300L,
                     stop.on.error = TRUE)$value,
    error = function(e) stop("Stop-loss numerical integration failed: ",
                             conditionMessage(e), call. = FALSE)
  )
  if (!is.finite(answer) || answer < 0) {
    stop("Stop-loss numerical integration did not converge to a finite premium.",
         call. = FALSE)
  }
  answer
}

.risk_stop_loss_method <- function(x) {
  if (.risk_exponential_payment(x)) return("closed_form")
  if (identical(x$dist, "custom") ||
      x$dist %in% c("burr", "loglogistic")) "numerical_survival" else
    "closed_form"
}

.risk_exponential_payment <- function(x) {
  inherits(x, "tidyact_payment_spec") &&
    identical(x$source$dist, "exponential")
}
