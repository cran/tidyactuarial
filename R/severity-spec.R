#' Specify an individual claim severity distribution
#'
#' @md
#' @description
#' Creates a lightweight specification of a nonnegative claim amount \eqn{X}.
#' A specification records the probabilistic model; it does not fit data,
#' draw samples, discretize amounts, or calculate aggregate losses.
#'
#' @param dist One of `"exponential"`, `"gamma"`, `"lognormal"`,
#'   `"weibull"`, `"pareto"`, `"burr"`, `"loglogistic"`, `"degenerate"`,
#'   `"discrete"`, or `"custom"`.
#' @param ... Named distribution parameters, as described in *Details*.
#' @param custom Named list of computational functions and support information
#'   when `dist = "custom"`; `NULL` for built-in families.
#' @param label Optional nonempty character scalar identifying the model.
#'
#' @details
#' All severities satisfy \eqn{X \ge 0}. The supported built-in models are:
#'
#' - **Exponential:** `rate > 0`, matching [stats::dexp()].
#' - **Gamma:** `shape > 0`, `rate > 0`, matching [stats::dgamma()].
#' - **Lognormal:** `meanlog` is finite and `sdlog > 0`, matching
#'   [stats::dlnorm()]. A point mass is specified with `"degenerate"`.
#' - **Weibull:** `shape > 0`, `scale > 0`, matching [stats::dweibull()].
#' - **Pareto II (Lomax):** `shape > 0`, `scale > 0`; the survival function is
#'   \eqn{P(X>x)=(1+x/\mathrm{scale})^{-\mathrm{shape}}} for \eqn{x \ge 0}.
#'   This is *not* the Pareto I distribution with positive lower endpoint.
#' - **Burr XII:** `shape1 > 0`, `shape2 > 0`, `scale > 0`, with survival
#'   function \eqn{P(X>x)=[1+(x/\mathrm{scale})^{\mathrm{shape1}}]^{-\mathrm{shape2}}}.
#' - **Log-logistic:** `shape > 0`, `scale > 0`, the Burr XII special case
#'   `shape1 = shape`, `shape2 = 1`.
#' - **Degenerate:** `value >= 0` is a finite amount with probability one.
#' - **Finite discrete:** `values` is a strictly increasing vector of distinct
#'   finite nonnegative amounts; `probs` contains strictly positive
#'   probabilities, with the same length and summing to one (to numerical
#'   tolerance). Amounts need not be integer multiples of a common unit.
#'
#' Unlike a frequency distribution, an individual claim severity may have an
#' atom at zero and need not be absolutely continuous. The `"degenerate"` and
#' `"discrete"` models are represented as discrete specifications rather than
#' being forced into a continuous family.
#'
#' For `dist = "custom"`, `custom` must be a named list with `support` and at
#' least one of `pdf` or `pmf`. The `support` list must contain `type`, `lower`,
#' and `upper`, with `type` in `"continuous"`, `"discrete"`, `"lattice"`,
#' or `"mixed"`; lattice support additionally requires a positive
#' `step`. The `support` component may also be a function of `params`
#' returning that list. A mixed specification requires `pdf`, `pmf`, and `cdf`.
#' Custom operations, when provided, follow these interfaces:
#' `pdf(x, params, log = FALSE)`, `pmf(x, params, log = FALSE)`,
#' `cdf(q, params, lower.tail = TRUE, log.p = FALSE)`,
#' `quantile(p, params, lower.tail = TRUE, log.p = FALSE)`,
#' `rng(n, params)`, `raw_moment(order, params)`,
#' `limited_moment(order, limit, params)`, `mgf(t, params)`, and
#' `laplace(t, params)`. The argument `params` is the named list collected
#' from `...`. As with `frequency_spec()`, the constructor checks the *shape*
#' of a custom definition but cannot certify its mathematical validity.
#'
#' The internal registry supplies probability functions and raw moments for
#' built-in families. Some advanced operations (for example limited moments,
#' the MGF of heavy-tailed distributions, and numerical Laplace transforms)
#' are intentionally left unavailable until their public interfaces are set.
#' In particular, no numerical approximation is performed by this constructor.
#'
#' @return An object with classes `"tidyact_severity_spec"` and
#'   `"tidyact_distribution_spec"`. The ordered fields are `schema_version`,
#'   `dist`, `params`, `support`, `data`, `custom`, `label`, and `metadata`.
#'   Its `data` field is `NULL`; empirical severities are handled separately.
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and
#' Nesbitt, C. J. (1997). *Actuarial Mathematics*, 2nd ed.
#' Society of Actuaries.
#'
#' @seealso [frequency_spec()]
#'
#' @examples
#' severity_spec("exponential", rate = 0.002)
#' severity_spec("gamma", shape = 2.5, rate = 0.001)
#' severity_spec("lognormal", meanlog = 8, sdlog = 0.8)
#' severity_spec("weibull", shape = 1.7, scale = 5000)
#'
#' # Pareto II and Burr XII use severity distributions starting at zero.
#' severity_spec("pareto", shape = 3, scale = 2000)
#' severity_spec("burr", shape1 = 1.5, shape2 = 2.5, scale = 1500)
#' severity_spec("loglogistic", shape = 2.5, scale = 1500)
#'
#' # A fixed benefit amount, including a point mass at zero.
#' severity_spec("degenerate", value = 1000)
#' severity_spec("degenerate", value = 0)
#'
#' # A claim amount distribution with a genuine atom at zero.
#' severity_spec(
#'   "discrete",
#'   values = c(0, 500, 2500),
#'   probs = c(0.1, 0.6, 0.3),
#'   label = "Illustrative claim payments"
#' )
#'
#' # An explicitly defined custom continuous model.
#' severity_spec(
#'   "custom",
#'   rate = 0.002,
#'   custom = list(
#'     pdf = function(x, params, log = FALSE) {
#'       stats::dexp(x, rate = params$rate, log = log)
#'     },
#'     support = list(type = "continuous", lower = 0, upper = Inf)
#'   )
#' )
#' @export
severity_spec <- function(dist, ..., custom = NULL, label = NULL) {
  dist <- .severity_validate_dist(dist)
  params <- list(...)
  .severity_validate_names(params, what = "Distribution parameters")
  .severity_validate_label(label)

  registry <- .severity_registry()
  if (!dist %in% names(registry)) {
    stop(
      "`dist` must be one of: ",
      paste(sprintf('"%s"', names(registry)), collapse = ", "), ".",
      call. = FALSE
    )
  }

  entry <- registry[[dist]]
  if (identical(dist, "custom")) {
    .severity_validate_custom(custom)
    support_definition <- if (is.function(custom$support)) {
      custom$support(params)
    } else {
      custom$support
    }
    support <- .severity_validate_support(support_definition)
    .severity_validate_custom_support(custom, support)
    custom$support <- NULL
  } else {
    if (!is.null(custom)) {
      stop("`custom` is allowed only when `dist = \"custom\"`.", call. = FALSE)
    }
    entry$validate(params)
    support <- .severity_validate_support(entry$support(params))
  }

  structure(
    list(
      schema_version = 1L,
      dist = dist,
      params = params,
      support = support,
      data = NULL,
      custom = custom,
      label = label,
      metadata = list()
    ),
    class = c("tidyact_severity_spec", "tidyact_distribution_spec")
  )
}


#' Print a severity distribution specification
#'
#' @param x A `tidyact_severity_spec` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.tidyact_severity_spec <- function(x, ...) {
  names_pretty <- c(
    exponential = "Exponential", gamma = "Gamma", lognormal = "Lognormal",
    weibull = "Weibull", pareto = "Pareto II (Lomax)", burr = "Burr XII",
    loglogistic = "Log-logistic", degenerate = "Degenerate",
    discrete = "Finite discrete", empirical = "Empirical", custom = "Custom"
  )
  cat("<tidyact_severity_spec>\n")
  cat("Distribution: ", unname(names_pretty[x$dist]), "\n", sep = "")
  if (length(x$params)) {
    cat("Parameters:\n")
    for (nm in names(x$params)) {
      value <- x$params[[nm]]
      if (is.atomic(value) && length(value) <= 5L) {
        value_text <- paste(format(value, trim = TRUE), collapse = ", ")
      } else if (is.atomic(value)) {
        value_text <- paste0(
          paste(format(utils::head(value, 4L), trim = TRUE), collapse = ", "), ", ..."
        )
      } else {
        value_text <- paste0("<", paste(class(value), collapse = "/"), ">")
      }
      cat("  ", nm, ": ", value_text, "\n", sep = "")
    }
  }
  upper <- if (is.infinite(x$support$upper)) "Inf" else
    format(x$support$upper, trim = TRUE)
  cat(
    "Support: ", x$support$type, " [",
    format(x$support$lower, trim = TRUE), ", ", upper, "]\n",
    sep = ""
  )
  if (!is.null(x$label)) cat("Label: ", x$label, "\n", sep = "")
  invisible(x)
}


# Registry -----------------------------------------------------------------

.severity_registry <- function() {
  continuous_support <- function(params) {
    list(type = "continuous", lower = 0, upper = Inf)
  }

  list(
    exponential = list(
      validate = function(params) {
        .severity_check_parameters(params, "rate")
        .severity_check_number(params$rate, "rate", lower = 0, strict = TRUE)
      },
      support = continuous_support,
      pdf = function(x, params, log = FALSE) stats::dexp(x, rate = params$rate, log = log),
      pmf = NULL,
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::pexp(q, rate = params$rate, lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qexp(prob, rate = params$rate, lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) stats::rexp(n, rate = params$rate),
      raw_moment = function(order, params) {
        .severity_check_order(order)
        exp(lgamma(order + 1) - order * log(params$rate))
      },
      limited_moment = NULL,
      mgf = function(t, params) {
        out <- rep(Inf, length(t))
        inside <- !is.na(t) & t < params$rate
        out[inside] <- params$rate / (params$rate - t[inside])
        out[is.na(t)] <- NA_real_
        out
      },
      laplace = function(t, params) {
        .severity_check_laplace_argument(t)
        params$rate / (params$rate + t)
      }
    ),
    gamma = list(
      validate = function(params) {
        .severity_check_parameters(params, c("shape", "rate"))
        .severity_check_number(params$shape, "shape", lower = 0, strict = TRUE)
        .severity_check_number(params$rate, "rate", lower = 0, strict = TRUE)
      },
      support = continuous_support,
      pdf = function(x, params, log = FALSE) {
        stats::dgamma(x, shape = params$shape, rate = params$rate, log = log)
      },
      pmf = NULL,
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::pgamma(q, shape = params$shape, rate = params$rate,
                      lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qgamma(prob, shape = params$shape, rate = params$rate,
                      lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) stats::rgamma(n, shape = params$shape, rate = params$rate),
      raw_moment = function(order, params) {
        .severity_check_order(order)
        exp(lgamma(params$shape + order) - lgamma(params$shape) - order * log(params$rate))
      },
      limited_moment = NULL,
      mgf = function(t, params) {
        out <- rep(Inf, length(t))
        inside <- !is.na(t) & t < params$rate
        out[inside] <- (params$rate / (params$rate - t[inside]))^params$shape
        out[is.na(t)] <- NA_real_
        out
      },
      laplace = function(t, params) {
        .severity_check_laplace_argument(t)
        (params$rate / (params$rate + t))^params$shape
      }
    ),
    lognormal = list(
      validate = function(params) {
        .severity_check_parameters(params, c("meanlog", "sdlog"))
        .severity_check_number(params$meanlog, "meanlog")
        .severity_check_number(params$sdlog, "sdlog", lower = 0, strict = TRUE)
      },
      support = continuous_support,
      pdf = function(x, params, log = FALSE) {
        stats::dlnorm(x, meanlog = params$meanlog, sdlog = params$sdlog, log = log)
      },
      pmf = NULL,
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::plnorm(q, meanlog = params$meanlog, sdlog = params$sdlog,
                      lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qlnorm(prob, meanlog = params$meanlog, sdlog = params$sdlog,
                      lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) stats::rlnorm(n, meanlog = params$meanlog, sdlog = params$sdlog),
      raw_moment = function(order, params) {
        .severity_check_order(order)
        exp(order * params$meanlog + 0.5 * order^2 * params$sdlog^2)
      },
      limited_moment = NULL,
      mgf = NULL,
      laplace = NULL
    ),
    weibull = list(
      validate = function(params) {
        .severity_check_parameters(params, c("shape", "scale"))
        .severity_check_number(params$shape, "shape", lower = 0, strict = TRUE)
        .severity_check_number(params$scale, "scale", lower = 0, strict = TRUE)
      },
      support = continuous_support,
      pdf = function(x, params, log = FALSE) {
        stats::dweibull(x, shape = params$shape, scale = params$scale, log = log)
      },
      pmf = NULL,
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::pweibull(q, shape = params$shape, scale = params$scale,
                        lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qweibull(prob, shape = params$shape, scale = params$scale,
                        lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) stats::rweibull(n, shape = params$shape, scale = params$scale),
      raw_moment = function(order, params) {
        .severity_check_order(order)
        exp(order * log(params$scale) + lgamma(1 + order / params$shape))
      },
      limited_moment = NULL,
      mgf = NULL,
      laplace = NULL
    ),
    pareto = list(
      validate = function(params) {
        .severity_check_parameters(params, c("shape", "scale"))
        .severity_check_number(params$shape, "shape", lower = 0, strict = TRUE)
        .severity_check_number(params$scale, "scale", lower = 0, strict = TRUE)
      },
      support = continuous_support,
      pdf = function(x, params, log = FALSE) {
        .severity_dpareto2(x, params$shape, params$scale, log = log)
      },
      pmf = NULL,
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        .severity_ppareto2(q, params$shape, params$scale,
                           lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        params$scale * expm1(-.severity_log_survival(prob, lower.tail, log.p) / params$shape)
      },
      rng = function(n, params) {
        params$scale * expm1(-log(stats::runif(n)) / params$shape)
      },
      raw_moment = function(order, params) {
        .severity_check_order(order)
        if (order == 0) return(1)
        if (order >= params$shape) return(Inf)
        exp(order * log(params$scale) + lgamma(1 + order) +
              lgamma(params$shape - order) - lgamma(params$shape))
      },
      limited_moment = NULL,
      mgf = NULL,
      laplace = NULL
    ),
    burr = list(
      validate = function(params) {
        .severity_check_parameters(params, c("shape1", "shape2", "scale"))
        .severity_check_number(params$shape1, "shape1", lower = 0, strict = TRUE)
        .severity_check_number(params$shape2, "shape2", lower = 0, strict = TRUE)
        .severity_check_number(params$scale, "scale", lower = 0, strict = TRUE)
      },
      support = continuous_support,
      pdf = function(x, params, log = FALSE) {
        .severity_dburr12(x, params$shape1, params$shape2, params$scale, log = log)
      },
      pmf = NULL,
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        .severity_pburr12(q, params$shape1, params$shape2, params$scale,
                          lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        .severity_qburr12(prob, params$shape1, params$shape2, params$scale,
                          lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) {
        .severity_qburr12(stats::runif(n), params$shape1, params$shape2, params$scale)
      },
      raw_moment = function(order, params) {
        .severity_check_order(order)
        if (order == 0) return(1)
        h <- order / params$shape1
        if (h >= params$shape2) return(Inf)
        exp(order * log(params$scale) + lgamma(1 + h) +
              lgamma(params$shape2 - h) - lgamma(params$shape2))
      },
      limited_moment = NULL,
      mgf = NULL,
      laplace = NULL
    ),
    loglogistic = list(
      validate = function(params) {
        .severity_check_parameters(params, c("shape", "scale"))
        .severity_check_number(params$shape, "shape", lower = 0, strict = TRUE)
        .severity_check_number(params$scale, "scale", lower = 0, strict = TRUE)
      },
      support = continuous_support,
      pdf = function(x, params, log = FALSE) {
        .severity_dburr12(x, params$shape, 1, params$scale, log = log)
      },
      pmf = NULL,
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        .severity_pburr12(q, params$shape, 1, params$scale,
                          lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        .severity_qburr12(prob, params$shape, 1, params$scale,
                          lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) .severity_qburr12(stats::runif(n), params$shape, 1, params$scale),
      raw_moment = function(order, params) {
        .severity_check_order(order)
        if (order == 0) return(1)
        h <- order / params$shape
        if (h >= 1) return(Inf)
        exp(order * log(params$scale) + lgamma(1 + h) + lgamma(1 - h))
      },
      limited_moment = NULL,
      mgf = NULL,
      laplace = NULL
    ),
    degenerate = list(
      validate = function(params) {
        .severity_check_parameters(params, "value")
        .severity_check_number(params$value, "value", lower = 0)
      },
      support = function(params) {
        list(type = "discrete", lower = params$value, upper = params$value)
      },
      pdf = NULL,
      pmf = function(x, params, log = FALSE) {
        out <- as.numeric(x == params$value)
        if (log) log(out) else out
      },
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        out <- as.numeric(q >= params$value)
        if (!lower.tail) out <- 1 - out
        if (log.p) log(out) else out
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        q <- .severity_valid_prob(prob, log.p)
        out <- rep(params$value, length(q))
        out[is.na(q)] <- NaN
        out
      },
      rng = function(n, params) rep(params$value, n),
      raw_moment = function(order, params) {
        .severity_check_order(order)
        params$value^order
      },
      limited_moment = function(order, limit, params) {
        .severity_check_order(order)
        .severity_check_limit(limit)
        pmin(params$value, limit)^order
      },
      mgf = function(t, params) exp(t * params$value),
      laplace = function(t, params) {
        .severity_check_laplace_argument(t)
        exp(-t * params$value)
      }
    ),
    discrete = list(
      validate = function(params) {
        .severity_check_parameters(params, c("values", "probs"))
        if (!is.numeric(params$values) || !length(params$values) ||
            anyNA(params$values) || any(!is.finite(params$values)) ||
            any(params$values < 0) || any(diff(params$values) <= 0)) {
          stop("`values` must contain distinct, strictly increasing finite nonnegative amounts.",
               call. = FALSE)
        }
        if (!is.numeric(params$probs) || length(params$probs) != length(params$values) ||
            anyNA(params$probs) || any(!is.finite(params$probs)) ||
            any(params$probs <= 0) || any(params$probs > 1) ||
            abs(sum(params$probs) - 1) > 1e-12) {
          stop("`probs` must contain strictly positive probabilities, match `values`, and sum to one.",
               call. = FALSE)
        }
        invisible(TRUE)
      },
      support = function(params) {
        list(type = "discrete", lower = min(params$values), upper = max(params$values))
      },
      pdf = NULL,
      pmf = function(x, params, log = FALSE) {
        idx <- match(x, params$values)
        out <- numeric(length(x))
        found <- !is.na(idx)
        out[found] <- params$probs[idx[found]]
        out[is.na(x)] <- NA_real_
        if (log) log(out) else out
      },
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        ans <- vapply(q, function(y) {
          if (is.na(y)) return(NA_real_)
          sum(params$probs[params$values <= y])
        }, numeric(1))
        ans <- pmin(1, ans)
        ans[!is.na(q) & q >= max(params$values)] <- 1
        if (!lower.tail) ans <- 1 - ans
        if (log.p) log(ans) else ans
      },
      quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
        pp <- .severity_valid_prob(prob, log.p)
        if (log.p) pp <- exp(pp)
        if (!lower.tail) pp <- 1 - pp

        cumulative <- cumsum(params$probs)
        cumulative[length(cumulative)] <- 1

        vapply(pp, function(y) {
          if (is.na(y)) return(NaN)
          if (y <= 0) return(params$values[1L])
          if (y >= 1) return(params$values[length(params$values)])

          # Protect exact jump points such as exp(log(0.7)) against
          # floating-point roundoff in cumulative probabilities.
          tol <- 8 * .Machine$double.eps * max(1, abs(y))
          i <- which(cumulative >= y - tol)[1L]

          if (is.na(i)) {
            params$values[length(params$values)]
          } else {
            params$values[i]
          }
        }, numeric(1))
      },
      rng = function(n, params) {
        params$values[base::sample.int(length(params$values), size = n,
                                   replace = TRUE, prob = params$probs)]
      },
      raw_moment = function(order, params) {
        .severity_check_order(order)
        sum(params$probs * params$values^order)
      },
      limited_moment = function(order, limit, params) {
        .severity_check_order(order)
        .severity_check_limit(limit)
        sum(params$probs * pmin(params$values, limit)^order)
      },
      mgf = function(t, params) {
        vapply(t, function(z) sum(params$probs * exp(z * params$values)), numeric(1))
      },
      laplace = function(t, params) {
        .severity_check_laplace_argument(t)
        vapply(t, function(z) sum(params$probs * exp(-z * params$values)), numeric(1))
      }
    ),
    empirical = .empirical_severity_registry_entry(),
    custom = list(
      validate = NULL, support = NULL, pdf = NULL, pmf = NULL,
      cdf = NULL, quantile = NULL, rng = NULL, raw_moment = NULL,
      limited_moment = NULL, mgf = NULL, laplace = NULL
    )
  )
}


# Distribution formulas ----------------------------------------------------

.severity_dpareto2 <- function(x, shape, scale, log = FALSE) {
  ans <- rep(if (log) -Inf else 0, length(x))
  ok <- !is.na(x) & x >= 0
  logval <- log(shape) - log(scale) - (shape + 1) * log1p(x[ok] / scale)
  ans[ok] <- if (log) logval else exp(logval)
  ans[is.na(x)] <- NA_real_
  ans
}

.severity_ppareto2 <- function(q, shape, scale,
                               lower.tail = TRUE, log.p = FALSE) {
  ans <- numeric(length(q))
  ok <- !is.na(q) & q >= 0
  logsurv <- -shape * log1p(q[ok] / scale)
  ans[ok] <- if (lower.tail) -expm1(logsurv) else exp(logsurv)
  if (!lower.tail) ans[!ok & !is.na(q)] <- 1
  ans[is.na(q)] <- NA_real_
  if (log.p) log(ans) else ans
}

.severity_log1pexp <- function(z) {
  ans <- numeric(length(z))
  pos <- !is.na(z) & z > 0
  nonpos <- !is.na(z) & !pos
  ans[pos] <- z[pos] + log1p(exp(-z[pos]))
  ans[nonpos] <- log1p(exp(z[nonpos]))
  ans[is.na(z)] <- NA_real_
  ans
}

.severity_dburr12 <- function(x, shape1, shape2, scale, log = FALSE) {
  ans <- rep(if (log) -Inf else 0, length(x))
  pos <- !is.na(x) & is.finite(x) & x > 0
  z <- log(x[pos]) - log(scale)
  logval <- log(shape1) + log(shape2) - log(scale) +
    (shape1 - 1) * z - (shape2 + 1) * .severity_log1pexp(shape1 * z)
  ans[pos] <- if (log) logval else exp(logval)
  zero <- !is.na(x) & x == 0
  if (shape1 < 1) {
    ans[zero] <- Inf
  } else if (shape1 == 1) {
    ans[zero] <- if (log) log(shape2) - log(scale) else shape2 / scale
  }
  ans[is.na(x)] <- NA_real_
  ans
}

.severity_pburr12 <- function(q, shape1, shape2, scale,
                              lower.tail = TRUE, log.p = FALSE) {
  ans <- numeric(length(q))
  ok <- !is.na(q) & q >= 0
  u <- shape1 * (log(q[ok]) - log(scale))
  logsurv <- -shape2 * .severity_log1pexp(u)
  ans[ok] <- if (lower.tail) -expm1(logsurv) else exp(logsurv)
  if (!lower.tail) ans[!ok & !is.na(q)] <- 1
  ans[is.na(q)] <- NA_real_
  if (log.p) log(ans) else ans
}

.severity_qburr12 <- function(prob, shape1, shape2, scale,
                              lower.tail = TRUE, log.p = FALSE) {
  logsurv <- .severity_log_survival(prob, lower.tail, log.p)
  scale * expm1(-logsurv / shape2)^(1 / shape1)
}

.severity_log_survival <- function(prob, lower.tail = TRUE, log.p = FALSE) {
  p <- .severity_valid_prob(prob, log.p)
  if (log.p) {
    ans <- if (lower.tail) log1p(-exp(p)) else p
  } else {
    ans <- if (lower.tail) log1p(-p) else log(p)
  }
  ans[is.na(p)] <- NaN
  ans
}


# Validation ---------------------------------------------------------------

.severity_validate_dist <- function(dist) {
  if (!is.character(dist) || length(dist) != 1L ||
      is.na(dist) || !nzchar(trimws(dist))) {
    stop("`dist` must be a nonempty character scalar.", call. = FALSE)
  }
  tolower(trimws(dist))
}

.severity_validate_label <- function(label) {
  if (is.null(label)) return(invisible(TRUE))
  if (!is.character(label) || length(label) != 1L ||
      is.na(label) || !nzchar(trimws(label))) {
    stop("`label` must be NULL or a nonempty character scalar.", call. = FALSE)
  }
  invisible(TRUE)
}

.severity_validate_names <- function(x, what) {
  if (!length(x)) return(invisible(TRUE))
  nm <- names(x)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop(what, " must have nonempty, distinct names.", call. = FALSE)
  }
  invisible(TRUE)
}

.severity_check_parameters <- function(params, required) {
  missing <- setdiff(required, names(params))
  extra <- setdiff(names(params), required)
  if (length(missing)) {
    stop("Missing required parameter(s): ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  if (length(extra)) {
    stop("Unused parameter(s): ", paste(extra, collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

.severity_check_number <- function(x, name, lower = -Inf, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("`", name, "` must be a finite numeric scalar.", call. = FALSE)
  }
  if (if (strict) x <= lower else x < lower) {
    stop("`", name, "` is outside its allowed parameter domain.", call. = FALSE)
  }
  invisible(TRUE)
}

.severity_check_order <- function(order) {
  if (!is.numeric(order) || length(order) != 1L || is.na(order) ||
      !is.finite(order) || order < 0 || order != floor(order)) {
    stop("`order` must be a nonnegative integer scalar.", call. = FALSE)
  }
  invisible(TRUE)
}

.severity_check_limit <- function(limit) {
  if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) ||
      limit < 0 || !(is.finite(limit) || limit == Inf)) {
    stop("`limit` must be a nonnegative numeric scalar (possibly Inf).",
         call. = FALSE)
  }
  invisible(TRUE)
}

.severity_check_laplace_argument <- function(t) {
  if (!is.numeric(t) || anyNA(t) || any(t < 0) || any(!is.finite(t))) {
    stop("Laplace argument `t` must contain finite nonnegative numbers.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.severity_valid_prob <- function(prob, log.p = FALSE) {
  if (!is.numeric(prob)) stop("Probabilities must be numeric.", call. = FALSE)
  invalid <- !is.na(prob) & if (log.p) prob > 0 else (prob < 0 | prob > 1)
  out <- prob
  out[invalid] <- NaN
  out
}

.severity_validate_support <- function(support) {
  if (!is.list(support)) {
    stop("`support` must be a list.", call. = FALSE)
  }
  .severity_validate_names(support, "Support fields")
  required <- c("type", "lower", "upper")
  if (!all(required %in% names(support))) {
    stop("`support` must contain `type`, `lower`, and `upper`.", call. = FALSE)
  }
  if (!is.character(support$type) || length(support$type) != 1L ||
      is.na(support$type) ||
      !support$type %in% c("continuous", "discrete", "lattice", "mixed")) {
    stop("Invalid `support$type`.", call. = FALSE)
  }
  .severity_check_number(support$lower, "support$lower", lower = 0)
  upper <- support$upper
  if (!is.numeric(upper) || length(upper) != 1L || is.na(upper) ||
      !(is.finite(upper) || upper == Inf) || upper < support$lower) {
    stop("`support$upper` must be at least `support$lower` or Inf.",
         call. = FALSE)
  }
  if (identical(support$type, "lattice")) {
    if (is.null(support$step)) {
      stop("Lattice support requires `step`.", call. = FALSE)
    }
    .severity_check_number(support$step, "support$step", lower = 0, strict = TRUE)
  } else if ("step" %in% names(support)) {
    stop("`support$step` is only allowed for lattice support.", call. = FALSE)
  }
  keep <- if (identical(support$type, "lattice")) {
    c("type", "lower", "upper", "step")
  } else {
    required
  }
  extra <- setdiff(names(support), keep)
  if (length(extra)) {
    stop("Unsupported field(s) in `support`: ", paste(extra, collapse = ", "),
         ".", call. = FALSE)
  }
  support[keep]
}

.severity_validate_custom <- function(custom) {
  if (!is.list(custom) || !length(custom)) {
    stop("`custom` must be a nonempty named list.", call. = FALSE)
  }
  .severity_validate_names(custom, "Custom components")
  allowed <- c("pdf", "pmf", "cdf", "quantile", "rng", "raw_moment",
               "limited_moment", "mgf", "laplace", "support")
  extra <- setdiff(names(custom), allowed)
  if (length(extra)) {
    stop("Unsupported component(s) in `custom`: ", paste(extra, collapse = ", "),
         ".", call. = FALSE)
  }
  if (is.null(custom$support) ||
      (!is.function(custom$pdf) && !is.function(custom$pmf))) {
    stop("`custom` requires `support` and at least one function: `pdf` or `pmf`.",
         call. = FALSE)
  }
  functions <- setdiff(names(custom), "support")
  if (any(!vapply(custom[functions], is.function, logical(1)))) {
    stop("Every computational component in `custom` must be a function.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.severity_validate_custom_support <- function(custom, support) {
  type <- support$type
  if (identical(type, "continuous") &&
      (!is.function(custom$pdf) || is.function(custom$pmf))) {
    stop("Continuous custom severity requires `pdf` without `pmf`; use `mixed` for atoms.",
         call. = FALSE)
  }
  if (type %in% c("discrete", "lattice") &&
      (!is.function(custom$pmf) || is.function(custom$pdf))) {
    stop("Discrete custom severity requires `pmf` without `pdf`.", call. = FALSE)
  }
  if (identical(type, "mixed") &&
      (!is.function(custom$pdf) || !is.function(custom$pmf) ||
       !is.function(custom$cdf))) {
    stop("Mixed custom severity requires `pdf`, `pmf`, and `cdf`.",
         call. = FALSE)
  }
  invisible(TRUE)
}
