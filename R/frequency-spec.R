#' Specify a claim frequency distribution
#'
#' Creates a lightweight specification of a claim count distribution for use in
#' actuarial risk models.
#'
#' @param dist Character scalar identifying the frequency distribution.
#'   Currently supported values are `"poisson"`, `"binomial"`,
#'   `"negative_binomial"`, `"geometric"`, and `"custom"`.
#' @param ... Named parameters defining the selected distribution.
#' @param custom Optional named list defining a custom count distribution.
#'   Used only when `dist = "custom"`. At minimum, the list must contain
#'   `pmf`, a probability mass function.
#' @param label Optional character scalar used to identify the specification.
#'
#' @details
#' A frequency specification represents a nonnegative integer-valued random
#' variable. The constructor stores the probabilistic model only; it does not
#' simulate observations, estimate parameters, or calculate aggregate losses.
#'
#' Built-in parameterizations follow `stats`: Poisson uses `lambda`, binomial
#' uses `size` and `prob`, negative binomial uses `size` and `prob`, and
#' geometric uses `prob` with support beginning at zero.
#'
#' @return An object with classes `"tidyact_frequency_spec"` and
#'   `"tidyact_distribution_spec"`.
#'
#' @examples
#' frequency_spec("poisson", lambda = 4)
#' frequency_spec("binomial", size = 20, prob = 0.15)
#' frequency_spec("negative_binomial", size = 2.5, prob = 0.4)
#' frequency_spec("geometric", prob = 0.25)
#'
#' @export
frequency_spec <- function(dist, ..., custom = NULL, label = NULL) {
  dist <- .frequency_validate_dist(dist)
  params <- list(...)

  .frequency_validate_named_params(params)
  .frequency_validate_label(label)

  registry <- .frequency_registry()

  if (!dist %in% names(registry)) {
    supported <- paste(paste0('"', names(registry), '"'), collapse = ", ")
    stop("`dist` must be one of: ", supported, ".", call. = FALSE)
  }

  entry <- registry[[dist]]

  if (identical(dist, "empirical")) {
    stop(
      "Use `empirical_frequency()` to specify an empirical count distribution.",
      call. = FALSE
    )
  }

  if (identical(dist, "custom")) {
    entry$validate(params = params, custom = custom)
    support <- entry$support(params = params, custom = custom)
    custom <- .frequency_normalize_custom(custom)
  } else {
    if (!is.null(custom)) {
      stop(
        "`custom` can only be supplied when `dist = \"custom\"`.",
        call. = FALSE
      )
    }

    entry$validate(params = params, custom = NULL)
    support <- entry$support(params = params, custom = NULL)
  }

  support <- .frequency_validate_support(support)

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
    class = c("tidyact_frequency_spec", "tidyact_distribution_spec")
  )
}


#' Print a frequency distribution specification
#'
#' @param x A `tidyact_frequency_spec` object.
#' @param ... Additional arguments. Currently ignored.
#'
#' @return The input object, invisibly.
#' @export
print.tidyact_frequency_spec <- function(x, ...) {
  cat("<tidyact_frequency_spec>\n")
  cat("Distribution: ", .frequency_pretty_name(x$dist), "\n", sep = "")

  if (length(x$params) > 0L) {
    cat("Parameters:\n")
    for (name in names(x$params)) {
      cat("  ", name, ": ", .frequency_format_value(x$params[[name]]), "\n", sep = "")
    }
  }

  cat("Support: ", .frequency_format_support(x$support), "\n", sep = "")

  if (!is.null(x$label)) {
    cat("Label: ", x$label, "\n", sep = "")
  }

  invisible(x)
}


# Internal registry ---------------------------------------------------------

.frequency_registry <- function() {
  list(
    poisson = list(
      validate = function(params, custom = NULL) {
        .frequency_validate_parameter_names(params, required = "lambda")
        .frequency_assert_numeric_scalar(
          params$lambda, "lambda", lower = 0, lower_closed = TRUE
        )
        invisible(TRUE)
      },
      support = function(params, custom = NULL) {
        list(type = "count", lower = 0, upper = Inf, step = 1)
      },
      pmf = function(x, params, log = FALSE) {
        stats::dpois(x = x, lambda = params$lambda, log = log)
      },
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::ppois(q = q, lambda = params$lambda,
                     lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(p, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qpois(p = p, lambda = params$lambda,
                     lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) {
        stats::rpois(n = n, lambda = params$lambda)
      },
      raw_moment = function(order, params) {
        .frequency_raw_from_factorial(
          order = order,
          factorial_moment = function(k) params$lambda^k
        )
      },
      factorial_moment = function(order, params) {
        .frequency_validate_moment_order(order)
        params$lambda^order
      },
      pgf = function(z, params) {
        exp(params$lambda * (z - 1))
      },
      panjer = function(params) {
        list(
          a = 0,
          b = params$lambda,
          p0 = stats::dpois(0, lambda = params$lambda),
          p1 = stats::dpois(1, lambda = params$lambda)
        )
      }
    ),

    binomial = list(
      validate = function(params, custom = NULL) {
        .frequency_validate_parameter_names(params, required = c("size", "prob"))
        .frequency_assert_whole_number(params$size, "size", lower = 0)
        .frequency_assert_numeric_scalar(
          params$prob, "prob", lower = 0, upper = 1,
          lower_closed = TRUE, upper_closed = TRUE
        )
        invisible(TRUE)
      },
      support = function(params, custom = NULL) {
        list(type = "count", lower = 0, upper = params$size, step = 1)
      },
      pmf = function(x, params, log = FALSE) {
        stats::dbinom(x = x, size = params$size, prob = params$prob, log = log)
      },
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::pbinom(q = q, size = params$size, prob = params$prob,
                      lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(p, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qbinom(p = p, size = params$size, prob = params$prob,
                      lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) {
        stats::rbinom(n = n, size = params$size, prob = params$prob)
      },
      raw_moment = function(order, params) {
        .frequency_raw_from_factorial(
          order = order,
          factorial_moment = function(k) {
            .frequency_falling_factorial(params$size, k) * params$prob^k
          }
        )
      },
      factorial_moment = function(order, params) {
        .frequency_validate_moment_order(order)
        .frequency_falling_factorial(params$size, order) * params$prob^order
      },
      pgf = function(z, params) {
        (1 - params$prob + params$prob * z)^params$size
      },
      panjer = function(params) {
        p <- params$prob
        m <- params$size
        if (p == 1) return(NULL)
        denominator <- 1 - p
        list(
          a = -p / denominator,
          b = (m + 1) * p / denominator,
          p0 = stats::dbinom(0, size = m, prob = p),
          p1 = stats::dbinom(1, size = m, prob = p)
        )
      }
    ),

    negative_binomial = list(
      validate = function(params, custom = NULL) {
        .frequency_validate_parameter_names(params, required = c("size", "prob"))
        .frequency_assert_numeric_scalar(
          params$size, "size", lower = 0, lower_closed = FALSE
        )
        .frequency_assert_numeric_scalar(
          params$prob, "prob", lower = 0, upper = 1,
          lower_closed = FALSE, upper_closed = TRUE
        )
        invisible(TRUE)
      },
      support = function(params, custom = NULL) {
        list(type = "count", lower = 0, upper = Inf, step = 1)
      },
      pmf = function(x, params, log = FALSE) {
        stats::dnbinom(x = x, size = params$size, prob = params$prob, log = log)
      },
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::pnbinom(q = q, size = params$size, prob = params$prob,
                       lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(p, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qnbinom(p = p, size = params$size, prob = params$prob,
                       lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) {
        stats::rnbinom(n = n, size = params$size, prob = params$prob)
      },
      raw_moment = function(order, params) {
        ratio <- (1 - params$prob) / params$prob
        .frequency_raw_from_factorial(
          order = order,
          factorial_moment = function(k) {
            .frequency_rising_factorial(params$size, k) * ratio^k
          }
        )
      },
      factorial_moment = function(order, params) {
        .frequency_validate_moment_order(order)
        ratio <- (1 - params$prob) / params$prob
        .frequency_rising_factorial(params$size, order) * ratio^order
      },
      pgf = function(z, params) {
        p <- params$prob
        q <- 1 - p
        (p / (1 - q * z))^params$size
      },
      panjer = function(params) {
        p <- params$prob
        q <- 1 - p
        r <- params$size
        list(
          a = q,
          b = (r - 1) * q,
          p0 = stats::dnbinom(0, size = r, prob = p),
          p1 = stats::dnbinom(1, size = r, prob = p)
        )
      }
    ),

    geometric = list(
      validate = function(params, custom = NULL) {
        .frequency_validate_parameter_names(params, required = "prob")
        .frequency_assert_numeric_scalar(
          params$prob, "prob", lower = 0, upper = 1,
          lower_closed = FALSE, upper_closed = TRUE
        )
        invisible(TRUE)
      },
      support = function(params, custom = NULL) {
        list(type = "count", lower = 0, upper = Inf, step = 1)
      },
      pmf = function(x, params, log = FALSE) {
        stats::dgeom(x = x, prob = params$prob, log = log)
      },
      cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
        stats::pgeom(q = q, prob = params$prob,
                     lower.tail = lower.tail, log.p = log.p)
      },
      quantile = function(p, params, lower.tail = TRUE, log.p = FALSE) {
        stats::qgeom(p = p, prob = params$prob,
                     lower.tail = lower.tail, log.p = log.p)
      },
      rng = function(n, params) {
        stats::rgeom(n = n, prob = params$prob)
      },
      raw_moment = function(order, params) {
        ratio <- (1 - params$prob) / params$prob
        .frequency_raw_from_factorial(
          order = order,
          factorial_moment = function(k) factorial(k) * ratio^k
        )
      },
      factorial_moment = function(order, params) {
        .frequency_validate_moment_order(order)
        ratio <- (1 - params$prob) / params$prob
        factorial(order) * ratio^order
      },
      pgf = function(z, params) {
        p <- params$prob
        p / (1 - (1 - p) * z)
      },
      panjer = function(params) {
        p <- params$prob
        q <- 1 - p
        list(
          a = q,
          b = 0,
          p0 = stats::dgeom(0, prob = p),
          p1 = stats::dgeom(1, prob = p)
        )
      }
    ),

    # Empirical frequencies are constructed by empirical_frequency().
    # Their computational operations live in empirical-frequency.R.
    empirical = .empirical_frequency_registry_entry(),

    custom = list(
      validate = function(params, custom = NULL) {
        .frequency_validate_custom(custom)
        invisible(TRUE)
      },
      support = function(params, custom = NULL) {
        if (is.null(custom$support)) {
          return(list(type = "count", lower = 0, upper = Inf, step = 1))
        }
        if (is.function(custom$support)) {
          return(custom$support(params))
        }
        custom$support
      },
      pmf = NULL,
      cdf = NULL,
      quantile = NULL,
      rng = NULL,
      raw_moment = NULL,
      factorial_moment = NULL,
      pgf = NULL,
      panjer = NULL
    )
  )
}


# Validation helpers -------------------------------------------------------

.frequency_validate_dist <- function(dist) {
  if (!is.character(dist) || length(dist) != 1L || is.na(dist) ||
      !nzchar(trimws(dist))) {
    stop("`dist` must be a single non-missing character string.", call. = FALSE)
  }
  tolower(trimws(dist))
}

.frequency_validate_named_params <- function(params) {
  if (length(params) == 0L) return(invisible(TRUE))
  param_names <- names(params)
  if (is.null(param_names) || any(is.na(param_names)) || any(param_names == "")) {
    stop("All distribution parameters supplied through `...` must be named.",
         call. = FALSE)
  }
  if (anyDuplicated(param_names)) {
    duplicated_names <- unique(param_names[duplicated(param_names)])
    stop(
      "Distribution parameter names must be unique. Duplicated: ",
      paste(duplicated_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.frequency_validate_label <- function(label) {
  if (is.null(label)) return(invisible(TRUE))
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(trimws(label))) {
    stop("`label` must be NULL or a single non-empty character string.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.frequency_validate_parameter_names <- function(params, required,
                                                optional = character()) {
  supplied <- names(params)
  allowed <- c(required, optional)
  missing <- setdiff(required, supplied)
  extra <- setdiff(supplied, allowed)

  if (length(missing) > 0L) {
    stop(
      "Missing required parameter", if (length(missing) > 1L) "s" else "",
      ": ", paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (length(extra) > 0L) {
    stop(
      "Unused parameter", if (length(extra) > 1L) "s" else "",
      " for this distribution: ", paste(extra, collapse = ", "), ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.frequency_assert_numeric_scalar <- function(x, name, lower = -Inf,
                                             upper = Inf,
                                             lower_closed = TRUE,
                                             upper_closed = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("`", name, "` must be a finite numeric scalar.", call. = FALSE)
  }

  lower_ok <- if (lower_closed) x >= lower else x > lower
  upper_ok <- if (upper_closed) x <= upper else x < upper

  if (!lower_ok || !upper_ok) {
    lower_symbol <- if (lower_closed) "[" else "("
    upper_symbol <- if (upper_closed) "]" else ")"
    stop(
      "`", name, "` must lie in ", lower_symbol, lower, ", ", upper,
      upper_symbol, ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.frequency_assert_whole_number <- function(x, name, lower = 0) {
  .frequency_assert_numeric_scalar(
    x = x, name = name, lower = lower, lower_closed = TRUE
  )
  tolerance <- sqrt(.Machine$double.eps)
  if (abs(x - round(x)) > tolerance) {
    stop("`", name, "` must be an integer-valued scalar.", call. = FALSE)
  }
  invisible(TRUE)
}

.frequency_validate_support <- function(support) {
  if (!is.list(support)) {
    stop("Frequency support must be supplied as a list.", call. = FALSE)
  }

  if (is.null(support$type)) support$type <- "count"
  required <- c("type", "lower", "upper", "step")
  missing <- setdiff(required, names(support))

  if (length(missing) > 0L) {
    stop("Frequency support is missing: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }

  if (!is.character(support$type) || length(support$type) != 1L ||
      is.na(support$type) || support$type != "count") {
    stop("Frequency support must have `type = \"count\"`.", call. = FALSE)
  }

  .frequency_assert_whole_number(support$lower, "support$lower", lower = 0)

  if (!is.numeric(support$upper) || length(support$upper) != 1L ||
      is.na(support$upper)) {
    stop("`support$upper` must be a numeric scalar.", call. = FALSE)
  }

  if (is.finite(support$upper)) {
    .frequency_assert_whole_number(
      support$upper, "support$upper", lower = support$lower
    )
  } else if (support$upper != Inf) {
    stop("`support$upper` must be finite or `Inf`.", call. = FALSE)
  }

  .frequency_assert_whole_number(support$step, "support$step", lower = 1)

  if (is.finite(support$upper)) {
    span <- (support$upper - support$lower) / support$step
    tolerance <- sqrt(.Machine$double.eps)
    if (abs(span - round(span)) > tolerance) {
      stop(
        "`support$upper` must lie on the support grid defined by ",
        "`support$lower` and `support$step`.",
        call. = FALSE
      )
    }
  }

  list(
    type = "count",
    lower = as.numeric(support$lower),
    upper = as.numeric(support$upper),
    step = as.numeric(support$step)
  )
}

.frequency_validate_custom <- function(custom) {
  if (!is.list(custom)) {
    stop("For `dist = \"custom\"`, `custom` must be a named list.",
         call. = FALSE)
  }

  if (is.null(names(custom)) || any(names(custom) == "") ||
      any(is.na(names(custom)))) {
    stop("`custom` must be a named list.", call. = FALSE)
  }

  if (anyDuplicated(names(custom))) {
    stop("Names in `custom` must be unique.", call. = FALSE)
  }

  allowed <- c(
    "pmf", "cdf", "quantile", "rng", "raw_moment", "factorial_moment",
    "pgf", "panjer", "support"
  )

  extra <- setdiff(names(custom), allowed)
  if (length(extra) > 0L) {
    stop(
      "Unsupported component", if (length(extra) > 1L) "s" else "",
      " in `custom`: ", paste(extra, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (is.null(custom$pmf) || !is.function(custom$pmf)) {
    stop("`custom$pmf` must be supplied as a function.", call. = FALSE)
  }

  function_components <- intersect(
    names(custom),
    c("pmf", "cdf", "quantile", "rng", "raw_moment",
      "factorial_moment", "pgf", "panjer")
  )

  for (name in function_components) {
    if (!is.function(custom[[name]])) {
      stop("`custom$", name, "` must be a function.", call. = FALSE)
    }
  }

  if (!is.null(custom$support) && !is.list(custom$support) &&
      !is.function(custom$support)) {
    stop("`custom$support` must be a list or a function of `params`.",
         call. = FALSE)
  }

  invisible(TRUE)
}

.frequency_normalize_custom <- function(custom) {
  custom$support <- NULL
  custom
}


# Moment helpers -----------------------------------------------------------

.frequency_validate_moment_order <- function(order) {
  if (!is.numeric(order) || length(order) != 1L || is.na(order) ||
      !is.finite(order) || order < 0) {
    stop("`order` must be a nonnegative integer.", call. = FALSE)
  }
  tolerance <- sqrt(.Machine$double.eps)
  if (abs(order - round(order)) > tolerance) {
    stop("`order` must be a nonnegative integer.", call. = FALSE)
  }
  invisible(TRUE)
}

.frequency_falling_factorial <- function(x, order) {
  .frequency_validate_moment_order(order)
  order <- as.integer(round(order))
  if (order == 0L) return(1)
  if (order > x) return(0)
  exp(lgamma(x + 1) - lgamma(x - order + 1))
}

.frequency_rising_factorial <- function(x, order) {
  .frequency_validate_moment_order(order)
  order <- as.integer(round(order))
  if (order == 0L) return(1)
  exp(lgamma(x + order) - lgamma(x))
}

.frequency_stirling2 <- function(n, k) {
  n <- as.integer(n)
  k <- as.integer(k)

  if (n == 0L && k == 0L) return(1)
  if (n == 0L || k == 0L || k > n) return(0)

  values <- matrix(0, nrow = n + 1L, ncol = k + 1L)
  values[1L, 1L] <- 1

  for (i in seq_len(n)) {
    max_j <- min(i, k)
    for (j in seq_len(max_j)) {
      previous_same <- values[i, j + 1L]
      previous_lower <- values[i, j]
      values[i + 1L, j + 1L] <- j * previous_same + previous_lower
    }
  }

  values[n + 1L, k + 1L]
}

.frequency_raw_from_factorial <- function(order, factorial_moment) {
  .frequency_validate_moment_order(order)
  order <- as.integer(round(order))
  if (order == 0L) return(1)

  result <- 0
  for (k in seq_len(order)) {
    result <- result +
      .frequency_stirling2(order, k) * factorial_moment(k)
  }
  result
}


# Print helpers ------------------------------------------------------------

.frequency_pretty_name <- function(dist) {
  switch(
    dist,
    poisson = "Poisson",
    binomial = "Binomial",
    negative_binomial = "Negative binomial",
    geometric = "Geometric",
    empirical = "Empirical",
    custom = "Custom",
    dist
  )
}

.frequency_format_value <- function(x) {
  if (is.null(x)) return("NULL")

  if (is.atomic(x) && length(x) <= 6L) {
    return(paste(format(x, trim = TRUE), collapse = ", "))
  }

  if (is.atomic(x)) {
    preview <- paste(format(x[seq_len(5L)], trim = TRUE), collapse = ", ")
    return(paste0(preview, ", ..."))
  }

  paste0("<", paste(class(x), collapse = "/"), ">")
}

.frequency_format_support <- function(support) {
  lower <- support$lower
  upper <- support$upper
  step <- support$step

  if (is.infinite(upper) && step == 1) {
    return(paste0(
      "{", format(lower, trim = TRUE), ", ",
      format(lower + 1, trim = TRUE), ", ",
      format(lower + 2, trim = TRUE), ", ...}"
    ))
  }

  if (is.infinite(upper)) {
    return(paste0(
      "{", format(lower, trim = TRUE), ", ",
      format(lower + step, trim = TRUE), ", ",
      format(lower + 2 * step, trim = TRUE), ", ...}"
    ))
  }

  if (lower == upper) {
    return(paste0("{", format(lower, trim = TRUE), "}"))
  }

  if (upper == lower + step) {
    return(paste0(
      "{", format(lower, trim = TRUE), ", ",
      format(upper, trim = TRUE), "}"
    ))
  }

  paste0(
    "{", format(lower, trim = TRUE), ", ",
    format(lower + step, trim = TRUE), ", ..., ",
    format(upper, trim = TRUE), "}"
  )
}
