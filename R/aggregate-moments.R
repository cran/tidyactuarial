#' Calculate moments of a collective aggregate loss
#'
#' @md
#' @description
#' Calculates analytical raw moments, central moments, or cumulants of a
#' collective risk model
#'
#' \deqn{S = \sum_{i=1}^{N} X_i,}
#'
#' using the moments of the claim count distribution and the individual
#' severity distribution stored in a `tidyact_collective_risk_spec` object.
#'
#' No simulation, numerical integration, discretization, or approximation is
#' performed by this function.
#'
#' @param x A `tidyact_collective_risk_spec` object created by
#'   [collective_risk_spec()].
#' @param order Nonempty numeric vector of nonnegative integer moment orders.
#'   Cumulant orders must be positive.
#' @param type One of `"raw"` (default), `"central"`, or `"cumulant"`.
#'
#' @details
#' Let
#'
#' \deqn{m_j = E[X^j]}
#'
#' denote the raw severity moments and
#'
#' \deqn{\nu_k = E[(N)_k]}
#'
#' the factorial moments of the claim count. For positive integer \eqn{r},
#' the raw aggregate moment is evaluated through the random-sum identity
#'
#' \deqn{
#' E[S^r] = \sum_{k=1}^{r} \nu_k
#' B_{r,k}(m_1,\ldots,m_{r-k+1}),
#' }
#'
#' where \eqn{B_{r,k}} is a partial exponential Bell polynomial.
#'
#' This formulation includes the familiar identities
#'
#' \deqn{E[S] = E[N]E[X]}
#'
#' and
#'
#' \deqn{
#' \operatorname{Var}(S)
#' = E[N]\operatorname{Var}(X)
#' + \operatorname{Var}(N)(E[X])^2.
#' }
#'
#' Central moments and cumulants are derived algebraically from aggregate raw
#' moments. For a compound Poisson model, cumulants use the direct identity
#'
#' \deqn{\kappa_r(S) = \lambda E[X^r]}
#'
#' whenever the corresponding severity moment is finite.
#'
#' The function respects the status information returned by
#' [distribution_moments()]. Results are classified as:
#'
#' - `"finite"`: a finite numerical result;
#' - `"infinite"`: a mathematically divergent raw or central moment;
#' - `"undefined"`: the requested quantity cannot be defined from the
#'   available moments;
#' - `"overflow"`: a theoretically finite calculation exceeded floating-point
#'   range.
#'
#' Because both frequency and severity are nonnegative in this module, a
#' finite first moment equal to zero identifies a degenerate-at-zero component.
#' If either \eqn{N=0} almost surely or \eqn{X=0} almost surely, then the
#' aggregate loss is degenerate at zero and all positive-order aggregate
#' moments and cumulants are zero.
#'
#' @return A [tibble::tibble()] with one row per requested order and columns:
#'   `type`, `order`, `value`, and `status`.
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and
#' Nesbitt, C. J. (1997). *Actuarial Mathematics*, 2nd ed.
#' Society of Actuaries.
#'
#' Klugman, S. A., Panjer, H. H., and Willmot, G. E. (2019).
#' *Loss Models: From Data to Decisions*, 5th ed. Wiley.
#'
#' @seealso [collective_risk_spec()], [distribution_moments()]
#'
#' @examples
#' frequency <- frequency_spec(
#'   "poisson",
#'   lambda = 4
#' )
#'
#' severity <- severity_spec(
#'   "gamma",
#'   shape = 2,
#'   rate = 0.001
#' )
#'
#' model <- collective_risk_spec(
#'   frequency,
#'   severity
#' )
#'
#' aggregate_moments(model)
#'
#' aggregate_moments(
#'   model,
#'   order = 1:4,
#'   type = "central"
#' )
#'
#' aggregate_moments(
#'   model,
#'   order = 1:4,
#'   type = "cumulant"
#' )
#'
#' @export
aggregate_moments <- function(
  x,
  order = 1:4,
  type = "raw"
) {
  if (!inherits(x, "tidyact_collective_risk_spec")) {
    stop(
      "`x` must be a `tidyact_collective_risk_spec` object.",
      call. = FALSE
    )
  }

  .aggregate_validate_order(order)

  allowed <- c(
    "raw",
    "central",
    "cumulant"
  )

  if (
    !is.character(type) ||
      length(type) != 1L ||
      is.na(type) ||
      !type %in% allowed
  ) {
    stop(
      "`type` must be one of: raw, central, cumulant.",
      call. = FALSE
    )
  }

  if (identical(type, "cumulant") && any(order == 0)) {
    stop(
      "Cumulant orders must be positive integers.",
      call. = FALSE
    )
  }

  order <- as.integer(order)
  max_order <- max(order)

  raw <- .aggregate_raw_moment_sequence(
    x,
    max_order = max_order
  )

  if (identical(type, "raw")) {
    results <- raw[order + 1L]
  } else if (identical(type, "central")) {
    central <- .aggregate_central_moment_sequence(
      raw,
      max_order = max_order
    )

    results <- central[order + 1L]
  } else {
    results <- .aggregate_cumulant_sequence(
      x,
      raw,
      max_order = max_order
    )[order]
  }

  tibble::tibble(
    type = rep(type, length(order)),
    order = order,
    value = vapply(
      results,
      function(z) z$value,
      numeric(1)
    ),
    status = vapply(
      results,
      function(z) z$status,
      character(1)
    )
  )
}


# Validation ---------------------------------------------------------------

.aggregate_validate_order <- function(order) {
  if (
    !is.numeric(order) ||
      length(order) == 0L ||
      anyNA(order) ||
      any(!is.finite(order)) ||
      any(order < 0) ||
      any(order != floor(order)) ||
      any(order > .Machine$integer.max)
  ) {
    stop(
      "`order` must be a nonempty vector of nonnegative integers.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# Internal result objects --------------------------------------------------

.aggregate_result <- function(
  value,
  status = "finite"
) {
  list(
    value = as.numeric(value),
    status = status
  )
}


.aggregate_undefined <- function() {
  .aggregate_result(
    NA_real_,
    "undefined"
  )
}


.aggregate_overflow <- function() {
  .aggregate_result(
    Inf,
    "overflow"
  )
}


.aggregate_infinite <- function() {
  .aggregate_result(
    Inf,
    "infinite"
  )
}


# Status-aware nonnegative arithmetic -------------------------------------

.aggregate_scale_nonnegative <- function(
  x,
  coefficient
) {
  if (!is.numeric(coefficient) || length(coefficient) != 1L ||
      is.na(coefficient) || coefficient < 0) {
    stop(
      "Internal aggregate-moment coefficients must be nonnegative.",
      call. = FALSE
    )
  }

  if (coefficient == 0) {
    return(.aggregate_result(0))
  }

  if (identical(x$status, "infinite")) {
    return(.aggregate_infinite())
  }

  if (identical(x$status, "undefined")) {
    return(.aggregate_undefined())
  }

  if (identical(x$status, "overflow")) {
    return(.aggregate_overflow())
  }

  value <- coefficient * x$value

  if (!is.finite(value)) {
    return(.aggregate_overflow())
  }

  .aggregate_result(value)
}


.aggregate_multiply_nonnegative <- function(
  x,
  y
) {
  if (
    identical(x$status, "finite") &&
      identical(x$value, 0)
  ) {
    return(.aggregate_result(0))
  }

  if (
    identical(y$status, "finite") &&
      identical(y$value, 0)
  ) {
    return(.aggregate_result(0))
  }

  if (
    identical(x$status, "infinite") ||
      identical(y$status, "infinite")
  ) {
    return(.aggregate_infinite())
  }

  if (
    identical(x$status, "undefined") ||
      identical(y$status, "undefined")
  ) {
    return(.aggregate_undefined())
  }

  if (
    identical(x$status, "overflow") ||
      identical(y$status, "overflow")
  ) {
    return(.aggregate_overflow())
  }

  value <- x$value * y$value

  if (!is.finite(value)) {
    return(.aggregate_overflow())
  }

  .aggregate_result(value)
}


.aggregate_sum_nonnegative <- function(parts) {
  if (!length(parts)) {
    return(.aggregate_result(0))
  }

  statuses <- vapply(
    parts,
    function(z) z$status,
    character(1)
  )

  if (any(statuses == "infinite")) {
    return(.aggregate_infinite())
  }

  if (any(statuses == "undefined")) {
    return(.aggregate_undefined())
  }

  if (any(statuses == "overflow")) {
    return(.aggregate_overflow())
  }

  value <- sum(
    vapply(
      parts,
      function(z) z$value,
      numeric(1)
    )
  )

  if (!is.finite(value)) {
    return(.aggregate_overflow())
  }

  .aggregate_result(value)
}


# Component moments --------------------------------------------------------

.aggregate_component_moments <- function(
  x,
  max_order
) {
  frequency <- distribution_moments(
    x$frequency,
    order = 0:max_order,
    type = "factorial"
  )

  severity <- distribution_moments(
    x$severity,
    order = 0:max_order,
    type = "raw"
  )

  list(
    frequency = lapply(
      seq_len(nrow(frequency)),
      function(i) {
        .aggregate_result(
          frequency$value[[i]],
          frequency$status[[i]]
        )
      }
    ),
    severity = lapply(
      seq_len(nrow(severity)),
      function(i) {
        .aggregate_result(
          severity$value[[i]],
          severity$status[[i]]
        )
      }
    )
  )
}


.aggregate_is_zero_component <- function(moment_sequence) {
  if (length(moment_sequence) < 2L) {
    return(FALSE)
  }

  first <- moment_sequence[[2L]]

  identical(first$status, "finite") &&
    identical(first$value, 0)
}


# Bell-polynomial machinery ------------------------------------------------

.aggregate_bell_table <- function(
  severity,
  max_order
) {
  bell <- vector(
    "list",
    max_order + 1L
  )

  for (n in 0:max_order) {
    bell[[n + 1L]] <- lapply(
      seq_len(n + 1L),
      function(i) .aggregate_result(0)
    )
  }

  bell[[1L]][[1L]] <- .aggregate_result(1)

  if (max_order == 0L) {
    return(bell)
  }

  for (n in seq_len(max_order)) {
    for (k in seq_len(n)) {
      parts <- vector(
        "list",
        n - k + 1L
      )

      for (j in seq_len(n - k + 1L)) {
        severity_j <- severity[[j + 1L]]
        previous <- bell[[n - j + 1L]][[k]]

        term <- .aggregate_multiply_nonnegative(
          severity_j,
          previous
        )

        parts[[j]] <- .aggregate_scale_nonnegative(
          term,
          choose(n - 1L, j - 1L)
        )
      }

      bell[[n + 1L]][[k + 1L]] <-
        .aggregate_sum_nonnegative(parts)
    }
  }

  bell
}


# Raw aggregate moments ----------------------------------------------------

.aggregate_raw_moment_sequence <- function(
  x,
  max_order
) {
  components <- .aggregate_component_moments(
    x,
    max_order
  )

  raw <- vector(
    "list",
    max_order + 1L
  )

  raw[[1L]] <- .aggregate_result(1)

  if (max_order == 0L) {
    return(raw)
  }

  if (
    .aggregate_is_zero_component(components$frequency) ||
      .aggregate_is_zero_component(components$severity)
  ) {
    for (n in seq_len(max_order)) {
      raw[[n + 1L]] <- .aggregate_result(0)
    }

    return(raw)
  }

  bell <- .aggregate_bell_table(
    components$severity,
    max_order
  )

  for (n in seq_len(max_order)) {
    terms <- vector(
      "list",
      n
    )

    for (k in seq_len(n)) {
      terms[[k]] <- .aggregate_multiply_nonnegative(
        components$frequency[[k + 1L]],
        bell[[n + 1L]][[k + 1L]]
      )
    }

    raw[[n + 1L]] <- .aggregate_sum_nonnegative(terms)
  }

  raw
}


# Central aggregate moments ------------------------------------------------

.aggregate_central_moment_sequence <- function(
  raw,
  max_order
) {
  central <- vector(
    "list",
    max_order + 1L
  )

  central[[1L]] <- .aggregate_result(1)

  if (max_order == 0L) {
    return(central)
  }

  mean <- raw[[2L]]

  if (identical(mean$status, "finite")) {
    central[[2L]] <- .aggregate_result(0)
  } else if (identical(mean$status, "overflow")) {
    central[[2L]] <- .aggregate_overflow()
  } else {
    central[[2L]] <- .aggregate_undefined()
  }

  if (max_order == 1L) {
    return(central)
  }

  if (!identical(mean$status, "finite")) {
    fill <- if (identical(mean$status, "overflow")) {
      .aggregate_overflow
    } else {
      .aggregate_undefined
    }

    for (n in 2:max_order) {
      central[[n + 1L]] <- fill()
    }

    return(central)
  }

  for (n in 2:max_order) {
    higher <- raw[[n + 1L]]

    if (identical(higher$status, "infinite")) {
      central[[n + 1L]] <- .aggregate_infinite()
      next
    }

    if (identical(higher$status, "undefined")) {
      central[[n + 1L]] <- .aggregate_undefined()
      next
    }

    if (identical(higher$status, "overflow")) {
      central[[n + 1L]] <- .aggregate_overflow()
      next
    }

    needed <- raw[seq_len(n + 1L)]
    statuses <- vapply(
      needed,
      function(z) z$status,
      character(1)
    )

    if (any(statuses == "undefined")) {
      central[[n + 1L]] <- .aggregate_undefined()
      next
    }

    if (any(statuses == "infinite")) {
      central[[n + 1L]] <- .aggregate_infinite()
      next
    }

    if (any(statuses == "overflow")) {
      central[[n + 1L]] <- .aggregate_overflow()
      next
    }

    k <- 0:n
    raw_values <- vapply(
      needed,
      function(z) z$value,
      numeric(1)
    )

    terms <- choose(n, k) *
      (-mean$value)^(n - k) *
      raw_values

    value <- sum(terms)

    if (!is.finite(value)) {
      central[[n + 1L]] <- .aggregate_overflow()
      next
    }

    if (n %% 2L == 0L && value < 0) {
      tolerance <- 1e-12 * max(
        1,
        sum(abs(terms))
      )

      if (abs(value) <= tolerance) {
        value <- 0
      } else {
        central[[n + 1L]] <- .aggregate_undefined()
        next
      }
    }

    central[[n + 1L]] <- .aggregate_result(value)
  }

  central
}


# Aggregate cumulants ------------------------------------------------------

.aggregate_cumulant_sequence <- function(
  x,
  raw,
  max_order
) {
  cumulants <- vector(
    "list",
    max_order
  )

  if (max_order == 0L) {
    return(cumulants)
  }

  if (
    identical(x$frequency$dist, "poisson") &&
      is.numeric(x$frequency$params$lambda) &&
      length(x$frequency$params$lambda) == 1L
  ) {
    lambda <- x$frequency$params$lambda

    if (lambda == 0) {
      for (n in seq_len(max_order)) {
        cumulants[[n]] <- .aggregate_result(0)
      }

      return(cumulants)
    }

    severity <- distribution_moments(
      x$severity,
      order = seq_len(max_order),
      type = "raw"
    )

    for (n in seq_len(max_order)) {
      status <- severity$status[[n]]
      value <- severity$value[[n]]

      if (identical(status, "finite")) {
        answer <- lambda * value

        cumulants[[n]] <- if (is.finite(answer)) {
          .aggregate_result(answer)
        } else {
          .aggregate_overflow()
        }
      } else if (identical(status, "overflow")) {
        cumulants[[n]] <- .aggregate_overflow()
      } else if (
        identical(status, "infinite") &&
          n == 1L
      ) {
        cumulants[[n]] <- .aggregate_infinite()
      } else {
        cumulants[[n]] <- .aggregate_undefined()
      }
    }

    return(cumulants)
  }

  numeric_cumulants <- numeric(max_order)

  for (n in seq_len(max_order)) {
    if (n == 1L) {
      cumulants[[n]] <- raw[[2L]]

      if (identical(cumulants[[n]]$status, "finite")) {
        numeric_cumulants[[n]] <- cumulants[[n]]$value
      }

      next
    }

    needed <- raw[seq_len(n + 1L)]
    statuses <- vapply(
      needed,
      function(z) z$status,
      character(1)
    )

    if (any(statuses == "undefined") || any(statuses == "infinite")) {
      cumulants[[n]] <- .aggregate_undefined()
      next
    }

    if (any(statuses == "overflow")) {
      cumulants[[n]] <- .aggregate_overflow()
      next
    }

    if (any(vapply(
      cumulants[seq_len(n - 1L)],
      is.null,
      logical(1)
    ))) {
      cumulants[[n]] <- .aggregate_undefined()
      next
    }

    previous_status <- vapply(
      cumulants[seq_len(n - 1L)],
      function(z) z$status,
      character(1)
    )

    if (any(previous_status != "finite")) {
      cumulants[[n]] <- if (any(previous_status == "overflow")) {
        .aggregate_overflow()
      } else {
        .aggregate_undefined()
      }

      next
    }

    moments <- vapply(
      needed,
      function(z) z$value,
      numeric(1)
    )

    kappa <- moments[[n + 1L]]

    for (j in seq_len(n - 1L)) {
      kappa <- kappa -
        choose(n - 1L, j - 1L) *
        numeric_cumulants[[j]] *
        moments[[n - j + 1L]]
    }

    if (!is.finite(kappa)) {
      cumulants[[n]] <- .aggregate_overflow()
      next
    }

    numeric_cumulants[[n]] <- kappa
    cumulants[[n]] <- .aggregate_result(kappa)
  }

  cumulants
}
