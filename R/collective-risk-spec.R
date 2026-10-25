#' Specify a collective risk model
#'
#' @md
#' @description
#' Creates a lightweight specification of the classical collective risk model
#'
#' \deqn{S = \sum_{i=1}^{N} X_i,}
#'
#' where \eqn{N} is a claim count random variable and \eqn{X_1, X_2, \ldots}
#' are identically distributed individual claim severities. Independence between
#' the frequency variable and the severities, and mutual independence of the
#' severities, are intrinsic assumptions of this specification.
#'
#' The constructor records the probabilistic model only. It does not calculate
#' aggregate moments, construct the distribution of \eqn{S}, discretize a
#' severity distribution, simulate losses, or apply an approximation.
#'
#' @param frequency A `tidyact_frequency_spec` object, typically created with
#'   [frequency_spec()] or [empirical_frequency()].
#' @param severity A `tidyact_severity_spec` object, typically created with
#'   [severity_spec()] or [empirical_severity()].
#' @param count_basis Either `NULL`, `"loss"`, or `"payment"`. The count basis
#'   indicates whether `frequency` counts underlying losses or positive
#'   payments. When `NULL`, the constructor uses basis information already
#'   carried by the component specifications when available; otherwise the
#'   default is `"loss"`.
#' @param label Optional nonempty character scalar identifying the model.
#'
#' @details
#' The collective model is represented by the random sum
#'
#' \deqn{S = X_1 + \cdots + X_N,}
#'
#' with the convention \eqn{S = 0} when \eqn{N = 0}.
#'
#' `collective_risk_spec()` deliberately keeps the model specification separate
#' from computational methods. In particular, no Panjer recursion, convolution,
#' FFT calculation, Monte Carlo simulation, numerical discretization, or normal
#' approximation is performed by this function.
#'
#' The `count_basis` field is part of the model because frequency and severity
#' must refer to the same event basis. A loss-basis specification counts
#' underlying losses, including losses that may later produce a zero payment
#' after coverage is applied. A payment-basis specification counts positive
#' payments. Conversion between these bases is a separate actuarial operation
#' and is not performed silently by this constructor.
#'
#' The constructor does not expose an `independence` argument. Independence is
#' part of the definition of this version of the collective risk model rather
#' than an optional computational setting.
#'
#' @return An object with classes `"tidyact_collective_risk_spec"` and
#'   `"tidyact_risk_spec"`. It contains the following fields:
#'
#' - `schema_version`: internal schema version.
#' - `frequency`: the frequency specification.
#' - `severity`: the severity specification.
#' - `count_basis`: either `"loss"` or `"payment"`.
#' - `label`: optional user label.
#' - `metadata`: reserved metadata list.
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and
#' Nesbitt, C. J. (1997). *Actuarial Mathematics*, 2nd ed.
#' Society of Actuaries.
#'
#' @seealso [frequency_spec()], [severity_spec()], [empirical_frequency()],
#'   [empirical_severity()], [distribution_moments()]
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
#' model
#'
#' # An explicit payment-basis model.
#' collective_risk_spec(
#'   frequency,
#'   severity,
#'   count_basis = "payment",
#'   label = "Illustrative payment model"
#' )
#'
#' # Empirical components can be used without fitting a parametric model.
#' empirical_n <- empirical_frequency(c(0, 0, 1, 1, 2, 3))
#' empirical_x <- empirical_severity(c(500, 750, 1000, 1500, 3000))
#'
#' collective_risk_spec(
#'   empirical_n,
#'   empirical_x
#' )
#'
#' @export
collective_risk_spec <- function(
  frequency,
  severity,
  count_basis = NULL,
  label = NULL
) {
  .collective_validate_frequency(frequency)
  .collective_validate_severity(severity)
  .collective_validate_label(label)

  if (is.null(count_basis)) {
    count_basis <- .collective_infer_count_basis(
      frequency = frequency,
      severity = severity
    )
  } else {
    count_basis <- .collective_validate_count_basis(count_basis)
  }

  structure(
    list(
      schema_version = 1L,
      frequency = frequency,
      severity = severity,
      count_basis = count_basis,
      label = label,
      metadata = list()
    ),
    class = c(
      "tidyact_collective_risk_spec",
      "tidyact_risk_spec"
    )
  )
}


#' Print a collective risk model specification
#'
#' @param x A `tidyact_collective_risk_spec` object.
#' @param ... Additional arguments. Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.tidyact_collective_risk_spec <- function(x, ...) {
  cat("<tidyact_collective_risk_spec>\n")
  cat("Model: S = sum(X_i, i = 1, ..., N)\n")
  cat(
    "Frequency: ",
    .collective_distribution_name(x$frequency),
    "\n",
    sep = ""
  )
  cat(
    "Severity: ",
    .collective_distribution_name(x$severity),
    "\n",
    sep = ""
  )
  cat(
    "Count basis: ",
    x$count_basis,
    "\n",
    sep = ""
  )

  if (!is.null(x$label)) {
    cat(
      "Label: ",
      x$label,
      "\n",
      sep = ""
    )
  }

  invisible(x)
}


# Validation ---------------------------------------------------------------

.collective_validate_frequency <- function(x) {
  if (!inherits(x, "tidyact_frequency_spec")) {
    stop(
      "`frequency` must be a `tidyact_frequency_spec` object.",
      call. = FALSE
    )
  }

  required <- c(
    "schema_version",
    "dist",
    "params",
    "support",
    "data",
    "custom",
    "label",
    "metadata"
  )

  if (!is.list(x) || !all(required %in% names(x))) {
    stop(
      "`frequency` is not a valid frequency specification.",
      call. = FALSE
    )
  }

  if (
    !is.list(x$support) ||
      !identical(x$support$type, "count")
  ) {
    stop(
      "`frequency` must have count support.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.collective_validate_severity <- function(x) {
  if (!inherits(x, "tidyact_severity_spec")) {
    stop(
      "`severity` must be a `tidyact_severity_spec` object.",
      call. = FALSE
    )
  }

  required <- c(
    "schema_version",
    "dist",
    "params",
    "support",
    "data",
    "custom",
    "label",
    "metadata"
  )

  if (!is.list(x) || !all(required %in% names(x))) {
    stop(
      "`severity` is not a valid severity specification.",
      call. = FALSE
    )
  }

  if (
    !is.list(x$support) ||
      is.null(x$support$type)
  ) {
    stop(
      "`severity` does not contain a valid support specification.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(x$support$lower) ||
      length(x$support$lower) != 1L ||
      is.na(x$support$lower) ||
      x$support$lower < 0
  ) {
    stop(
      "`severity` must have nonnegative support.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.collective_validate_count_basis <- function(count_basis) {
  if (
    !is.character(count_basis) ||
      length(count_basis) != 1L ||
      is.na(count_basis)
  ) {
    stop(
      "`count_basis` must be NULL, \"loss\", or \"payment\".",
      call. = FALSE
    )
  }

  count_basis <- tolower(trimws(count_basis))

  if (!count_basis %in% c("loss", "payment")) {
    stop(
      "`count_basis` must be NULL, \"loss\", or \"payment\".",
      call. = FALSE
    )
  }

  count_basis
}


.collective_validate_label <- function(label) {
  if (is.null(label)) {
    return(invisible(TRUE))
  }

  if (
    !is.character(label) ||
      length(label) != 1L ||
      is.na(label) ||
      !nzchar(trimws(label))
  ) {
    stop(
      "`label` must be NULL or a single nonempty character string.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# Count-basis inference ----------------------------------------------------

.collective_infer_count_basis <- function(frequency, severity) {
  frequency_basis <- .collective_basis_from_metadata(frequency)
  severity_basis <- .collective_basis_from_metadata(severity)

  known <- c(
    frequency_basis,
    severity_basis
  )

  known <- known[
    !is.na(known) & nzchar(known)
  ]

  if (length(known) == 0L) {
    return("loss")
  }

  known <- unique(known)

  if (length(known) > 1L) {
    stop(
      "The component specifications contain conflicting count-basis metadata. ",
      "Supply `count_basis` explicitly after verifying the intended model.",
      call. = FALSE
    )
  }

  .collective_validate_count_basis(known)
}


.collective_basis_from_metadata <- function(x) {
  if (!is.list(x$metadata)) {
    return(NA_character_)
  }

  basis <- x$metadata$count_basis

  if (!is.null(basis)) {
    if (
      is.character(basis) &&
        length(basis) == 1L &&
        !is.na(basis)
    ) {
      basis <- tolower(trimws(basis))

      if (basis %in% c("loss", "payment")) {
        return(basis)
      }
    }
  }

  payment_basis <- x$metadata$basis

  if (
    inherits(x, "tidyact_payment_spec") &&
      is.character(payment_basis) &&
      length(payment_basis) == 1L &&
      !is.na(payment_basis)
  ) {
    payment_basis <- tolower(trimws(payment_basis))

    if (identical(payment_basis, "per_payment")) {
      return("payment")
    }

    if (identical(payment_basis, "per_loss")) {
      return("loss")
    }
  }

  NA_character_
}


# Print helpers ------------------------------------------------------------

.collective_distribution_name <- function(x) {
  dist <- x$dist

  name <- switch(
    dist,
    poisson = "Poisson",
    binomial = "Binomial",
    negative_binomial = "Negative binomial",
    geometric = "Geometric",
    exponential = "Exponential",
    gamma = "Gamma",
    lognormal = "Lognormal",
    weibull = "Weibull",
    pareto = "Pareto II (Lomax)",
    burr = "Burr XII",
    loglogistic = "Log-logistic",
    degenerate = "Degenerate",
    discrete = "Finite discrete",
    empirical = "Empirical",
    custom = "Custom",
    dist
  )

  if (!is.null(x$label)) {
    paste0(
      name,
      " [",
      x$label,
      "]"
    )
  } else {
    name
  }
}
