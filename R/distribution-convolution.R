#' Specify a convolution of probability distributions
#'
#' @md
#' @description
#' Creates a lazy specification of the convolution of two independent
#' distributions or of an iid sum of a single distribution.
#'
#' If \eqn{X} and \eqn{Y} are independent, the pairwise form represents
#'
#' \deqn{Z = X + Y.}
#'
#' The iid form represents
#'
#' \deqn{Z_n = X_1 + \cdots + X_n,}
#'
#' where \eqn{X_1,\ldots,X_n} are independent copies of the distribution
#' specified by `x`.
#'
#' The function records the convolution model but does not evaluate a density,
#' mass function, distribution function, or numerical convolution.
#'
#' @param x A `tidyact_distribution_spec` object.
#' @param y Optional second `tidyact_distribution_spec` object. Supply `y` for
#'   a pairwise convolution `x + y`.
#' @param n Optional nonnegative integer scalar. Supply `n` for the iid
#'   convolution of `n` copies of `x`.
#' @param label Optional nonempty character scalar identifying the convolution.
#'
#' @details
#' Exactly one of `y` and `n` must be supplied.
#'
#' The pairwise form
#'
#' ```r
#' distribution_convolution(x, y)
#' ```
#'
#' represents the distribution of the sum of two independent random variables.
#'
#' The iid form
#'
#' ```r
#' distribution_convolution(x, n = 5)
#' ```
#'
#' represents the five-fold convolution of `x`.
#'
#' The cases `n = 0` and `n = 1` follow the usual convolution conventions:
#'
#' - `n = 0` represents the degenerate distribution at zero;
#' - `n = 1` returns the original distribution specification. If `label` is
#'   supplied, only the label of the returned copy is changed.
#'
#' The operation is deliberately *lazy*. It does not choose a computational
#' algorithm, discretize a continuous distribution, truncate an infinite
#' support, or apply FFT or Panjer recursion. Those decisions belong to later
#' computational functions.
#'
#' Support bounds are propagated algebraically. For example, if `x` has
#' support bounded below by \eqn{a} and `y` by \eqn{b}, the convolution is
#' bounded below by \eqn{a+b}. Infinite upper bounds are preserved.
#'
#' Support type is inferred conservatively:
#'
#' - count + count remains count;
#' - iid lattice sums remain lattice;
#' - a continuous component convolved with a discrete component is continuous;
#' - any convolution involving a mixed distribution is classified as mixed;
#' - otherwise discrete supports are retained as discrete.
#'
#' Independence is part of the meaning of a convolution specification and is
#' therefore not stored as a separate flag.
#'
#' @return
#' For a genuine convolution, an object with classes
#' `"tidyact_convolution_spec"` and `"tidyact_distribution_spec"`.
#'
#' The object follows the common distribution specification schema:
#'
#' \describe{
#'   \item{`schema_version`}{Internal schema version.}
#'   \item{`dist`}{The identifier `"convolution"`.}
#'   \item{`params`}{The lazy convolution definition.}
#'   \item{`support`}{Inferred support information.}
#'   \item{`data`}{Always `NULL`.}
#'   \item{`custom`}{Always `NULL`.}
#'   \item{`label`}{Optional user label.}
#'   \item{`metadata`}{Reserved metadata list.}
#' }
#'
#' When `n = 1`, the original distribution specification is returned instead
#' of wrapping it in a convolution object.
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and
#' Nesbitt, C. J. (1997). *Actuarial Mathematics*, 2nd ed.
#' Society of Actuaries.
#'
#' Klugman, S. A., Panjer, H. H., and Willmot, G. E. (2019).
#' *Loss Models: From Data to Decisions*, 5th ed. Wiley.
#'
#' @seealso [frequency_spec()], [severity_spec()], [collective_risk_spec()]
#'
#' @examples
#' x <- severity_spec(
#'   "gamma",
#'   shape = 2,
#'   rate = 0.001
#' )
#'
#' y <- severity_spec(
#'   "lognormal",
#'   meanlog = 7,
#'   sdlog = 0.6
#' )
#'
#' # Pairwise convolution: X + Y
#' distribution_convolution(x, y)
#'
#' # Four iid claim amounts: X1 + ... + X4
#' distribution_convolution(
#'   x,
#'   n = 4
#' )
#'
#' # The zero-fold convolution is degenerate at zero.
#' distribution_convolution(
#'   x,
#'   n = 0
#' )
#'
#' # The one-fold convolution is the original specification.
#' identical(
#'   distribution_convolution(x, n = 1),
#'   x
#' )
#'
#' # Count distributions can also be convolved.
#' n_claims <- frequency_spec(
#'   "poisson",
#'   lambda = 3
#' )
#'
#' distribution_convolution(
#'   n_claims,
#'   n = 2
#' )
#'
#' @export
distribution_convolution <- function(
  x,
  y = NULL,
  n = NULL,
  label = NULL
) {
  .convolution_validate_distribution(x, "x")
  .convolution_validate_label(label)

  has_y <- !is.null(y)
  has_n <- !is.null(n)

  if (identical(has_y, has_n)) {
    stop(
      "Supply exactly one of `y` or `n`.",
      call. = FALSE
    )
  }

  if (has_y) {
    .convolution_validate_distribution(y, "y")

    support <- .convolution_pair_support(
      x$support,
      y$support
    )

    params <- list(
      mode = "pair",
      components = list(
        x,
        y
      )
    )
  } else {
    n <- .convolution_validate_n(n)

    if (n == 1L) {
      if (is.null(label)) {
        return(x)
      }

      out <- x
      out$label <- label
      return(out)
    }

    support <- .convolution_iid_support(
      x$support,
      n
    )

    params <- list(
      mode = "iid",
      component = x,
      n = n
    )
  }

  structure(
    list(
      schema_version = 1L,
      dist = "convolution",
      params = params,
      support = support,
      data = NULL,
      custom = NULL,
      label = label,
      metadata = list()
    ),
    class = c(
      "tidyact_convolution_spec",
      "tidyact_distribution_spec"
    )
  )
}


#' Print a convolution specification
#'
#' @param x A `tidyact_convolution_spec` object.
#' @param ... Additional arguments. Currently ignored.
#'
#' @return The input object, invisibly.
#'
#' @export
print.tidyact_convolution_spec <- function(x, ...) {
  cat("<tidyact_convolution_spec>\n")

  mode <- x$params$mode

  if (identical(mode, "pair")) {
    components <- x$params$components

    cat("Mode: pairwise\n")
    cat(
      "Components: ",
      .convolution_component_name(components[[1L]]),
      " + ",
      .convolution_component_name(components[[2L]]),
      "\n",
      sep = ""
    )
  } else {
    cat("Mode: iid\n")
    cat(
      "Component: ",
      .convolution_component_name(x$params$component),
      "\n",
      sep = ""
    )
    cat(
      "Copies: ",
      x$params$n,
      "\n",
      sep = ""
    )
  }

  cat(
    "Support: ",
    .convolution_format_support(x$support),
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

.convolution_validate_distribution <- function(x, name) {
  if (!inherits(x, "tidyact_distribution_spec")) {
    stop(
      "`",
      name,
      "` must be a `tidyact_distribution_spec` object.",
      call. = FALSE
    )
  }

  if (
    !is.list(x$support) ||
      is.null(x$support$type) ||
      is.null(x$support$lower) ||
      is.null(x$support$upper)
  ) {
    stop(
      "`",
      name,
      "` contains an invalid support specification.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.convolution_validate_label <- function(label) {
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
      "`label` must be NULL or a single non-empty character string.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.convolution_validate_n <- function(n) {
  if (
    !is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      !is.finite(n) ||
      n < 0 ||
      n != floor(n) ||
      n > .Machine$integer.max
  ) {
    stop(
      "`n` must be a nonnegative integer scalar.",
      call. = FALSE
    )
  }

  as.integer(n)
}


# Support propagation ------------------------------------------------------

.convolution_pair_support <- function(x_support, y_support) {
  lower <- x_support$lower + y_support$lower

  upper <- if (
    is.infinite(x_support$upper) ||
      is.infinite(y_support$upper)
  ) {
    Inf
  } else {
    x_support$upper + y_support$upper
  }

  support_type <- .convolution_pair_support_type(
    x_support,
    y_support
  )

  out <- list(
    type = support_type,
    lower = lower,
    upper = upper
  )

  if (identical(support_type, "count")) {
    out$step <- 1
  }

  if (identical(support_type, "lattice")) {
    out$step <- .convolution_common_step(
      x_support,
      y_support
    )
  }

  out
}


.convolution_iid_support <- function(support, n) {
  if (n == 0L) {
    if (identical(support$type, "count")) {
      return(
        list(
          type = "count",
          lower = 0,
          upper = 0,
          step = 1
        )
      )
    }

    return(
      list(
        type = "discrete",
        lower = 0,
        upper = 0
      )
    )
  }

  upper <- if (is.infinite(support$upper)) {
    Inf
  } else {
    n * support$upper
  }

  out <- list(
    type = support$type,
    lower = n * support$lower,
    upper = upper
  )

  if (
    support$type %in% c("count", "lattice") &&
      !is.null(support$step)
  ) {
    out$step <- support$step
  }

  out
}


.convolution_pair_support_type <- function(
  x_support,
  y_support
) {
  x_type <- x_support$type
  y_type <- y_support$type

  if (
    identical(x_type, "mixed") ||
      identical(y_type, "mixed")
  ) {
    return("mixed")
  }

  if (
    identical(x_type, "continuous") ||
      identical(y_type, "continuous")
  ) {
    return("continuous")
  }

  if (
    identical(x_type, "count") &&
      identical(y_type, "count")
  ) {
    return("count")
  }

  discrete_types <- c(
    "count",
    "discrete",
    "lattice"
  )

  if (
    x_type %in% discrete_types &&
      y_type %in% discrete_types
  ) {
    common_step <- .convolution_common_step(
      x_support,
      y_support
    )

    if (!is.null(common_step)) {
      return("lattice")
    }

    return("discrete")
  }

  "mixed"
}


.convolution_common_step <- function(
  x_support,
  y_support
) {
  x_step <- .convolution_support_step(x_support)
  y_step <- .convolution_support_step(y_support)

  if (
    is.null(x_step) ||
      is.null(y_step)
  ) {
    return(NULL)
  }

  tolerance <- sqrt(.Machine$double.eps) *
    max(
      1,
      abs(x_step),
      abs(y_step)
    )

  if (abs(x_step - y_step) <= tolerance) {
    return(
      (x_step + y_step) / 2
    )
  }

  NULL
}


.convolution_support_step <- function(support) {
  if (identical(support$type, "count")) {
    return(
      if (is.null(support$step)) 1 else support$step
    )
  }

  if (
    identical(support$type, "lattice") &&
      !is.null(support$step)
  ) {
    return(support$step)
  }

  NULL
}


# Print helpers ------------------------------------------------------------

.convolution_component_name <- function(x) {
  if (!is.null(x$label)) {
    return(
      paste0(
        '"',
        x$label,
        '"'
      )
    )
  }

  if (
    !is.null(x$dist) &&
      is.character(x$dist) &&
      length(x$dist) == 1L
  ) {
    return(x$dist)
  }

  class(x)[1L]
}


.convolution_format_support <- function(support) {
  lower <- format(
    support$lower,
    trim = TRUE
  )

  upper <- if (is.infinite(support$upper)) {
    "Inf"
  } else {
    format(
      support$upper,
      trim = TRUE
    )
  }

  if (identical(support$type, "count")) {
    if (
      is.finite(support$upper) &&
        support$lower == support$upper
    ) {
      return(
        paste0(
          "count {",
          lower,
          "}"
        )
      )
    }

    return(
      paste0(
        "count [",
        lower,
        ", ",
        upper,
        "]"
      )
    )
  }

  if (identical(support$type, "lattice")) {
    return(
      paste0(
        "lattice [",
        lower,
        ", ",
        upper,
        "], step = ",
        format(
          support$step,
          trim = TRUE
        )
      )
    )
  }

  paste0(
    support$type,
    " [",
    lower,
    ", ",
    upper,
    "]"
  )
}
