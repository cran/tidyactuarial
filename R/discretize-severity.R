#' Discretize a continuous claim severity distribution
#'
#' @md
#' @description
#' Converts a continuous nonnegative severity specification into an explicit
#' lattice distribution. The transformation is intended for actuarial
#' algorithms that require discrete severities, such as Panjer recursion or
#' FFT-based aggregate-loss methods.
#'
#' Discretization is always explicit: this function is never called silently by
#' aggregate-loss functions.
#'
#' @param x A continuous `tidyact_severity_spec` object.
#' @param step Positive finite numeric scalar giving the lattice spacing.
#' @param to Optional positive finite upper grid point. It must be an integer
#'   multiple of `step`. If `NULL`, the upper point is selected from the
#'   `(1 - tail_prob)` quantile and rounded upward to the next grid point.
#' @param method Discretization method. One of `"unbiased"`, `"upper"`,
#'   `"lower"`, or `"rounding"`.
#' @param tail_prob Probability used to select `to` when `to = NULL`.
#'   Must satisfy `0 < tail_prob < 1`.
#' @param label Optional nonempty character scalar. If `NULL`, the label of
#'   `x` is retained.
#'
#' @details
#' Let \eqn{F} be the severity CDF, \eqn{h} the lattice step, and
#' \eqn{b = \mathrm{to}}. The lattice is
#'
#' \deqn{0, h, 2h, \ldots, b.}
#'
#' The four methods follow standard actuarial discretization conventions.
#'
#' ## Upper method
#'
#' For lattice points below `to`,
#'
#' \deqn{
#' p_x = F(x+h)-F(x).
#' }
#'
#' This assigns each interval to its lower lattice endpoint. Probability above
#' `to` is assigned explicitly to the upper grid point.
#'
#' ## Lower method
#'
#' \deqn{
#' p_0 = F(0), \qquad
#' p_x = F(x)-F(x-h).
#' }
#'
#' This assigns each interval to its upper lattice endpoint. Probability above
#' `to` is added to the upper grid point.
#'
#' ## Rounding method
#'
#' Probability is allocated according to the midpoint intervals around each
#' lattice point. Thus the first cell extends from zero to `step / 2`, interior
#' cells have width `step`, and the last cell contains everything above
#' `to - step / 2`.
#'
#' ## Unbiased method
#'
#' The unbiased method performs local first-moment matching. If
#' \eqn{L(d) = E[\min(X,d)]}, the finite-grid probabilities before the final
#' tail adjustment are
#'
#' \deqn{
#' p_0 = 1 - F(0) - \frac{L(h)-L(0)}{h},
#' }
#'
#' \deqn{
#' p_x =
#' \frac{2L(x)-L(x-h)-L(x+h)}{h},
#' \qquad 0 < x < b,
#' }
#'
#' and
#'
#' \deqn{
#' p_b =
#' \frac{L(b)-L(b-h)}{h} - 1 + F(b).
#' }
#'
#' The residual tail probability above `to` is then assigned to `to`.
#' Consequently, the resulting distribution is a proper probability
#' distribution and its first moment matches \eqn{E[\min(X,b)]}, rather than
#' silently discarding probability above the truncation point.
#'
#' For the exponential, gamma, lognormal, Weibull, and Pareto II families,
#' the first limited moment used by the unbiased method is evaluated
#' analytically. For other continuous specifications with a CDF, including
#' Burr XII, log-logistic, and eligible custom models, the limited mean is
#' obtained by numerical integration of the survival function. The source used
#' is recorded in `metadata$discretization$limited_mean_source`.
#'
#' Only continuous severity specifications are accepted in this initial
#' implementation. Finite discrete and empirical severities already have
#' probability masses and should not be silently re-discretized. Support for
#' mixed severities can be added later without changing this public interface.
#'
#' The returned object uses the existing `"discrete"` severity operations but
#' has `support$type = "lattice"`. Zero-probability grid points need not be
#' stored explicitly in `params$values`; the lattice information is retained
#' in `support`.
#'
#' @return
#' A `tidyact_severity_spec` and `tidyact_distribution_spec` object with
#' `dist = "discrete"`, lattice support, and discretization information in
#' `metadata$discretization`.
#'
#' The metadata records:
#'
#' - `method`;
#' - `step`;
#' - `to`;
#' - whether `to` was supplied or selected from `tail_prob`;
#' - the requested `tail_prob` when applicable;
#' - the actual probability above `to`;
#' - the explicit tail treatment;
#' - the original distribution identifier; and
#' - the source of the limited mean for the unbiased method.
#'
#' @references
#' Klugman, S. A., Panjer, H. H., and Willmot, G. E. (2019).
#' *Loss Models: From Data to Decisions*, 5th ed. Wiley.
#'
#' @seealso [severity_spec()], [distribution_convolution()]
#'
#' @examples
#' x <- severity_spec(
#'   "gamma",
#'   shape = 2,
#'   rate = 0.001
#' )
#'
#' # Automatic upper grid point from a small tail probability.
#' dx <- discretize_severity(
#'   x,
#'   step = 100
#' )
#' dx
#'
#' # Inspect the resulting lattice probabilities.
#' tibble::tibble(
#'   amount = dx$params$values,
#'   probability = dx$params$probs
#' )
#'
#' # Use a fixed truncation point and another discretization rule.
#' discretize_severity(
#'   x,
#'   step = 250,
#'   to = 10000,
#'   method = "rounding"
#' )
#'
#' # The discretization choices remain auditable.
#' dx$metadata$discretization
#'
#' @export
discretize_severity <- function(
  x,
  step,
  to = NULL,
  method = c(
    "unbiased",
    "upper",
    "lower",
    "rounding"
  ),
  tail_prob = 1e-8,
  label = NULL
) {
  .discretize_validate_severity(x)
  step <- .discretize_validate_step(step)
  method <- .discretize_validate_method(method)
  .discretize_validate_tail_prob(tail_prob)
  .discretize_validate_label(label)

  if (!identical(x$support$type, "continuous")) {
    stop(
      "`x` must have continuous support. ",
      "Discrete, empirical, lattice, and mixed severities are not ",
      "re-discretized by this implementation.",
      call. = FALSE
    )
  }

  cdf <- .discretize_get_operation(
    x,
    "cdf"
  )

  if (!is.function(cdf)) {
    stop(
      "The severity specification must provide a CDF for discretization.",
      call. = FALSE
    )
  }

  cutoff_source <- "supplied"
  requested_tail_prob <- NA_real_

  if (is.null(to)) {
    quantile <- .discretize_get_operation(
      x,
      "quantile"
    )

    if (!is.function(quantile)) {
      stop(
        "`to` must be supplied when the severity specification does not ",
        "provide a quantile function.",
        call. = FALSE
      )
    }

    target_probability <- 1 - tail_prob

    raw_to <- .discretize_eval_quantile(
      quantile,
      target_probability,
      x$params
    )

    if (
      length(raw_to) != 1L ||
        is.na(raw_to) ||
        !is.finite(raw_to) ||
        raw_to < 0
    ) {
      stop(
        "The requested tail quantile is not a finite nonnegative value. ",
        "Supply `to` explicitly.",
        call. = FALSE
      )
    }

    grid_index <- ceiling(
      raw_to / step -
        8 * .Machine$double.eps
    )

    grid_index <- max(
      1,
      grid_index
    )

    to <- grid_index * step
    cutoff_source <- "tail_prob"
    requested_tail_prob <- tail_prob
  } else {
    to <- .discretize_validate_to(
      to,
      step
    )
  }

  to <- .discretize_validate_to(
    to,
    step
  )

  n_steps <- as.numeric(
    round(to / step)
  )

  if (
    !is.finite(n_steps) ||
      n_steps < 1 ||
      n_steps > .Machine$integer.max - 1
  ) {
    stop(
      "The combination of `to` and `step` creates an unsupported grid size.",
      call. = FALSE
    )
  }

  n_steps <- as.integer(n_steps)

  grid <- (0:n_steps) * step
  grid[length(grid)] <- to

  F0 <- .discretize_eval_cdf(
    cdf,
    0,
    x$params
  )

  if (
    !is.finite(F0) ||
      abs(F0) >
        1e-10
  ) {
    stop(
      "A continuous nonnegative severity must satisfy F(0) = 0 for these ",
      "discretization methods.",
      call. = FALSE
    )
  }

  F_to <- .discretize_eval_cdf(
    cdf,
    to,
    x$params
  )

  F_to <- .discretize_clean_probability(
    F_to,
    "F(to)"
  )

  actual_tail <- 1 - F_to

  limited_mean_source <- NA_character_

  probabilities <- switch(
    method,

    upper = {
      F_grid <- .discretize_eval_cdf(
        cdf,
        grid,
        x$params
      )

      F_grid <- .discretize_clean_cdf(
        F_grid
      )

      c(
        diff(F_grid),
        actual_tail
      )
    },

    lower = {
      F_grid <- .discretize_eval_cdf(
        cdf,
        grid,
        x$params
      )

      F_grid <- .discretize_clean_cdf(
        F_grid
      )

      out <- c(
        F_grid[1L],
        diff(F_grid)
      )

      out[length(out)] <- out[length(out)] +
        actual_tail

      out
    },

    rounding = {
      midpoints <- (
        seq_len(n_steps) -
          0.5
      ) * step

      F_mid <- .discretize_eval_cdf(
        cdf,
        midpoints,
        x$params
      )

      F_mid <- .discretize_clean_cdf(
        F_mid
      )

      c(
        F_mid[1L],
        if (length(F_mid) > 1L) diff(F_mid) else numeric(),
        1 - F_mid[length(F_mid)]
      )
    },

    unbiased = {
      limited <- .discretize_limited_mean(
        x,
        grid,
        cdf
      )

      limited_mean_source <- attr(
        limited,
        "source",
        exact = TRUE
      )

      if (
        length(limited) != length(grid) ||
          anyNA(limited) ||
          any(!is.finite(limited))
      ) {
        stop(
          "The first limited moment could not be evaluated on the ",
          "discretization grid.",
          call. = FALSE
        )
      }

      out <- numeric(
        length(grid)
      )

      out[1L] <- (
        limited[1L] -
          limited[2L]
      ) / step +
        1 -
        F0

      if (length(grid) > 2L) {
        interior <- 2:(length(grid) - 1L)

        out[interior] <- (
          2 * limited[interior] -
            limited[interior - 1L] -
            limited[interior + 1L]
        ) / step
      }

      out[length(out)] <- (
        limited[length(limited)] -
          limited[length(limited) - 1L]
      ) / step -
        1 +
        F_to

      out[length(out)] <- out[length(out)] +
        actual_tail

      out
    }
  )

  probabilities <- .discretize_clean_masses(
    probabilities
  )

  if (length(probabilities) != length(grid)) {
    stop(
      "Internal discretization error: probability and grid lengths differ.",
      call. = FALSE
    )
  }

  positive <- probabilities > 0

  if (!any(positive)) {
    stop(
      "Discretization produced no positive probability masses.",
      call. = FALSE
    )
  }

  values <- grid[positive]
  probs <- probabilities[positive]

  out_label <- if (is.null(label)) {
    x$label
  } else {
    label
  }

  structure(
    list(
      schema_version = 1L,
      dist = "discrete",
      params = list(
        values = values,
        probs = probs
      ),
      support = list(
        type = "lattice",
        lower = 0,
        upper = to,
        step = step
      ),
      data = NULL,
      custom = NULL,
      label = out_label,
      metadata = list(
        discretization = list(
          method = method,
          step = step,
          to = to,
          cutoff_source = cutoff_source,
          requested_tail_prob = requested_tail_prob,
          actual_tail_probability = actual_tail,
          tail_handling = "assigned_to_upper_grid_point",
          source_dist = x$dist,
          limited_mean_source = limited_mean_source
        )
      )
    ),
    class = c(
      "tidyact_severity_spec",
      "tidyact_distribution_spec"
    )
  )
}


# Validation ---------------------------------------------------------------

.discretize_validate_severity <- function(x) {
  if (!inherits(
    x,
    "tidyact_severity_spec"
  )) {
    stop(
      "`x` must be a `tidyact_severity_spec` object.",
      call. = FALSE
    )
  }

  if (
    !inherits(
      x,
      "tidyact_distribution_spec"
    ) ||
      !is.list(x$support)
  ) {
    stop(
      "`x` contains an invalid severity specification.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.discretize_validate_step <- function(step) {
  if (
    !is.numeric(step) ||
      length(step) != 1L ||
      is.na(step) ||
      !is.finite(step) ||
      step <= 0
  ) {
    stop(
      "`step` must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }

  as.numeric(step)
}


.discretize_validate_to <- function(
  to,
  step
) {
  if (
    !is.numeric(to) ||
      length(to) != 1L ||
      is.na(to) ||
      !is.finite(to) ||
      to <= 0
  ) {
    stop(
      "`to` must be a positive finite numeric scalar.",
      call. = FALSE
    )
  }

  ratio <- to / step

  tolerance <- sqrt(
    .Machine$double.eps
  ) * max(
    1,
    abs(ratio)
  )

  if (
    abs(
      ratio -
        round(ratio)
    ) >
      tolerance
  ) {
    stop(
      "`to` must be an integer multiple of `step`.",
      call. = FALSE
    )
  }

  round(ratio) * step
}


.discretize_validate_method <- function(method) {
  if (
    !is.character(method) ||
      length(method) < 1L ||
      anyNA(method)
  ) {
    stop(
      "`method` must be one of \"unbiased\", \"upper\", \"lower\", ",
      "or \"rounding\".",
      call. = FALSE
    )
  }

  match.arg(
    method,
    c(
      "unbiased",
      "upper",
      "lower",
      "rounding"
    )
  )
}


.discretize_validate_tail_prob <- function(tail_prob) {
  if (
    !is.numeric(tail_prob) ||
      length(tail_prob) != 1L ||
      is.na(tail_prob) ||
      !is.finite(tail_prob) ||
      tail_prob <= 0 ||
      tail_prob >= 1
  ) {
    stop(
      "`tail_prob` must be a finite numeric scalar strictly between 0 and 1.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.discretize_validate_label <- function(label) {
  if (is.null(label)) {
    return(
      invisible(TRUE)
    )
  }

  if (
    !is.character(label) ||
      length(label) != 1L ||
      is.na(label) ||
      !nzchar(
        trimws(label)
      )
  ) {
    stop(
      "`label` must be NULL or a single non-empty character string.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# Distribution operations --------------------------------------------------

.discretize_get_operation <- function(
  x,
  operation
) {
  if (identical(
    x$dist,
    "custom"
  )) {
    return(
      x$custom[[operation]]
    )
  }

  registry <- .severity_registry()

  if (!x$dist %in% names(registry)) {
    return(NULL)
  }

  registry[[x$dist]][[operation]]
}


.discretize_eval_cdf <- function(
  cdf,
  q,
  params
) {
  out <- cdf(
    q,
    params,
    lower.tail = TRUE,
    log.p = FALSE
  )

  if (
    !is.numeric(out) ||
      length(out) != length(q)
  ) {
    stop(
      "The severity CDF must return one numeric value per evaluation point.",
      call. = FALSE
    )
  }

  out
}


.discretize_eval_quantile <- function(
  quantile,
  prob,
  params
) {
  quantile(
    prob,
    params,
    lower.tail = TRUE,
    log.p = FALSE
  )
}


.discretize_clean_probability <- function(
  p,
  name
) {
  if (
    !is.numeric(p) ||
      length(p) != 1L ||
      is.na(p) ||
      !is.finite(p)
  ) {
    stop(
      "`",
      name,
      "` is not a finite probability.",
      call. = FALSE
    )
  }

  tolerance <- 64 *
    .Machine$double.eps

  if (
    p < -tolerance ||
      p > 1 + tolerance
  ) {
    stop(
      "`",
      name,
      "` lies outside [0, 1].",
      call. = FALSE
    )
  }

  min(
    1,
    max(
      0,
      p
    )
  )
}


.discretize_clean_cdf <- function(Fx) {
  if (
    !is.numeric(Fx) ||
      anyNA(Fx) ||
      any(!is.finite(Fx))
  ) {
    stop(
      "The severity CDF returned non-finite values on the discretization grid.",
      call. = FALSE
    )
  }

  tolerance <- 64 *
    .Machine$double.eps

  if (
    any(Fx < -tolerance) ||
      any(Fx > 1 + tolerance)
  ) {
    stop(
      "The severity CDF returned values outside [0, 1].",
      call. = FALSE
    )
  }

  Fx <- pmin(
    1,
    pmax(
      0,
      Fx
    )
  )

  if (
    length(Fx) > 1L &&
      any(
        diff(Fx) <
          -sqrt(.Machine$double.eps)
      )
  ) {
    stop(
      "The severity CDF is not nondecreasing on the discretization grid.",
      call. = FALSE
    )
  }

  cummax(Fx)
}


# Limited means ------------------------------------------------------------

.discretize_limited_mean <- function(
  x,
  limit,
  cdf
) {
  if (identical(
    x$dist,
    "custom"
  )) {
    custom_limited <- x$custom$limited_moment

    if (is.function(custom_limited)) {
      out <- vapply(
        limit,
        function(d) {
          custom_limited(
            1,
            d,
            x$params
          )
        },
        numeric(1)
      )

      attr(
        out,
        "source"
      ) <- "custom_analytic"

      return(out)
    }

    out <- .discretize_limited_mean_numeric(
      cdf,
      limit,
      x$params
    )

    attr(
      out,
      "source"
    ) <- "numerical_cdf_integration"

    return(out)
  }

  analytical <- .discretize_builtin_limited_mean(
    x,
    limit
  )

  if (!is.null(analytical)) {
    attr(
      analytical,
      "source"
    ) <- "analytic"

    return(analytical)
  }

  entry <- .severity_registry()[[
    x$dist
  ]]

  if (
    !is.null(entry) &&
      is.function(
        entry$limited_moment
      )
  ) {
    out <- vapply(
      limit,
      function(d) {
        entry$limited_moment(
          1,
          d,
          x$params
        )
      },
      numeric(1)
    )

    attr(
      out,
      "source"
    ) <- "registry_analytic"

    return(out)
  }

  out <- .discretize_limited_mean_numeric(
    cdf,
    limit,
    x$params
  )

  attr(
    out,
    "source"
  ) <- "numerical_cdf_integration"

  out
}


.discretize_builtin_limited_mean <- function(
  x,
  limit
) {
  params <- x$params

  switch(
    x$dist,

    exponential = {
      rate <- params$rate

      -expm1(
        -rate * limit
      ) / rate
    },

    gamma = {
      shape <- params$shape
      rate <- params$rate

      shape / rate *
        stats::pgamma(
          limit,
          shape = shape + 1,
          rate = rate
        ) +
        limit *
        stats::pgamma(
          limit,
          shape = shape,
          rate = rate,
          lower.tail = FALSE
        )
    },

    lognormal = {
      meanlog <- params$meanlog
      sdlog <- params$sdlog

      out <- numeric(
        length(limit)
      )

      positive <- limit > 0

      if (any(positive)) {
        d <- limit[positive]

        truncated_first <- exp(
          meanlog +
            sdlog^2 / 2
        ) *
          stats::pnorm(
            (
              log(d) -
                meanlog -
                sdlog^2
            ) / sdlog
          )

        out[positive] <- truncated_first +
          d *
          stats::plnorm(
            d,
            meanlog = meanlog,
            sdlog = sdlog,
            lower.tail = FALSE
          )
      }

      out
    },

    weibull = {
      shape <- params$shape
      scale <- params$scale

      z <- (
        limit / scale
      )^shape

      scale *
        gamma(
          1 +
            1 / shape
        ) *
        stats::pgamma(
          z,
          shape = 1 +
            1 / shape
        ) +
        limit *
        exp(
          -z
        )
    },

    pareto = {
      shape <- params$shape
      scale <- params$scale

      if (
        abs(
          shape -
            1
        ) <=
          sqrt(.Machine$double.eps)
      ) {
        scale *
          log1p(
            limit / scale
          )
      } else {
        scale *
          (
            expm1(
              (
                1 -
                  shape
              ) *
                log1p(
                  limit / scale
                )
            )
          ) /
          (
            1 -
              shape
          )
      }
    },

    NULL
  )
}


.discretize_limited_mean_numeric <- function(
  cdf,
  limit,
  params
) {
  vapply(
    limit,
    function(d) {
      if (d == 0) {
        return(0)
      }

      integrand <- function(z) {
        Fz <- .discretize_eval_cdf(
          cdf,
          z,
          params
        )

        1 -
          Fz
      }

      result <- tryCatch(
        stats::integrate(
          integrand,
          lower = 0,
          upper = d,
          subdivisions = 200L,
          rel.tol = 1e-9,
          stop.on.error = TRUE
        ),
        error = function(e) {
          stop(
            "Numerical evaluation of the first limited moment failed: ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )

      result$value
    },
    numeric(1)
  )
}


# Probability cleanup ------------------------------------------------------

.discretize_clean_masses <- function(probabilities) {
  if (
    !is.numeric(probabilities) ||
      !length(probabilities) ||
      anyNA(probabilities) ||
      any(!is.finite(probabilities))
  ) {
    stop(
      "Discretization produced invalid probability masses.",
      call. = FALSE
    )
  }

  tolerance <- 1e-10

  small_negative <- probabilities < 0 &
    probabilities >= -tolerance

  probabilities[small_negative] <- 0

  if (
    any(
      probabilities <
        -tolerance
    )
  ) {
    stop(
      "Discretization produced materially negative probability masses. ",
      "Try a finer `step` or another discretization method.",
      call. = FALSE
    )
  }

  total <- sum(
    probabilities
  )

  if (
    !is.finite(total) ||
      total <= 0
  ) {
    stop(
      "Discretized probabilities do not have a positive finite total.",
      call. = FALSE
    )
  }

  residual <- 1 -
    total

  if (
    abs(residual) >
      1e-8
  ) {
    stop(
      "Discretized probabilities do not sum to one within numerical tolerance.",
      call. = FALSE
    )
  }

  probabilities[length(probabilities)] <-
    probabilities[length(probabilities)] +
    residual

  if (
    probabilities[length(probabilities)] <
      -tolerance
  ) {
    stop(
      "Numerical correction would create a negative final probability mass.",
      call. = FALSE
    )
  }

  if (
    probabilities[length(probabilities)] < 0
  ) {
    probabilities[length(probabilities)] <- 0
  }

  probabilities
}
