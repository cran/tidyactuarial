#' Compute the aggregate loss distribution
#'
#' @md
#' @description
#' Evaluates the distribution of the collective aggregate loss
#'
#' \deqn{S = \sum_{i=1}^{N} X_i}
#'
#' from a `tidyact_collective_risk_spec` object.
#'
#' The function separates the actuarial model from the computational
#' algorithm. It never fits a distribution, never simulates the aggregate
#' loss, and never discretizes a continuous severity implicitly.
#'
#' @param x A `tidyact_collective_risk_spec` object.
#' @param method Computational method. One of `"auto"`, `"analytic"`,
#'   `"convolution"`, `"panjer"`, or `"fft"`.
#' @param tol Positive finite probability tolerance smaller than one. It is
#'   used when an infinite claim-count distribution must be truncated
#'   numerically or when Panjer recursion is allowed to stop after capturing
#'   at least `1 - tol` of the aggregate probability.
#' @param max_loss Optional nonnegative finite numeric scalar giving the
#'   largest aggregate-loss value to retain in the returned grid. Probability
#'   above the retained grid is not renormalized; it is reported separately as
#'   `tail_probability`, together with any omitted claim-count probability.
#' @param ... Reserved for future computational controls. Additional arguments
#'   are currently rejected.
#'
#' @details
#' `aggregate_distribution()` operates only on an already specified collective
#' model. The probabilistic model is never changed in order to make an
#' algorithm applicable.
#'
#' In particular, if the severity distribution is continuous and no supported
#' analytic form is available, the function stops and asks the user to apply
#' [discretize_severity()] explicitly before using Panjer recursion or FFT.
#'
#' ## Automatic method selection
#'
#' With `method = "auto"`, the function uses the following hierarchy:
#'
#' 1. an implemented analytic representation for deterministic aggregate
#'    structures;
#' 2. Panjer recursion when the severity is on a lattice and the frequency
#'    distribution supplies valid Panjer parameters;
#' 3. exact finite convolution when both frequency and severity have finite
#'    discrete support;
#' 4. FFT when the severity is on a lattice and the frequency distribution
#'    supplies the probability operations required to truncate its count
#'    support.
#'
#' If none of these conditions is satisfied, an informative error is returned.
#'
#' ## Analytic method
#'
#' The initial analytic engine handles aggregate losses that are degenerate at
#' zero and deterministic severities. If \eqn{X = c} almost surely, then
#'
#' \deqn{S = cN}
#'
#' and the aggregate probabilities are obtained directly from the count
#' distribution. Infinite count support is truncated according to `tol` unless
#' `max_loss` supplies a tighter explicit cutoff.
#'
#' Additional analytic compound families may be added later without changing
#' the public interface.
#'
#' ## Convolution method
#'
#' The convolution engine is exact for finite claim-count support and finite
#' discrete severity support. It evaluates
#'
#' \deqn{
#' P(S \in A)
#' =
#' \sum_n P(N=n) P(X_1+\cdots+X_n \in A)
#' }
#'
#' by repeated finite convolution. If `max_loss` is supplied, the full finite
#' distribution is evaluated first and probability above the requested cutoff
#' is reported as tail probability.
#'
#' ## Panjer method
#'
#' For a frequency distribution in the \eqn{(a,b,0)} class,
#'
#' \deqn{
#' p_n = \left(a + \frac{b}{n}\right)p_{n-1},
#' }
#'
#' and a lattice severity with probabilities \eqn{f_j}, the aggregate
#' probabilities satisfy
#'
#' \deqn{
#' g_0 = G_N(f_0)
#' }
#'
#' and, for \eqn{k \ge 1},
#'
#' \deqn{
#' g_k =
#' \frac{
#' \sum_{j=1}^{k}
#' \left(a + \frac{bj}{k}\right)
#' f_j g_{k-j}
#' }{
#' 1-a f_0
#' }.
#' }
#'
#' Only severity probabilities available on the lattice are used. No
#' discretization is performed by this engine.
#'
#' ## FFT method
#'
#' FFT evaluates a finite polynomial approximation to the compound
#' distribution. For an infinite claim-count distribution, the count support is
#' first truncated at a quantile whose omitted probability is no greater than
#' `tol`. This omitted count probability is carried explicitly into
#' `tail_probability`; it is not wrapped into the retained FFT grid.
#'
#' The FFT engine therefore uses the truncated count polynomial
#'
#' \deqn{
#' \sum_{n=0}^{n_{\max}} p_n z^n
#' }
#'
#' rather than the full PGF. This avoids silently treating circular aliasing as
#' aggregate probability.
#'
#' ## Returned probabilities
#'
#' The returned `grid` contains:
#'
#' - `loss`: aggregate loss amount;
#' - `probability`: probability mass assigned to that amount;
#' - `cdf`: cumulative probability captured on the retained grid.
#'
#' The grid is intentionally not renormalized when a positive tail remains.
#' Consequently,
#'
#' \deqn{
#' \sum p_s + p_{\mathrm{omitted}} = 1
#' }
#'
#' up to numerical tolerance.
#' With loss-grid truncation alone, omitted mass is above the grid. FFT count
#' truncation also omits events with large claim counts but potentially small
#' aggregate losses (especially when zero severities are possible).
#' Thus `tail_probability` is total unrepresented mass, not necessarily
#' `P(S > max(grid$loss))`. VaR checks whether omitted count mass could change
#' its selected quantile; TVaR and stop-loss reject unresolved mass.
#'
#' @return
#' An object with classes `"tidyact_aggregate_distribution"` and
#' `"tidyact_risk_distribution"`. It contains:
#'
#' \describe{
#'   \item{`schema_version`}{Internal schema version.}
#'   \item{`method`}{The computational method actually used.}
#'   \item{`grid`}{A tibble with `loss`, `probability`, and `cdf`.}
#'   \item{`support`}{Support information for the retained grid.}
#'   \item{`tail_probability`}{Probability not represented on the retained
#'     grid.}
#'   \item{`model`}{The source collective risk specification.}
#'   \item{`label`}{The model label, if any.}
#'   \item{`metadata`}{Audit information about method selection, truncation,
#'     exactness, and numerical controls.}
#' }
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and
#' Nesbitt, C. J. (1997). *Actuarial Mathematics*, 2nd ed.
#' Society of Actuaries.
#'
#' Klugman, S. A., Panjer, H. H., and Willmot, G. E. (2019).
#' *Loss Models: From Data to Decisions*, 5th ed. Wiley.
#'
#' @seealso [collective_risk_spec()], [discretize_severity()],
#'   [aggregate_moments()]
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
#' # Continuous severities are never discretized implicitly.
#' \dontrun{
#' aggregate_distribution(model)
#' }
#'
#' lattice_severity <- discretize_severity(
#'   severity,
#'   step = 100
#' )
#'
#' lattice_model <- collective_risk_spec(
#'   frequency,
#'   lattice_severity
#' )
#'
#' aggregate_distribution(
#'   lattice_model,
#'   method = "panjer",
#'   tol = 1e-8
#' )
#'
#' aggregate_distribution(
#'   lattice_model,
#'   method = "fft",
#'   tol = 1e-8
#' )
#'
#' @export
aggregate_distribution <- function(
  x,
  method = c(
    "auto",
    "analytic",
    "convolution",
    "panjer",
    "fft"
  ),
  tol = 1e-8,
  max_loss = NULL,
  ...
) {
  .aggregate_distribution_validate_model(x)
  method <- match.arg(method)
  tol <- .aggregate_distribution_validate_tol(tol)
  max_loss <- .aggregate_distribution_validate_max_loss(max_loss)
  .aggregate_distribution_check_dots(...)

  requested_method <- method

  if (identical(method, "auto")) {
    method <- .aggregate_distribution_select_method(x)
  }

  result <- switch(
    method,
    analytic = .aggregate_distribution_analytic(
      x = x,
      tol = tol,
      max_loss = max_loss
    ),
    convolution = .aggregate_distribution_convolution(
      x = x,
      tol = tol,
      max_loss = max_loss
    ),
    panjer = .aggregate_distribution_panjer(
      x = x,
      tol = tol,
      max_loss = max_loss
    ),
    fft = .aggregate_distribution_fft(
      x = x,
      tol = tol,
      max_loss = max_loss
    )
  )

  result$metadata$requested_method <- requested_method
  result$metadata$selected_method <- method
  result$metadata$tol <- tol
  result$metadata$requested_max_loss <- max_loss

  result
}


#' Print an aggregate loss distribution
#'
#' @param x A `tidyact_aggregate_distribution` object.
#' @param ... Additional arguments. Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.tidyact_aggregate_distribution <- function(x, ...) {
  cat("<tidyact_aggregate_distribution>\n")
  cat(
    "Method: ",
    x$method,
    "\n",
    sep = ""
  )
  cat(
    "Grid points: ",
    nrow(x$grid),
    "\n",
    sep = ""
  )

  if (nrow(x$grid) > 0L) {
    cat(
      "Retained loss range: [",
      format(x$grid$loss[1L], trim = TRUE),
      ", ",
      format(x$grid$loss[nrow(x$grid)], trim = TRUE),
      "]\n",
      sep = ""
    )
  }

  cat(
    "Captured probability: ",
    format(
      sum(x$grid$probability),
      digits = 10,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )

  cat(
    "Tail probability: ",
    format(
      x$tail_probability,
      digits = 10,
      trim = TRUE
    ),
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

.aggregate_distribution_validate_model <- function(x) {
  if (!inherits(
    x,
    "tidyact_collective_risk_spec"
  )) {
    stop(
      "`x` must be a `tidyact_collective_risk_spec` object.",
      call. = FALSE
    )
  }

  if (
    !inherits(
      x$frequency,
      "tidyact_frequency_spec"
    ) ||
      !inherits(
        x$severity,
        "tidyact_severity_spec"
      )
  ) {
    stop(
      "`x` contains invalid frequency or severity components.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.aggregate_distribution_validate_tol <- function(tol) {
  if (
    !is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol <= 0 ||
      tol >= 1
  ) {
    stop(
      "`tol` must be a finite numeric scalar strictly between 0 and 1.",
      call. = FALSE
    )
  }

  as.numeric(tol)
}


.aggregate_distribution_validate_max_loss <- function(max_loss) {
  if (is.null(max_loss)) {
    return(NULL)
  }

  if (
    !is.numeric(max_loss) ||
      length(max_loss) != 1L ||
      is.na(max_loss) ||
      !is.finite(max_loss) ||
      max_loss < 0
  ) {
    stop(
      "`max_loss` must be NULL or a nonnegative finite numeric scalar.",
      call. = FALSE
    )
  }

  as.numeric(max_loss)
}


.aggregate_distribution_check_dots <- function(...) {
  dots <- list(...)

  if (length(dots) > 0L) {
    stop(
      "Additional arguments in `...` are reserved for future computational ",
      "controls and are not currently supported.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# Automatic method selection ----------------------------------------------

.aggregate_distribution_select_method <- function(x) {
  if (
    .aggregate_distribution_frequency_zero(x$frequency) ||
      .aggregate_distribution_severity_deterministic(x$severity)
  ) {
    return("analytic")
  }

  lattice <- .aggregate_distribution_lattice_severity(
    x$severity,
    error = FALSE
  )

  if (!is.null(lattice)) {
    panjer <- .aggregate_distribution_frequency_panjer(
      x$frequency
    )

    if (!is.null(panjer)) {
      return("panjer")
    }
  }

  if (
    .aggregate_distribution_frequency_finite(x$frequency) &&
      !is.null(
        .aggregate_distribution_severity_masses(
          x$severity,
          error = FALSE
        )
      )
  ) {
    return("convolution")
  }

  if (!is.null(lattice)) {
    return("fft")
  }

  if (identical(
    x$severity$support$type,
    "continuous"
  )) {
    stop(
      "No supported analytic aggregate distribution is available for this ",
      "continuous severity. Apply `discretize_severity()` explicitly before ",
      "using Panjer recursion or FFT.",
      call. = FALSE
    )
  }

  stop(
    "No available aggregate-distribution method is compatible with the ",
    "specified frequency and severity distributions.",
    call. = FALSE
  )
}


# Result constructor -------------------------------------------------------

.aggregate_distribution_result <- function(
  x,
  method,
  loss,
  probability,
  tail_probability,
  support_type = "discrete",
  step = NULL,
  exact = FALSE,
  metadata = list()
) {
  cleaned <- .aggregate_distribution_clean_result(
    loss = loss,
    probability = probability,
    tail_probability = tail_probability
  )

  loss <- cleaned$loss
  probability <- cleaned$probability
  tail_probability <- cleaned$tail_probability

  grid <- tibble::tibble(
    loss = loss,
    probability = probability,
    cdf = cumsum(probability)
  )

  support <- list(
    type = support_type,
    lower = if (length(loss)) min(loss) else 0,
    upper = if (length(loss)) max(loss) else 0
  )

  if (!is.null(step)) {
    support$step <- step
  }

  structure(
    list(
      schema_version = 1L,
      method = method,
      grid = grid,
      support = support,
      tail_probability = tail_probability,
      model = x,
      label = x$label,
      metadata = c(
        list(
          exact = exact
        ),
        metadata
      )
    ),
    class = c(
      "tidyact_aggregate_distribution",
      "tidyact_risk_distribution"
    )
  )
}


.aggregate_distribution_clean_result <- function(
  loss,
  probability,
  tail_probability
) {
  if (
    !is.numeric(loss) ||
      !is.numeric(probability) ||
      length(loss) != length(probability) ||
      length(loss) == 0L ||
      anyNA(loss) ||
      anyNA(probability) ||
      any(!is.finite(loss)) ||
      any(!is.finite(probability)) ||
      any(loss < 0)
  ) {
    stop(
      "Internal aggregate-distribution error: invalid probability grid.",
      call. = FALSE
    )
  }

  tolerance <- 1e-10

  probability[
    probability < 0 &
      probability >= -tolerance
  ] <- 0

  if (any(probability < -tolerance)) {
    stop(
      "The aggregate algorithm produced materially negative probabilities.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(tail_probability) ||
      length(tail_probability) != 1L ||
      is.na(tail_probability) ||
      !is.finite(tail_probability)
  ) {
    stop(
      "Internal aggregate-distribution error: invalid tail probability.",
      call. = FALSE
    )
  }

  if (
    tail_probability < -tolerance ||
      tail_probability > 1 + tolerance
  ) {
    stop(
      "Internal aggregate-distribution error: tail probability lies outside ",
      "[0, 1].",
      call. = FALSE
    )
  }

  tail_probability <- min(
    1,
    max(
      0,
      tail_probability
    )
  )

  ord <- order(loss)
  loss <- loss[ord]
  probability <- probability[ord]

  merged <- .aggregate_distribution_merge_masses(
    loss,
    probability
  )

  loss <- merged$values
  probability <- merged$probs

  total <- sum(probability) +
    tail_probability

  residual <- 1 - total

  if (abs(residual) <= 1e-8) {
    tail_probability <- tail_probability +
      residual
  } else {
    stop(
      "Aggregate probabilities plus tail probability do not sum to one ",
      "within numerical tolerance.",
      call. = FALSE
    )
  }

  if (tail_probability < 0 &&
      tail_probability >= -1e-10) {
    tail_probability <- 0
  }

  list(
    loss = loss,
    probability = probability,
    tail_probability = tail_probability
  )
}


# Distribution callbacks --------------------------------------------------

.aggregate_distribution_frequency_entry <- function(frequency) {
  registry <- .frequency_registry()

  if (!frequency$dist %in% names(registry)) {
    return(NULL)
  }

  registry[[frequency$dist]]
}


.aggregate_distribution_frequency_callback <- function(
  frequency,
  operation
) {
  if (identical(
    frequency$dist,
    "custom"
  )) {
    return(
      frequency$custom[[operation]]
    )
  }

  entry <- .aggregate_distribution_frequency_entry(
    frequency
  )

  if (is.null(entry)) {
    return(NULL)
  }

  entry[[operation]]
}


.aggregate_distribution_frequency_pmf <- function(
  frequency,
  n
) {
  pmf <- .aggregate_distribution_frequency_callback(
    frequency,
    "pmf"
  )

  if (!is.function(pmf)) {
    stop(
      "The frequency specification does not provide a PMF required by this ",
      "aggregate-distribution method.",
      call. = FALSE
    )
  }

  out <- pmf(
    n,
    frequency$params,
    log = FALSE
  )

  if (
    !is.numeric(out) ||
      length(out) != length(n) ||
      anyNA(out) ||
      any(!is.finite(out))
  ) {
    stop(
      "The frequency PMF returned invalid values.",
      call. = FALSE
    )
  }

  out
}


.aggregate_distribution_frequency_cdf <- function(
  frequency,
  n
) {
  cdf <- .aggregate_distribution_frequency_callback(
    frequency,
    "cdf"
  )

  if (is.function(cdf)) {
    out <- cdf(
      n,
      frequency$params,
      lower.tail = TRUE,
      log.p = FALSE
    )

    if (
      is.numeric(out) &&
        length(out) == length(n) &&
        !anyNA(out) &&
        all(is.finite(out))
    ) {
      return(out)
    }
  }

  vapply(
    n,
    function(k) {
      sum(
        .aggregate_distribution_frequency_pmf(
          frequency,
          0:floor(k)
        )
      )
    },
    numeric(1)
  )
}


.aggregate_distribution_frequency_quantile <- function(
  frequency,
  p
) {
  quantile <- .aggregate_distribution_frequency_callback(
    frequency,
    "quantile"
  )

  if (!is.function(quantile)) {
    stop(
      "The frequency specification does not provide a quantile function ",
      "required to truncate its infinite support.",
      call. = FALSE
    )
  }

  out <- quantile(
    p,
    frequency$params,
    lower.tail = TRUE,
    log.p = FALSE
  )

  if (
    !is.numeric(out) ||
      length(out) != length(p) ||
      anyNA(out) ||
      any(!is.finite(out)) ||
      any(out < 0)
  ) {
    stop(
      "The frequency quantile function returned invalid count values.",
      call. = FALSE
    )
  }

  as.integer(
    ceiling(out)
  )
}


.aggregate_distribution_frequency_pgf <- function(
  frequency,
  z
) {
  pgf <- .aggregate_distribution_frequency_callback(
    frequency,
    "pgf"
  )

  if (is.function(pgf)) {
    return(
      pgf(
        z,
        frequency$params
      )
    )
  }

  if (
    .aggregate_distribution_frequency_finite(
      frequency
    )
  ) {
    n_max <- as.integer(
      frequency$support$upper
    )

    probs <- .aggregate_distribution_frequency_pmf(
      frequency,
      0:n_max
    )

    return(
      vapply(
        z,
        function(w) {
          sum(
            probs *
              w^(0:n_max)
          )
        },
        if (is.complex(z)) 0 + 0i else 0
      )
    )
  }

  stop(
    "The frequency specification does not provide a PGF required by Panjer ",
    "recursion.",
    call. = FALSE
  )
}


.aggregate_distribution_frequency_panjer <- function(frequency) {
  panjer <- .aggregate_distribution_frequency_callback(
    frequency,
    "panjer"
  )

  if (!is.function(panjer)) {
    return(NULL)
  }

  out <- panjer(
    frequency$params
  )

  if (is.null(out)) {
    return(NULL)
  }

  required <- c(
    "a",
    "b",
    "p0",
    "p1"
  )

  if (
    !is.list(out) ||
      !all(required %in% names(out)) ||
      any(
        !vapply(
          out[required],
          function(z) {
            is.numeric(z) &&
              length(z) == 1L &&
              !is.na(z) &&
              is.finite(z)
          },
          logical(1)
        )
      )
  ) {
    stop(
      "The frequency Panjer callback must return finite numeric `a`, `b`, ",
      "`p0`, and `p1` values.",
      call. = FALSE
    )
  }

  out[required]
}


.aggregate_distribution_frequency_finite <- function(frequency) {
  is.list(frequency$support) &&
    is.numeric(frequency$support$upper) &&
    length(frequency$support$upper) == 1L &&
    !is.na(frequency$support$upper) &&
    is.finite(frequency$support$upper)
}


.aggregate_distribution_frequency_zero <- function(frequency) {
  p0 <- tryCatch(
    .aggregate_distribution_frequency_pmf(
      frequency,
      0
    ),
    error = function(e) NA_real_
  )

  is.numeric(p0) &&
    length(p0) == 1L &&
    is.finite(p0) &&
    abs(p0 - 1) <= 1e-12
}


# Severity masses and lattice ---------------------------------------------

.aggregate_distribution_severity_masses <- function(
  severity,
  error = TRUE
) {
  out <- NULL

  if (identical(
    severity$dist,
    "degenerate"
  )) {
    out <- list(
      values = severity$params$value,
      probs = 1
    )
  }

  if (
    is.null(out) &&
      severity$dist %in% c(
        "discrete",
        "empirical"
      )
  ) {
    out <- list(
      values = severity$params$values,
      probs = severity$params$probs
    )
  }

  if (
    is.null(out) &&
      severity$support$type %in% c(
        "lattice",
        "discrete"
      ) &&
      is.finite(severity$support$upper) &&
      !is.null(severity$support$step)
  ) {
    pmf <- if (identical(
      severity$dist,
      "custom"
    )) {
      severity$custom$pmf
    } else {
      entry <- .severity_registry()[[
        severity$dist
      ]]

      if (is.null(entry)) NULL else entry$pmf
    }

    if (is.function(pmf)) {
      grid <- seq(
        severity$support$lower,
        severity$support$upper,
        by = severity$support$step
      )

      probs <- pmf(
        grid,
        severity$params,
        log = FALSE
      )

      out <- list(
        values = grid,
        probs = probs
      )
    }
  }

  if (is.null(out)) {
    if (error) {
      stop(
        "This method requires a finite discrete severity distribution.",
        call. = FALSE
      )
    }

    return(NULL)
  }

  if (
    !is.numeric(out$values) ||
      !is.numeric(out$probs) ||
      length(out$values) != length(out$probs) ||
      length(out$values) == 0L ||
      anyNA(out$values) ||
      anyNA(out$probs) ||
      any(!is.finite(out$values)) ||
      any(!is.finite(out$probs)) ||
      any(out$values < 0) ||
      any(out$probs < 0)
  ) {
    stop(
      "The severity specification contains invalid discrete probability ",
      "masses.",
      call. = FALSE
    )
  }

  positive <- out$probs > 0

  out$values <- out$values[positive]
  out$probs <- out$probs[positive]

  total <- sum(out$probs)

  if (
    !is.finite(total) ||
      abs(total - 1) > 1e-10
  ) {
    stop(
      "Discrete severity probabilities must sum to one.",
      call. = FALSE
    )
  }

  out$probs <- out$probs / total

  ord <- order(out$values)

  out$values <- out$values[ord]
  out$probs <- out$probs[ord]

  out
}


.aggregate_distribution_severity_deterministic <- function(severity) {
  masses <- .aggregate_distribution_severity_masses(
    severity,
    error = FALSE
  )

  !is.null(masses) &&
    length(masses$values) == 1L &&
    abs(masses$probs - 1) <= 1e-12
}


.aggregate_distribution_lattice_severity <- function(
  severity,
  error = TRUE
) {
  masses <- .aggregate_distribution_severity_masses(
    severity,
    error = FALSE
  )

  if (is.null(masses)) {
    if (error) {
      stop(
        "Panjer recursion and FFT require a finite lattice severity. ",
        "Use `discretize_severity()` explicitly for continuous severities.",
        call. = FALSE
      )
    }

    return(NULL)
  }

  values <- masses$values
  probs <- masses$probs

  step <- NULL

  if (
    identical(
      severity$support$type,
      "lattice"
    ) &&
      !is.null(
        severity$support$step
      )
  ) {
    step <- severity$support$step
  } else {
    step <- .aggregate_distribution_infer_decimal_lattice(
      values
    )
  }

  if (
    is.null(step) ||
      !is.numeric(step) ||
      length(step) != 1L ||
      !is.finite(step) ||
      step <= 0
  ) {
    if (error) {
      stop(
        "The severity amounts do not define a supported arithmetic lattice. ",
        "Use `discretize_severity()` to create an explicit lattice.",
        call. = FALSE
      )
    }

    return(NULL)
  }

  indices <- values / step

  tolerance <- sqrt(
    .Machine$double.eps
  ) *
    pmax(
      1,
      abs(indices)
    )

  if (
    any(
      abs(
        indices -
          round(indices)
      ) >
        tolerance
    )
  ) {
    if (error) {
      stop(
        "Severity values are not integer multiples of the lattice step.",
        call. = FALSE
      )
    }

    return(NULL)
  }

  indices <- as.integer(
    round(indices)
  )

  if (any(indices < 0L)) {
    stop(
      "Severity lattice indices must be nonnegative.",
      call. = FALSE
    )
  }

  max_index <- max(indices)

  if (max_index > 2000000L) {
    stop(
      "The severity lattice is too large for the current aggregate engines.",
      call. = FALSE
    )
  }

  f <- numeric(
    max_index + 1L
  )

  for (i in seq_along(indices)) {
    f[indices[i] + 1L] <-
      f[indices[i] + 1L] +
      probs[i]
  }

  list(
    step = step,
    f = f,
    max_index = max_index,
    values = values,
    probs = probs
  )
}


.aggregate_distribution_infer_decimal_lattice <- function(values) {
  if (length(values) == 1L) {
    if (values == 0) {
      return(1)
    }

    return(values)
  }

  # Work in units of the support magnitude so monetary units cannot turn
  # strictly positive amounts into zero lattice indices.
  magnitude <- max(values)
  unit <- if (magnitude > 0) 10^floor(log10(magnitude)) else 1
  if (!is.finite(unit) || unit == 0) return(NULL)
  values <- values / unit

  for (digits in 0:10) {
    scale <- 10^digits
    scaled <- values * scale
    rounded <- round(scaled)

    tolerance <- 1e-8 *
      pmax(
        1,
        abs(scaled)
      )

    if (
      all(
        abs(
          scaled -
            rounded
        ) <= tolerance
      ) && !any(values > 0 & rounded == 0) &&
        max(abs(rounded)) <
          2^52
    ) {
      integers <- abs(
        rounded
      )

      integers <- integers[
        integers > 0
      ]

      if (length(integers) == 0L) {
        return(1 / scale)
      }

      gcd_value <- integers[1L]

      if (length(integers) > 1L) {
        for (z in integers[-1L]) {
          gcd_value <-
            .aggregate_distribution_numeric_gcd(
              gcd_value,
              z
            )
        }
      }

      if (gcd_value > 0) {
        return(
          (gcd_value / scale) * unit
        )
      }
    }
  }

  NULL
}


.aggregate_distribution_numeric_gcd <- function(a, b) {
  a <- abs(round(a))
  b <- abs(round(b))

  while (b > 0) {
    r <- a %% b
    a <- b
    b <- r
  }

  a
}


# Analytic engine ----------------------------------------------------------

.aggregate_distribution_analytic <- function(
  x,
  tol,
  max_loss
) {
  frequency <- x$frequency
  severity <- x$severity

  if (
    .aggregate_distribution_frequency_zero(
      frequency
    )
  ) {
    return(
      .aggregate_distribution_result(
        x = x,
        method = "analytic",
        loss = 0,
        probability = 1,
        tail_probability = 0,
        support_type = "discrete",
        exact = TRUE,
        metadata = list(
          analytic_case = "zero_frequency"
        )
      )
    )
  }

  masses <- .aggregate_distribution_severity_masses(
    severity,
    error = FALSE
  )

  if (
    is.null(masses) ||
      length(masses$values) != 1L
  ) {
    stop(
      "The current analytic engine is available only for a zero aggregate ",
      "or deterministic severity.",
      call. = FALSE
    )
  }

  amount <- masses$values[1L]

  if (amount == 0) {
    return(
      .aggregate_distribution_result(
        x = x,
        method = "analytic",
        loss = 0,
        probability = 1,
        tail_probability = 0,
        support_type = "discrete",
        exact = TRUE,
        metadata = list(
          analytic_case = "zero_severity"
        )
      )
    )
  }

  finite_frequency <-
    .aggregate_distribution_frequency_finite(
      frequency
    )

  if (!is.null(max_loss)) {
    n_max <- floor(
      max_loss / amount
    )

    n_max <- max(
      0,
      as.integer(n_max)
    )

    if (finite_frequency) {
      n_max <- min(
        n_max,
        as.integer(
          frequency$support$upper
        )
      )
    }
  } else if (finite_frequency) {
    n_max <- as.integer(
      frequency$support$upper
    )
  } else {
    n_max <-
      .aggregate_distribution_frequency_quantile(
        frequency,
        1 - tol
      )[1L]
  }

  counts <- 0:n_max

  probs <-
    .aggregate_distribution_frequency_pmf(
      frequency,
      counts
    )

  tail <- max(
    0,
    1 -
      sum(probs)
  )

  exact <- finite_frequency &&
    n_max >=
      frequency$support$upper

  .aggregate_distribution_result(
    x = x,
    method = "analytic",
    loss = amount *
      counts,
    probability = probs,
    tail_probability = tail,
    support_type = "lattice",
    step = amount,
    exact = exact,
    metadata = list(
      analytic_case = "deterministic_severity",
      count_truncation = n_max
    )
  )
}


# Exact finite convolution engine -----------------------------------------

.aggregate_distribution_convolution <- function(
  x,
  tol,
  max_loss
) {
  frequency <- x$frequency
  severity <- x$severity

  severity_masses <-
    .aggregate_distribution_severity_masses(
      severity,
      error = TRUE
    )

  if (
    !.aggregate_distribution_frequency_finite(
      frequency
    )
  ) {
    stop(
      "The convolution engine currently requires finite claim-count support. ",
      "Use Panjer recursion or FFT for supported infinite-count models.",
      call. = FALSE
    )
  }

  n_max <- as.integer(
    frequency$support$upper
  )

  counts <- 0:n_max

  count_probs <-
    .aggregate_distribution_frequency_pmf(
      frequency,
      counts
    )

  current_values <- 0
  current_probs <- 1

  result_values <- numeric()
  result_probs <- numeric()

  for (n in counts) {
    weight <- count_probs[n + 1L]

    if (weight > 0) {
      mixed <-
        .aggregate_distribution_combine_masses(
          result_values,
          result_probs,
          current_values,
          weight *
            current_probs
        )

      result_values <- mixed$values
      result_probs <- mixed$probs
    }

    if (n < n_max) {
      convolved <-
        .aggregate_distribution_convolve_masses(
          current_values,
          current_probs,
          severity_masses$values,
          severity_masses$probs
        )

      current_values <- convolved$values
      current_probs <- convolved$probs
    }
  }

  tail <- 0

  if (!is.null(max_loss)) {
    keep <- result_values <=
      max_loss +
        1e-12 *
          max(
            max_loss
          )

    tail <- sum(
      result_probs[!keep]
    )

    result_values <- result_values[keep]
    result_probs <- result_probs[keep]
  }

  lattice <-
    .aggregate_distribution_lattice_severity(
      severity,
      error = FALSE
    )

  support_type <- if (is.null(lattice)) {
    "discrete"
  } else {
    "lattice"
  }

  step <- if (is.null(lattice)) {
    NULL
  } else {
    lattice$step
  }

  .aggregate_distribution_result(
    x = x,
    method = "convolution",
    loss = result_values,
    probability = result_probs,
    tail_probability = tail,
    support_type = support_type,
    step = step,
    exact = tail == 0,
    metadata = list(
      count_support_upper = n_max,
      convolution_type = "finite_exact"
    )
  )
}


.aggregate_distribution_convolve_masses <- function(
  x_values,
  x_probs,
  y_values,
  y_probs
) {
  terms <- length(x_values) *
    length(y_values)

  if (terms > 5000000) {
    stop(
      "Finite convolution would create more than 5,000,000 intermediate ",
      "terms. Use a lattice method such as Panjer recursion or FFT.",
      call. = FALSE
    )
  }

  values <- as.vector(
    outer(
      x_values,
      y_values,
      "+"
    )
  )

  probs <- as.vector(
    outer(
      x_probs,
      y_probs,
      "*"
    )
  )

  .aggregate_distribution_merge_masses(
    values,
    probs
  )
}


.aggregate_distribution_combine_masses <- function(
  x_values,
  x_probs,
  y_values,
  y_probs
) {
  if (length(x_values) == 0L) {
    return(
      .aggregate_distribution_merge_masses(
        y_values,
        y_probs
      )
    )
  }

  .aggregate_distribution_merge_masses(
    c(
      x_values,
      y_values
    ),
    c(
      x_probs,
      y_probs
    )
  )
}


.aggregate_distribution_merge_masses <- function(
  values,
  probs
) {
  if (length(values) == 0L) {
    return(
      list(
        values = numeric(),
        probs = numeric()
      )
    )
  }

  keep <- probs > 0

  values <- values[keep]
  probs <- probs[keep]

  ord <- order(values)

  values <- values[ord]
  probs <- probs[ord]

  if (length(values) == 1L) {
    return(
      list(
        values = values,
        probs = probs
      )
    )
  }

  out_values <- numeric()
  out_probs <- numeric()

  current_value <- values[1L]
  current_prob <- probs[1L]

  for (i in 2:length(values)) {
    tolerance <-
      32 *
        .Machine$double.eps *
        max(
          abs(current_value),
          abs(values[i])
        )

    if (
      abs(
        values[i] -
          current_value
      ) <= tolerance
    ) {
      # Values within numerical tolerance represent the same support point.
      # Keep the first canonical value instead of averaging coordinates:
      # weighted averaging can introduce floating-point drift in exact
      # arithmetic lattices (for example, 30 becoming 30.000000000000004).
      current_prob <-
        current_prob +
        probs[i]
    } else {
      out_values <- c(
        out_values,
        current_value
      )

      out_probs <- c(
        out_probs,
        current_prob
      )

      current_value <- values[i]
      current_prob <- probs[i]
    }
  }

  list(
    values = c(
      out_values,
      current_value
    ),
    probs = c(
      out_probs,
      current_prob
    )
  )
}


# Panjer engine ------------------------------------------------------------

.aggregate_distribution_panjer <- function(
  x,
  tol,
  max_loss
) {
  frequency <- x$frequency
  severity <- x$severity

  lattice <-
    .aggregate_distribution_lattice_severity(
      severity,
      error = TRUE
    )

  panjer <-
    .aggregate_distribution_frequency_panjer(
      frequency
    )

  if (is.null(panjer)) {
    stop(
      "The frequency distribution does not provide an applicable Panjer ",
      "(a, b, 0) representation.",
      call. = FALSE
    )
  }

  f <- lattice$f
  m <- lattice$max_index
  step <- lattice$step

  f0 <- f[1L]

  denominator <- 1 -
    panjer$a *
      f0

  if (
    !is.finite(denominator) ||
      abs(denominator) <=
        sqrt(.Machine$double.eps)
  ) {
    stop(
      "Panjer recursion is numerically singular because `1 - a * f0` is ",
      "too close to zero.",
      call. = FALSE
    )
  }

  g0 <-
    .aggregate_distribution_frequency_pgf(
      frequency,
      f0
    )

  g0 <- Re(g0[1L])

  if (
    !is.finite(g0) ||
      g0 < -1e-12 ||
      g0 > 1 + 1e-12
  ) {
    stop(
      "The frequency PGF returned an invalid initial Panjer probability.",
      call. = FALSE
    )
  }

  g0 <- min(
    1,
    max(
      0,
      g0
    )
  )

  finite_frequency <-
    .aggregate_distribution_frequency_finite(
      frequency
    )

  finite_upper <- NULL

  if (finite_frequency) {
    finite_upper <-
      as.integer(
        frequency$support$upper
      ) *
      m
  }

  if (!is.null(max_loss)) {
    k_limit <- floor(
      max_loss / step
    )

    k_limit <- max(
      0,
      as.integer(k_limit)
    )

    if (!is.null(finite_upper)) {
      k_limit <- min(
        k_limit,
        finite_upper
      )
    }

    stop_when_tol <- FALSE
  } else if (!is.null(finite_upper)) {
    k_limit <- finite_upper
    stop_when_tol <- FALSE
  } else {
    k_limit <- 1000000L
    stop_when_tol <- TRUE
  }

  g <- numeric(
    k_limit + 1L
  )

  g[1L] <- g0
  cumulative <- g0
  last_k <- 0L

  if (k_limit > 0L) {
    for (k in seq_len(k_limit)) {
      j_max <- min(
        k,
        m
      )

      if (j_max == 0L) {
        value <- 0
      } else {
        j <- seq_len(j_max)

        value <- sum(
          (
            panjer$a +
              panjer$b *
                j /
                k
          ) *
            f[j + 1L] *
            g[k - j + 1L]
        ) /
          denominator
      }

      if (
        value < 0 &&
          value >= -1e-12
      ) {
        value <- 0
      }

      if (
        !is.finite(value) ||
          value < -1e-12
      ) {
        stop(
          "Panjer recursion produced an invalid aggregate probability.",
          call. = FALSE
        )
      }

      g[k + 1L] <- value
      cumulative <- cumulative +
        value

      last_k <- k

      if (
        stop_when_tol &&
          cumulative >=
            1 - tol
      ) {
        break
      }
    }
  }

  if (
    stop_when_tol &&
      cumulative <
        1 - tol
  ) {
    stop(
      "Panjer recursion did not reach the requested probability tolerance ",
      "within 1,000,000 lattice points. Supply `max_loss`, increase `tol`, ",
      "or reconsider the lattice step.",
      call. = FALSE
    )
  }

  g <- g[
    seq_len(
      last_k +
        1L
    )
  ]

  losses <- (
    seq_along(g) -
      1L
  ) *
    step

  tail <- max(
    0,
    1 -
      sum(g)
  )

  exact <- !is.null(finite_upper) &&
    last_k >=
      finite_upper &&
    is.null(max_loss)

  .aggregate_distribution_result(
    x = x,
    method = "panjer",
    loss = losses,
    probability = g,
    tail_probability = tail,
    support_type = "lattice",
    step = step,
    exact = exact,
    metadata = list(
      panjer_a = panjer$a,
      panjer_b = panjer$b,
      severity_grid_points = length(f),
      recursion_last_index = last_k,
      stopped_by_tolerance = stop_when_tol &&
        cumulative >=
          1 - tol
    )
  )
}


# FFT engine ---------------------------------------------------------------

.aggregate_distribution_fft <- function(
  x,
  tol,
  max_loss
) {
  frequency <- x$frequency
  severity <- x$severity

  lattice <-
    .aggregate_distribution_lattice_severity(
      severity,
      error = TRUE
    )

  f <- lattice$f
  m <- lattice$max_index
  step <- lattice$step

  if (m == 0L) {
    return(
      .aggregate_distribution_analytic(
        x,
        tol = tol,
        max_loss = max_loss
      )
    )
  }

  finite_frequency <-
    .aggregate_distribution_frequency_finite(
      frequency
    )

  if (finite_frequency) {
    n_max <- as.integer(
      frequency$support$upper
    )
  } else {
    n_max <-
      .aggregate_distribution_frequency_quantile(
        frequency,
        1 - tol
      )[1L]
  }

  n <- 0:n_max

  count_probs <-
    .aggregate_distribution_frequency_pmf(
      frequency,
      n
    )

  represented_count_mass <-
    sum(
      count_probs
    )

  count_tail <- max(
    0,
    1 -
      represented_count_mass
  )

  degree <- n_max *
    m

  if (
    degree >
      4000000
  ) {
    stop(
      "The required FFT polynomial exceeds 4,000,000 aggregate lattice ",
      "indices. Increase the severity step, increase `tol`, or use Panjer ",
      "recursion when available.",
      call. = FALSE
    )
  }

  fft_length <- 1L

  while (
    fft_length <
      degree +
        1L
  ) {
    fft_length <-
      fft_length *
      2L
  }

  if (
    fft_length >
      4194304L
  ) {
    stop(
      "The required FFT length exceeds 2^22. Increase the severity step or ",
      "increase `tol`.",
      call. = FALSE
    )
  }

  f_pad <- numeric(
    fft_length
  )

  f_pad[
    seq_along(f)
  ] <- f

  severity_values <- stats::fft(
    f_pad
  )

  aggregate_values <- rep(
    count_probs[n_max + 1L],
    fft_length
  ) +
    0i

  if (n_max > 0L) {
    for (
      k in seq(
        from = n_max - 1L,
        to = 0L,
        by = -1L
      )
    ) {
      aggregate_values <-
        aggregate_values *
        severity_values +
        count_probs[k + 1L]
    }
  }

  coefficients <-
    stats::fft(
      aggregate_values,
      inverse = TRUE
    ) /
    fft_length

  imaginary_error <-
    max(
      abs(
        Im(coefficients)
      )
    )

  probs <- Re(
    coefficients[
      seq_len(
        degree +
          1L
      )
    ]
  )

  probs[
    probs < 0 &
      probs >= -1e-10
  ] <- 0

  if (any(probs < -1e-10)) {
    stop(
      "FFT inversion produced materially negative aggregate probabilities.",
      call. = FALSE
    )
  }

  numerical_residual <-
    represented_count_mass -
    sum(probs)

  if (
    abs(
      numerical_residual
    ) <=
      1e-8
  ) {
    largest <- which.max(
      probs
    )

    probs[largest] <-
      probs[largest] +
      numerical_residual
  } else {
    stop(
      "FFT probabilities do not recover the represented claim-count mass ",
      "within numerical tolerance.",
      call. = FALSE
    )
  }

  losses <- (
    seq_along(probs) -
      1L
  ) *
    step

  tail <- count_tail

  if (!is.null(max_loss)) {
    keep <- losses <=
      max_loss +
        1e-12 *
          max(
            max_loss
          )

    tail <- tail +
      sum(
        probs[!keep]
      )

    losses <- losses[keep]
    probs <- probs[keep]
  }

  .aggregate_distribution_result(
    x = x,
    method = "fft",
    loss = losses,
    probability = probs,
    tail_probability = tail,
    support_type = "lattice",
    step = step,
    exact = finite_frequency &&
      is.null(max_loss),
    metadata = list(
      count_truncation = n_max,
      represented_count_probability =
        represented_count_mass,
      count_tail_probability = count_tail,
      fft_length = fft_length,
      polynomial_degree = degree,
      maximum_imaginary_error =
        imaginary_error
    )
  )
}
