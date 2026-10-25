#' Summarise Monte Carlo simulation output
#'
#' Computes tidy summary statistics from simulated actuarial values, using a
#' column-oriented interface consistent with the Monte Carlo functions in
#' \code{tidyactuarial}.
#'
#' This function is intentionally generic. It can summarise simulated present
#' values from annuities, insurances, premiums, reserves, losses, or any other
#' numeric actuarial indicator stored in a tidy simulation table.
#'
#' @param .data A data.frame or tibble containing simulation output.
#' @param col_value Character string. Name of the numeric column to summarise.
#'   Defaults to \code{"present_value"}.
#' @param by Optional character vector of grouping columns. If supplied, the
#'   summary is computed separately within each group. If \code{by = NULL} and
#'   \code{.data} is already grouped with \code{dplyr::group_by()}, the current
#'   grouping structure is used.
#' @param probs Numeric vector of probabilities for quantiles. Use
#'   \code{numeric()} to omit these columns.
#' @param var_probs Numeric vector of probabilities for VaR and TVaR. Use
#'   \code{numeric()} to omit these columns.
#' @param na_rm Logical scalar. If \code{TRUE}, missing values are removed
#'   before computing statistics.
#' @param ... Transitional compatibility for older calls using
#'   \code{data}, \code{value_col}, and \code{group_cols}.
#'
#' @return A tibble with summary statistics.
#'
#' @details
#' This function follows the compact column-naming convention used throughout
#' the Monte Carlo layer of \code{tidyactuarial}: column arguments are prefixed
#' with \code{col_}. Thus, \code{col_value} identifies the simulated actuarial
#' value to summarise.
#'
#' The function computes the number of simulations used, number of missing
#' values, mean, variance, standard deviation, standard error, minimum, maximum,
#' selected quantiles, empirical VaR, and empirical TVaR.
#'
#' Empirical VaR uses the lower quantile (the inverse empirical CDF; type 1),
#' without interpolation between observed loss amounts. In contrast, the
#' descriptive `q...` columns retain the default interpolated quantiles from
#' [stats::quantile()].
#'
#' TVaR uses the mean of the upper `1 - p` fraction of the empirical loss
#' distribution, with the appropriate fraction of any probability mass at VaR.
#' Equivalently, for `0 <= p < 1` and `v = VaR_p(X)`,
#' `TVaR_p(X) = v + E[(X - v)_+] / (1 - p)`.
#' For `p = 1`, the finite empirical upper-endpoint limit (`max(x)`) is
#' reported. At `p = 0`, TVaR is the empirical mean. Negative simulated
#' values are permitted because `summary_mc()` is a generic numeric summary.
#' TVaR is evaluated from ordered empirical probability masses, avoiding
#' subtraction of extreme losses. A fractional boundary atom is retained.
#' Variance, standard deviation and error are evaluated independently in
#' scaled form when squaring loses range: variance may be zero or infinite
#' while the deviation or error remains representable. Undefined moments
#' for nonfinite observations remain missing; extended upper tails can be
#' infinite. A sample smaller than two has missing variance, deviation and error.
#'
#' The column \code{se_mean} is the usual \eqn{s/\sqrt{n}} estimate under
#' independent, identically distributed observations. It is not a general
#' Monte Carlo error estimate for antithetic, clustered, weighted or otherwise
#' dependent samples. A method label alone does not identify independent
#' sampling units after filtering, grouping or reordering, so pairs are not
#' inferred automatically. For known independent antithetic pairs, summarise
#' their pair means instead. An odd unpaired observation or unequal weights
#' requires an error calculation appropriate to that sampling design.
#'
#' Quantile labels retain the existing rounding to thousandths of probability;
#' probabilities within a selector must yield distinct labels. Empty grouped
#' inputs retain the statistical columns with zero output rows. An ungrouped
#' empty effective sample retains one row of counts and missing statistics.
#'
#' If \code{na_rm = TRUE}, missing values in \code{col_value} are removed before
#' computing statistics. If \code{na_rm = FALSE} and missing values are present,
#' most numerical statistics are returned as \code{NA_real_}, avoiding accidental
#' silent deletion of missing reserve or loss scenarios.
#'
#' @family monte-carlo
#'
#' @examples
#' sim <- tibble::tibble(
#'   sim_id = 1:5,
#'   t = c(0, 0, 1, 1, 1),
#'   L_t = c(10, 12, 8, 15, 11)
#' )
#'
#' summary_mc(sim, col_value = "L_t")
#' summary_mc(sim, col_value = "L_t", by = "t")
#'
#' # Transitional compatibility with older argument names
#' summary_mc(sim, value_col = "L_t", group_cols = "t")
#'
#' @export
summary_mc <- function(
    .data = NULL,
    col_value = "present_value",
    by = NULL,
    probs = c(0.025, 0.5, 0.975),
    var_probs = c(0.95, 0.99),
    na_rm = TRUE,
    ...
) {
  dots <- list(...)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c("data", "value_col", "group_cols")
  bad_dots <- setdiff(names(dots), allowed_old)

  if (length(bad_dots) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad_dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(dots$data)) {
    if (!is.null(.data)) {
      stop("Provide only one of `.data` or deprecated `data`.", call. = FALSE)
    }

    .data <- dots$data
  }

  if (!is.null(dots$value_col)) {
    if (!identical(col_value, "present_value")) {
      stop("Provide only one of `col_value` or deprecated `value_col`.",
           call. = FALSE)
    }

    col_value <- dots$value_col
  }

  if (!is.null(dots$group_cols)) {
    if (!is.null(by)) {
      stop("Provide only one of `by` or deprecated `group_cols`.",
           call. = FALSE)
    }

    by <- dots$group_cols
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (!inherits(.data, "data.frame")) {
    stop("`.data` must be a data.frame or tibble.", call. = FALSE)
  }

  if (!is.character(col_value) ||
      length(col_value) != 1L ||
      is.na(col_value) ||
      !nzchar(col_value)) {
    stop("`col_value` must be a single non-empty string.", call. = FALSE)
  }

  if (!col_value %in% names(.data)) {
    stop("`col_value` was not found in `.data`.", call. = FALSE)
  }

  if (!is.numeric(.data[[col_value]])) {
    stop("`col_value` must identify a numeric column.", call. = FALSE)
  }

  if (!is.null(by)) {
    if (!is.character(by) || anyNA(by) || any(!nzchar(by))) {
      stop("`by` must be NULL or a character vector of column names.",
           call. = FALSE)
    }

    missing_groups <- setdiff(by, names(.data))

    if (length(missing_groups) > 0L) {
      stop(
        "The following grouping columns were not found in `.data`: ",
        paste(sprintf("`%s`", missing_groups), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(probs) ||
      anyNA(probs) ||
      any(!is.finite(probs)) ||
      any(probs < 0) ||
      any(probs > 1)) {
    stop("`probs` must contain finite probabilities between 0 and 1.",
         call. = FALSE)
  }

  if (!is.numeric(var_probs) ||
      anyNA(var_probs) ||
      any(!is.finite(var_probs)) ||
      any(var_probs < 0) ||
      any(var_probs > 1)) {
    stop("`var_probs` must contain finite probabilities between 0 and 1.",
         call. = FALSE)
  }

  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("`na_rm` must be TRUE or FALSE.", call. = FALSE)
  }

  # If no explicit grouping is supplied, respect existing dplyr groups.
  active_by <- by

  if (is.null(active_by)) {
    active_by <- dplyr::group_vars(.data)
  }

  if (length(active_by) == 0L) {
    active_by <- NULL
  }

  # -------------------------------------------------------------------------
  # Helpers
  # -------------------------------------------------------------------------

  make_prob_name <- function(prefix, p) {
    if (!length(p)) return(character())
    paste0(prefix, sprintf("%03.0f", p * 1000))
  }

  safe_numeric_stat <- function(x, fun, min_n = 1L) {
    if (length(x) < min_n) {
      return(NA_real_)
    }

    out <- fun(x)

    if (length(out) != 1L || is.na(out)) {
      return(NA_real_)
    }

    as.numeric(out)
  }

  empirical_tail_mean <- function(ordered, p, mean_value = mean(ordered)) {
    n <- length(ordered)
    if (p == 0) return(mean_value)
    if (p == 1 || ordered[1L] == ordered[n]) return(ordered[1L])
    finite <- all(is.finite(ordered))
    if (finite && max(abs(ordered)) < .Machine$double.xmin) {
      # A power-of-two lift keeps subnormal masses until after averaging.
      return(empirical_tail_mean(ordered * 2^1023, p) / 2^1023)
    }
    if (finite && n * p < 1) {
      # Remove a fraction of the lowest atom without rounding 1-p to 1
      # before forming the excluded monetary mass.
      return((mean_value - p * ordered[n]) / (1 - p))
    }
    mass <- n * (1 - p)
    if (mass <= 1) return(ordered[1L])
    whole <- floor(mass)
    fraction <- mass - whole
    values <- ordered[seq_len(whole)]
    weights <- rep(1 / mass, whole)
    if (fraction > 0) {
      values <- c(values, ordered[whole + 1])
      weights <- c(weights, fraction / mass)
    }
    out <- sum(values * weights)
    if (is.nan(out)) NA_real_ else out
  }

  summarise_one <- function(df) {
    x_raw <- df[[col_value]]
    n_total <- length(x_raw)
    n_missing <- sum(is.na(x_raw))

    if (isTRUE(na_rm)) {
      x <- x_raw[!is.na(x_raw)]
    } else {
      x <- x_raw
    }

    n <- length(x)
    has_missing <- anyNA(x)

    if (n == 0L || has_missing) {
      q_values <- rep(NA_real_, length(probs))
      names(q_values) <- make_prob_name("q", probs)

      var_values <- rep(NA_real_, length(var_probs))
      names(var_values) <- make_prob_name("VaR_", var_probs)

      tvar_values <- rep(NA_real_, length(var_probs))
      names(tvar_values) <- make_prob_name("TVaR_", var_probs)

      return(dplyr::bind_cols(
        tibble::tibble(
          n_sim = n,
          n_total = n_total,
          n_missing = n_missing,
          mean = NA_real_,
          variance = NA_real_,
          sd = NA_real_,
          se_mean = NA_real_,
          min = NA_real_,
          max = NA_real_
        ),
        tibble::as_tibble_row(as.list(q_values)),
        tibble::as_tibble_row(as.list(var_values)),
        tibble::as_tibble_row(as.list(tvar_values))
      ))
    }

    variance_value <- safe_numeric_stat(x, stats::var, min_n = 2L)
    sd_value <- sqrt(variance_value)
    se_value <- if (n > 1L) sd_value / sqrt(n) else NA_real_
    if (n > 1L && all(is.finite(x)) && any(x != x[1L]) &&
        (!is.finite(variance_value) || variance_value < .Machine$double.xmin)) {
      scale <- max(abs(x))
      centered <- x - mean(x)
      normalized <- if (all(is.finite(centered))) centered / scale else x / scale
      relative_sd <- stats::sd(normalized)
      sd_value <- scale * relative_sd
      se_value <- scale * (relative_sd / sqrt(n))
      variance_value <- .cash_flow_value(scale, log(scale) + 2 * log(relative_sd))
    }
    base <- tibble::tibble(
      n_sim = n,
      n_total = n_total,
      n_missing = n_missing,
      mean = safe_numeric_stat(x, mean),
      variance = variance_value,
      sd = sd_value,
      se_mean = se_value,
      min = safe_numeric_stat(x, min),
      max = safe_numeric_stat(x, max)
    )

    q_values <- as.numeric(stats::quantile(
      x,
      probs = probs,
      na.rm = FALSE,
      names = FALSE
    ))

    names(q_values) <- make_prob_name("q", probs)

    # Actuarial VaR: the lower quantile of the empirical distribution.
    # Descriptive quantiles above retain the usual type-7 interpolation.
    var_values <- as.numeric(stats::quantile(
      x,
      probs = var_probs,
      na.rm = FALSE,
      names = FALSE,
      type = 1L
    ))

    names(var_values) <- make_prob_name("VaR_", var_probs)

    # Empirical TVaR (Expected Shortfall): the required fraction of the
    # atom at VaR is included, rather than all observations equal to VaR.
    # At p = 1 the value is its bounded empirical right-endpoint limit.
    ordered <- if (length(var_probs)) sort(x, decreasing = TRUE) else numeric()
    tvar_values <- vapply(var_probs, function(p)
      empirical_tail_mean(ordered, p, base$mean), numeric(1L))

    names(tvar_values) <- make_prob_name("TVaR_", var_probs)

    dplyr::bind_cols(
      base,
      tibble::as_tibble_row(as.list(q_values)),
      tibble::as_tibble_row(as.list(var_values)),
      tibble::as_tibble_row(as.list(tvar_values))
    )
  }

  # -------------------------------------------------------------------------
  # Computation
  # -------------------------------------------------------------------------

  work <- tibble::as_tibble(.data)

  if (is.null(active_by)) {
    return(summarise_one(work))
  }

  grouped <- work |>
    dplyr::group_by(dplyr::across(dplyr::all_of(active_by))) |>
    dplyr::group_split()

  key_tbl <- work |>
    dplyr::group_by(dplyr::across(dplyr::all_of(active_by))) |>
    dplyr::group_keys()

  res <- if (length(grouped)) {
    dplyr::bind_rows(lapply(grouped, summarise_one))
  } else {
    summarise_one(work[0, , drop = FALSE])[0, , drop = FALSE]
  }

  dplyr::bind_cols(key_tbl, res)
}
