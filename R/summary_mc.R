#' Summarise Monte Carlo simulation output
#'
#' Computes tidy summary statistics from simulated actuarial values.
#'
#' This function is intentionally generic. It can summarise simulated present
#' values from annuities, insurances, premiums, reserves, losses, or any other
#' numeric actuarial indicator stored in a tidy simulation table.
#'
#' @param .data A data.frame or tibble containing simulation output.
#' @param value_col Character string. Name of the numeric column to summarise.
#'   Defaults to `"present_value"`.
#' @param group_cols Optional character vector of grouping columns.
#' @param by Optional character vector of grouping columns. This is a convenient
#'   alias for `group_cols`, useful in examples and pipelines.
#' @param probs Numeric vector of probabilities for quantiles.
#' @param var_probs Numeric vector of probabilities for VaR and TVaR.
#' @param na_rm Logical. If `TRUE`, missing values are removed.
#'
#' @return A tibble with summary statistics.
#'
#' @details
#' The function computes the number of simulations, mean, variance, standard
#' deviation, standard error, minimum, maximum, selected quantiles, VaR, and
#' TVaR.
#'
#' TVaR is computed empirically as the mean of simulated values greater than
#' or equal to the corresponding empirical VaR.
#'
#' @family simulation
#'
#' @examples
#' sim <- tibble::tibble(
#'   simulation_id = 1:5,
#'   duration = c(0, 0, 1, 1, 1),
#'   present_value = c(10, 12, 8, 15, 11)
#' )
#'
#' summary_mc(sim)
#' summary_mc(sim, by = "duration")
#'
#' @export
summary_mc <- function(
    .data,
    value_col = "present_value",
    group_cols = NULL,
    by = NULL,
    probs = c(0.025, 0.5, 0.975),
    var_probs = c(0.95, 0.99),
    na_rm = TRUE
) {
  if (!inherits(.data, "data.frame")) {
    stop("`.data` must be a data.frame or tibble.", call. = FALSE)
  }

  if (!is.null(group_cols) && !is.null(by)) {
    stop("Use only one of `group_cols` or `by`, not both.", call. = FALSE)
  }

  if (is.null(group_cols) && !is.null(by)) {
    group_cols <- by
  }

  if (!is.character(value_col) || length(value_col) != 1L ||
      is.na(value_col) || !nzchar(value_col)) {
    stop("`value_col` must be a single non-empty string.", call. = FALSE)
  }

  if (!value_col %in% names(.data)) {
    stop("`value_col` was not found in `.data`.", call. = FALSE)
  }

  if (!is.null(group_cols)) {
    if (!is.character(group_cols) || anyNA(group_cols) ||
        any(!nzchar(group_cols))) {
      stop("`group_cols`/`by` must be NULL or a character vector of column names.", call. = FALSE)
    }

    missing_groups <- setdiff(group_cols, names(.data))

    if (length(missing_groups) > 0L) {
      stop(
        "The following grouping columns were not found in `.data`: ",
        paste(missing_groups, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(probs) || anyNA(probs) || any(probs < 0) || any(probs > 1)) {
    stop("`probs` must contain probabilities between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(var_probs) || anyNA(var_probs) ||
      any(var_probs < 0) || any(var_probs > 1)) {
    stop("`var_probs` must contain probabilities between 0 and 1.", call. = FALSE)
  }

  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("`na_rm` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(.data[[value_col]])) {
    stop("`value_col` must identify a numeric column.", call. = FALSE)
  }

  summarise_one <- function(df) {
    x <- df[[value_col]]

    if (isTRUE(na_rm)) {
      x <- x[!is.na(x)]
    }

    n <- length(x)

    base <- tibble::tibble(
      n_sim = n,
      mean = if (n > 0L) mean(x) else NA_real_,
      variance = if (n > 1L) stats::var(x) else NA_real_,
      sd = if (n > 1L) stats::sd(x) else NA_real_,
      se_mean = if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_,
      min = if (n > 0L) min(x) else NA_real_,
      max = if (n > 0L) max(x) else NA_real_
    )

    q_values <- if (n > 0L) {
      as.numeric(stats::quantile(
        x,
        probs = probs,
        na.rm = na_rm,
        names = FALSE
      ))
    } else {
      rep(NA_real_, length(probs))
    }

    names(q_values) <- paste0("q", sprintf("%03.0f", probs * 1000))

    var_values <- if (n > 0L) {
      as.numeric(stats::quantile(
        x,
        probs = var_probs,
        na.rm = na_rm,
        names = FALSE
      ))
    } else {
      rep(NA_real_, length(var_probs))
    }

    names(var_values) <- paste0("VaR_", sprintf("%03.0f", var_probs * 1000))

    tvar_values <- vapply(seq_along(var_probs), function(j) {
      threshold <- var_values[[j]]

      if (is.na(threshold)) {
        return(NA_real_)
      }

      tail_x <- x[x >= threshold]

      if (length(tail_x) == 0L) {
        return(NA_real_)
      }

      mean(tail_x)
    }, numeric(1L))

    names(tvar_values) <- paste0("TVaR_", sprintf("%03.0f", var_probs * 1000))

    dplyr::bind_cols(
      base,
      tibble::as_tibble_row(as.list(q_values)),
      tibble::as_tibble_row(as.list(var_values)),
      tibble::as_tibble_row(as.list(tvar_values))
    )
  }

  work <- tibble::as_tibble(.data)

  if (is.null(group_cols)) {
    return(summarise_one(work))
  }

  grouped <- work |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_split()

  key_tbl <- work |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_keys()

  res <- lapply(grouped, summarise_one) |>
    dplyr::bind_rows()

  dplyr::bind_cols(key_tbl, res)
}

