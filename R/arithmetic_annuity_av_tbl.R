#' Accumulated value of an arithmetic progression annuity
#'
#' Computes the accumulated value of an annuity whose payments follow an
#' arithmetic progression.
#'
#' The first payment is \code{amount}, and each subsequent payment changes by
#' \code{step}. The annuity may be either:
#' \itemize{
#'   \item annuity-immediate (\code{timing = "immediate"} or \code{"vencida"}),
#'   \item annuity-due (\code{timing = "due"} or \code{"anticipada"}).
#' }
#'
#' This is a tibble-first mutate-style function: each input row is one case.
#'
#' Assumptions:
#' \itemize{
#'   \item Time is discrete.
#'   \item \code{i} is the effective interest rate per payment period.
#'   \item \code{n} is the number of payments.
#'   \item All payments must remain nonnegative: the final payment
#'         \code{amount + (n - 1) * step} must be >= 0.
#'   \item \code{perpetuity = TRUE} is not supported because the accumulated value
#'         of a perpetuity is not finite.
#' }
#'
#' @param .data A data.frame or tibble. If \code{NULL}, inputs must be supplied
#'   directly as scalars or equal-length vectors.
#' @param amount Numeric first payment when \code{.data = NULL}.
#' @param step Numeric arithmetic increment per payment period when \code{.data = NULL}.
#' @param n Number of payments when \code{.data = NULL}.
#' @param i Effective interest rate per payment period when \code{.data = NULL}.
#' @param timing Payment timing. Accepted values are \code{"immediate"}, \code{"due"},
#'   \code{"vencida"}, and \code{"anticipada"}.
#' @param perpetuity Logical. Must be \code{FALSE}. Included only for interface
#'   consistency with the present-value companion function.
#' @param col_amount Name of the first-payment column.
#' @param col_step Name of the arithmetic-step column.
#' @param col_n Name of the number-of-payments column.
#' @param col_i Name of the interest-rate column.
#' @param .out Name of the output column containing accumulated value.
#' @param .keep One of \code{"all"}, \code{"used"}, or \code{"none"}.
#' @param .na NA handling policy: \code{"propagate"}, \code{"error"}, or \code{"drop"}.
#'
#' @return A tibble with a new numeric column named by \code{.out}.
#'
#' @details
#' For an annuity-immediate with \code{n} payments, first payment \code{P}, step \code{Q},
#' and effective rate \code{i} per period, the accumulated value at time \code{n} is:
#' \deqn{AV = P \, s_{\overline{n|}} + Q \, \frac{s_{\overline{n|}} - n}{i}}{AV = P s_n| + Q s_n| - ni}
#' where \eqn{s_{\overline{n|}} = \frac{(1+i)^n - 1}{i}}.
#'
#' For an annuity-due the result is \eqn{(1+i)} times the immediate value.
#'
#' When \eqn{i = 0}, the accumulated value simplifies to:
#' \deqn{AV = n P + Q \frac{n(n-1)}{2}.}{AV = n P + Q (n(n-1))/(2).}
#'
#' @seealso \code{\link{arithmetic_annuity_pv_tbl}}, \code{\link{a_angle}}, \code{\link{s_angle}}
#'
#' @family annuities
#'
#' @examples
#' arithmetic_annuity_av_tbl(
#'   amount = 100,
#'   step = 5,
#'   n = 5,
#'   i = 0.05,
#'   timing = "immediate"
#' )
#'
#' cases <- tibble::tibble(
#'   amount = c(100, 200, 50),
#'   step   = c(5, 10, 0),
#'   n      = c(5, 4, 10),
#'   i      = c(0.05, 0.04, 0.03)
#' )
#'
#' arithmetic_annuity_av_tbl(
#'   cases,
#'   timing = "due",
#'   .out = "av_due"
#' )
#'
#' @references
#' Marcel B. Finan, *A Basic Course in the Theory of Interest and
#' Derivatives Markets: A Preparation for the Actuarial Exam FM/2*.
#'
#' Kellison, S. G. *The Theory of Interest*.
#'
#' @export
arithmetic_annuity_av_tbl <- function(
    .data = NULL,
    amount = NULL,
    step = NULL,
    n = NULL,
    i = NULL,
    timing = c("immediate", "due", "vencida", "anticipada"),
    perpetuity = FALSE,
    col_amount = "amount",
    col_step = "step",
    col_n = "n",
    col_i = "i",
    .out = "av",
    .keep = c("all", "used", "none"),
    .na = c("propagate", "error", "drop")
) {
  .keep <- match.arg(.keep)
  .na <- match.arg(.na)

  abort <- function(msg) rlang::abort(msg)

  timing <- .normalize_annuity_timing(timing, abort)

  if (!is.logical(perpetuity) || length(perpetuity) != 1L || is.na(perpetuity)) {
    abort("`perpetuity` must be TRUE or FALSE.")
  }

  if (isTRUE(perpetuity)) {
    abort("`perpetuity = TRUE` is not supported here because the accumulated value of a perpetuity is not finite.")
  }

  if (!is.character(.out) || length(.out) != 1L || is.na(.out) || !nzchar(.out)) {
    abort("`.out` must be a single non-empty string.")
  }

  zero_tol <- sqrt(.Machine$double.eps)

  if (is.null(.data)) {
    required_inputs <- list(amount = amount, step = step, n = n, i = i)

    if (any(vapply(required_inputs, is.null, logical(1)))) {
      abort("When `.data = NULL`, supply `amount`, `step`, `n`, and `i`.")
    }

    lens <- vapply(required_inputs, length, integer(1))
    n_rows <- max(lens)

    if (!all(lens %in% c(1L, n_rows))) {
      abort("When `.data = NULL`, inputs must be scalars or have a common length.")
    }

    recycle1 <- function(x, n_rows) {
      if (length(x) == 1L) rep(x, n_rows) else x
    }

    data_in <- tibble::tibble(
      amount = recycle1(amount, n_rows),
      step   = recycle1(step, n_rows),
      n      = recycle1(n, n_rows),
      i      = recycle1(i, n_rows)
    )

    amount_name <- "amount"
    step_name <- "step"
    n_name <- "n"
    i_name <- "i"
    original_names <- names(data_in)
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)
    original_names <- names(data_in)

    required_cols <- c(col_amount, col_step, col_n, col_i)
    missing_cols <- setdiff(required_cols, names(data_in))

    if (length(missing_cols) > 0L) {
      abort(
        paste0(
          "Missing required column(s): ",
          paste(sprintf("`%s`", missing_cols), collapse = ", "),
          ". Pass the correct names via the corresponding `col_*` arguments."
        )
      )
    }

    amount_name <- col_amount
    step_name <- col_step
    n_name <- col_n
    i_name <- col_i
  }

  used_names <- unique(c(amount_name, step_name, n_name, i_name))

  data_in$.row_id <- seq_len(nrow(data_in))
  data_in$.amount <- data_in[[amount_name]]
  data_in$.step <- data_in[[step_name]]
  data_in$.n <- data_in[[n_name]]
  data_in$.i <- data_in[[i_name]]

  if (!is.numeric(data_in$.amount)) {
    abort(paste0("`", amount_name, "` must be numeric."))
  }
  if (!is.numeric(data_in$.step)) {
    abort(paste0("`", step_name, "` must be numeric."))
  }
  if (!is.numeric(data_in$.n)) {
    abort(paste0("`", n_name, "` must be numeric/integer."))
  }
  if (!is.numeric(data_in$.i)) {
    abort(paste0("`", i_name, "` must be numeric."))
  }

  data_in$.bad_na <- is.na(data_in$.amount) | is.na(data_in$.step) |
    is.na(data_in$.n) | is.na(data_in$.i)

  if (.na == "error" && any(data_in$.bad_na)) {
    bad_rows <- data_in$.row_id[data_in$.bad_na]
    abort(
      paste0(
        "Missing values found in required inputs at row(s): ",
        paste(bad_rows, collapse = ", "),
        "."
      )
    )
  }

  if (.na == "drop" && any(data_in$.bad_na)) {
    bad_rows <- data_in$.row_id[data_in$.bad_na]
    warning(
      paste0(
        "Dropping row(s) with missing values in required inputs: ",
        paste(bad_rows, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
    data_work <- data_in[!data_in$.bad_na, , drop = FALSE]
  } else {
    data_work <- data_in
  }

  data_valid <- data_work[!data_work$.bad_na, , drop = FALSE]

  if (nrow(data_valid) > 0L) {
    bad_amount <- !is.finite(data_valid$.amount)
    if (any(bad_amount)) {
      bad_rows <- data_valid$.row_id[bad_amount]
      abort(
        paste0(
          "`", amount_name, "` must be finite. Problem at row(s): ",
          paste(bad_rows, collapse = ", "),
          "."
        )
      )
    }

    bad_step <- !is.finite(data_valid$.step)
    if (any(bad_step)) {
      bad_rows <- data_valid$.row_id[bad_step]
      abort(
        paste0(
          "`", step_name, "` must be finite. Problem at row(s): ",
          paste(bad_rows, collapse = ", "),
          "."
        )
      )
    }

    bad_n <- !is.finite(data_valid$.n) |
      data_valid$.n < 0 |
      abs(data_valid$.n - round(data_valid$.n)) > 1e-10

    if (any(bad_n)) {
      bad_rows <- data_valid$.row_id[bad_n]
      abort(
        paste0(
          "`", n_name, "` must be a finite integer >= 0. Problem at row(s): ",
          paste(bad_rows, collapse = ", "),
          "."
        )
      )
    }

    bad_i <- !is.finite(data_valid$.i) | data_valid$.i <= -1
    if (any(bad_i)) {
      bad_rows <- data_valid$.row_id[bad_i]
      abort(
        paste0(
          "`", i_name, "` must be finite and > -1. Problem at row(s): ",
          paste(bad_rows, collapse = ", "),
          "."
        )
      )
    }

    n_int <- as.integer(round(data_valid$.n))
    last_payment <- data_valid$.amount + pmax(n_int - 1L, 0L) * data_valid$.step
    bad_final_payment <- last_payment < 0

    if (any(bad_final_payment)) {
      bad_rows <- data_valid$.row_id[bad_final_payment]
      abort(
        paste0(
          "The implied final payment is negative. Arithmetic-payment annuities must remain nonnegative. Problem at row(s): ",
          paste(bad_rows, collapse = ", "),
          "."
        )
      )
    }
  }

  compute_av <- function(amount, step, n, i, timing) {
    n <- as.integer(round(n))

    if (n == 0L) {
      return(0)
    }

    if (abs(i) <= zero_tol) {
      return(n * amount + step * n * (n - 1) / 2)
    }

    s_n <- ((1 + i)^n - 1) / i
    av_immediate <- amount * s_n + step * (s_n - n) / i

    if (timing == "immediate") {
      return(av_immediate)
    } else {
      return((1 + i) * av_immediate)
    }
  }

  if (nrow(data_work) == 0L) {
    if (.keep == "all") {
      out <- data_work[, original_names, drop = FALSE]
    } else if (.keep == "used") {
      out <- data_work[, used_names, drop = FALSE]
    } else {
      out <- tibble::tibble()
    }

    out[[.out]] <- numeric(0)
    return(tibble::as_tibble(out))
  }

  # --- Vectorized computation ---
  out_val <- rep(NA_real_, nrow(data_work))
  computable <- !data_work$.bad_na

  if (any(computable)) {
    out_val[computable] <- vapply(
      which(computable),
      function(r) {
        compute_av(
          amount = data_work$.amount[[r]],
          step = data_work$.step[[r]],
          n = data_work$.n[[r]],
          i = data_work$.i[[r]],
          timing = timing
        )
      },
      numeric(1)
    )
  }

  if (.keep == "all") {
    out <- data_work[, original_names, drop = FALSE]
  } else if (.keep == "used") {
    out <- data_work[, used_names, drop = FALSE]
  } else {
    out <- tibble::tibble()
  }

  out[[.out]] <- out_val
  tibble::as_tibble(out)
}


# --- Internal utility: normalize timing strings (shared across annuity functions) ---
# This could be moved to a utils.R file to avoid duplication with _pv_tbl.
.normalize_annuity_timing <- function(timing, abort = stop) {
  if (!is.character(timing) || length(timing) != 1L || is.na(timing)) {
    abort("`timing` must be a single non-missing string.")
  }

  x <- tolower(trimws(timing))

  if (x %in% c("immediate", "vencida")) {
    return("immediate")
  }

  if (x %in% c("due", "anticipada")) {
    return("due")
  }

  abort("`timing` must be one of: \"immediate\", \"due\", \"vencida\", \"anticipada\".")
}
