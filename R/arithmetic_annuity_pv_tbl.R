#' Present value of an arithmetic progression annuity
#'
#' Computes the present value of an annuity whose payments follow an
#' arithmetic progression.
#'
#' The first payment is \code{amount}, and each subsequent payment changes by
#' \code{step}. The annuity may be either:
#' \itemize{
#'   \item annuity-immediate (\code{timing = "immediate"} or \code{"vencida"}),
#'   \item annuity-due (\code{timing = "due"} or \code{"anticipada"}),
#'   \item finite (\code{perpetuity = FALSE}),
#'   \item or perpetual (\code{perpetuity = TRUE}).
#' }
#'
#' This is a tibble-first mutate-style function: each input row is one case.
#'
#' Assumptions:
#' \itemize{
#'   \item Time is discrete.
#'   \item \code{i} is the effective interest rate per payment period.
#'   \item If \code{perpetuity = FALSE}, \code{n} is the number of payments.
#'   \item If \code{perpetuity = TRUE}, \code{n} is ignored.
#'   \item All payments must remain nonnegative: for finite annuities the
#'         final payment \code{amount + (n - 1) * step} must be >= 0; for
#'         perpetuities, \code{step >= 0} is required.
#' }
#'
#' @param .data A data.frame or tibble. If \code{NULL}, inputs must be supplied
#'   directly as scalars or equal-length vectors.
#' @param amount Numeric first payment when \code{.data = NULL}.
#' @param step Numeric arithmetic increment per payment period when \code{.data = NULL}.
#' @param n Number of payments when \code{.data = NULL} and \code{perpetuity = FALSE}.
#' @param i Effective interest rate per payment period when \code{.data = NULL}.
#' @param timing Payment timing. Accepted values are \code{"immediate"}, \code{"due"},
#'   \code{"vencida"}, and \code{"anticipada"}.
#' @param perpetuity Logical; if \code{TRUE}, computes a perpetuity instead of a
#'   finite annuity.
#' @param col_amount Name of the first-payment column.
#' @param col_step Name of the arithmetic-step column.
#' @param col_n Name of the number-of-payments column.
#' @param col_i Name of the interest-rate column.
#' @param .out Name of the output column containing present value.
#' @param .keep One of \code{"all"}, \code{"used"}, or \code{"none"}.
#' @param .na NA handling policy: \code{"propagate"}, \code{"error"}, or \code{"drop"}.
#'
#' @return A tibble with a new numeric column named by \code{.out}.
#'
#' @details
#' For a finite annuity-immediate with \code{n} payments, first payment \code{P},
#' arithmetic step \code{Q}, and effective rate \code{i} per period:
#' \deqn{PV = P \, a_{\overline{n|}} + Q \, \frac{a_{\overline{n|}} - n v^n}{i}}{PV = P a_n| + Q a_n| - n v^ni}
#' where \eqn{a_{\overline{n|}} = (1 - v^n)/i} and \eqn{v = 1/(1+i)}.
#'
#' For an annuity-due, replace \eqn{a_{\overline{n|}}} by
#' \eqn{\ddot{a}_{\overline{n|}} = (1+i) a_{\overline{n|}}} and divide the step
#' component by \eqn{d = i/(1+i)} instead of \eqn{i}.
#'
#' For a perpetuity-immediate (\eqn{i > 0}):
#' \deqn{PV = \frac{P}{i} + \frac{Q}{i^2}.}{PV = (P)/(i) + (Q)/(i^2).}
#'
#' For a perpetuity-due:
#' \deqn{PV = \frac{P}{d} + \frac{Q}{i \, d} = (1+i)\left(\frac{P}{i} + \frac{Q}{i^2}\right).}{PV = (P)/(d) + (Q)/(i d) = (1+i)((P)/(i) + (Q)/(i^2)).}
#'
#' When \eqn{i = 0} (finite case only), the present value simplifies to:
#' \deqn{PV = n P + Q \frac{n(n-1)}{2}.}{PV = n P + Q (n(n-1))/(2).}
#'
#' @seealso \code{\link{arithmetic_annuity_av_tbl}}, \code{\link{a_angle}}, \code{\link{s_angle}}
#'
#' @family annuities
#'
#' @examples
#' # Simple example: finite annuity-immediate
#' arithmetic_annuity_pv_tbl(
#'   amount = 100,
#'   step = 5,
#'   n = 5,
#'   i = 0.05,
#'   timing = "immediate"
#' )
#'
#' # Medium example: finite annuity-due for multiple rows
#' cases <- tibble::tibble(
#'   amount = c(100, 200, 50),
#'   step   = c(5, 10, 0),
#'   n      = c(5, 4, 10),
#'   i      = c(0.05, 0.04, 0.03)
#' )
#'
#' arithmetic_annuity_pv_tbl(
#'   cases,
#'   timing = "due",
#'   perpetuity = FALSE,
#'   .out = "pv_due"
#' )
#'
#' # Perpetuity-due
#' arithmetic_annuity_pv_tbl(
#'   amount = 1000,
#'   step = 40,
#'   i = 0.10,
#'   timing = "due",
#'   perpetuity = TRUE
#' )
#'
#' @references
#' Marcel B. Finan, *A Basic Course in the Theory of Interest and
#' Derivatives Markets: A Preparation for the Actuarial Exam FM/2*.
#'
#' Kellison, S. G. *The Theory of Interest*.
#'
#' @export
arithmetic_annuity_pv_tbl <- function(
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
    .out = "pv",
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

  if (!is.character(.out) || length(.out) != 1L || is.na(.out) || !nzchar(.out)) {
    abort("`.out` must be a single non-empty string.")
  }

  zero_tol <- sqrt(.Machine$double.eps)

  if (is.null(.data)) {
    required_inputs <- list(amount = amount, step = step, i = i)
    if (!perpetuity) {
      required_inputs$n <- n
    }

    if (any(vapply(required_inputs, is.null, logical(1)))) {
      abort("When `.data = NULL`, supply `amount`, `step`, `i`, and `n` if `perpetuity = FALSE`.")
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
      i      = recycle1(i, n_rows)
    )

    if (!perpetuity) {
      data_in$n <- recycle1(n, n_rows)
    }

    amount_name <- "amount"
    step_name <- "step"
    i_name <- "i"
    n_name <- "n"
    original_names <- names(data_in)
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)
    original_names <- names(data_in)

    required_cols <- c(col_amount, col_step, col_i)
    if (!perpetuity) {
      required_cols <- c(required_cols, col_n)
    }

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
    i_name <- col_i
    n_name <- col_n
  }

  used_names <- c(amount_name, step_name, i_name)
  if (!perpetuity) {
    used_names <- c(used_names, n_name)
  }
  used_names <- unique(used_names)

  data_in$.row_id <- seq_len(nrow(data_in))
  data_in$.amount <- data_in[[amount_name]]
  data_in$.step <- data_in[[step_name]]
  data_in$.i <- data_in[[i_name]]
  data_in$.n <- if (perpetuity) NA_real_ else data_in[[n_name]]

  if (!is.numeric(data_in$.amount)) {
    abort(paste0("`", amount_name, "` must be numeric."))
  }
  if (!is.numeric(data_in$.step)) {
    abort(paste0("`", step_name, "` must be numeric."))
  }
  if (!is.numeric(data_in$.i)) {
    abort(paste0("`", i_name, "` must be numeric."))
  }
  if (!perpetuity && !is.numeric(data_in$.n)) {
    abort(paste0("`", n_name, "` must be numeric/integer."))
  }

  data_in$.bad_na <- is.na(data_in$.amount) | is.na(data_in$.step) | is.na(data_in$.i)
  if (!perpetuity) {
    data_in$.bad_na <- data_in$.bad_na | is.na(data_in$.n)
  }

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

    if (!perpetuity) {
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
    } else {
      bad_perp_i <- data_valid$.i <= 0
      if (any(bad_perp_i)) {
        bad_rows <- data_valid$.row_id[bad_perp_i]
        abort(
          paste0(
            "Perpetuities require `", i_name, "` > 0. Problem at row(s): ",
            paste(bad_rows, collapse = ", "),
            "."
          )
        )
      }

      bad_perp_step <- data_valid$.step < 0
      if (any(bad_perp_step)) {
        bad_rows <- data_valid$.row_id[bad_perp_step]
        abort(
          paste0(
            "For `perpetuity = TRUE`, `", step_name, "` must be >= 0. Problem at row(s): ",
            paste(bad_rows, collapse = ", "),
            "."
          )
        )
      }
    }
  }

  compute_pv <- function(amount, step, n, i, timing, perpetuity) {
    if (perpetuity) {
      pv_imm <- amount / i + step / (i^2)

      if (timing == "immediate") {
        return(pv_imm)
      } else {
        return((1 + i) * pv_imm)
      }
    }

    n <- as.integer(round(n))

    if (n == 0L) {
      return(0)
    }

    if (abs(i) <= zero_tol) {
      return(n * amount + step * n * (n - 1) / 2)
    }

    v <- 1 / (1 + i)
    a_n <- (1 - v^n) / i
    step_component <- (a_n - n * v^n) / i

    if (timing == "immediate") {
      return(amount * a_n + step * step_component)
    } else {
      # Due = (1 + i) * immediate
      return((1 + i) * (amount * a_n + step * step_component))
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
        compute_pv(
          amount = data_work$.amount[[r]],
          step = data_work$.step[[r]],
          n = data_work$.n[[r]],
          i = data_work$.i[[r]],
          timing = timing,
          perpetuity = perpetuity
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
