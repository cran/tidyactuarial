#' Compute Monte Carlo net premiums for life contingencies
#'
#' Estimates net premiums from simulated present values of benefits and premium
#' annuities, while distinguishing the annualized premium from the amount paid
#' at each premium date.
#'
#' @param .data A data frame or tibble containing simulated present values.
#'   Usually obtained after applying \code{\link{mc_insurance}} and
#'   \code{\link{mc_annuity}} to the same simulated lifetime sample.
#' @param col_Z Character scalar naming the simulated present value of benefits.
#'   Default is \code{"pv_benefit"}.
#' @param col_Y Character scalar naming the simulated present value of the
#'   premium annuity. Default is \code{"pv_annuity"}.
#' @param col_P Character scalar naming the compatibility output containing the
#'   coefficient that directly multiplies \code{col_Y}. Default is \code{"P"}.
#'   For a normalized \eqn{k}-thly premium annuity constructed with
#'   \code{payment = 1 / k}, this coefficient is the annualized premium.
#' @param by Optional character vector of grouping columns. If supplied, the
#'   premium is estimated separately within each group. If \code{by = NULL}
#'   and \code{.data} is already grouped, the current grouping is used.
#' @param na_rm Logical scalar. If \code{TRUE}, simulations with a missing
#'   benefit or annuity present value are removed as complete pairs. The same
#'   simulation rows are therefore used in both means.
#' @param k Optional positive integer number of premium payments per year.
#'   By default it is inferred from \code{payments_per_year} or \code{k} in
#'   \code{.data}; if neither exists, \code{k = 1} is used.
#' @param annuity_payment Optional positive numeric scalar giving the amount
#'   paid at each date when \code{col_Y} was constructed. It is normally
#'   inferred from \code{payment_per_payment}, \code{payment}, or
#'   \code{payment_annualized} in the output of \code{\link{mc_annuity}}.
#' @param tol Nonnegative numeric tolerance for consistency checks.
#' @param ... Transitional compatibility for older calls using \code{data},
#'   \code{benefit_col}, \code{annuity_col}, and \code{premium_col}.
#'
#' @details
#' Let \eqn{Z} denote the simulated present value of benefits. Suppose the
#' simulated premium annuity is
#' \deqn{
#' Y_c = c\sum_j v^{t_j},
#' }
#' where \eqn{c} is the amount used at each of the \eqn{k} premium dates per
#' year. The Monte Carlo coefficient that directly multiplies this annuity is
#' \deqn{
#' \widehat q = \frac{\overline Z}{\overline{Y_c}}.
#' }
#'
#' The corresponding premium amount at each payment date is
#' \deqn{
#' \widehat P_{\mathrm{per\ payment}} = \widehat q c,
#' }
#' and the annualized premium is
#' \deqn{
#' \widehat P^{(k)} = k\widehat q c.
#' }
#'
#' The recommended premium-annuity normalization uses
#' \code{payment = 1 / k} in \code{\link{mc_annuity}}. In that case,
#' \deqn{
#' \widehat q = \widehat P^{(k)},
#' \qquad
#' \widehat P_{\mathrm{per\ payment}}
#' = \frac{\widehat P^{(k)}}{k}.
#' }
#'
#' If instead \code{payment = 1}, then \eqn{\widehat q} is the amount paid at
#' each premium date, while the annualized premium is \eqn{k\widehat q}.
#'
#' When \code{na_rm = TRUE}, missing values are removed jointly: a simulation
#' contributes only if both \code{col_Z} and \code{col_Y} are observed. This
#' prevents the numerator and denominator from being estimated from different
#' simulated samples.
#'
#' @return A tibble containing the original simulations and standardized
#' premium columns:
#' \describe{
#'   \item{premium_annualized}{Estimated annualized premium
#'   \eqn{\widehat P^{(k)}}.}
#'   \item{premium_per_payment}{Amount collected at each premium date.}
#'   \item{payments_per_year}{Premium payment frequency \code{k}.}
#'   \item{premium_annuity_scale}{Coefficient \eqn{\widehat q} multiplying
#'   \code{col_Y}.}
#'   \item{mc_equivalence_residual}{Difference
#'   \eqn{\overline Z-\widehat q\overline{Y_c}} within the estimation group.}
#' }
#'
#' The column selected by \code{col_P} stores
#' \code{premium_annuity_scale} for backward compatibility. Unless
#' \code{col_P = "premium"}, a compatibility column \code{premium} is also
#' created with the same scale coefficient.
#'
#' @seealso
#' \code{\link{mc_insurance}}, \code{\link{mc_annuity}},
#' \code{\link{mc_loss}}, \code{\link{mc_reserve}},
#' \code{\link{premium_x}}, \code{\link{premium_xy}}
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). \emph{Actuarial Mathematics}. Second Edition.
#' Society of Actuaries.
#'
#' @family monte-carlo
#'
#' @examples
#' simulated <- tibble::tibble(
#'   pv_benefit = c(800, 1200),
#'   pv_annuity = c(8, 12),
#'   payment_per_payment = 1 / 12,
#'   payment_annualized = 1,
#'   payments_per_year = 12
#' )
#'
#' simulated |>
#'   mc_premium()
#'
#' # Here P is annualized because the annuity uses payment = 1 / 12.
#'
#' # With payment = 1 at each monthly date, P is the monthly amount,
#' # while premium_annualized is twelve times P.
#' simulated_unit_payments <- tibble::tibble(
#'   pv_benefit = c(800, 1200),
#'   pv_annuity = c(96, 144),
#'   payment_per_payment = 1,
#'   payment_annualized = 12,
#'   payments_per_year = 12
#' )
#'
#' simulated_unit_payments |>
#'   mc_premium()
#'
#' @export
mc_premium <- function(
    .data = NULL,
    col_Z = "pv_benefit",
    col_Y = "pv_annuity",
    col_P = "P",
    by = NULL,
    na_rm = TRUE,
    k = NULL,
    annuity_payment = NULL,
    tol = 1e-10,
    ...
) {
  data_missing <- missing(.data)
  col_Z_missing <- missing(col_Z)
  col_Y_missing <- missing(col_Y)
  col_P_missing <- missing(col_P)

  dots <- list(...)

  dot_has <- function(name) {
    name %in% names(dots)
  }

  dot_get <- function(name) {
    dots[[name]]
  }

  allowed_old <- c(
    "data",
    "benefit_col",
    "annuity_col",
    "premium_col"
  )

  bad_dots <- setdiff(names(dots), allowed_old)

  if (length(bad_dots) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad_dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (dot_has("data")) {
    if (!data_missing) {
      stop("Provide only one of `.data` or deprecated `data`.", call. = FALSE)
    }

    .data <- dot_get("data")
  }

  if (dot_has("benefit_col")) {
    if (!col_Z_missing) {
      stop(
        "Provide only one of `col_Z` or deprecated `benefit_col`.",
        call. = FALSE
      )
    }

    col_Z <- dot_get("benefit_col")
  }

  if (dot_has("annuity_col")) {
    if (!col_Y_missing) {
      stop(
        "Provide only one of `col_Y` or deprecated `annuity_col`.",
        call. = FALSE
      )
    }

    col_Y <- dot_get("annuity_col")
  }

  if (dot_has("premium_col")) {
    if (!col_P_missing) {
      stop(
        "Provide only one of `col_P` or deprecated `premium_col`.",
        call. = FALSE
      )
    }

    col_P <- dot_get("premium_col")
  }

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_column(.data, col_Z, "col_Z")
  .mc_assert_numeric_column(.data, col_Y, "col_Y")
  .mc_assert_character_scalar(col_P, "col_P")

  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("`na_rm` must be a logical scalar.", call. = FALSE)
  }

  if (!is.null(by)) {
    if (!is.character(by) || anyNA(by)) {
      stop("`by` must be `NULL` or a character vector.", call. = FALSE)
    }

    if (!all(by %in% names(.data))) {
      stop("All columns supplied in `by` must exist in `.data`.", call. = FALSE)
    }
  }

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  reserved_names <- c(
    "premium_annualized",
    "premium_per_payment",
    "payments_per_year",
    "premium_annuity_scale",
    "mc_equivalence_residual"
  )

  if (col_P %in% reserved_names) {
    stop(
      "`col_P` must not use a standardized premium column name. ",
      "Use the default `P` or another compatibility name.",
      call. = FALSE
    )
  }

  Z_all <- .data[[col_Z]]
  Y_all <- .data[[col_Y]]

  if (any(is.infinite(Z_all))) {
    stop("`col_Z` must not contain infinite values.", call. = FALSE)
  }

  if (any(is.infinite(Y_all))) {
    stop("`col_Y` must not contain infinite values.", call. = FALSE)
  }

  if (any(Z_all < 0, na.rm = TRUE)) {
    stop("`col_Z` must contain nonnegative present values.", call. = FALSE)
  }

  if (any(Y_all < 0, na.rm = TRUE)) {
    stop("`col_Y` must contain nonnegative present values.", call. = FALSE)
  }

  active_by <- by

  if (is.null(active_by)) {
    active_by <- dplyr::group_vars(.data)
  }

  out <- dplyr::ungroup(.data)
  n_rows <- nrow(out)

  if (length(active_by) == 0L) {
    group_ids <- rep(1L, n_rows)
  } else {
    group_ids <- out |>
      dplyr::group_by(
        dplyr::across(dplyr::all_of(active_by)),
        .drop = FALSE
      ) |>
      dplyr::group_indices()
  }

  group_rows <- split(seq_len(n_rows), group_ids)

  premium_annuity_scale <- rep(NA_real_, n_rows)
  premium_per_payment <- rep(NA_real_, n_rows)
  premium_annualized <- rep(NA_real_, n_rows)
  payments_per_year <- rep(NA_integer_, n_rows)
  mc_equivalence_residual <- rep(NA_real_, n_rows)

  constant_in_group <- function(data, column) {
    values <- data[[column]]

    if (!is.numeric(values)) {
      stop("`", column, "` must be numeric.", call. = FALSE)
    }

    if (anyNA(values) || any(!is.finite(values))) {
      stop(
        "`", column, "` must contain finite nonmissing metadata within ",
        "each premium group.",
        call. = FALSE
      )
    }

    observed <- unique(values)

    if (length(observed) != 1L) {
      stop(
        "`", column, "` must be constant within each premium group.",
        call. = FALSE
      )
    }

    observed[[1L]]
  }

  validate_k <- function(value) {
    if (!is.numeric(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !is.finite(value) ||
        value < 1 ||
        abs(value - round(value)) > tol) {
      stop("`k` must be a single positive integer.", call. = FALSE)
    }

    as.integer(round(value))
  }

  if (!is.null(k)) {
    k <- validate_k(k)
  }

  if (!is.null(annuity_payment)) {
    if (!is.numeric(annuity_payment) ||
        length(annuity_payment) != 1L ||
        is.na(annuity_payment) ||
        !is.finite(annuity_payment) ||
        annuity_payment <= 0) {
      stop(
        "`annuity_payment` must be a single positive finite number.",
        call. = FALSE
      )
    }
  }

  for (rows in group_rows) {
    data_group <- out[rows, , drop = FALSE]

    k_from_data <- NULL

    if ("payments_per_year" %in% names(data_group)) {
      k_from_data <- constant_in_group(
        data_group,
        "payments_per_year"
      )
    }

    if ("k" %in% names(data_group)) {
      k_legacy <- constant_in_group(data_group, "k")

      if (!is.null(k_from_data) &&
          abs(k_from_data - k_legacy) > tol) {
        stop(
          "`payments_per_year` and `k` are inconsistent within a ",
          "premium group.",
          call. = FALSE
        )
      }

      if (is.null(k_from_data)) {
        k_from_data <- k_legacy
      }
    }

    k_group <- if (is.null(k)) {
      if (is.null(k_from_data)) 1L else validate_k(k_from_data)
    } else {
      if (!is.null(k_from_data) &&
          abs(k - k_from_data) > tol) {
        stop(
          "`k` does not agree with the payment frequency stored in ",
          "`.data`.",
          call. = FALSE
        )
      }

      k
    }

    payment_from_data <- NULL

    if ("payment_per_payment" %in% names(data_group)) {
      payment_from_data <- constant_in_group(
        data_group,
        "payment_per_payment"
      )
    } else if ("payment" %in% names(data_group)) {
      payment_from_data <- constant_in_group(
        data_group,
        "payment"
      )
    } else if ("payment_annualized" %in% names(data_group)) {
      payment_from_data <-
        constant_in_group(
          data_group,
          "payment_annualized"
        ) / k_group
    }

    payment_group <- if (is.null(annuity_payment)) {
      if (!is.null(payment_from_data)) {
        payment_from_data
      } else if (k_group == 1L) {
        1
      } else {
        stop(
          "The payment amount used to construct `col_Y` cannot be ",
          "inferred. Supply `annuity_payment`, or retain the payment ",
          "metadata returned by `mc_annuity()`.",
          call. = FALSE
        )
      }
    } else {
      if (!is.null(payment_from_data) &&
          abs(annuity_payment - payment_from_data) >
            tol * max(1, abs(payment_from_data))) {
        stop(
          "`annuity_payment` does not agree with the payment amount ",
          "stored in `.data`.",
          call. = FALSE
        )
      }

      annuity_payment
    }

    if (!is.numeric(payment_group) ||
        length(payment_group) != 1L ||
        is.na(payment_group) ||
        !is.finite(payment_group) ||
        payment_group <= 0) {
      stop(
        "The annuity payment amount must be a single positive finite ",
        "number within each premium group.",
        call. = FALSE
      )
    }

    if ("payment_annualized" %in% names(data_group)) {
      annualized_basis <- constant_in_group(
        data_group,
        "payment_annualized"
      )

      if (abs(annualized_basis - k_group * payment_group) >
          tol * max(1, abs(annualized_basis))) {
        stop(
          "`payment_annualized` is inconsistent with the payment ",
          "frequency and per-payment amount.",
          call. = FALSE
        )
      }
    }

    Z <- data_group[[col_Z]]
    Y <- data_group[[col_Y]]

    if (isTRUE(na_rm)) {
      keep <- !is.na(Z) & !is.na(Y)

      if (!any(keep)) {
        stop(
          "No complete simulated benefit-annuity pairs are available ",
          "within a premium group.",
          call. = FALSE
        )
      }

      Z <- Z[keep]
      Y <- Y[keep]
    } else if (anyNA(Z) || anyNA(Y)) {
      stop(
        "Missing simulated values are present. Use `na_rm = TRUE` to ",
        "remove complete benefit-annuity pairs.",
        call. = FALSE
      )
    }

    mean_Z <- mean(Z)
    mean_Y <- mean(Y)

    if (!is.finite(mean_Z)) {
      stop(
        "The mean simulated benefit present value could not be computed.",
        call. = FALSE
      )
    }

    if (!is.finite(mean_Y) || mean_Y <= 0) {
      stop(
        "The mean simulated premium annuity present value must be ",
        "positive.",
        call. = FALSE
      )
    }

    scale_group <- mean_Z / mean_Y
    per_payment_group <- scale_group * payment_group
    annualized_group <- k_group * per_payment_group
    residual_group <- mean_Z - scale_group * mean_Y

    premium_annuity_scale[rows] <- scale_group
    premium_per_payment[rows] <- per_payment_group
    premium_annualized[rows] <- annualized_group
    payments_per_year[rows] <- k_group
    mc_equivalence_residual[rows] <- residual_group
  }

  out[[col_P]] <- premium_annuity_scale

  out <- out |>
    dplyr::mutate(
      premium_annualized = premium_annualized,
      premium_per_payment = premium_per_payment,
      payments_per_year = payments_per_year,
      premium_annuity_scale = premium_annuity_scale,
      mc_equivalence_residual = mc_equivalence_residual
    )

  if (!identical(col_P, "premium")) {
    if ("premium" %in% names(out)) {
      existing <- out[["premium"]]

      if (!is.numeric(existing) ||
          anyNA(existing) ||
          any(!is.finite(existing)) ||
          any(
            abs(existing - premium_annuity_scale) >
              tol * pmax(1, abs(premium_annuity_scale))
          )) {
        stop(
          "The existing `premium` column is inconsistent with the ",
          "estimated premium-annuity scale.",
          call. = FALSE
        )
      }
    } else {
      out <- out |>
        dplyr::mutate(
          premium = premium_annuity_scale
        )
    }
  }

  out
}
