#' Compute Monte Carlo loss random variables for life contingencies
#'
#' Constructs simulated actuarial losses from present values of benefits and
#' premium-payment annuities while keeping annualized premiums distinct from
#' amounts paid at each premium date.
#'
#' @param .data A data frame or tibble containing simulated present values of
#'   benefits and premium annuities.
#' @param col_Z Character scalar naming the simulated present value of benefits.
#'   Default is \code{"pv_benefit"}.
#' @param col_Y Character scalar naming the simulated present value of the
#'   premium annuity. Default is \code{"pv_annuity"}.
#' @param col_P Character scalar naming the premium input column. When omitted,
#'   the function searches, in order, for \code{premium_annualized},
#'   \code{premium_per_payment}, \code{P}, and \code{premium}.
#' @param col_L Character scalar naming the simulated loss output. Default is
#'   \code{"L"}.
#' @param P Optional nonnegative numeric scalar supplied directly instead of
#'   reading a premium column.
#' @param premium_unit Unit of \code{P} or the selected premium column:
#'   \code{"annualized"}, \code{"per_payment"}, or
#'   \code{"annuity_scale"}. The latter is the coefficient that directly
#'   multiplies \code{col_Y}. With \code{"auto"}, explicit modern column names
#'   determine the unit, while \code{P}, \code{premium}, and a direct
#'   \code{P} argument retain their historical annuity-scale interpretation.
#' @param k Optional positive integer number of premium payments per year.
#'   By default it is inferred from \code{payments_per_year} or \code{k} in
#'   \code{.data}; if neither exists, \code{k = 1} is used.
#' @param annuity_payment Optional positive numeric scalar giving the payment
#'   amount used to construct \code{col_Y}. It is normally inferred from
#'   \code{payment_per_payment}, \code{payment}, or
#'   \code{payment_annualized} in the output of \code{\link{mc_annuity}}.
#' @param tol Nonnegative numeric tolerance for consistency checks.
#' @param ... Transitional compatibility for older calls using \code{data},
#'   \code{benefit_col}, \code{annuity_col}, \code{premium_col},
#'   \code{loss_col}, and \code{premium}.
#'
#' @details
#' Let \eqn{Z} be the present value random variable of benefits. Suppose
#' \code{col_Y} contains
#' \deqn{
#' Y_c = c\sum_j v^{t_j},
#' }
#' where \eqn{c} is the amount used in \code{\link{mc_annuity}} at each of the
#' \eqn{k} premium dates per year.
#'
#' If \eqn{P^{(k)}} is the annualized premium, the actual installment is
#' \eqn{P^{(k)}/k}, and the simulated present value of premiums is
#' \deqn{
#' \Pi =
#' \frac{P^{(k)}}{k c}Y_c.
#' }
#' Hence the loss at issue is
#' \deqn{
#' L = Z - \Pi.
#' }
#'
#' The recommended normalized premium annuity uses
#' \code{payment = 1 / k} in \code{\link{mc_annuity}}. Then \eqn{c=1/k},
#' \eqn{Y_c=Y^{(k)}}, and the expression simplifies to
#' \deqn{
#' L = Z - P^{(k)}Y^{(k)}.
#' }
#'
#' This convention agrees with \code{\link{premium_x}},
#' \code{\link{premium_xy}}, \code{\link{reserve_x}}, and
#' \code{\link{reserve_xy}}: the premium is annualized, while the amount
#' collected at each date is the annualized premium divided by \code{k}.
#'
#' For backward compatibility, legacy columns \code{P} and \code{premium}
#' are interpreted as coefficients that multiply the supplied annuity present
#' value directly. If such a coefficient is \eqn{q}, then
#' \deqn{
#' P_{\mathrm{per\ payment}} = q c,
#' \qquad
#' P^{(k)} = k q c.
#' }
#'
#' @return A tibble containing the original simulation and standardized columns:
#' \describe{
#'   \item{premium_annualized}{Annualized premium \eqn{P^{(k)}}.}
#'   \item{premium_per_payment}{Amount collected at each premium date.}
#'   \item{payments_per_year}{Premium payment frequency \code{k}.}
#'   \item{premium_annuity_scale}{Coefficient multiplying \code{col_Y}.}
#'   \item{pv_premiums}{Simulated present value of premium income.}
#'   \item{L}{Simulated loss, or the name supplied through \code{col_L}.}
#' }
#'
#' A compatibility column \code{loss} is synchronized with \code{col_L}.
#'
#' @seealso
#' \code{\link{mc_annuity}}, \code{\link{mc_insurance}},
#' \code{\link{mc_premium}}, \code{\link{mc_reserve}},
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
#'   pv_benefit = c(900, 700, 1100),
#'   pv_annuity = c(8, 7, 9),
#'   payment_per_payment = 1 / 12,
#'   payment_annualized = 1,
#'   payments_per_year = 12,
#'   premium_annualized = 120
#' )
#'
#' simulated |>
#'   mc_loss()
#'
#' # A premium supplied directly as an annualized amount
#' simulated |>
#'   mc_loss(
#'     P = 120,
#'     premium_unit = "annualized"
#'   )
#'
#' # Historical coefficient multiplying the annuity PV directly
#' tibble::tibble(
#'   pv_benefit = c(900, 700),
#'   pv_annuity = c(8, 7),
#'   P = 100
#' ) |>
#'   mc_loss()
#'
#' @export
mc_loss <- function(
    .data = NULL,
    col_Z = "pv_benefit",
    col_Y = "pv_annuity",
    col_P = "P",
    col_L = "L",
    P = NULL,
    premium_unit = c(
      "auto",
      "annualized",
      "per_payment",
      "annuity_scale"
    ),
    k = NULL,
    annuity_payment = NULL,
    tol = 1e-10,
    ...
) {
  data_missing <- missing(.data)
  col_Z_missing <- missing(col_Z)
  col_Y_missing <- missing(col_Y)
  col_P_missing <- missing(col_P)
  col_L_missing <- missing(col_L)
  P_missing <- missing(P)

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
    "premium_col",
    "loss_col",
    "premium"
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
    data_missing <- FALSE
  }

  if (dot_has("benefit_col")) {
    if (!col_Z_missing) {
      stop(
        "Provide only one of `col_Z` or deprecated `benefit_col`.",
        call. = FALSE
      )
    }

    col_Z <- dot_get("benefit_col")
    col_Z_missing <- FALSE
  }

  if (dot_has("annuity_col")) {
    if (!col_Y_missing) {
      stop(
        "Provide only one of `col_Y` or deprecated `annuity_col`.",
        call. = FALSE
      )
    }

    col_Y <- dot_get("annuity_col")
    col_Y_missing <- FALSE
  }

  if (dot_has("premium_col")) {
    if (!col_P_missing) {
      stop(
        "Provide only one of `col_P` or deprecated `premium_col`.",
        call. = FALSE
      )
    }

    col_P <- dot_get("premium_col")
    col_P_missing <- FALSE
  }

  if (dot_has("loss_col")) {
    if (!col_L_missing) {
      stop(
        "Provide only one of `col_L` or deprecated `loss_col`.",
        call. = FALSE
      )
    }

    col_L <- dot_get("loss_col")
    col_L_missing <- FALSE
  }

  if (dot_has("premium")) {
    if (!P_missing) {
      stop(
        "Provide only one of `P` or deprecated `premium`.",
        call. = FALSE
      )
    }

    P <- dot_get("premium")
    P_missing <- FALSE
  }

  premium_unit <- match.arg(premium_unit)

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_column(.data, col_Z, "col_Z")
  .mc_assert_numeric_column(.data, col_Y, "col_Y")
  .mc_assert_character_scalar(col_P, "col_P")
  .mc_assert_character_scalar(col_L, "col_L")

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  n_rows <- nrow(.data)
  Z <- .data[[col_Z]]
  Y <- .data[[col_Y]]

  constant_column <- function(column, label) {
    values <- .data[[column]]

    if (!is.numeric(values)) {
      stop("`", column, "` must be numeric.", call. = FALSE)
    }

    observed <- unique(values[!is.na(values)])

    if (length(observed) != 1L || !is.finite(observed[[1L]])) {
      stop(
        "`", column, "` must contain one constant finite value for the ",
        "simulated contract.",
        call. = FALSE
      )
    }

    observed[[1L]]
  }

  k_from_data <- NULL

  if ("payments_per_year" %in% names(.data)) {
    k_from_data <- constant_column(
      "payments_per_year",
      "payments_per_year"
    )
  }

  if ("k" %in% names(.data)) {
    k_legacy <- constant_column("k", "k")

    if (!is.null(k_from_data) && abs(k_from_data - k_legacy) > tol) {
      stop(
        "`payments_per_year` and `k` are inconsistent in `.data`.",
        call. = FALSE
      )
    }

    if (is.null(k_from_data)) {
      k_from_data <- k_legacy
    }
  }

  if (is.null(k)) {
    k <- if (is.null(k_from_data)) 1L else k_from_data
  } else if (!is.null(k_from_data) &&
             is.numeric(k) &&
             length(k) == 1L &&
             is.finite(k) &&
             abs(k - k_from_data) > tol) {
    stop(
      "`k` does not agree with the payment frequency stored in `.data`.",
      call. = FALSE
    )
  }

  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k < 1 ||
      abs(k - round(k)) > tol) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  k <- as.integer(round(k))

  payment_from_data <- NULL

  if ("payment_per_payment" %in% names(.data)) {
    payment_from_data <- constant_column(
      "payment_per_payment",
      "payment_per_payment"
    )
  } else if ("payment" %in% names(.data)) {
    payment_from_data <- constant_column("payment", "payment")
  } else if ("payment_annualized" %in% names(.data)) {
    payment_from_data <-
      constant_column("payment_annualized", "payment_annualized") / k
  }

  if (is.null(annuity_payment)) {
    if (!is.null(payment_from_data)) {
      annuity_payment <- payment_from_data
    } else if (k == 1L) {
      annuity_payment <- 1
    } else {
      stop(
        "The payment amount used to construct `col_Y` cannot be inferred. ",
        "Supply `annuity_payment`, or retain the payment metadata returned ",
        "by `mc_annuity()`.",
        call. = FALSE
      )
    }
  } else if (!is.null(payment_from_data) &&
             is.numeric(annuity_payment) &&
             length(annuity_payment) == 1L &&
             is.finite(annuity_payment) &&
             abs(annuity_payment - payment_from_data) >
               tol * max(1, abs(payment_from_data))) {
    stop(
      "`annuity_payment` does not agree with the payment amount stored ",
      "in `.data`.",
      call. = FALSE
    )
  }

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

  if ("payment_annualized" %in% names(.data)) {
    annualized_basis <- constant_column(
      "payment_annualized",
      "payment_annualized"
    )

    if (abs(annualized_basis - k * annuity_payment) >
        tol * max(1, abs(annualized_basis))) {
      stop(
        "`payment_annualized` is inconsistent with ",
        "`k * annuity_payment`.",
        call. = FALSE
      )
    }
  }

  premium_column <- NULL

  if (!is.null(P)) {
    .mc_assert_numeric_scalar(P, "P", min = 0)
    premium_input <- rep(P, n_rows)
    premium_column <- col_P
  } else {
    if (col_P_missing) {
      candidates <- c(
        "premium_annualized",
        "premium_per_payment",
        "P",
        "premium"
      )

      premium_column <- candidates[
        candidates %in% names(.data)
      ][1L]
    } else {
      premium_column <- col_P
    }

    if (length(premium_column) == 0L ||
        is.na(premium_column) ||
        !premium_column %in% names(.data)) {
      stop(
        "No premium was found. Supply `P`, or provide one of ",
        "`premium_annualized`, `premium_per_payment`, `P`, or `premium`.",
        call. = FALSE
      )
    }

    .mc_assert_numeric_column(.data, premium_column, "col_P")
    premium_input <- .data[[premium_column]]

    if (any(
      is.na(premium_input) |
        !is.finite(premium_input) |
        premium_input < 0
    )) {
      stop(
        "The premium column must contain nonnegative finite values.",
        call. = FALSE
      )
    }
  }

  resolved_unit <- premium_unit

  if (identical(resolved_unit, "auto")) {
    if (identical(premium_column, "premium_annualized")) {
      resolved_unit <- "annualized"
    } else if (identical(premium_column, "premium_per_payment")) {
      resolved_unit <- "per_payment"
    } else {
      resolved_unit <- "annuity_scale"
    }
  }

  if (identical(resolved_unit, "annualized")) {
    premium_annualized <- premium_input
    premium_per_payment <- premium_annualized / k
    premium_annuity_scale <-
      premium_per_payment / annuity_payment
  } else if (identical(resolved_unit, "per_payment")) {
    premium_per_payment <- premium_input
    premium_annualized <- k * premium_per_payment
    premium_annuity_scale <-
      premium_per_payment / annuity_payment
  } else {
    premium_annuity_scale <- premium_input
    premium_per_payment <-
      premium_annuity_scale * annuity_payment
    premium_annualized <-
      k * premium_per_payment
  }

  pv_premiums <- premium_annuity_scale * Y
  L <- Z - pv_premiums

  reserved_loss_names <- c(
    col_Z,
    col_Y,
    premium_column,
    "premium_annualized",
    "premium_per_payment",
    "payments_per_year",
    "premium_annuity_scale",
    "pv_premiums"
  )

  if (col_L %in% reserved_loss_names) {
    stop(
      "`col_L` conflicts with an input or standardized premium column.",
      call. = FALSE
    )
  }

  out <- .data |>
    dplyr::mutate(
      "{premium_column}" := premium_input,
      premium_annualized = premium_annualized,
      premium_per_payment = premium_per_payment,
      payments_per_year = k,
      premium_annuity_scale = premium_annuity_scale,
      pv_premiums = pv_premiums,
      "{col_L}" := L
    )

  if (!identical(col_L, "loss")) {
    out <- out |>
      dplyr::mutate(
        loss = L
      )
  }

  out
}
