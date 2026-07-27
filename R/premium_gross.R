#' Gross premium under a simple expense-loaded equivalence principle
#'
#' Converts a net-premium result from \code{\link{premium_x}} or
#' \code{\link{premium_xy}} into an expense-loaded gross premium.
#'
#' The function preserves the package-wide distinction between an annualized
#' premium and the amount paid at each of the \code{k} payment dates.
#'
#' @param prem A one-row data frame or tibble. The preferred input is
#'   \code{premium_x(..., output = "summary")} or
#'   \code{premium_xy(..., output = "summary")}. It must contain an annualized
#'   net premium and the APV of the normalized premium annuity.
#' @param alpha Nonnegative numeric scalar. Initial acquisition expense as a
#'   multiple of one gross premium installment. The expense at issue is
#'   \eqn{\alpha G^{(k)} / k}.
#' @param beta Numeric scalar in \eqn{[0,1)}. Proportional expense charged on
#'   every gross premium installment.
#' @param gamma Nonnegative numeric scalar. Fixed monetary expense incurred at
#'   every premium-payment date. Its APV is
#'   \eqn{\gamma k a^{(k)}}, because \code{a_premiums} is normalized with
#'   payments of \eqn{1/k}.
#' @param k Optional positive integer payment frequency. By default it is
#'   inferred from \code{payments_per_year} or \code{k} in \code{prem}; if
#'   neither is present, \code{k = 1} is used.
#' @param output Output level. \code{"value"} returns the annualized gross
#'   premium, \code{"summary"} returns a compact one-row tibble, and
#'   \code{"audit"} returns the extended-equivalence components in long
#'   format. \code{"table"} is accepted as a compatibility alias for
#'   \code{"summary"}.
#' @param tidy Deprecated compatibility argument. \code{TRUE} maps to
#'   \code{output = "summary"} and \code{FALSE} to
#'   \code{output = "value"}.
#' @param check Logical scalar. If \code{TRUE}, validates the inputs.
#' @param tol Nonnegative numeric tolerance used for consistency checks.
#' @param ... Transitional compatibility. The deprecated argument
#'   \code{payments_per_year} is mapped to \code{k}.
#'
#' @return
#' For \code{output = "value"}, a numeric scalar containing the annualized
#' gross premium \eqn{G^{(k)}}.
#'
#' For \code{output = "summary"}, a one-row tibble with six columns:
#' annualized gross premium, gross premium per payment, annualized net premium,
#' annualized loading, payment frequency, and equivalence residual.
#'
#' For \code{output = "audit"}, a long-format tibble with the APVs of gross
#' premiums, benefits, and each expense component.
#'
#' @details
#' Let \eqn{a^{(k)}} denote the APV of the premium annuity returned by
#' \code{\link{premium_x}} or \code{\link{premium_xy}}. That annuity has
#' \eqn{k} payments per year, each of amount \eqn{1/k}; hence it represents an
#' annual payment rate of 1.
#'
#' If \eqn{P^{(k)}} is the annualized net premium and \eqn{G^{(k)}} is the
#' annualized gross premium, the extended equivalence equation implemented is
#' \deqn{
#' G^{(k)}a^{(k)}
#' =
#' P^{(k)}a^{(k)}
#' +
#' \alpha\frac{G^{(k)}}{k}
#' +
#' \beta G^{(k)}a^{(k)}
#' +
#' \gamma k a^{(k)}.
#' }
#'
#' Therefore,
#' \deqn{
#' G^{(k)}
#' =
#' \frac{P^{(k)} + k\gamma}
#' {(1-\beta)-\alpha/(k a^{(k)})}.
#' }
#'
#' The actual premium installment is
#' \deqn{
#' G_{\mathrm{per\ payment}} = \frac{G^{(k)}}{k}.
#' }
#'
#' For \code{k = 1}, this reduces to the previous annual formula:
#' \deqn{
#' G =
#' \frac{P_{\mathrm{net}}+\gamma}
#' {(1-\beta)-\alpha/a}.
#' }
#'
#' This intentionally simple model assumes that \code{gamma} is incurred at
#' the same dates and under the same contingency as premium payments. Expenses
#' with a different timing or contingency require their own APV and are not
#' represented by \code{gamma}.
#'
#' Preferred column names in \code{prem} are
#' \code{premium_annualized}, \code{premium_per_payment},
#' \code{payments_per_year}, and \code{apv_premium_annuity}. Legacy annual
#' tables using \code{P}, \code{premium}, \code{P_net},
#' \code{a_premiums}, or \code{apv_premiums} remain supported.
#'
#' Legacy premium columns are accepted automatically only when \code{k = 1},
#' because their unit is ambiguous for subannual premiums.
#'
#' @seealso \code{\link{premium_x}}, \code{\link{premium_xy}},
#'   \code{\link{annuity_x}}, \code{\link{annuity_xy}}
#'
#' @family life-contingencies
#'
#' @examples
#' net <- tibble::tibble(
#'   premium_annualized = 1200,
#'   premium_per_payment = 100,
#'   payments_per_year = 12,
#'   apv_premium_annuity = 10
#' )
#'
#' premium_gross(
#'   net,
#'   alpha = 0.5,
#'   beta = 0.05,
#'   gamma = 20,
#'   output = "summary"
#' )
#'
#' # Pipeline from a net-premium calculation
#' \dontrun{
#' life_contract(
#'   lt = lt,
#'   lives = "single",
#'   x = 40,
#'   i = 0.05
#' ) |>
#'   add_insurance(
#'     type = "term",
#'     benefit = 100000,
#'     n = 20
#'   ) |>
#'   add_premium_schedule(
#'     k = 12,
#'     n_prem = 10
#'   ) |>
#'   premium_x(output = "summary") |>
#'   premium_gross(
#'     alpha = 0.5,
#'     beta = 0.05,
#'     gamma = 20,
#'     output = "summary"
#'   )
#' }
#'
#' @export
premium_gross <- function(
    prem,
    alpha = 0,
    beta = 0,
    gamma = 0,
    k = NULL,
    output = c("value", "summary", "audit", "table"),
    tidy = NULL,
    check = TRUE,
    tol = 1e-10,
    ...
) {
  k_missing <- missing(k)
  output_missing <- missing(output)
  dots <- list(...)

  allowed_old <- "payments_per_year"
  bad_dots <- setdiff(names(dots), allowed_old)

  if (length(bad_dots) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad_dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if ("payments_per_year" %in% names(dots)) {
    if (!k_missing && !is.null(k)) {
      stop(
        "Provide only one of `k` or deprecated `payments_per_year`.",
        call. = FALSE
      )
    }

    k <- dots[["payments_per_year"]]
    k_missing <- FALSE
  }

  if (!is.null(tidy)) {
    if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
      stop("`tidy` must be `NULL` or a logical scalar.", call. = FALSE)
    }

    if (!output_missing) {
      stop(
        "Provide only one of `output` or deprecated `tidy`.",
        call. = FALSE
      )
    }

    output <- if (isTRUE(tidy)) "summary" else "value"
  }

  output <- match.arg(output)

  if (identical(output, "table")) {
    output <- "summary"
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
  }

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  if (!inherits(prem, "data.frame") || nrow(prem) != 1L) {
    stop(
      "`prem` must be a one-row premium summary from `premium_x()`, ",
      "`premium_xy()`, or a compatible table.",
      call. = FALSE
    )
  }

  k_column <- if ("payments_per_year" %in% names(prem)) {
    "payments_per_year"
  } else if ("k" %in% names(prem)) {
    "k"
  } else {
    NA_character_
  }

  k_from_prem <- if (!is.na(k_column)) {
    prem[[k_column]][[1L]]
  } else {
    NULL
  }

  if (is.null(k)) {
    k <- if (is.null(k_from_prem)) 1L else k_from_prem
  } else if (!is.null(k_from_prem) &&
             is.numeric(k) &&
             is.numeric(k_from_prem) &&
             length(k) == 1L &&
             length(k_from_prem) == 1L &&
             is.finite(k) &&
             is.finite(k_from_prem) &&
             abs(k - k_from_prem) > tol) {
    stop(
      "`k` does not agree with the payment frequency stored in `prem`.",
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

  apv_column <- if ("apv_premium_annuity" %in% names(prem)) {
    "apv_premium_annuity"
  } else if ("a_premiums" %in% names(prem)) {
    "a_premiums"
  } else if ("apv_premiums" %in% names(prem)) {
    "apv_premiums"
  } else {
    NA_character_
  }

  if (is.na(apv_column)) {
    stop(
      "The APV of the premium annuity was not found. Expected ",
      "`apv_premium_annuity`, `a_premiums`, or `apv_premiums`.",
      call. = FALSE
    )
  }

  a_premiums <- prem[[apv_column]][[1L]]

  annualized_column <- if ("premium_annualized" %in% names(prem)) {
    "premium_annualized"
  } else if ("P_net_annualized" %in% names(prem)) {
    "P_net_annualized"
  } else {
    NA_character_
  }

  per_payment_column <- if ("premium_per_payment" %in% names(prem)) {
    "premium_per_payment"
  } else if ("P_net_per_payment" %in% names(prem)) {
    "P_net_per_payment"
  } else {
    NA_character_
  }

  legacy_column <- if ("P_net" %in% names(prem)) {
    "P_net"
  } else if ("P" %in% names(prem)) {
    "P"
  } else if ("premium" %in% names(prem)) {
    "premium"
  } else {
    NA_character_
  }

  if (!is.na(annualized_column)) {
    P_net_annualized <- prem[[annualized_column]][[1L]]
  } else if (!is.na(per_payment_column)) {
    P_net_annualized <- k * prem[[per_payment_column]][[1L]]
  } else if (!is.na(legacy_column)) {
    if (k != 1L) {
      stop(
        "Legacy premium columns are ambiguous when `k > 1`. Supply ",
        "`premium_annualized` or `premium_per_payment` explicitly.",
        call. = FALSE
      )
    }

    P_net_annualized <- prem[[legacy_column]][[1L]]
  } else {
    stop(
      "The net premium was not found. Expected `premium_annualized`, ",
      "`premium_per_payment`, or a supported legacy annual column.",
      call. = FALSE
    )
  }

  if (!is.na(annualized_column) && !is.na(per_payment_column)) {
    P_from_payment <- k * prem[[per_payment_column]][[1L]]

    if (is.numeric(P_net_annualized) &&
        is.numeric(P_from_payment) &&
        length(P_net_annualized) == 1L &&
        length(P_from_payment) == 1L &&
        is.finite(P_net_annualized) &&
        is.finite(P_from_payment) &&
        abs(P_net_annualized - P_from_payment) >
          tol * max(1, abs(P_net_annualized))) {
      stop(
        "`premium_annualized` is inconsistent with ",
        "`premium_per_payment * payments_per_year`.",
        call. = FALSE
      )
    }
  }

  if (isTRUE(check)) {
    if (!is.numeric(P_net_annualized) ||
        length(P_net_annualized) != 1L ||
        is.na(P_net_annualized) ||
        !is.finite(P_net_annualized) ||
        P_net_annualized < 0) {
      stop(
        "The annualized net premium must be a single nonnegative ",
        "finite number.",
        call. = FALSE
      )
    }

    if (!is.numeric(a_premiums) ||
        length(a_premiums) != 1L ||
        is.na(a_premiums) ||
        !is.finite(a_premiums) ||
        a_premiums <= 0) {
      stop(
        "The APV of the premium annuity must be a single positive ",
        "finite number.",
        call. = FALSE
      )
    }

    if (!is.numeric(alpha) ||
        length(alpha) != 1L ||
        is.na(alpha) ||
        !is.finite(alpha) ||
        alpha < 0) {
      stop("`alpha` must be a single nonnegative finite number.", call. = FALSE)
    }

    if (!is.numeric(beta) ||
        length(beta) != 1L ||
        is.na(beta) ||
        !is.finite(beta) ||
        beta < 0 ||
        beta >= 1) {
      stop("`beta` must be a single finite number in [0, 1).", call. = FALSE)
    }

    if (!is.numeric(gamma) ||
        length(gamma) != 1L ||
        is.na(gamma) ||
        !is.finite(gamma) ||
        gamma < 0) {
      stop("`gamma` must be a single nonnegative finite number.", call. = FALSE)
    }
  }

  denominator <- (1 - beta) - alpha / (k * a_premiums)

  if (!is.finite(denominator) || denominator <= 0) {
    stop(
      "No level gross premium exists: the expense structure makes the ",
      "extended-equivalence denominator nonpositive.",
      call. = FALSE
    )
  }

  gross_premium_annualized <-
    (P_net_annualized + k * gamma) / denominator

  gross_premium_per_payment <- gross_premium_annualized / k
  net_premium_per_payment <- P_net_annualized / k

  apv_benefits <- P_net_annualized * a_premiums
  apv_gross_premiums <- gross_premium_annualized * a_premiums

  apv_initial_expense <-
    alpha * gross_premium_per_payment

  apv_proportional_expense <-
    beta * apv_gross_premiums

  apv_fixed_per_payment_expense <-
    gamma * k * a_premiums

  apv_total_expenses <-
    apv_initial_expense +
    apv_proportional_expense +
    apv_fixed_per_payment_expense

  equivalence_residual <-
    apv_gross_premiums -
    apv_benefits -
    apv_total_expenses

  loading_annualized <-
    gross_premium_annualized - P_net_annualized

  loading_per_payment <-
    gross_premium_per_payment - net_premium_per_payment

  if (identical(output, "value")) {
    return(gross_premium_annualized)
  }

  if (identical(output, "summary")) {
    return(
      tibble::tibble(
        gross_premium_annualized = gross_premium_annualized,
        gross_premium_per_payment = gross_premium_per_payment,
        net_premium_annualized = P_net_annualized,
        loading_annualized = loading_annualized,
        payments_per_year = k,
        equivalence_residual = equivalence_residual
      )
    )
  }

  tibble::tibble(
    component = c(
      "apv_gross_premiums",
      "apv_benefits",
      "apv_initial_expense",
      "apv_proportional_expense",
      "apv_fixed_per_payment_expense",
      "apv_total_expenses",
      "gross_premium_annualized",
      "gross_premium_per_payment",
      "net_premium_annualized",
      "net_premium_per_payment",
      "loading_annualized",
      "loading_per_payment",
      "payments_per_year",
      "equivalence_residual"
    ),
    value = c(
      apv_gross_premiums,
      apv_benefits,
      apv_initial_expense,
      apv_proportional_expense,
      apv_fixed_per_payment_expense,
      apv_total_expenses,
      gross_premium_annualized,
      gross_premium_per_payment,
      P_net_annualized,
      net_premium_per_payment,
      loading_annualized,
      loading_per_payment,
      k,
      equivalence_residual
    ),
    unit = c(
      "currency",
      "currency",
      "currency",
      "currency",
      "currency",
      "currency",
      "currency per year",
      "currency per payment",
      "currency per year",
      "currency per payment",
      "currency per year",
      "currency per payment",
      "payments per year",
      "currency"
    )
  )
}
