#' Net premium for two-life insurance by the equivalence principle
#'
#' Computes the net benefit premium for a joint-life or last-survivor
#' insurance contract.
#'
#' For premiums payable \eqn{k} times per year, the function distinguishes
#' between the annualized premium \eqn{P^{(k)}} and the amount paid at each
#' installment, \eqn{P^{(k)} / k}.
#'
#' @param lt A life table, a list of two life tables
#'   \code{list(lt_x, lt_y)}, or a two-life contract created with
#'   \code{\link{life_contract}}.
#' @param x Integer actuarial age for the first life. Optional for a
#'   \code{life_contract}.
#' @param y Integer actuarial age for the second life. Optional for a
#'   \code{life_contract}.
#' @param i Numeric scalar. Annual interest-rate input. Optional for a
#'   \code{life_contract}.
#' @param i_type Interest-rate type: \code{"effective"},
#'   \code{"nominal_interest"}, \code{"nominal_discount"}, or \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates.
#' @param type Insurance type: \code{"whole"}, \code{"term"},
#'   \code{"endowment"}, or \code{"pure_endowment"}.
#' @param benefit Nonnegative insurance benefit.
#' @param n Insurance term in years after deferment. Use \code{Inf} for
#'   whole-life insurance.
#' @param h Nonnegative integer deferment period in years.
#' @param k Positive integer. Number of premium payments per year.
#' @param frac Fractional-age assumption: \code{"UDD"}, \code{"CF"},
#'   \code{"CML"}, or \code{"Balducci"}.
#' @param timing Premium timing: \code{"due"} or \code{"immediate"}.
#' @param premium_start Start of premiums: \code{"issue"} or
#'   \code{"deferred"}.
#' @param n_prem Premium-paying term in years. A fractional value is permitted
#'   when \code{n_prem * k} is an integer.
#' @param status Two-life status: \code{"joint"} or \code{"last"}. For a
#'   \code{life_contract}, the value is inferred from \code{lives} unless
#'   supplied explicitly.
#' @param woolhouse Woolhouse approximation order for the premium annuity:
#'   \code{"none"}, \code{"first"}, or \code{"second"}.
#' @param output Output level. \code{"value"} returns the annualized premium;
#'   \code{"summary"} returns a compact one-row result; \code{"audit"} returns
#'   the equivalence components in long format. \code{"table"} is accepted as
#'   a deprecated alias for \code{"summary"}.
#' @param tidy Deprecated compatibility argument. \code{TRUE} maps to
#'   \code{output = "summary"} and \code{FALSE} to
#'   \code{output = "value"}.
#' @param check Logical scalar. If \code{TRUE}, validates inputs.
#' @param tol Numeric tolerance for integer-grid checks.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age_x}, \code{age_y}, \code{rate},
#'   \code{rate_type}, \code{insurance_type}, \code{term_years},
#'   \code{deferment_years}, \code{payments_per_year},
#'   \code{premium_timing}, \code{premium_term_years}, and \code{cohort}.
#'
#' @return
#' For \code{output = "value"}, a numeric scalar containing the annualized
#' premium \eqn{P^{(k)}}.
#'
#' For \code{output = "summary"}, a one-row tibble with six columns:
#' annualized premium, premium per payment, payment frequency, APV of
#' benefits, APV of the premium annuity, and equivalence residual.
#'
#' For \code{output = "audit"}, a long-format tibble.
#'
#' @details
#' Let \eqn{Z} denote the present-value random variable of the two-life
#' insurance benefit and let \eqn{Y^{(k)}} denote the present value of the
#' contingent premium annuity normalized to an annual payment rate of 1.
#' The loss at issue is
#' \deqn{
#' L_0 = Z - P^{(k)}Y^{(k)}.
#' }
#'
#' The equivalence principle gives
#' \deqn{
#' P^{(k)}
#' =
#' \frac{\operatorname{APV}(\text{benefits})}
#'      {\operatorname{APV}(\text{premium annuity})}.
#' }
#'
#' Since \code{\link{annuity_xy}} assigns amount \eqn{1/k} to each k-thly
#' payment, this quotient is the annualized premium. The installment is
#' \deqn{
#' P_{\text{per payment}} = \frac{P^{(k)}}{k}.
#' }
#'
#' Standard whole-life, term, and endowment benefits are valued through
#' \code{\link{insurance_xy}}. A pure endowment is obtained as the difference
#' between the corresponding endowment and term insurance values.
#'
#' A contract assembled with pipes may be valued directly:
#'
#' \preformatted{
#' life_contract(...) |>
#'   add_insurance(...) |>
#'   add_premium_schedule(...) |>
#'   premium_xy(output = "summary")
#' }
#'
#' Explicit arguments supplied to \code{premium_xy()} override values stored
#' in the contract components.
#'
#' @seealso \code{\link{premium_x}}, \code{\link{insurance_xy}},
#'   \code{\link{annuity_xy}}, \code{\link{reserve_xy}},
#'   \code{\link{life_contract}}, \code{\link{add_insurance}},
#'   \code{\link{add_premium_schedule}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 40:100,
#'   lx = round(100000 * exp(-0.012 * (0:60)^1.35))
#' )
#' lt$lx[nrow(lt)] <- 0
#'
#' premium_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   type = "term",
#'   status = "joint",
#'   benefit = 100000,
#'   n = 20,
#'   k = 12,
#'   n_prem = 10,
#'   output = "summary"
#' )
#'
#' life_contract(
#'   lt = lt,
#'   lives = "joint",
#'   x = 60,
#'   y = 62,
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
#'   premium_xy(output = "summary")
#'
#' @export
premium_xy <- function(
    lt,
    x = NULL,
    y = NULL,
    i = NULL,
    i_type = "effective",
    m = 1L,
    type = c("whole", "term", "endowment", "pure_endowment"),
    benefit = 1,
    n = Inf,
    h = 0L,
    k = 1L,
    frac = c("UDD", "CF", "CML", "Balducci"),
    timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    n_prem = NULL,
    status = c("joint", "last"),
    woolhouse = c("none", "first", "second"),
    output = c("value", "summary", "audit", "table"),
    tidy = NULL,
    check = TRUE,
    tol = 1e-10,
    ...
) {
  lt_missing <- missing(lt)
  x_missing <- missing(x)
  y_missing <- missing(y)
  i_missing <- missing(i)
  i_type_missing <- missing(i_type)
  m_missing <- missing(m)
  type_missing <- missing(type)
  benefit_missing <- missing(benefit)
  n_missing <- missing(n)
  h_missing <- missing(h)
  k_missing <- missing(k)
  frac_missing <- missing(frac)
  timing_missing <- missing(timing)
  premium_start_missing <- missing(premium_start)
  n_prem_missing <- missing(n_prem)
  status_missing <- missing(status)
  woolhouse_missing <- missing(woolhouse)
  output_missing <- missing(output)

  dots <- list(...)

  allowed_old <- c(
    "mortality_table",
    "age_x",
    "age_y",
    "rate",
    "rate_type",
    "insurance_type",
    "term_years",
    "deferment_years",
    "payments_per_year",
    "premium_timing",
    "premium_term_years",
    "cohort"
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

  if ("mortality_table" %in% names(dots)) {
    if (!lt_missing) {
      stop(
        "Provide only one of `lt` or deprecated `mortality_table`.",
        call. = FALSE
      )
    }

    lt <- dots[["mortality_table"]]
    lt_missing <- FALSE
  }

  if ("age_x" %in% names(dots)) {
    if (!x_missing && !is.null(x)) {
      stop("Provide only one of `x` or deprecated `age_x`.", call. = FALSE)
    }

    x <- dots[["age_x"]]
    x_missing <- FALSE
  }

  if ("age_y" %in% names(dots)) {
    if (!y_missing && !is.null(y)) {
      stop("Provide only one of `y` or deprecated `age_y`.", call. = FALSE)
    }

    y <- dots[["age_y"]]
    y_missing <- FALSE
  }

  if ("rate" %in% names(dots)) {
    if (!i_missing && !is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dots[["rate"]]
    i_missing <- FALSE
  }

  if ("rate_type" %in% names(dots)) {
    if (!i_type_missing) {
      stop(
        "Provide only one of `i_type` or deprecated `rate_type`.",
        call. = FALSE
      )
    }

    i_type <- dots[["rate_type"]]
    i_type_missing <- FALSE
  }

  if ("insurance_type" %in% names(dots)) {
    if (!type_missing) {
      stop(
        "Provide only one of `type` or deprecated `insurance_type`.",
        call. = FALSE
      )
    }

    type <- dots[["insurance_type"]]
    type_missing <- FALSE
  }

  if ("term_years" %in% names(dots)) {
    if (!n_missing) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }

    n <- dots[["term_years"]]
    n_missing <- FALSE
  }

  if ("deferment_years" %in% names(dots)) {
    if (!h_missing) {
      stop(
        "Provide only one of `h` or deprecated `deferment_years`.",
        call. = FALSE
      )
    }

    h <- dots[["deferment_years"]]
    h_missing <- FALSE
  }

  if ("payments_per_year" %in% names(dots)) {
    if (!k_missing) {
      stop(
        "Provide only one of `k` or deprecated `payments_per_year`.",
        call. = FALSE
      )
    }

    k <- dots[["payments_per_year"]]
    k_missing <- FALSE
  }

  if ("premium_timing" %in% names(dots)) {
    if (!timing_missing) {
      stop(
        "Provide only one of `timing` or deprecated `premium_timing`.",
        call. = FALSE
      )
    }

    timing <- dots[["premium_timing"]]
    timing_missing <- FALSE
  }

  if ("premium_term_years" %in% names(dots)) {
    if (!n_prem_missing) {
      stop(
        "Provide only one of `n_prem` or deprecated `premium_term_years`.",
        call. = FALSE
      )
    }

    n_prem <- dots[["premium_term_years"]]
    n_prem_missing <- FALSE
  }

  if ("cohort" %in% names(dots)) {
    if (!status_missing) {
      stop(
        "Provide only one of `status` or deprecated `cohort`.",
        call. = FALSE
      )
    }

    old_cohort <- match.arg(dots[["cohort"]], c("first", "last"))
    status <- if (identical(old_cohort, "first")) "joint" else "last"
    status_missing <- FALSE
  }

  if (!is.null(tidy)) {
    if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
      stop("`tidy` must be `NULL` or a logical scalar.", call. = FALSE)
    }

    if (!output_missing) {
      stop("Provide only one of `output` or deprecated `tidy`.", call. = FALSE)
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

  insurance_spec <- NULL
  premium_spec <- NULL

  if (!lt_missing &&
      exists(".as_life_contract", mode = "function") &&
      .as_life_contract(lt)) {
    contract <- lt

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop(
        "`premium_xy()` requires a two-life `life_contract()` object.",
        call. = FALSE
      )
    }

    lt <- contract$lt
    insurance_spec <- contract$insurance
    premium_spec <- contract$premium_schedule

    if (x_missing || is.null(x)) {
      x <- contract$x
      x_missing <- FALSE
    }

    if (y_missing || is.null(y)) {
      y <- contract$y
      y_missing <- FALSE
    }

    if (i_missing || is.null(i)) {
      i <- contract$i
      i_missing <- FALSE
    }

    if (i_type_missing || is.null(i_type)) {
      i_type <- contract$i_type
      i_type_missing <- FALSE
    }

    if (m_missing || is.null(m)) {
      m <- contract$m
      m_missing <- FALSE
    }

    if (status_missing) {
      status <- if (identical(contract$lives, "joint")) {
        "joint"
      } else {
        "last"
      }
      status_missing <- FALSE
    }

    if (!is.null(insurance_spec)) {
      if (identical(insurance_spec$type, "variable_k")) {
        stop(
          "`premium_xy()` does not support `type = 'variable_k'`.",
          call. = FALSE
        )
      }

      if (type_missing) {
        type <- insurance_spec$type
        type_missing <- FALSE
      }

      if (benefit_missing) {
        benefit <- insurance_spec$benefit
        benefit_missing <- FALSE
      }

      if (n_missing) {
        n <- insurance_spec$n
        n_missing <- FALSE
      }

      if (h_missing) {
        h <- insurance_spec$h
        h_missing <- FALSE
      }

      if (status_missing && !is.null(insurance_spec$status)) {
        status <- insurance_spec$status
        status_missing <- FALSE
      }
    }

    if (!is.null(premium_spec)) {
      if (k_missing) {
        k <- premium_spec$k
        k_missing <- FALSE
      }

      if (timing_missing) {
        timing <- premium_spec$timing
        timing_missing <- FALSE
      }

      if (premium_start_missing) {
        premium_start <- premium_spec$premium_start
        premium_start_missing <- FALSE
      }

      if (n_prem_missing) {
        n_prem <- premium_spec$n_prem
        n_prem_missing <- FALSE
      }

      if (woolhouse_missing) {
        woolhouse <- premium_spec$woolhouse
        woolhouse_missing <- FALSE
      }
    }
  }

  type <- match.arg(type)
  timing <- match.arg(timing)
  premium_start <- match.arg(premium_start)
  status <- match.arg(status)
  woolhouse <- match.arg(woolhouse)

  if (frac_missing) {
    frac_benefit <- if (!is.null(insurance_spec$frac)) {
      insurance_spec$frac
    } else if (!is.null(premium_spec$frac)) {
      premium_spec$frac
    } else {
      "UDD"
    }

    frac_premium <- if (!is.null(premium_spec$frac)) {
      premium_spec$frac
    } else {
      frac_benefit
    }
  } else {
    frac_common <- match.arg(frac)
    frac_benefit <- frac_common
    frac_premium <- frac_common
  }

  valid_frac <- c("UDD", "CF", "CML", "Balducci")

  frac_benefit <- match.arg(frac_benefit, valid_frac)
  frac_premium <- match.arg(frac_premium, valid_frac)

  if (identical(frac_benefit, "CML")) {
    frac_benefit <- "CF"
  }

  if (identical(frac_premium, "CML")) {
    frac_premium <- "CF"
  }

  if (isTRUE(check)) {
    if (lt_missing) {
      stop("`lt` must be provided.", call. = FALSE)
    }

    if (x_missing || is.null(x)) {
      stop("`x` must be provided.", call. = FALSE)
    }

    if (y_missing || is.null(y)) {
      stop("`y` must be provided.", call. = FALSE)
    }

    if (i_missing || is.null(i)) {
      stop("`i` must be provided.", call. = FALSE)
    }

    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        abs(x - round(x)) > tol) {
      stop("`x` must be a single integer age.", call. = FALSE)
    }

    if (!is.numeric(y) ||
        length(y) != 1L ||
        is.na(y) ||
        !is.finite(y) ||
        abs(y - round(y)) > tol) {
      stop("`y` must be a single integer age.", call. = FALSE)
    }

    if (!is.numeric(i) ||
        length(i) != 1L ||
        is.na(i) ||
        !is.finite(i)) {
      stop("`i` must be a single finite numeric value.", call. = FALSE)
    }

    valid_i_type <- c(
      "effective",
      "nominal_interest",
      "nominal_discount",
      "force"
    )

    if (!is.character(i_type) ||
        length(i_type) != 1L ||
        is.na(i_type) ||
        !i_type %in% valid_i_type) {
      stop(
        "`i_type` must be one of: ",
        paste(sprintf("'%s'", valid_i_type), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.numeric(m) ||
        length(m) != 1L ||
        is.na(m) ||
        !is.finite(m) ||
        m < 1 ||
        abs(m - round(m)) > tol) {
      stop("`m` must be a single positive integer.", call. = FALSE)
    }

    if (!is.numeric(benefit) ||
        length(benefit) != 1L ||
        is.na(benefit) ||
        !is.finite(benefit) ||
        benefit < 0) {
      stop(
        "`benefit` must be a single finite nonnegative number.",
        call. = FALSE
      )
    }

    if (!is.numeric(n) ||
        length(n) != 1L ||
        is.na(n) ||
        n <= 0 ||
        (!is.infinite(n) &&
         (!is.finite(n) || abs(n - round(n)) > tol))) {
      stop("`n` must be `Inf` or a single positive integer.", call. = FALSE)
    }

    if (type %in% c("term", "endowment", "pure_endowment") &&
        is.infinite(n)) {
      stop(
        "`n` must be finite for term, endowment, and pure endowment insurance.",
        call. = FALSE
      )
    }

    if (!is.numeric(h) ||
        length(h) != 1L ||
        is.na(h) ||
        !is.finite(h) ||
        h < 0 ||
        abs(h - round(h)) > tol) {
      stop("`h` must be a single nonnegative integer.", call. = FALSE)
    }

    if (!is.numeric(k) ||
        length(k) != 1L ||
        is.na(k) ||
        !is.finite(k) ||
        k < 1 ||
        abs(k - round(k)) > tol) {
      stop("`k` must be a single positive integer.", call. = FALSE)
    }

    if (!is.null(n_prem) &&
        (!is.numeric(n_prem) ||
         length(n_prem) != 1L ||
         is.na(n_prem) ||
         n_prem <= 0 ||
         (!is.infinite(n_prem) &&
          (!is.finite(n_prem) ||
           abs(n_prem * k - round(n_prem * k)) > tol)))) {
      stop(
        "`n_prem` must be NULL, Inf, or a positive value satisfying ",
        "`n_prem * k` integer.",
        call. = FALSE
      )
    }
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (!is.null(n_prem) && !is.infinite(n_prem)) {
    n_prem <- round(n_prem * k) / k
  }

  if (is.null(n_prem)) {
    n_prem <- if (identical(type, "whole")) Inf else n
  }

  premium_start_time <- if (identical(premium_start, "issue")) 0 else h

  if (!identical(type, "whole")) {
    premium_end_time <- if (is.infinite(n_prem)) {
      Inf
    } else {
      premium_start_time + n_prem
    }

    coverage_end_time <- h + n

    if (premium_end_time > coverage_end_time + tol) {
      stop(
        "Premium payments must not extend beyond the end of coverage. ",
        "Premiums end at time ", premium_end_time,
        ", while coverage ends at time ", coverage_end_time, ".",
        call. = FALSE
      )
    }
  }

  benefit_args <- list(
    lt = lt,
    x = x,
    y = y,
    i = i,
    i_type = i_type,
    m = m,
    status = status,
    n = n,
    h = h,
    benefit = benefit,
    frac = frac_benefit,
    tidy = FALSE
  )

  if (identical(type, "pure_endowment")) {
    endowment_args <- benefit_args
    endowment_args$type <- "endowment"

    term_args <- benefit_args
    term_args$type <- "term"

    apv_benefits <- do.call(insurance_xy, endowment_args) -
      do.call(insurance_xy, term_args)

    if (apv_benefits < 0 && abs(apv_benefits) <= tol) {
      apv_benefits <- 0
    }
  } else {
    benefit_args$type <- type
    apv_benefits <- do.call(insurance_xy, benefit_args)
  }

  h_prem <- if (identical(premium_start, "issue")) 0 else h

  apv_premium_annuity <- annuity_xy(
    lt = lt,
    x = x,
    y = y,
    i = i,
    i_type = i_type,
    m = m,
    status = status,
    n = n_prem,
    h = h_prem,
    k = k,
    timing = timing,
    woolhouse = woolhouse,
    frac = frac_premium,
    tidy = FALSE
  )

  if (!is.finite(apv_premium_annuity) ||
      apv_premium_annuity <= 0) {
    stop(
      "APV of the premium annuity is nonpositive or not finite.",
      call. = FALSE
    )
  }

  premium_annualized <- apv_benefits / apv_premium_annuity
  premium_per_payment <- premium_annualized / k

  equivalence_residual <- premium_annualized *
    apv_premium_annuity -
    apv_benefits

  if (identical(output, "value")) {
    return(premium_annualized)
  }

  if (identical(output, "summary")) {
    return(
      tibble::tibble(
        premium_annualized = premium_annualized,
        premium_per_payment = premium_per_payment,
        payments_per_year = k,
        apv_benefits = apv_benefits,
        apv_premium_annuity = apv_premium_annuity,
        equivalence_residual = equivalence_residual
      )
    )
  }

  tibble::tibble(
    component = c(
      "apv_benefits",
      "apv_premium_annuity",
      "premium_annualized",
      "premium_per_payment",
      "payments_per_year",
      "equivalence_residual"
    ),
    value = c(
      apv_benefits,
      apv_premium_annuity,
      premium_annualized,
      premium_per_payment,
      k,
      equivalence_residual
    ),
    unit = c(
      "currency",
      "annuity factor",
      "currency per year",
      "currency per payment",
      "payments per year",
      "currency"
    )
  )
}
