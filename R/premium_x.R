#' Net premium for life insurance by the equivalence principle
#'
#' Computes the net benefit premium of a life insurance contract using the
#' equivalence principle:
#' \deqn{P = \frac{\text{APV of benefits}}{\text{APV of premium annuity}}.}
#'
#' The premium returned corresponds to one premium payment. For example, when
#' \code{payments_per_year = 12}, the returned value is the monthly premium.
#'
#' @param mortality_table A life table data frame containing at least columns
#'   \code{x} and \code{lx}.
#' @param age Integer actuarial age at issue.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Compounding frequency for nominal rates. Ignored
#'   for \code{rate_type = "effective"} and \code{rate_type = "force"}.
#' @param insurance_type Type of insurance contract. One of \code{"whole"},
#'   \code{"term"}, \code{"endowment"}, or \code{"variable_k"}.
#' @param benefit Benefit amount. For standard products, a single nonnegative
#'   numeric value. For \code{insurance_type = "variable_k"}, a numeric vector
#'   or a function of time may be supplied and is passed to
#'   \code{\link{insurance_variable_k}}.
#' @param term_years Term of the insurance contract in years. Use \code{Inf}
#'   for whole-life insurance. Required as a finite integer for term and
#'   endowment insurance.
#' @param deferral_years Integer deferral period in years.
#' @param payments_per_year Positive integer. Number of premium payments per
#'   year.
#' @param frac Fractional-age assumption used only for
#'   \code{insurance_type = "variable_k"}. One of \code{"UDD"}, \code{"CF"},
#'   \code{"CML"}, or \code{"Balducci"}.
#' @param premium_timing Timing of premium payments. Use \code{"due"} for
#'   payments in advance or \code{"immediate"} for payments in arrears.
#' @param premium_start Start of premium payments. Use \code{"issue"} for
#'   premiums starting at issue, or \code{"deferred"} for premiums starting
#'   after \code{deferral_years}.
#' @param premium_term_years Optional premium-paying term in years, counted from
#'   \code{premium_start}. If \code{NULL}, it defaults to whole life for whole
#'   life insurance and to \code{term_years} for temporary products. For finite
#'   term contracts, premiums must not extend beyond the end of coverage.
#' @param woolhouse Woolhouse order for the premium annuity when
#'   \code{payments_per_year > 1}. One of \code{"none"}, \code{"first"}, or
#'   \code{"second"}.
#' @param output Character string. Use \code{"value"} to return a numeric
#'   premium, or \code{"table"} to return a one-row tibble with details.
#' @param check Logical. If \code{TRUE}, performs input validation.
#'
#' @return
#' If \code{output = "value"}, a numeric scalar with the net premium per
#' payment.
#'
#' If \code{output = "table"}, a one-row tibble with the main inputs, APV of
#' benefits, APV of premiums, premium per payment, and annualized premium.
#'
#' @details
#' The benefit premium is the level payment satisfying the equivalence
#' principle: the APV of premiums equals the APV of benefits at issue.
#'
#' For standard products, the APV of benefits is computed with
#' \code{\link{insurance_x}}. The APV of the premium annuity is computed with
#' \code{\link{annuity_x}}, supporting k-thly payments and Woolhouse
#' approximations.
#'
#' @seealso \code{\link{insurance_x}}, \code{\link{annuity_x}},
#'   \code{\link{premium_xy}}, \code{\link{premium_gross}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Whole life insurance, annual premium
#' premium_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.05,
#'   insurance_type = "whole",
#'   benefit = 100000
#' )
#'
#' # Verify manually: P = A / a-double-dot
#' A <- insurance_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.05,
#'   insurance_type = "whole",
#'   benefit = 100000
#' )
#'
#' ad <- annuity_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.05,
#'   timing = "due"
#' )
#'
#' A / ad
#'
#' # Five-year term insurance
#' premium_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.05,
#'   insurance_type = "term",
#'   term_years = 5,
#'   benefit = 100000
#' )
#'
#' # Table output
#' premium_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.05,
#'   insurance_type = "term",
#'   term_years = 5,
#'   benefit = 100000,
#'   output = "table"
#' )
#'
#' # Monthly premiums paid for a shorter period than the coverage term
#' premium_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.05,
#'   insurance_type = "term",
#'   term_years = 5,
#'   benefit = 100000,
#'   payments_per_year = 12,
#'   premium_term_years = 3
#' )
#'
#' @export
premium_x <- function(
    mortality_table,
    age,
    rate,
    rate_type = "effective",
    m = 1L,
    insurance_type = c("whole", "term", "endowment", "variable_k"),
    benefit = 1,
    term_years = Inf,
    deferral_years = 0L,
    payments_per_year = 1L,
    frac = c("UDD", "CF", "CML", "Balducci"),
    premium_timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    premium_term_years = NULL,
    woolhouse = c("none", "first", "second"),
    output = c("value", "table"),
    check = TRUE
) {
  insurance_type <- match.arg(insurance_type)
  frac <- match.arg(frac)
  premium_timing <- match.arg(premium_timing)
  premium_start <- match.arg(premium_start)
  woolhouse <- match.arg(woolhouse)
  output <- match.arg(output)

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (.as_life_contract(mortality_table)) {
    contract <- mortality_table

    if (!identical(contract$lives, "single")) {
      stop(
        "`premium_x()` currently supports only single-life `life_contract()` objects.",
        call. = FALSE
      )
    }

    mortality_table <- contract$mortality_table

    if (missing(age) || is.null(age)) {
      age <- contract$age
    }

    if (missing(rate) || is.null(rate)) {
      rate <- contract$rate
    }

    if (missing(rate_type) || is.null(rate_type)) {
      rate_type <- contract$rate_type
    }

    if (missing(m) || is.null(m)) {
      m <- contract$m
    }
  }

  if (isTRUE(check)) {
    if (!is.data.frame(mortality_table)) {
      stop("`mortality_table` must be a data.frame or tibble.", call. = FALSE)
    }

    if (!all(c("x", "lx") %in% names(mortality_table))) {
      stop("`mortality_table` must contain columns `x` and `lx`.", call. = FALSE)
    }

    if (!is.numeric(age) ||
        length(age) != 1L ||
        is.na(age) ||
        !is.finite(age) ||
        abs(age - round(age)) > 1e-10) {
      stop("`age` must be a single integer age.", call. = FALSE)
    }

    if (!is.numeric(rate) ||
        length(rate) != 1L ||
        is.na(rate) ||
        !is.finite(rate)) {
      stop("`rate` must be a single finite numeric value.", call. = FALSE)
    }

    if (!is.character(rate_type) ||
        length(rate_type) != 1L ||
        is.na(rate_type)) {
      stop("`rate_type` must be a single character string.", call. = FALSE)
    }

    if (!is.numeric(m) ||
        length(m) != 1L ||
        is.na(m) ||
        !is.finite(m) ||
        m < 1 ||
        abs(m - round(m)) > 1e-10) {
      stop("`m` must be a single positive integer.", call. = FALSE)
    }

    if (!is.numeric(deferral_years) ||
        length(deferral_years) != 1L ||
        is.na(deferral_years) ||
        !is.finite(deferral_years) ||
        deferral_years < 0 ||
        abs(deferral_years - round(deferral_years)) > 1e-10) {
      stop("`deferral_years` must be a single nonnegative integer.", call. = FALSE)
    }

    if (!is.numeric(payments_per_year) ||
        length(payments_per_year) != 1L ||
        is.na(payments_per_year) ||
        !is.finite(payments_per_year) ||
        payments_per_year < 1 ||
        abs(payments_per_year - round(payments_per_year)) > 1e-10) {
      stop("`payments_per_year` must be a single positive integer.", call. = FALSE)
    }

    if (!is.numeric(term_years) ||
        length(term_years) != 1L ||
        is.na(term_years) ||
        term_years <= 0 ||
        (!is.infinite(term_years) &&
         (!is.finite(term_years) ||
          abs(term_years - round(term_years)) > 1e-10))) {
      stop("`term_years` must be `Inf` or a single positive integer.", call. = FALSE)
    }

    if (insurance_type %in% c("term", "endowment") && is.infinite(term_years)) {
      stop(
        "`term_years` must be finite for term and endowment insurance.",
        call. = FALSE
      )
    }

    if (!is.null(premium_term_years) &&
        (!is.numeric(premium_term_years) ||
         length(premium_term_years) != 1L ||
         is.na(premium_term_years) ||
         premium_term_years <= 0 ||
         (!is.infinite(premium_term_years) &&
          (!is.finite(premium_term_years) ||
           abs(premium_term_years - round(premium_term_years)) > 1e-10)))) {
      stop(
        "`premium_term_years` must be NULL, Inf, or a single positive integer.",
        call. = FALSE
      )
    }

    if (insurance_type != "variable_k") {
      if (!is.numeric(benefit) ||
          length(benefit) != 1L ||
          is.na(benefit) ||
          !is.finite(benefit) ||
          benefit < 0) {
        stop(
          "For standard products, `benefit` must be a single finite nonnegative number.",
          call. = FALSE
        )
      }
    }
  }

  age <- as.integer(round(age))
  m <- as.integer(round(m))
  deferral_years <- as.integer(round(deferral_years))
  payments_per_year <- as.integer(round(payments_per_year))

  if (!is.infinite(term_years)) {
    term_years <- as.integer(round(term_years))
  }

  if (!is.null(premium_term_years) && !is.infinite(premium_term_years)) {
    premium_term_years <- as.integer(round(premium_term_years))
  }

  # -------------------------------------------------------------------------
  # 1. APV of benefits, valued at issue
  # -------------------------------------------------------------------------

  if (insurance_type == "variable_k") {
    i_effective <- standardize_interest(
      type = rate_type,
      rate = rate,
      m = m
    )

    n_old <- if (is.infinite(term_years)) NULL else term_years

    apv_benefits <- insurance_variable_k(
      lt = mortality_table,
      x = age,
      i = i_effective,
      benefit = benefit,
      n = n_old,
      m = deferral_years,
      k = payments_per_year,
      frac = frac
    )

    if (is.null(premium_term_years)) {
      if (!is.infinite(term_years)) {
        premium_term_years <- term_years
      } else if (!is.function(benefit)) {
        premium_term_years <- as.integer(length(benefit) / payments_per_year)
      } else {
        stop(
          "Cannot infer `premium_term_years` for variable benefits supplied as a function.",
          call. = FALSE
        )
      }
    }
  } else {
    apv_benefits <- insurance_x(
      mortality_table = mortality_table,
      age = age,
      rate = rate,
      rate_type = rate_type,
      m = m,
      term_years = term_years,
      deferral_years = deferral_years,
      insurance_type = insurance_type,
      benefit = benefit,
      output = "value"
    )

    if (is.null(premium_term_years)) {
      premium_term_years <- if (insurance_type == "whole") Inf else term_years
    }
  }

  # -------------------------------------------------------------------------
  # 1b. Premium-paying horizon validation
  # -------------------------------------------------------------------------

  if (insurance_type %in% c("term", "endowment", "variable_k") &&
      !is.infinite(term_years)) {
    if (is.infinite(premium_term_years)) {
      stop(
        "`premium_term_years` cannot be `Inf` for a finite-term insurance contract.",
        call. = FALSE
      )
    }

    premium_start_time <- if (premium_start == "issue") {
      0L
    } else {
      deferral_years
    }

    premium_end_time <- premium_start_time + premium_term_years
    coverage_end_time <- deferral_years + term_years

    if (premium_end_time > coverage_end_time) {
      stop(
        "For finite-term insurance contracts, premium payments must not extend ",
        "beyond the end of coverage. With the current inputs, premiums end at ",
        "time ", premium_end_time, ", while coverage ends at time ",
        coverage_end_time, ".",
        call. = FALSE
      )
    }
  }


  # -------------------------------------------------------------------------
  # 2. APV of premium annuity, valued at issue
  # -------------------------------------------------------------------------

  premium_deferral_years <- if (premium_start == "issue") {
    0L
  } else {
    deferral_years
  }

  apv_premiums <- annuity_x(
    mortality_table = mortality_table,
    age = age,
    rate = rate,
    rate_type = rate_type,
    m = m,
    term_years = premium_term_years,
    deferral_years = premium_deferral_years,
    payments_per_year = payments_per_year,
    timing = premium_timing,
    woolhouse = woolhouse,
    output = "value"
  )

  if (!is.finite(apv_premiums) || apv_premiums <= 0) {
    stop("APV of premium annuity is nonpositive or not finite.", call. = FALSE)
  }

  premium <- apv_benefits / apv_premiums

  if (output == "value") {
    return(premium)
  }

  tibble::tibble(
    age = age,
    rate = rate,
    rate_type = rate_type,
    m = m,
    insurance_type = insurance_type,
    benefit = if (is.function(benefit)) NA_real_ else suppressWarnings(as.numeric(benefit)[1]),
    term_years = term_years,
    deferral_years = deferral_years,
    payments_per_year = payments_per_year,
    frac = frac,
    premium_timing = premium_timing,
    premium_start = premium_start,
    premium_term_years = premium_term_years,
    woolhouse = woolhouse,
    premium = premium,
    premium_annual = payments_per_year * premium,
    apv_benefits = apv_benefits,
    apv_premiums = apv_premiums
  )
}
