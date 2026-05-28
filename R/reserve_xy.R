#' Benefit reserve schedule for two-life insurance
#'
#' Computes terminal benefit reserves at selected policy durations for a
#' fully discrete two-life insurance contract, assuming independent future
#' lifetimes.
#'
#' @param mortality_table Either a single life table used for both lives, or a
#'   list of two life tables `list(table_x, table_y)`. Each table must contain
#'   columns `x` and `lx`. A `tidyact_life_contract` object created by
#'   [life_contract()] is also accepted.
#' @param age_x Integer actuarial age for the first life at issue.
#' @param age_y Integer actuarial age for the second life at issue.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are `"effective"`, `"nominal_interest"`, `"nominal_discount"`, and `"force"`.
#' @param m Positive integer. Compounding frequency for nominal rates.
#' @param insurance_type Insurance type. One of `"whole"`, `"term"`, or
#'   `"endowment"`.
#' @param cohort Status definition. Use `"first"` for joint-life status or
#'   `"last"` for last-survivor status.
#' @param term_years Insurance term in years. Required for term and endowment
#'   insurance. For whole life insurance, use `Inf` or omit it.
#' @param benefit Numeric benefit amount.
#' @param premium Net premium per payment. If `NULL`, it is computed internally
#'   using [premium_xy()].
#' @param premium_term_years Premium-paying term in years. If `NULL`, premiums
#'   are payable for the full contract term.
#' @param payments_per_year Number of premium payments per year. Default is 1.
#' @param premium_timing Timing of premium payments. Use `"due"` or
#'   `"immediate"`. The recursive method currently requires annual due
#'   premiums.
#' @param at Integer vector of policy durations at which to compute reserves.
#'   If `NULL`, reserves are computed for all integer durations.
#' @param method Computation method. Use `"prospective"` or `"recursive"`.
#' @param output Character string. Use `"table"` for a tibble schedule or
#'   `"value"` for a named numeric vector.
#'
#' @details
#' The prospective reserve at duration `k` is computed as:
#'
#' \deqn{
#' {}_kV = APV(\text{future benefits}) -
#' P \cdot APV(\text{future premiums}).
#' }
#'
#' Future benefit values are computed with [insurance_xy()] at shifted ages
#' `age_x + k` and `age_y + k`. Future premium values are computed with
#' [annuity_xy()] on the same two-life status.
#'
#' The recursive method is implemented for annual due premiums and follows the
#' usual one-year recursion for the selected two-life status.
#'
#' This function assumes independent future lifetimes.
#'
#' @return
#' If `output = "table"`, a tibble with reserve schedule details.
#' If `output = "value"`, a named numeric vector.
#'
#' @seealso [reserve_x()], [premium_xy()], [insurance_xy()], [annuity_xy()],
#'   [t_pxy()]
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:70,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000,
#'          86000, 81000, 75000, 68000, 60000)
#' )
#'
#' reserve_xy(
#'   mortality_table = lt,
#'   age_x = 60,
#'   age_y = 62,
#'   rate = 0.06,
#'   insurance_type = "term",
#'   cohort = "first",
#'   term_years = 4
#' )
#'
#' reserve_xy(
#'   mortality_table = lt,
#'   age_x = 60,
#'   age_y = 62,
#'   rate = 0.06,
#'   insurance_type = "endowment",
#'   cohort = "last",
#'   term_years = 5,
#'   benefit = 100000
#' )
#'
#' @export
reserve_xy <- function(
    mortality_table,
    age_x = NULL,
    age_y = NULL,
    rate = NULL,
    rate_type = NULL,
    m = NULL,
    insurance_type = c("whole", "term", "endowment"),
    cohort = c("first", "last"),
    term_years = Inf,
    benefit = 1,
    premium = NULL,
    premium_term_years = NULL,
    payments_per_year = 1L,
    premium_timing = c("due", "immediate"),
    at = NULL,
    method = c("prospective", "recursive"),
    output = c("table", "value")
) {
  insurance_type <- match.arg(insurance_type)
  cohort <- match.arg(cohort)
  premium_timing <- match.arg(premium_timing)
  method <- match.arg(method)
  output <- match.arg(output)

  # -------------------------------------------------------------------------
  # Resolve life_contract input
  # -------------------------------------------------------------------------

  if (inherits(mortality_table, "tidyact_life_contract")) {
    contract <- mortality_table

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop(
        "`reserve_xy()` requires a two-life `life_contract()` object.",
        call. = FALSE
      )
    }

    mortality_table <- contract$mortality_table

    if (is.null(age_x)) age_x <- contract$age_x
    if (is.null(age_y)) age_y <- contract$age_y
    if (is.null(rate)) rate <- contract$rate
    if (is.null(rate_type)) rate_type <- contract$rate_type
    if (is.null(m)) m <- contract$m

    if (contract$lives == "last_survivor" && missing(cohort)) {
      cohort <- "last"
    }
  }

  if (is.null(rate_type)) rate_type <- "effective"
  if (is.null(m)) m <- 1L

  # -------------------------------------------------------------------------
  # Resolve life tables
  # -------------------------------------------------------------------------

  if (is.data.frame(mortality_table)) {
    if (!all(c("x", "lx") %in% names(mortality_table))) {
      stop("`mortality_table` must contain columns `x` and `lx`.", call. = FALSE)
    }

    table_use <- mortality_table
    table_x <- mortality_table
    table_y <- mortality_table
  } else if (
    is.list(mortality_table) &&
    length(mortality_table) == 2L &&
    all(vapply(mortality_table, is.data.frame, logical(1L)))
  ) {
    if (!all(c("x", "lx") %in% names(mortality_table[[1L]]))) {
      stop("The first life table must contain columns `x` and `lx`.", call. = FALSE)
    }

    if (!all(c("x", "lx") %in% names(mortality_table[[2L]]))) {
      stop("The second life table must contain columns `x` and `lx`.", call. = FALSE)
    }

    table_use <- mortality_table
    table_x <- mortality_table[[1L]]
    table_y <- mortality_table[[2L]]
  } else {
    stop(
      "`mortality_table` must be one life table, a list of two life tables, ",
      "or a `tidyact_life_contract` object.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (is.null(rate) ||
      !is.numeric(rate) ||
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

  if (is.null(age_x) ||
      !is.numeric(age_x) ||
      length(age_x) != 1L ||
      is.na(age_x) ||
      !is.finite(age_x) ||
      abs(age_x - round(age_x)) > 1e-10) {
    stop("`age_x` must be a single integer age.", call. = FALSE)
  }

  if (is.null(age_y) ||
      !is.numeric(age_y) ||
      length(age_y) != 1L ||
      is.na(age_y) ||
      !is.finite(age_y) ||
      abs(age_y - round(age_y)) > 1e-10) {
    stop("`age_y` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      is.na(m) ||
      !is.finite(m) ||
      m < 1 ||
      abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit) ||
      benefit <= 0) {
    stop("`benefit` must be a single positive finite number.", call. = FALSE)
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
      term_years < 0 ||
      (!is.infinite(term_years) &&
       (!is.finite(term_years) ||
        abs(term_years - round(term_years)) > 1e-10))) {
    stop("`term_years` must be `Inf` or a single nonnegative integer.", call. = FALSE)
  }

  if (!is.null(premium)) {
    if (!is.numeric(premium) ||
        length(premium) != 1L ||
        is.na(premium) ||
        !is.finite(premium) ||
        premium < 0) {
      stop("`premium` must be NULL or a single nonnegative finite number.", call. = FALSE)
    }
  }

  age_x <- as.integer(round(age_x))
  age_y <- as.integer(round(age_y))
  m <- as.integer(round(m))
  payments_per_year <- as.integer(round(payments_per_year))

  if (!is.infinite(term_years)) {
    term_years <- as.integer(round(term_years))
  }

  if (insurance_type %in% c("term", "endowment") && is.infinite(term_years)) {
    stop(
      "`term_years` must be finite for term and endowment insurance.",
      call. = FALSE
    )
  }

  if (method == "recursive" &&
      (payments_per_year != 1L || premium_timing != "due")) {
    stop(
      "The recursive method currently requires annual due premiums: ",
      "`payments_per_year = 1` and `premium_timing = 'due'`.",
      call. = FALSE
    )
  }

  i_effective <- standardize_interest(
    type = rate_type,
    rate = rate,
    m = m
  )

  if (i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Contract horizon
  # -------------------------------------------------------------------------

  omega_x <- max(table_x$x, na.rm = TRUE)
  omega_y <- max(table_y$x, na.rm = TRUE)

  horizon_x <- max(0L, omega_x - age_x)
  horizon_y <- max(0L, omega_y - age_y)

  max_term <- if (cohort == "first") {
    min(horizon_x, horizon_y)
  } else {
    max(horizon_x, horizon_y)
  }

  term_used <- if (insurance_type == "whole") {
    if (is.infinite(term_years)) {
      max_term
    } else {
      term_years
    }
  } else {
    term_years
  }

  if (term_used > max_term) {
    stop(
      "`term_years` exceeds the maximum term implied by the life tables for this status.",
      call. = FALSE
    )
  }

  term_used <- as.integer(round(term_used))

  if (is.null(premium_term_years)) {
    premium_term_years <- term_used
  } else {
    if (!is.numeric(premium_term_years) ||
        length(premium_term_years) != 1L ||
        is.na(premium_term_years) ||
        !is.finite(premium_term_years) ||
        premium_term_years < 1 ||
        abs(premium_term_years - round(premium_term_years)) > 1e-10) {
      stop(
        "`premium_term_years` must be NULL or a single positive integer.",
        call. = FALSE
      )
    }

    premium_term_years <- as.integer(round(premium_term_years))
  }

  if (premium_term_years > term_used) {
    stop("`premium_term_years` must not exceed the contract term.", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Net premium
  # -------------------------------------------------------------------------

  if (is.null(premium)) {
    premium <- premium_xy(
      mortality_table = table_use,
      age_x = age_x,
      age_y = age_y,
      rate = rate,
      rate_type = rate_type,
      m = m,
      insurance_type = insurance_type,
      cohort = cohort,
      term_years = if (insurance_type == "whole") Inf else term_used,
      benefit = benefit,
      premium_term_years = premium_term_years,
      payments_per_year = payments_per_year,
      premium_timing = premium_timing,
      output = "value"
    )
  }

  # -------------------------------------------------------------------------
  # Durations
  # -------------------------------------------------------------------------

  if (is.null(at)) {
    duration_vec <- 0:term_used
  } else {
    if (!is.numeric(at) ||
        any(is.na(at)) ||
        any(!is.finite(at)) ||
        any(abs(at - round(at)) > 1e-10)) {
      stop("`at` must be NULL or an integer-valued numeric vector.", call. = FALSE)
    }

    duration_vec <- sort(unique(as.integer(round(at))))

    if (any(duration_vec < 0) || any(duration_vec > term_used)) {
      stop("`at` durations must be between 0 and the contract term.", call. = FALSE)
    }
  }

  status <- if (cohort == "first") "joint" else "last"

  # -------------------------------------------------------------------------
  # Prospective method
  # -------------------------------------------------------------------------

  if (method == "prospective") {
    reserves <- vapply(duration_vec, function(duration) {
      if (duration == 0L) return(0)
      if (insurance_type == "endowment" && duration == term_used) return(benefit)
      if (insurance_type == "term" && duration >= term_used) return(0)
      if (duration >= term_used) return(0)

      remaining_term <- term_used - duration

      apv_benefits <- insurance_xy(
        mortality_table = table_use,
        age_x = age_x + duration,
        age_y = age_y + duration,
        rate = rate,
        rate_type = rate_type,
        m = m,
        insurance_type = insurance_type,
        cohort = cohort,
        term_years = if (insurance_type == "whole") Inf else remaining_term,
        benefit = benefit,
        output = "value"
      )

      remaining_premium_term <- max(0L, premium_term_years - duration)

      if (remaining_premium_term == 0L) {
        apv_premiums <- 0
      } else {
        apv_premiums <- premium * annuity_xy(
          mortality_table = table_use,
          age_x = age_x + duration,
          age_y = age_y + duration,
          rate = rate,
          rate_type = rate_type,
          m = m,
          cohort = cohort,
          term_years = remaining_premium_term,
          payments_per_year = payments_per_year,
          timing = premium_timing,
          output = "value"
        )
      }

      apv_benefits - apv_premiums
    }, numeric(1L))
  } else {
    # -----------------------------------------------------------------------
    # Recursive method
    # -----------------------------------------------------------------------

    reserves_full <- numeric(term_used + 1L)
    reserves_full[[1L]] <- 0

    for (duration in 0:(term_used - 1L)) {
      premium_due <- if (duration < premium_term_years) premium else 0
      benefit_due <- benefit

      p_status <- t_pxy(
        lt = table_use,
        x = age_x + duration,
        y = age_y + duration,
        t = 1,
        frac = "UDD",
        status = status
      )

      if (is.na(p_status) || p_status <= 0) {
        reserves_full[(duration + 2L):(term_used + 1L)] <- 0
        break
      }

      q_status <- 1 - p_status

      reserves_full[[duration + 2L]] <- (
        (reserves_full[[duration + 1L]] + premium_due) *
          (1 + i_effective) -
          benefit_due * q_status
      ) / p_status
    }

    if (insurance_type == "endowment") {
      reserves_full[[term_used + 1L]] <- benefit
    }

    if (insurance_type == "term") {
      reserves_full[[term_used + 1L]] <- 0
    }

    reserves <- reserves_full[duration_vec + 1L]
  }

  # -------------------------------------------------------------------------
  # Output
  # -------------------------------------------------------------------------

  if (output == "value") {
    names(reserves) <- paste0("duration_", duration_vec)
    return(reserves)
  }

  premium_paid <- vapply(duration_vec, function(duration) {
    if (duration >= term_used) {
      return(0)
    }

    if (duration < premium_term_years) {
      premium
    } else {
      0
    }
  }, numeric(1L))

  benefit_due <- vapply(duration_vec, function(duration) {
    if (duration >= term_used) {
      return(0)
    }

    benefit
  }, numeric(1L))

  tibble::tibble(
    duration = duration_vec,
    age_x = age_x + duration_vec,
    age_y = age_y + duration_vec,
    reserve = reserves,
    premium_paid = premium_paid,
    benefit_due = benefit_due,
    insurance_type = insurance_type,
    cohort = cohort,
    benefit = benefit,
    premium = premium,
    premium_term_years = premium_term_years,
    payments_per_year = payments_per_year,
    premium_timing = premium_timing,
    contract_term = term_used,
    method = method,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective
  )
}
