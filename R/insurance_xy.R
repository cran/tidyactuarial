#' Actuarial present value of a two-life insurance
#'
#' Computes the actuarial present value of a discrete two-life insurance with
#' benefit payable at the end of the year of the triggering death, assuming
#' independent future lifetimes.
#'
#' Supported contracts:
#' \itemize{
#'   \item whole life insurance,
#'   \item term insurance,
#'   \item endowment insurance.
#' }
#'
#' The cohort determines the two-life status:
#' \itemize{
#'   \item `cohort = "first"`: joint-life status; insurance is triggered by
#'     the first death.
#'   \item `cohort = "last"`: last-survivor status; insurance is triggered by
#'     the second death.
#' }
#'
#' @param mortality_table Either a single life table used for both lives, or a
#'   list of two life tables `list(table_x, table_y)`. Each table must contain
#'   columns `x` and `lx`. A `tidyact_life_contract` object created by
#'   [life_contract()] is also accepted.
#' @param age_x Integer actuarial age for the first life.
#' @param age_y Integer actuarial age for the second life.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are `"effective"`, `"nominal_interest"`, `"nominal_discount"`, and `"force"`.
#' @param m Positive integer. Compounding frequency for nominal rates.
#' @param insurance_type Character string. One of `"whole"`, `"term"`, or
#'   `"endowment"`.
#' @param cohort Character string. Use `"first"` for joint-life insurance or
#'   `"last"` for last-survivor insurance.
#' @param term_years Term in years. Required for term and endowment insurance.
#'   Use `Inf` or omit it for whole life insurance.
#' @param deferment_years Integer deferment period in years.
#' @param benefit Numeric scalar. Insurance benefit.
#' @param output Character string. Use `"value"` for a numeric APV or
#'   `"table"` for a one-row tibble.
#'
#' @details
#' This function uses the standard identities that express insurance values
#' through two-life annuity-due values.
#'
#' For a whole life contract:
#' \deqn{A = 1 - d \ddot{a}}
#'
#' For an n-year term insurance:
#' \deqn{A^1_{:\overline{n}|} =
#' 1 - d\ddot{a}_{:\overline{n}|} - v^n\,{}_np}
#'
#' For an endowment insurance:
#' \deqn{A_{:\overline{n}|} =
#' 1 - d\ddot{a}_{:\overline{n}|}}
#'
#' The function assumes independent future lifetimes.
#'
#' @return
#' If `output = "value"`, a numeric scalar.
#' If `output = "table"`, a one-row tibble.
#'
#' @seealso [annuity_xy()], [premium_xy()], [insurance_x()], [t_pxy()]
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:110,
#'   lx = seq(100000, 0, length.out = 51)
#' )
#'
#' insurance_xy(
#'   mortality_table = lt,
#'   age_x = 60,
#'   age_y = 62,
#'   rate = 0.06,
#'   insurance_type = "whole",
#'   cohort = "first"
#' )
#'
#' lt |>
#'   life_contract(lives = "joint", age_x = 60, age_y = 62, rate = 0.06) |>
#'   insurance_xy(
#'     insurance_type = "term",
#'     term_years = 4,
#'     cohort = "first"
#'   )
#'
#' @export
insurance_xy <- function(
    mortality_table,
    age_x = NULL,
    age_y = NULL,
    rate = NULL,
    rate_type = NULL,
    m = NULL,
    insurance_type = c("whole", "term", "endowment"),
    cohort = c("first", "last"),
    term_years = Inf,
    deferment_years = 0L,
    benefit = 1,
    output = c("value", "table")
) {
  insurance_type <- match.arg(insurance_type)
  cohort <- match.arg(cohort)
  output <- match.arg(output)

  # -------------------------------------------------------------------------
  # Resolve life_contract input
  # -------------------------------------------------------------------------

  if (inherits(mortality_table, "tidyact_life_contract")) {
    contract <- mortality_table

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop(
        "`insurance_xy()` requires a two-life `life_contract()` object.",
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
  # Resolve life table input for t_pxy() and annuity_xy()
  # -------------------------------------------------------------------------

  if (is.data.frame(mortality_table)) {
    if (!all(c("x", "lx") %in% names(mortality_table))) {
      stop("`mortality_table` must contain columns `x` and `lx`.", call. = FALSE)
    }

    table_use <- mortality_table
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

  if (!is.numeric(deferment_years) ||
      length(deferment_years) != 1L ||
      is.na(deferment_years) ||
      !is.finite(deferment_years) ||
      deferment_years < 0 ||
      abs(deferment_years - round(deferment_years)) > 1e-10) {
    stop("`deferment_years` must be a single nonnegative integer.", call. = FALSE)
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

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit) ||
      benefit < 0) {
    stop("`benefit` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (insurance_type %in% c("term", "endowment") && is.infinite(term_years)) {
    stop(
      "`term_years` must be finite for term and endowment insurance.",
      call. = FALSE
    )
  }

  age_x <- as.integer(round(age_x))
  age_y <- as.integer(round(age_y))
  m <- as.integer(round(m))
  deferment_years <- as.integer(round(deferment_years))

  if (!is.infinite(term_years)) {
    term_years <- as.integer(round(term_years))
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

  v_fun <- function(tt) (1 + i_effective)^(-tt)
  d <- i_effective / (1 + i_effective)

  status <- if (cohort == "first") "joint" else "last"

  # -------------------------------------------------------------------------
  # Deferment
  # -------------------------------------------------------------------------

  p_def <- t_pxy(
    lt = table_use,
    x = age_x,
    y = age_y,
    t = deferment_years,
    frac = "UDD",
    status = status
  )

  if (is.na(p_def)) {
    stop("Cannot compute survival to deferment under the life table.", call. = FALSE)
  }

  deferment_factor <- v_fun(deferment_years) * p_def

  start_age_x <- age_x + deferment_years
  start_age_y <- age_y + deferment_years

  # -------------------------------------------------------------------------
  # APV computation
  # -------------------------------------------------------------------------

  if (insurance_type == "whole") {
    annuity_due <- annuity_xy(
      mortality_table = table_use,
      age_x = start_age_x,
      age_y = start_age_y,
      rate = rate,
      rate_type = rate_type,
      m = m,
      cohort = cohort,
      term_years = Inf,
      deferment_years = 0L,
      payments_per_year = 1L,
      timing = "due",
      woolhouse = "none",
      output = "value"
    )

    unit_apv <- deferment_factor * (1 - d * annuity_due)
  } else if (term_years == 0L) {
    unit_apv <- 0
  } else {
    annuity_due_term <- annuity_xy(
      mortality_table = table_use,
      age_x = start_age_x,
      age_y = start_age_y,
      rate = rate,
      rate_type = rate_type,
      m = m,
      cohort = cohort,
      term_years = term_years,
      deferment_years = 0L,
      payments_per_year = 1L,
      timing = "due",
      woolhouse = "none",
      output = "value"
    )

    if (insurance_type == "endowment") {
      unit_apv <- deferment_factor * (1 - d * annuity_due_term)
    } else {
      p_n <- t_pxy(
        lt = table_use,
        x = start_age_x,
        y = start_age_y,
        t = term_years,
        frac = "UDD",
        status = status
      )

      if (is.na(p_n)) {
        stop("Cannot compute n-year survival at deferred ages.", call. = FALSE)
      }

      unit_apv <- deferment_factor *
        (1 - d * annuity_due_term - v_fun(term_years) * p_n)
    }
  }

  result <- benefit * unit_apv

  if (output == "value") {
    return(result)
  }

  tibble::tibble(
    age_x = age_x,
    age_y = age_y,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    insurance_type = insurance_type,
    cohort = cohort,
    term_years = if (insurance_type == "whole") Inf else term_years,
    deferment_years = deferment_years,
    benefit = benefit,
    unit_apv = unit_apv,
    apv = result
  )
}
