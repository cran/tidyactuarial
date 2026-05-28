#' Actuarial present value of a life insurance
#'
#' Computes the actuarial present value of a discrete life insurance using a
#' life table.
#'
#' The benefit is paid at the end of the year of death for whole-life and term
#' insurance. For endowment insurance, the same benefit is paid either at death
#' within the term or at the end of the term if the life survives.
#'
#' Supported contracts:
#' \itemize{
#'   \item \code{"whole"}: whole-life insurance.
#'   \item \code{"term"}: n-year term insurance.
#'   \item \code{"endowment"}: n-year endowment insurance.
#' }
#'
#' @param mortality_table A life table as produced by \code{\link{lifetable}}.
#'   It must contain columns \code{x} and \code{lx}.
#' @param age Integer actuarial age.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Compounding frequency for nominal rates. Ignored
#'   for \code{rate_type = "effective"} and \code{rate_type = "force"}.
#' @param term_years Integer term in years. Required for
#'   \code{insurance_type = "term"} and \code{insurance_type = "endowment"}.
#'   Use \code{Inf} only for whole-life insurance.
#' @param deferral_years Integer deferral period in years.
#' @param insurance_type Character string. One of \code{"whole"},
#'   \code{"term"}, or \code{"endowment"}.
#' @param benefit Numeric scalar. Benefit amount.
#' @param output Character string. Use \code{"value"} to return a numeric APV
#'   or \code{"table"} to return a one-row tibble with intermediate quantities.
#'
#' @return
#' If \code{output = "value"}, a numeric scalar containing the actuarial present
#' value.
#'
#' If \code{output = "table"}, a one-row tibble with the main input values,
#' equivalent interest rate, deferral factor, pure endowment factor, annuity-due
#' value used in the identity, and APV.
#'
#' @details
#' This function uses standard annual discrete identities.
#'
#' For whole-life insurance,
#' \deqn{A_x = 1 - d\,\ddot{a}_x,}
#' where \eqn{d = i/(1+i)}.
#'
#' For n-year term insurance,
#' \deqn{A^1_{x:\overline{n}|} =
#' 1 - d\,\ddot{a}_{x:\overline{n}|} - v^n\,{}_np_x.}
#'
#' For n-year endowment insurance,
#' \deqn{A_{x:\overline{n}|} =
#' 1 - d\,\ddot{a}_{x:\overline{n}|}.}
#'
#' Deferral is handled by multiplying the value at age
#' \code{age + deferral_years} by \eqn{v^h {}_hp_x}, where
#' \eqn{h =} \code{deferral_years}.
#'
#' @seealso \code{\link{annuity_x}}, \code{\link{premium_x}},
#'   \code{\link{reserve_x}}, \code{\link{t_Ex}}, \code{\link{insurance_xy}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:65,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000)
#' )
#'
#' # Whole-life insurance
#' insurance_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   insurance_type = "whole"
#' )
#'
#' # 5-year term insurance
#' insurance_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   term_years = 5,
#'   insurance_type = "term"
#' )
#'
#' # 5-year endowment insurance
#' insurance_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   term_years = 5,
#'   insurance_type = "endowment"
#' )
#'
#' # Deferred whole-life insurance
#' insurance_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   deferral_years = 2,
#'   insurance_type = "whole"
#' )
#'
#' # Table output
#' insurance_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   term_years = 5,
#'   insurance_type = "term",
#'   output = "table"
#' )
#'
#' @export
insurance_x <- function(
    mortality_table,
    age,
    rate,
    rate_type = "effective",
    m = 1L,
    term_years = Inf,
    deferral_years = 0L,
    insurance_type = c("whole", "term", "endowment"),
    benefit = 1,
    output = c("value", "table")
) {
  insurance_type <- match.arg(insurance_type)
  output <- match.arg(output)

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (.as_life_contract(mortality_table)) {
    contract <- mortality_table

    if (!identical(contract$lives, "single")) {
      stop(
        "`insurance_x()` currently supports only single-life `life_contract()` objects.",
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

  # -------------------------------------------------------------------------
  # Basic validation
  # -------------------------------------------------------------------------

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

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit)) {
    stop("`benefit` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(term_years) ||
      length(term_years) != 1L ||
      is.na(term_years) ||
      term_years < 0 ||
      (!is.infinite(term_years) &&
       (!is.finite(term_years) ||
        abs(term_years - round(term_years)) > 1e-10))) {
    stop(
      "`term_years` must be `Inf` or a single nonnegative integer.",
      call. = FALSE
    )
  }

  if (insurance_type %in% c("term", "endowment") && is.infinite(term_years)) {
    stop(
      "`term_years` must be finite for term and endowment insurance.",
      call. = FALSE
    )
  }

  age <- as.integer(round(age))
  m <- as.integer(round(m))
  deferral_years <- as.integer(round(deferral_years))

  if (!is.infinite(term_years)) {
    term_years <- as.integer(round(term_years))
  }

  # -------------------------------------------------------------------------
  # Interest conversion
  # -------------------------------------------------------------------------

  i_effective <- standardize_interest(
    type = rate_type,
    rate = rate,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  v_fun <- function(tt) {
    (1 + i_effective)^(-tt)
  }

  d_effective <- i_effective / (1 + i_effective)

  # -------------------------------------------------------------------------
  # Life table preparation
  # -------------------------------------------------------------------------

  lt <- mortality_table[order(mortality_table$x), , drop = FALSE]

  if (!is.numeric(lt$x)) {
    stop("Column `x` in `mortality_table` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(lt$lx)) {
    stop("Column `lx` in `mortality_table` must be numeric.", call. = FALSE)
  }

  if (any(is.na(lt$x)) || any(!is.finite(lt$x))) {
    stop("Column `x` must contain finite non-missing values.", call. = FALSE)
  }

  if (any(abs(lt$x - round(lt$x)) > 1e-10)) {
    stop("Column `x` must contain integer ages.", call. = FALSE)
  }

  if (anyDuplicated(lt$x)) {
    stop("Life table ages in column `x` must be unique.", call. = FALSE)
  }

  if (any(is.na(lt$lx)) || any(!is.finite(lt$lx)) || any(lt$lx < 0)) {
    stop("Column `lx` must contain finite nonnegative values.", call. = FALSE)
  }

  ages <- as.integer(round(lt$x))
  lx <- as.numeric(lt$lx)
  omega <- max(ages)

  get_lx <- function(current_age) {
    idx <- match(current_age, ages)

    if (!is.na(idx)) {
      return(lx[[idx]])
    }

    if (current_age == omega + 1L) {
      return(0)
    }

    NA_real_
  }

  t_p_int <- function(current_age, tt) {
    if (tt == 0) {
      return(1)
    }

    l0 <- get_lx(current_age)
    l1 <- get_lx(current_age + tt)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    l1 / l0
  }

  # -------------------------------------------------------------------------
  # Deferral
  # -------------------------------------------------------------------------

  deferment_factor <- v_fun(deferral_years) *
    t_p_int(age, deferral_years)

  if (is.na(deferment_factor)) {
    stop(
      "The deferral age `age + deferral_years` is outside the life table ",
      "or `lx(age)` is zero.",
      call. = FALSE
    )
  }

  start_age <- age + deferral_years

  # -------------------------------------------------------------------------
  # Insurance value at the deferred starting age
  # -------------------------------------------------------------------------

  pure_endowment_factor <- NA_real_
  annuity_due_value <- NA_real_

  if (insurance_type == "whole") {
    annuity_due_value <- annuity_x(
      mortality_table = lt,
      age = start_age,
      rate = rate,
      rate_type = rate_type,
      m = m,
      term_years = Inf,
      deferral_years = 0L,
      payments_per_year = 1L,
      timing = "due",
      woolhouse = "none",
      output = "value"
    )

    value_at_start <- 1 - d_effective * annuity_due_value
  } else if (term_years == 0L) {
    value_at_start <- 0
    annuity_due_value <- 0
    pure_endowment_factor <- 1
  } else {
    annuity_due_value <- annuity_x(
      mortality_table = lt,
      age = start_age,
      rate = rate,
      rate_type = rate_type,
      m = m,
      term_years = term_years,
      deferral_years = 0L,
      payments_per_year = 1L,
      timing = "due",
      woolhouse = "none",
      output = "value"
    )

    n_p_start <- t_p_int(start_age, term_years)

    if (is.na(n_p_start)) {
      stop(
        "The life table does not support the survival probability needed ",
        "for this term insurance calculation.",
        call. = FALSE
      )
    }

    pure_endowment_factor <- v_fun(term_years) * n_p_start

    if (insurance_type == "endowment") {
      value_at_start <- 1 - d_effective * annuity_due_value
    } else {
      value_at_start <- 1 - d_effective * annuity_due_value -
        pure_endowment_factor
    }
  }

  result <- benefit * deferment_factor * value_at_start

  if (output == "value") {
    return(result)
  }

  tibble::tibble(
    age = age,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    d_effective = d_effective,
    term_years = term_years,
    deferral_years = deferral_years,
    start_age = start_age,
    insurance_type = insurance_type,
    benefit = benefit,
    deferment_factor = deferment_factor,
    annuity_due_value = annuity_due_value,
    pure_endowment_factor = pure_endowment_factor,
    value_at_start = value_at_start,
    apv = result
  )
}
