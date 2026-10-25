#' Actuarial present value of a life insurance
#'
#' Computes the actuarial present value of a discrete single-life insurance
#' using a life table and compact actuarial notation.
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
#' @param lt A life table as produced by \code{\link{lifetable}}. It must
#'   contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age at issue, or a vector of scenario ages.
#' @param i Numeric annual interest-rate input, scalar or by scenario.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param n Integer insurance term in years. Required for
#'   \code{type = "term"} and \code{type = "endowment"}. Use \code{Inf} only
#'   for whole-life insurance.
#' @param h Integer deferment period in years.
#' @param type Character string. One of \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param benefit Numeric benefit amount, scalar or by scenario. A vector
#'   represents separate contracts, not a benefit schedule within one contract.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric APV. If
#'   \code{TRUE}, returns one tibble row per scenario with intermediate quantities.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{term_years}, \code{deferral_years}, \code{insurance_type}, and
#'   \code{output}.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of actuarial present values,
#' of length one for a single scenario.
#'
#' If \code{tidy = TRUE}, a tibble with one row per scenario, the main input values,
#' equivalent interest rate, deferral factor, pure endowment factor,
#' annuity-due value used in the standard identities, and APV.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age at
#' issue, \code{i} is the interest-rate input, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' \code{n} is the insurance term, and \code{h} is the deferment period.
#'
#' The function computes APVs directly from \code{lx}. For a deferred
#' insurance, the value at age \code{x + h} is multiplied by
#' \deqn{v^h\,{}_hp_x.}
#'
#' For whole-life insurance, the death benefit APV at the deferred starting age
#' is computed over the available life-table horizon. For term and endowment
#' insurance, the death benefit is computed over the first \code{n} years.
#' Endowment insurance additionally includes the pure endowment component
#' \deqn{v^n\,{}_np_x.}
#'
#' \strong{Scenarios.} \code{x}, \code{i}, \code{i_type}, \code{m},
#' \code{n}, \code{h}, \code{type}, and \code{benefit} may have a common
#' positive length; only length-one arguments recycle. The table and \code{tidy}
#' control are shared. Scenarios are paired by position, not expanded to a
#' Cartesian product. Missing or invalid cases are rejected with the scenario
#' number. Use the numeric result in \code{dplyr::mutate()}, or request a tibble
#' for intermediate values. Explicit interest arguments override contract defaults.
#'
#' \strong{Valuation and scope.} For annual deaths in interval \eqn{[r,r+1)},
#' the unconditional payment probability is
#' \eqn{(l_{x+r}-l_{x+r+1})/l_x}. The benefit is discounted at \eqn{r+1},
#' not at the exact moment of death. Decrements are formed before division,
#' and log ratios and original interest rates retain representable APVs when
#' intermediate conditional values or discount factors overflow or underflow.
#' Intermediate diagnostics may therefore be zero or infinite while APV is finite.
#' Tables require an annual age grid, nonincreasing finite \code{lx}, and a
#' positive count at issue. For compatibility the terminal count one year beyond
#' the table is zero, even for an open input table. Extinct tails contribute zero.
#' A zero-term endowment pays the benefit at time \code{h} on survival;
#' zero-term death insurance pays nothing.
#'
#' Use \code{\link{insurance_variable_k}} for benefits varying by subperiod,
#' including level benefits paid at the end of the subperiod of death, and
#' \code{\link{t_Ex}} for pure endowments. Continuous death-time insurance,
#' fractional issue ages and time-varying interest curves require additional
#' model specifications and are not implied by this annual contract.
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
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   type = "whole"
#' )
#'
#' # 5-year term insurance
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 5,
#'   type = "term"
#' )
#'
#' # 5-year endowment insurance
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 5,
#'   type = "endowment"
#' )
#'
#' # Deferred whole-life insurance
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   h = 2,
#'   type = "whole"
#' )
#'
#' # Tidy output
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 5,
#'   type = "term",
#'   tidy = TRUE
#' )
#'
#' # Parallel single-life scenarios with a shared mortality table
#' insurance_x(lt, x = c(60, 61), i = 0.04, n = c(2, 3),
#'             type = "term", benefit = c(10000, 20000), tidy = TRUE)
#'
#' @export
insurance_x <- function(
    lt,
    x,
    i,
    i_type = "effective",
    m = 1L,
    n = Inf,
    h = 0L,
    type = c("whole", "term", "endowment"),
    benefit = 1,
    tidy = FALSE,
    ...
) {
  dots <- list(...)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "mortality_table",
    "age",
    "rate",
    "rate_type",
    "term_years",
    "deferral_years",
    "insurance_type",
    "output"
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

  if (!is.null(dots$mortality_table)) {
    if (!missing(lt)) {
      stop("Provide only one of `lt` or deprecated `mortality_table`.", call. = FALSE)
    }
    lt <- dots$mortality_table
  }

  if (!is.null(dots$age)) {
    if (!missing(x)) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }
    x <- dots$age
  }

  if (!is.null(dots$rate)) {
    if (!missing(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }
    i <- dots$rate
  }

  if (!is.null(dots$rate_type)) {
    if (!identical(i_type, "effective")) {
      stop("Provide only one of `i_type` or deprecated `rate_type`.", call. = FALSE)
    }
    i_type <- dots$rate_type
  }

  if (!is.null(dots$term_years)) {
    if (!is.infinite(n)) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }
    n <- dots$term_years
  }

  if (!is.null(dots$deferral_years)) {
    if (!identical(h, 0L) && !identical(h, 0)) {
      stop("Provide only one of `h` or deprecated `deferral_years`.", call. = FALSE)
    }
    h <- dots$deferral_years
  }

  if (!is.null(dots$insurance_type)) {
    type <- dots$insurance_type
  }

  if (!is.null(dots$output)) {
    if (!identical(tidy, FALSE)) {
      stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
    }

    output <- match.arg(dots$output, c("value", "table"))
    tidy <- identical(output, "table")
  }

  if (missing(type) && is.null(dots$insurance_type)) type <- "whole"

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (!missing(lt) &&
      exists(".as_life_contract", mode = "function") &&
      .as_life_contract(lt)) {
    contract <- lt

    if (!identical(contract$lives, "single")) {
      stop(
        "`insurance_x()` currently supports only single-life `life_contract()` objects.",
        call. = FALSE
      )
    }

    lt <- contract$mortality_table

    if (missing(x) || is.null(x)) {
      x <- contract$age %||% contract$x
    }

    if (missing(i) || is.null(i)) {
      i <- contract$rate %||% contract$i
    }

    if ((missing(i_type) || is.null(i_type)) && is.null(dots$rate_type)) {
      i_type <- contract$rate_type %||% contract$i_type %||% i_type
    }

    if (missing(m) || is.null(m)) {
      m <- contract$m %||% m
    }
  }

  # -------------------------------------------------------------------------
  # Basic validation
  # -------------------------------------------------------------------------

  if (missing(lt)) {
    stop("`lt` must be provided.", call. = FALSE)
  }

  if (missing(x)) {
    stop("`x` must be provided.", call. = FALSE)
  }

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  if (!is.data.frame(lt)) {
    stop("`lt` must be a data.frame or tibble.", call. = FALSE)
  }

  if (!all(c("x", "lx") %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
  }

  scenarios <- .single_life_scenarios(insurance_x, lt,
    list(x = x, i = i, i_type = i_type, m = m, n = n, h = h,
         type = type, benefit = benefit), tidy)
  if (!is.null(scenarios)) return(scenarios)
  type <- match.arg(type, c("whole", "term", "endowment"))

  if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      abs(x - round(x)) > 1e-10 || abs(x) > .Machine$integer.max) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(i) ||
      length(i) != 1L ||
      is.na(i) ||
      !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) ||
      length(i_type) != 1L ||
      is.na(i_type)) {
    stop("`i_type` must be a single character string.", call. = FALSE)
  }

  valid_i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  if (!i_type %in% valid_i_type) {
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
      abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(h) ||
      length(h) != 1L ||
      is.na(h) ||
      !is.finite(h) ||
      h < 0 ||
      abs(h - round(h)) > 1e-10 || h > .Machine$integer.max) {
    stop("`h` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit)) {
    stop("`benefit` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n < 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) || abs(n - round(n)) > 1e-10))) {
    stop(
      "`n` must be `Inf` or a single nonnegative integer.",
      call. = FALSE
    )
  }

  if (type %in% c("term", "endowment") && is.infinite(n)) {
    stop("`n` must be finite for term and endowment insurance.", call. = FALSE)
  }

  x <- as.integer(round(x))
  m <- if (m <= .Machine$integer.max) as.integer(round(m)) else round(m)
  h <- as.integer(round(h))

  if (!is.infinite(n)) {
    n <- if (n <= .Machine$integer.max) as.integer(round(n)) else round(n)
  }

  # -------------------------------------------------------------------------
  # Interest conversion
  # -------------------------------------------------------------------------

  basis <- .interest_basis(i_type, i, m)
  i_effective <- basis$effective
  d_effective <- -expm1(-basis$growth)

  # -------------------------------------------------------------------------
  # Life table preparation
  # -------------------------------------------------------------------------

  lt <- lt[order(lt$x), , drop = FALSE]

  if (!is.numeric(lt$x)) {
    stop("Column `x` in `lt` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(lt$lx)) {
    stop("Column `lx` in `lt` must be numeric.", call. = FALSE)
  }

  if (any(is.na(lt$x)) || any(!is.finite(lt$x))) {
    stop("Column `x` must contain finite non-missing values.", call. = FALSE)
  }

  if (!nrow(lt) || any(abs(lt$x) > .Machine$integer.max) ||
      any(abs(lt$x - round(lt$x)) > 1e-10)) {
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
  if (any(diff(lx) > 0)) {
    stop("Column `lx` must be nonincreasing.", call. = FALSE)
  }
  if (any(diff(ages) != 1)) {
    stop("Life table ages must form an annual age grid.", call. = FALSE)
  }
  omega <- max(ages)

  get_lx <- function(current_age) {
    idx <- match(current_age, ages)

    if (!is.na(idx)) {
      return(lx[[idx]])
    }

    if (current_age == as.double(omega) + 1) {
      return(0)
    }

    NA_real_
  }

  log_p_int <- function(current_age, tt) {
    if (tt == 0L) {
      return(0)
    }

    l0 <- get_lx(current_age)
    l1 <- get_lx(as.double(current_age) + tt)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    .single_life_log_ratio(l1, l0)
  }

  annuity_due_from <- function(current_age, years) {
    years <- as.integer(round(years))

    if (years <= 0L) {
      return(0)
    }

    sum(vapply(0:(years - 1L), function(r) {
      px_r <- log_p_int(current_age, r)

      if (is.na(px_r)) {
        stop("The life table does not support the requested annuity APV.", call. = FALSE)
      }

      if (px_r == -Inf) 0 else exp(px_r - basis$growth * r)
    }, numeric(1L)))
  }

  death_apv_from <- function(current_age, years, origin = current_age,
                            offset = 0, amount = 1) {
    years <- as.integer(round(years))

    if (years <= 0L) {
      return(0)
    }

    sum(vapply(0:(years - 1L), function(r) {
      deaths <- get_lx(as.double(current_age) + r) -
        get_lx(as.double(current_age) + r + 1)
      if (is.na(deaths)) {
        stop("The life table does not support the requested insurance APV.", call. = FALSE)
      }
      if (deaths == 0 || amount == 0) return(0)

      .cash_flow_value(amount, .single_life_log_ratio(deaths, get_lx(origin)) -
                         basis$growth * (offset + r + 1))
    }, numeric(1L)))
  }

  # -------------------------------------------------------------------------
  # Deferral
  # -------------------------------------------------------------------------

  if (is.na(get_lx(x)) || get_lx(x) <= 0) {
    stop("The issue age `x` must have positive `lx` in the life table.", call. = FALSE)
  }
  log_survival_h <- log_p_int(x, h)
  deferment_factor <- if (isTRUE(log_survival_h == -Inf)) 0 else
    exp(log_survival_h - basis$growth * h)

  if (is.na(deferment_factor)) {
    stop(
      "The deferral age `x + h` is outside the life table or `lx(x)` is zero.",
      call. = FALSE
    )
  }

  start_x <- as.double(x) + h
  if (abs(start_x) <= .Machine$integer.max) start_x <- as.integer(start_x)

  # If deferment reaches a point where survival is zero, the APV is zero.
  if (log_survival_h == -Inf) {
    result <- 0

    if (!tidy) {
      return(result)
    }

    return(tibble::tibble(
      x = x,
      age = x,
      i = i,
      rate = i,
      i_type = i_type,
      rate_type = i_type,
      m = m,
      i_effective = i_effective,
      d_effective = d_effective,
      n = n,
      term_years = n,
      h = h,
      deferral_years = h,
      start_x = start_x,
      start_age = start_x,
      type = type,
      insurance_type = type,
      benefit = benefit,
      deferment_factor = deferment_factor,
      annuity_due_value = NA_real_,
      pure_endowment_factor = NA_real_,
      value_at_start = 0,
      apv = result
    ))
  }

  if (start_x > omega) {
    stop(
      "The deferred starting age `x + h` is outside the available life table.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Insurance value at the deferred starting age
  # -------------------------------------------------------------------------

  pure_endowment_factor <- NA_real_

  contract_years <- if (type == "whole") {
    as.double(omega) - start_x + 1
  } else {
    n
  }

  if (contract_years < 0L) {
    stop("The requested contract term is invalid.", call. = FALSE)
  }

  if (type %in% c("term", "endowment") &&
      as.double(start_x) + contract_years > as.double(omega) + 1) {
    stop(
      "The requested term exceeds the available life table horizon.",
      call. = FALSE
    )
  }

  if (type == "whole" && contract_years == 0L) {
    value_at_start <- 0
    annuity_due_value <- 0
  } else if (type %in% c("term", "endowment") && contract_years == 0L) {
    # Approved C10: an empty endowment matures at the deferred start.
    value_at_start <- if (type == "endowment") 1 else 0
    annuity_due_value <- 0
    pure_endowment_factor <- 1
  } else {
    death_apv <- death_apv_from(
      current_age = start_x,
      years = contract_years
    )

    annuity_due_value <- annuity_due_from(
      current_age = start_x,
      years = contract_years
    )

    if (type == "endowment") {
      n_p_start <- log_p_int(start_x, contract_years)

      if (is.na(n_p_start)) {
        stop(
          "The life table does not support the survival probability needed ",
          "for this endowment insurance calculation.",
          call. = FALSE
        )
      }

      pure_endowment_factor <- if (n_p_start == -Inf) 0 else
        exp(n_p_start - basis$growth * contract_years)
      value_at_start <- death_apv + pure_endowment_factor
    } else {
      value_at_start <- death_apv
    }
  }

  result <- death_apv_from(start_x, contract_years, origin = x,
                           offset = h, amount = benefit)
  if (type == "endowment") {
    log_survival <- log_p_int(x, h + contract_years)
    if (log_survival != -Inf && benefit != 0) {
      result <- result + .cash_flow_value(benefit,
        log_survival - basis$growth * (h + contract_years))
    }
  }

  if (!tidy) {
    return(result)
  }

  tibble::tibble(
    x = x,
    age = x,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective,
    d_effective = d_effective,
    n = n,
    term_years = n,
    h = h,
    deferral_years = h,
    start_x = start_x,
    start_age = start_x,
    type = type,
    insurance_type = type,
    benefit = benefit,
    deferment_factor = deferment_factor,
    annuity_due_value = annuity_due_value,
    pure_endowment_factor = pure_endowment_factor,
    value_at_start = value_at_start,
    apv = result
  )
}
