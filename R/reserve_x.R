#' Benefit reserve schedule for single-life insurance
#'
#' Computes terminal benefit reserves at selected policy durations for a
#' fully discrete single-life insurance contract.
#'
#' The function supports whole-life, term, and endowment insurance. Reserves
#' may be computed prospectively or recursively. Premiums are assumed payable
#' annually in advance. Limited-payment policies are supported through
#' \code{premium_term_years}.
#'
#' @param mortality_table A life table data frame with columns \code{x} and
#'   \code{lx}.
#' @param age Integer actuarial age at issue.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Compounding frequency for nominal rates. Ignored
#'   for \code{rate_type = "effective"} and \code{rate_type = "force"}.
#' @param insurance_type Character string. One of \code{"whole"},
#'   \code{"term"}, or \code{"endowment"}.
#' @param term_years Insurance term in years. Use \code{Inf} for whole-life
#'   insurance. For term and endowment insurance, this must be finite.
#' @param benefit Numeric scalar. Insurance benefit amount.
#' @param premium Optional numeric scalar. Premium per annual payment. If
#'   \code{NULL}, the net premium is computed internally using
#'   \code{\link{premium_x}}.
#' @param premium_term_years Optional premium-paying term in years. If
#'   \code{NULL}, premiums are payable for the full contract duration.
#' @param durations Optional integer vector of policy durations at which to
#'   compute reserves. If \code{NULL}, reserves are computed for all integer
#'   durations from issue to the contract horizon.
#' @param method Character string. Either \code{"prospective"} or
#'   \code{"recursive"}.
#' @param output Character string. Use \code{"table"} to return a reserve
#'   schedule, or \code{"value"} to return a named numeric vector.
#'
#' @return
#' If \code{output = "table"}, a tibble with one row per selected duration.
#' If \code{output = "value"}, a named numeric vector of reserves.
#'
#' @details
#' The prospective reserve is computed as
#' \deqn{{}_kV = APV_k(\text{future benefits}) - P\,APV_k(\text{future premiums}).}
#'
#' The recursive method uses the annual fully discrete recursion
#' \deqn{{}_{k+1}V = \frac{({}_kV + P_k)(1+i) - b_{k+1}q_{x+k}}{p_{x+k}}.}
#'
#' @seealso \code{\link{premium_x}}, \code{\link{insurance_x}},
#'   \code{\link{annuity_x}}, \code{\link{t_px}}
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
#' reserve_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   insurance_type = "whole"
#' )
#'
#' reserve_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   insurance_type = "endowment",
#'   term_years = 5,
#'   benefit = 100000
#' )
#'
#' @export
reserve_x <- function(
    mortality_table,
    age,
    rate,
    rate_type = "effective",
    m = 1L,
    insurance_type = c("whole", "term", "endowment"),
    term_years = Inf,
    benefit = 1,
    premium = NULL,
    premium_term_years = NULL,
    durations = NULL,
    method = c("prospective", "recursive"),
    output = c("table", "value")
) {
  insurance_type <- match.arg(insurance_type)
  method <- match.arg(method)
  output <- match.arg(output)

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (.as_life_contract(mortality_table)) {
    contract <- mortality_table

    if (!identical(contract$lives, "single")) {
      stop(
        "`reserve_x()` currently supports only single-life `life_contract()` objects.",
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

  if (!is.numeric(age) || length(age) != 1L || is.na(age) ||
      !is.finite(age) || abs(age - round(age)) > 1e-10) {
    stop("`age` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) ||
      !is.finite(rate)) {
    stop("`rate` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(rate_type) || length(rate_type) != 1L || is.na(rate_type)) {
    stop("`rate_type` must be a single character string.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m < 1 || abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) || length(benefit) != 1L || is.na(benefit) ||
      !is.finite(benefit) || benefit <= 0) {
    stop("`benefit` must be a single positive finite number.", call. = FALSE)
  }

  if (!is.null(premium) &&
      (!is.numeric(premium) || length(premium) != 1L || is.na(premium) ||
       !is.finite(premium))) {
    stop("`premium` must be NULL or a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(term_years) || length(term_years) != 1L || is.na(term_years) ||
      term_years <= 0 ||
      (!is.infinite(term_years) &&
       (!is.finite(term_years) || abs(term_years - round(term_years)) > 1e-10))) {
    stop("`term_years` must be `Inf` or a single positive integer.", call. = FALSE)
  }

  if (insurance_type %in% c("term", "endowment") && is.infinite(term_years)) {
    stop("`term_years` must be finite for term and endowment insurance.", call. = FALSE)
  }

  if (!is.null(premium_term_years) &&
      (!is.numeric(premium_term_years) || length(premium_term_years) != 1L ||
       is.na(premium_term_years) || premium_term_years <= 0 ||
       (!is.infinite(premium_term_years) &&
        (!is.finite(premium_term_years) ||
         abs(premium_term_years - round(premium_term_years)) > 1e-10)))) {
    stop(
      "`premium_term_years` must be NULL, Inf, or a single positive integer.",
      call. = FALSE
    )
  }

  age <- as.integer(round(age))
  m <- as.integer(round(m))

  if (!is.infinite(term_years)) {
    term_years <- as.integer(round(term_years))
  }

  if (!is.null(premium_term_years) && !is.infinite(premium_term_years)) {
    premium_term_years <- as.integer(round(premium_term_years))
  }

  # -------------------------------------------------------------------------
  # Interest and life table preparation
  # -------------------------------------------------------------------------

  i_effective <- standardize_interest(type = rate_type, rate = rate, m = m)

  if (!is.numeric(i_effective) || length(i_effective) != 1L ||
      is.na(i_effective) || !is.finite(i_effective) || i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  lt <- mortality_table[order(mortality_table$x), , drop = FALSE]

  if (!is.numeric(lt$x) || !is.numeric(lt$lx)) {
    stop("Columns `x` and `lx` in `mortality_table` must be numeric.", call. = FALSE)
  }

  if (any(is.na(lt$x)) || any(!is.finite(lt$x)) ||
      any(abs(lt$x - round(lt$x)) > 1e-10)) {
    stop("Column `x` must contain finite integer ages.", call. = FALSE)
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

  one_year_qx <- function(current_age) {
    l0 <- get_lx(current_age)
    l1 <- get_lx(current_age + 1L)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    1 - l1 / l0
  }

  # -------------------------------------------------------------------------
  # Contract and premium-paying horizons
  # -------------------------------------------------------------------------

  contract_horizon <- if (insurance_type == "whole") {
    omega - age
  } else {
    term_years
  }

  if (!is.finite(contract_horizon) || contract_horizon < 1) {
    stop(
      "The life table does not provide a positive contract horizon from `age`.",
      call. = FALSE
    )
  }

  if (age + contract_horizon > omega) {
    stop(
      "The requested contract horizon exceeds the available life table ages.",
      call. = FALSE
    )
  }

  if (is.null(premium_term_years)) {
    premium_horizon <- contract_horizon
  } else if (is.infinite(premium_term_years)) {
    premium_horizon <- contract_horizon
  } else {
    premium_horizon <- min(premium_term_years, contract_horizon)
  }

  premium_horizon <- as.integer(round(premium_horizon))

  # -------------------------------------------------------------------------
  # Premium
  # -------------------------------------------------------------------------

  premium_was_computed <- is.null(premium)

  if (is.null(premium)) {
    premium <- premium_x(
      mortality_table = mortality_table,
      age = age,
      rate = rate,
      rate_type = rate_type,
      m = m,
      insurance_type = insurance_type,
      benefit = benefit,
      term_years = if (insurance_type == "whole") Inf else term_years,
      premium_term_years = premium_horizon,
      payments_per_year = 1L,
      premium_timing = "due",
      output = "value"
    )
  }

  # -------------------------------------------------------------------------
  # Durations
  # -------------------------------------------------------------------------

  if (is.null(durations)) {
    duration_vec <- 0:contract_horizon
  } else {
    if (!is.numeric(durations) || any(is.na(durations)) ||
        any(!is.finite(durations)) ||
        any(abs(durations - round(durations)) > 1e-10)) {
      stop("`durations` must be an integer vector.", call. = FALSE)
    }

    duration_vec <- sort(unique(as.integer(round(durations))))

    if (any(duration_vec < 0) || any(duration_vec > contract_horizon)) {
      stop(
        "`durations` must lie between 0 and the contract horizon.",
        call. = FALSE
      )
    }
  }

  # -------------------------------------------------------------------------
  # Prospective method
  # -------------------------------------------------------------------------

  if (method == "prospective") {
    reserves <- vapply(duration_vec, function(k) {
      if (k >= contract_horizon) {
        if (insurance_type == "endowment" && k == contract_horizon) {
          return(benefit)
        }

        return(0)
      }

      remaining_contract <- contract_horizon - k

      apv_benefits <- if (insurance_type == "whole") {
        insurance_x(
          mortality_table = mortality_table,
          age = age + k,
          rate = rate,
          rate_type = rate_type,
          m = m,
          insurance_type = "whole",
          benefit = benefit,
          output = "value"
        )
      } else {
        insurance_x(
          mortality_table = mortality_table,
          age = age + k,
          rate = rate,
          rate_type = rate_type,
          m = m,
          insurance_type = insurance_type,
          term_years = remaining_contract,
          benefit = benefit,
          output = "value"
        )
      }

      remaining_premiums <- max(0L, premium_horizon - k)

      apv_premiums <- if (remaining_premiums == 0L) {
        0
      } else {
        premium * annuity_x(
          mortality_table = mortality_table,
          age = age + k,
          rate = rate,
          rate_type = rate_type,
          m = m,
          term_years = remaining_premiums,
          payments_per_year = 1L,
          timing = "due",
          woolhouse = "none",
          output = "value"
        )
      }

      apv_benefits - apv_premiums
    }, numeric(1L))

    if (premium_was_computed && 0L %in% duration_vec) {
      reserves[duration_vec == 0L] <- 0
    }
  } else {
    # -----------------------------------------------------------------------
    # Recursive method
    # -----------------------------------------------------------------------

    full_reserves <- numeric(contract_horizon + 1L)
    full_reserves[[1]] <- 0

    for (kk in 0:(contract_horizon - 1L)) {
      premium_paid <- if (kk < premium_horizon) premium else 0
      qx_k <- one_year_qx(age + kk)

      if (is.na(qx_k)) {
        stop("The life table does not support recursive reserve computation.", call. = FALSE)
      }

      px_k <- 1 - qx_k

      if (px_k <= 0) {
        if (kk + 2L <= length(full_reserves)) {
          full_reserves[(kk + 2L):length(full_reserves)] <- 0
        }
        break
      }

      full_reserves[[kk + 2L]] <-
        ((full_reserves[[kk + 1L]] + premium_paid) * (1 + i_effective) -
           benefit * qx_k) / px_k
    }

    if (insurance_type == "endowment") {
      full_reserves[[contract_horizon + 1L]] <- benefit
    }

    if (insurance_type == "term") {
      full_reserves[[contract_horizon + 1L]] <- 0
    }

    reserves <- full_reserves[duration_vec + 1L]
  }

  if (output == "value") {
    names(reserves) <- paste0("duration=", duration_vec)
    return(reserves)
  }

  premium_paid <- vapply(duration_vec, function(k) {
    if (k >= contract_horizon) {
      return(0)
    }

    if (k < premium_horizon) premium else 0
  }, numeric(1L))

  benefit_due <- vapply(duration_vec, function(k) {
    if (k >= contract_horizon) {
      return(0)
    }

    benefit
  }, numeric(1L))

  tibble::tibble(
    duration = duration_vec,
    age = age + duration_vec,
    reserve = reserves,
    premium_paid = premium_paid,
    benefit_due = benefit_due,
    insurance_type = insurance_type,
    benefit = benefit,
    premium = premium,
    premium_term_years = premium_horizon,
    contract_horizon = contract_horizon,
    method = method,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective
  )
}
