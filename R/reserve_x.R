#' Benefit reserve schedule for single-life insurance
#'
#' Computes terminal benefit reserves at selected policy durations for a fully
#' discrete single-life insurance contract, using compact actuarial notation.
#'
#' The function supports whole-life, term, and endowment insurance. Reserves may
#' be computed prospectively or recursively. Premiums are assumed payable
#' annually in advance. Limited-payment policies are supported through
#' \code{n_prem}.
#'
#' @param lt A life table data frame with columns \code{x} and \code{lx}.
#' @param x Integer actuarial age at issue.
#' @param i Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param type Character string. One of \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param n Insurance term in years. Use \code{Inf} for whole-life insurance.
#'   For term and endowment insurance, this must be finite.
#' @param benefit Numeric scalar. Insurance benefit amount.
#' @param P Optional numeric scalar. Premium per annual payment. If \code{NULL},
#'   the net premium is computed internally by the equivalence principle.
#' @param n_prem Optional premium-paying term in years. If \code{NULL}, premiums
#'   are payable for the full contract duration.
#' @param t Optional integer vector of policy durations at which to compute
#'   reserves. If \code{NULL}, reserves are computed for all integer durations
#'   from issue to the contract horizon.
#' @param method Character string. Either \code{"prospective"} or
#'   \code{"recursive"}.
#' @param tidy Logical scalar. If \code{TRUE}, returns a reserve schedule as a
#'   tibble. If \code{FALSE}, returns a named numeric vector.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{insurance_type}, \code{term_years}, \code{premium},
#'   \code{premium_term_years}, \code{durations}, and \code{output}.
#'
#' @return
#' If \code{tidy = TRUE}, a tibble with one row per selected duration.
#' If \code{tidy = FALSE}, a named numeric vector of reserves.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age at
#' issue, \code{i} is the interest-rate input, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' \code{n} is the contract term, \code{P} is the annual premium, and \code{t}
#' is the policy duration.
#'
#' The prospective reserve is computed as
#' \deqn{{}_tV_x =
#' APV_t(\text{future benefits}) -
#' P\,APV_t(\text{future premiums}).}
#'
#' The recursive method uses the annual fully discrete recursion
#' \deqn{{}_{k+1}V =
#' \frac{({}_kV + P_k)(1+i) - b_{k+1}q_{x+k}}{p_{x+k}}.}
#'
#' When \code{P = NULL}, the net premium is computed directly from the life
#' table by applying the equivalence principle at issue.
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
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   type = "whole"
#' )
#'
#' reserve_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   type = "endowment",
#'   n = 5,
#'   benefit = 100000
#' )
#'
#' reserve_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   type = "term",
#'   n = 5,
#'   benefit = 100000,
#'   t = c(0, 1, 2, 3, 4, 5),
#'   tidy = FALSE
#' )
#'
#' @export
reserve_x <- function(
    lt,
    x,
    i,
    i_type = "effective",
    m = 1L,
    type = c("whole", "term", "endowment"),
    n = Inf,
    benefit = 1,
    P = NULL,
    n_prem = NULL,
    t = NULL,
    method = c("prospective", "recursive"),
    tidy = TRUE,
    ...
) {
  dots <- list(...)

  # Use exact matching for deprecated arguments.
  # This avoids partial matching problems such as `premium` matching
  # `premium_term_years`.
  dot_has <- function(nm) {
    nm %in% names(dots)
  }

  dot_get <- function(nm) {
    dots[[nm]]
  }

  # local infix helper without importing anything
  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  # Use exact matching for deprecated arguments.
  # This avoids partial matching problems such as `premium` matching
  # `premium_term_years`.
  dot_has <- function(nm) {
    nm %in% names(dots)
  }

  dot_get <- function(nm) {
    dots[[nm]]
  }

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "mortality_table",
    "age",
    "rate",
    "rate_type",
    "insurance_type",
    "term_years",
    "premium",
    "premium_term_years",
    "durations",
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

  if (dot_has("mortality_table")) {
    if (!missing(lt)) {
      stop("Provide only one of `lt` or deprecated `mortality_table`.", call. = FALSE)
    }
    lt <- dot_get("mortality_table")
  }

  if (dot_has("age")) {
    if (!missing(x)) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }
    x <- dot_get("age")
  }

  if (dot_has("rate")) {
    if (!missing(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }
    i <- dot_get("rate")
  }

  if (dot_has("rate_type")) {
    if (!identical(i_type, "effective")) {
      stop("Provide only one of `i_type` or deprecated `rate_type`.", call. = FALSE)
    }
    i_type <- dot_get("rate_type")
  }

  if (dot_has("insurance_type")) {
    type <- dot_get("insurance_type")
  }

  if (dot_has("term_years")) {
    if (!is.infinite(n)) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }
    n <- dot_get("term_years")
  }

  if (dot_has("premium")) {
    if (!is.null(P)) {
      stop("Provide only one of `P` or deprecated `premium`.", call. = FALSE)
    }
    P <- dot_get("premium")
  }

  if (dot_has("premium_term_years")) {
    if (!is.null(n_prem)) {
      stop("Provide only one of `n_prem` or deprecated `premium_term_years`.", call. = FALSE)
    }
    n_prem <- dot_get("premium_term_years")
  }

  if (dot_has("durations")) {
    if (!is.null(t)) {
      stop("Provide only one of `t` or deprecated `durations`.", call. = FALSE)
    }
    t <- dot_get("durations")
  }

  if (dot_has("output")) {
    if (!identical(tidy, TRUE)) {
      stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
    }

    output <- match.arg(dot_get("output"), c("table", "value"))
    tidy <- identical(output, "table")
  }

  type <- match.arg(type)
  method <- match.arg(method)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (exists(".as_life_contract", mode = "function") && .as_life_contract(lt)) {
    contract <- lt

    if (!identical(contract$lives, "single")) {
      stop(
        "`reserve_x()` currently supports only single-life `life_contract()` objects.",
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

    if (is.null(i_type) || identical(i_type, "effective")) {
      i_type <- contract$rate_type %||% contract$i_type %||% i_type
    }

    if (is.null(m)) {
      m <- contract$m
    }
  }

  # local infix helper without importing anything
  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
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

  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || abs(x - round(x)) > 1e-10) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(i) || length(i) != 1L || is.na(i) ||
      !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
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

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m < 1 || abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) || length(benefit) != 1L || is.na(benefit) ||
      !is.finite(benefit) || benefit <= 0) {
    stop("`benefit` must be a single positive finite number.", call. = FALSE)
  }

  if (!is.null(P) &&
      (!is.numeric(P) || length(P) != 1L || is.na(P) ||
       !is.finite(P))) {
    stop("`P` must be NULL or a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n <= 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) || abs(n - round(n)) > 1e-10))) {
    stop("`n` must be `Inf` or a single positive integer.", call. = FALSE)
  }

  if (type %in% c("term", "endowment") && is.infinite(n)) {
    stop("`n` must be finite for term and endowment insurance.", call. = FALSE)
  }

  if (!is.null(n_prem) &&
      (!is.numeric(n_prem) || length(n_prem) != 1L ||
       is.na(n_prem) || n_prem <= 0 ||
       (!is.infinite(n_prem) &&
        (!is.finite(n_prem) ||
         abs(n_prem - round(n_prem)) > 1e-10)))) {
    stop(
      "`n_prem` must be NULL, Inf, or a single positive integer.",
      call. = FALSE
    )
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (!is.null(n_prem) && !is.infinite(n_prem)) {
    n_prem <- as.integer(round(n_prem))
  }

  # -------------------------------------------------------------------------
  # Interest and life table preparation
  # -------------------------------------------------------------------------

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (!is.numeric(i_effective) || length(i_effective) != 1L ||
      is.na(i_effective) || !is.finite(i_effective) || i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  v <- 1 / (1 + i_effective)

  lt <- lt[order(lt$x), , drop = FALSE]

  if (!is.numeric(lt$x) || !is.numeric(lt$lx)) {
    stop("Columns `x` and `lx` in `lt` must be numeric.", call. = FALSE)
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

  survival <- function(current_age, years) {
    if (years == 0L) {
      return(1)
    }

    l0 <- get_lx(current_age)
    l1 <- get_lx(current_age + years)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    l1 / l0
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

  contract_horizon <- if (type == "whole") {
    omega - x
  } else {
    n
  }

  if (!is.finite(contract_horizon) || contract_horizon < 1) {
    stop(
      "The life table does not provide a positive contract horizon from `x`.",
      call. = FALSE
    )
  }

  if (x + contract_horizon > omega) {
    stop(
      "The requested contract horizon exceeds the available life table ages.",
      call. = FALSE
    )
  }

  contract_horizon <- as.integer(round(contract_horizon))

  if (is.null(n_prem)) {
    premium_horizon <- contract_horizon
  } else if (is.infinite(n_prem)) {
    premium_horizon <- contract_horizon
  } else {
    premium_horizon <- min(n_prem, contract_horizon)
  }

  premium_horizon <- as.integer(round(premium_horizon))

  # -------------------------------------------------------------------------
  # Internal APV helpers
  # -------------------------------------------------------------------------

  apv_benefits_from <- function(current_age, remaining_contract, current_type) {
    remaining_contract <- as.integer(round(remaining_contract))

    if (remaining_contract <= 0L) {
      return(0)
    }

    death_apv <- sum(vapply(0:(remaining_contract - 1L), function(r) {
      px_r <- survival(current_age, r)
      qx_r <- one_year_qx(current_age + r)

      if (is.na(px_r) || is.na(qx_r)) {
        stop("The life table does not support the requested benefit APV.", call. = FALSE)
      }

      benefit * v^(r + 1L) * px_r * qx_r
    }, numeric(1L)))

    if (current_type == "endowment") {
      pure_endowment <- survival(current_age, remaining_contract)

      if (is.na(pure_endowment)) {
        stop("The life table does not support the requested endowment APV.", call. = FALSE)
      }

      death_apv + benefit * v^remaining_contract * pure_endowment
    } else {
      death_apv
    }
  }

  apv_premiums_from <- function(current_age, remaining_premiums) {
    remaining_premiums <- as.integer(round(remaining_premiums))

    if (remaining_premiums <= 0L) {
      return(0)
    }

    sum(vapply(0:(remaining_premiums - 1L), function(r) {
      px_r <- survival(current_age, r)

      if (is.na(px_r)) {
        stop("The life table does not support the requested premium APV.", call. = FALSE)
      }

      v^r * px_r
    }, numeric(1L)))
  }

  # -------------------------------------------------------------------------
  # Premium
  # -------------------------------------------------------------------------

  premium_was_computed <- is.null(P)

  if (is.null(P)) {
    apv_b0 <- apv_benefits_from(
      current_age = x,
      remaining_contract = contract_horizon,
      current_type = type
    )

    apv_p0 <- apv_premiums_from(
      current_age = x,
      remaining_premiums = premium_horizon
    )

    if (!is.finite(apv_p0) || apv_p0 <= 0) {
      stop("The premium annuity APV is zero; the net premium is undefined.",
           call. = FALSE)
    }

    P <- apv_b0 / apv_p0
  }

  # -------------------------------------------------------------------------
  # Durations
  # -------------------------------------------------------------------------

  if (is.null(t)) {
    t_vec <- 0:contract_horizon
  } else {
    if (!is.numeric(t) || any(is.na(t)) ||
        any(!is.finite(t)) ||
        any(abs(t - round(t)) > 1e-10)) {
      stop("`t` must be an integer vector of policy durations.", call. = FALSE)
    }

    t_vec <- sort(unique(as.integer(round(t))))

    if (any(t_vec < 0) || any(t_vec > contract_horizon)) {
      stop(
        "`t` must lie between 0 and the contract horizon.",
        call. = FALSE
      )
    }
  }

  # -------------------------------------------------------------------------
  # Prospective method
  # -------------------------------------------------------------------------

  if (method == "prospective") {
    reserves <- vapply(t_vec, function(tt) {
      if (tt >= contract_horizon) {
        if (type == "endowment" && tt == contract_horizon) {
          return(benefit)
        }

        return(0)
      }

      remaining_contract <- contract_horizon - tt
      current_age <- x + tt

      apv_benefits <- apv_benefits_from(
        current_age = current_age,
        remaining_contract = remaining_contract,
        current_type = type
      )

      remaining_premiums <- max(0L, premium_horizon - tt)

      apv_premiums <- if (remaining_premiums == 0L) {
        0
      } else {
        P * apv_premiums_from(
          current_age = current_age,
          remaining_premiums = remaining_premiums
        )
      }

      apv_benefits - apv_premiums
    }, numeric(1L))

    if (premium_was_computed && 0L %in% t_vec) {
      reserves[t_vec == 0L] <- 0
    }
  } else {
    # -----------------------------------------------------------------------
    # Recursive method
    # -----------------------------------------------------------------------

    full_reserves <- numeric(contract_horizon + 1L)
    full_reserves[[1L]] <- 0

    for (kk in 0:(contract_horizon - 1L)) {
      premium_paid <- if (kk < premium_horizon) P else 0
      qx_k <- one_year_qx(x + kk)

      if (is.na(qx_k)) {
        stop("The life table does not support recursive reserve computation.",
             call. = FALSE)
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

    if (type == "endowment") {
      full_reserves[[contract_horizon + 1L]] <- benefit
    }

    if (type == "term") {
      full_reserves[[contract_horizon + 1L]] <- 0
    }

    reserves <- full_reserves[t_vec + 1L]
  }

  if (!tidy) {
    names(reserves) <- paste0("t=", t_vec)
    return(reserves)
  }

  P_paid <- vapply(t_vec, function(tt) {
    if (tt >= contract_horizon) {
      return(0)
    }

    if (tt < premium_horizon) P else 0
  }, numeric(1L))

  benefit_due <- vapply(t_vec, function(tt) {
    if (tt >= contract_horizon) {
      return(0)
    }

    benefit
  }, numeric(1L))

  n_out <- if (identical(type, "whole")) {
    Inf
  } else {
    n
  }

  tibble::tibble(
    t = t_vec,
    duration = t_vec,
    x_t = x + t_vec,
    age = x + t_vec,
    V = reserves,
    reserve = reserves,
    P_paid = P_paid,
    premium_paid = P_paid,
    benefit_due = benefit_due,
    type = type,
    insurance_type = type,
    benefit = benefit,
    P = P,
    premium = P,
    n_prem = premium_horizon,
    premium_term_years = premium_horizon,
    n = n_out,
    term_years = n_out,
    contract_horizon = contract_horizon,
    method = method,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective
  )
}
