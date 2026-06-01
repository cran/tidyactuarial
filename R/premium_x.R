#' Net premium for single-life insurance by the equivalence principle
#'
#' Computes the net benefit premium of a single-life insurance contract using
#' the equivalence principle:
#' \deqn{
#' P = \frac{APV(\text{benefits})}{APV(\text{premium annuity})}.
#' }
#'
#' The premium returned corresponds to one premium payment. For example, when
#' \code{k = 12}, the returned value is the monthly premium.
#'
#' @param lt A life table data frame containing at least columns \code{x} and
#'   \code{lx}.
#' @param x Integer actuarial age at issue.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal interest-rate
#'   inputs. Ignored for \code{i_type = "effective"} and
#'   \code{i_type = "force"}.
#' @param type Character string. Type of insurance contract. One of
#'   \code{"whole"}, \code{"term"}, \code{"endowment"}, or
#'   \code{"variable_k"}.
#' @param benefit Benefit amount. For standard products, a single nonnegative
#'   numeric value. For \code{type = "variable_k"}, a numeric vector or a
#'   function of time may be supplied and is passed to
#'   \code{\link{insurance_variable_k}}.
#' @param n Insurance term in years. Use \code{Inf} for whole-life insurance.
#'   Required as a finite value for term and endowment insurance. For
#'   \code{type = "variable_k"}, \code{n = Inf} allows the term to be inferred
#'   from a numeric \code{benefit} vector.
#' @param h Integer deferment period in years.
#' @param k Positive integer. Number of premium payments per year.
#' @param frac Fractional-age assumption used only for
#'   \code{type = "variable_k"}. One of \code{"UDD"}, \code{"CF"},
#'   \code{"CML"}, or \code{"Balducci"}.
#' @param timing Timing of premium payments. Use \code{"due"} for payments in
#'   advance or \code{"immediate"} for payments in arrears.
#' @param premium_start Start of premium payments. Use \code{"issue"} for
#'   premiums starting at issue, or \code{"deferred"} for premiums starting
#'   after \code{h}.
#' @param n_prem Optional premium-paying term in years, counted from
#'   \code{premium_start}. If \code{NULL}, it defaults to whole life for
#'   whole-life insurance and to \code{n} for temporary products. For finite
#'   term contracts, premiums must not extend beyond the end of coverage.
#' @param woolhouse Woolhouse order for the premium annuity when \code{k > 1}.
#'   One of \code{"none"}, \code{"first"}, or \code{"second"}.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric premium. If
#'   \code{TRUE}, returns a one-row tibble with details.
#' @param check Logical. If \code{TRUE}, performs input validation.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{insurance_type}, \code{term_years}, \code{deferral_years},
#'   \code{payments_per_year}, \code{premium_timing},
#'   \code{premium_term_years}, and \code{output}.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric scalar with the net premium per payment.
#'
#' If \code{tidy = TRUE}, a one-row tibble with the main inputs, APV of
#' benefits, APV of premiums, premium per payment, and annualized premium.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age at
#' issue, \code{i} is the interest-rate input, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' \code{n} is the insurance term, \code{h} is the deferment period, and
#' \code{k} is the premium payment frequency.
#'
#' The benefit premium is the level payment satisfying the equivalence
#' principle: the APV of premiums equals the APV of benefits at issue.
#'
#' For standard products, the APV of benefits is computed with
#' \code{\link{insurance_x}}. For \code{type = "variable_k"}, it is computed
#' with \code{\link{insurance_variable_k}}. The APV of the premium annuity is
#' computed with \code{\link{annuity_x}}, supporting k-thly payments and
#' Woolhouse approximations.
#'
#' @seealso \code{\link{insurance_x}}, \code{\link{insurance_variable_k}},
#'   \code{\link{annuity_x}}, \code{\link{premium_xy}},
#'   \code{\link{premium_gross}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Whole-life insurance, annual premium
#' premium_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   type = "whole",
#'   benefit = 100000
#' )
#'
#' # Verify manually: P = A / a-double-dot
#' A <- insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   type = "whole",
#'   benefit = 100000
#' )
#'
#' ad <- annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   timing = "due"
#' )
#'
#' A / ad
#'
#' # Five-year term insurance
#' premium_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   type = "term",
#'   n = 5,
#'   benefit = 100000
#' )
#'
#' # Tidy output
#' premium_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   type = "term",
#'   n = 5,
#'   benefit = 100000,
#'   tidy = TRUE
#' )
#'
#' # Monthly premiums paid for a shorter period than the coverage term
#' premium_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   type = "term",
#'   n = 5,
#'   benefit = 100000,
#'   k = 12,
#'   n_prem = 3
#' )
#'
#' @export
premium_x <- function(
    lt,
    x,
    i,
    i_type = "effective",
    m = 1L,
    type = c("whole", "term", "endowment", "variable_k"),
    benefit = 1,
    n = Inf,
    h = 0L,
    k = 1L,
    frac = c("UDD", "CF", "CML", "Balducci"),
    timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    n_prem = NULL,
    woolhouse = c("none", "first", "second"),
    tidy = FALSE,
    check = TRUE,
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
    "insurance_type",
    "term_years",
    "deferral_years",
    "payments_per_year",
    "premium_timing",
    "premium_term_years",
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

  if (!is.null(dots$insurance_type)) {
    type <- dots$insurance_type
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

  if (!is.null(dots$payments_per_year)) {
    if (!identical(k, 1L) && !identical(k, 1)) {
      stop("Provide only one of `k` or deprecated `payments_per_year`.", call. = FALSE)
    }
    k <- dots$payments_per_year
  }

  if (!is.null(dots$premium_timing)) {
    timing <- dots$premium_timing
  }

  if (!is.null(dots$premium_term_years)) {
    if (!is.null(n_prem)) {
      stop("Provide only one of `n_prem` or deprecated `premium_term_years`.", call. = FALSE)
    }
    n_prem <- dots$premium_term_years
  }

  if (!is.null(dots$output)) {
    if (!identical(tidy, FALSE)) {
      stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
    }

    output <- match.arg(dots$output, c("value", "table"))
    tidy <- identical(output, "table")
  }

  type <- match.arg(type)
  frac <- match.arg(frac)
  timing <- match.arg(timing)
  premium_start <- match.arg(premium_start)
  woolhouse <- match.arg(woolhouse)

  if (frac == "CML") {
    frac <- "CF"
  }

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
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
        "`premium_x()` currently supports only single-life `life_contract()` objects.",
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

    if (identical(i_type, "effective")) {
      i_type <- contract$rate_type %||% contract$i_type %||% i_type
    }

    if (identical(m, 1L) || identical(m, 1)) {
      m <- contract$m %||% m
    }
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (isTRUE(check)) {
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

    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        abs(x - round(x)) > 1e-10) {
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
        abs(h - round(h)) > 1e-10) {
      stop("`h` must be a single nonnegative integer.", call. = FALSE)
    }

    if (!is.numeric(k) ||
        length(k) != 1L ||
        is.na(k) ||
        !is.finite(k) ||
        k < 1 ||
        abs(k - round(k)) > 1e-10) {
      stop("`k` must be a single positive integer.", call. = FALSE)
    }

    if (!is.numeric(n) ||
        length(n) != 1L ||
        is.na(n) ||
        n <= 0 ||
        (!is.infinite(n) &&
         (!is.finite(n) ||
          abs(n - round(n)) > 1e-10))) {
      stop("`n` must be `Inf` or a single positive integer.", call. = FALSE)
    }

    if (type %in% c("term", "endowment") && is.infinite(n)) {
      stop(
        "`n` must be finite for term and endowment insurance.",
        call. = FALSE
      )
    }

    if (type == "variable_k" && is.infinite(n) && is.function(benefit)) {
      stop(
        "For `type = 'variable_k'` with a functional benefit, `n` must be finite.",
        call. = FALSE
      )
    }

    if (!is.null(n_prem) &&
        (!is.numeric(n_prem) ||
         length(n_prem) != 1L ||
         is.na(n_prem) ||
         n_prem <= 0 ||
         (!is.infinite(n_prem) &&
          (!is.finite(n_prem) ||
           abs(n_prem * k - round(n_prem * k)) > 1e-10)))) {
      stop(
        "`n_prem` must be NULL, Inf, or a positive value satisfying `n_prem * k` integer.",
        call. = FALSE
      )
    }

    if (type != "variable_k") {
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

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (!is.null(n_prem) && !is.infinite(n_prem)) {
    # n_prem may be fractional when k > 1, but must lie on the k-thly grid.
    n_prem <- round(n_prem * k) / k
  }

  # -------------------------------------------------------------------------
  # 1. APV of benefits, valued at issue
  # -------------------------------------------------------------------------

  if (type == "variable_k") {
    n_variable <- if (is.infinite(n)) NULL else n

    apv_benefits <- insurance_variable_k(
      lt = lt,
      x = x,
      i = i,
      i_type = i_type,
      m = m,
      benefit = benefit,
      n = n_variable,
      h = h,
      k = k,
      frac = frac,
      tidy = FALSE,
      check = check
    )

    if (is.null(n_prem)) {
      if (!is.infinite(n)) {
        n_prem <- n
      } else if (!is.function(benefit)) {
        n_prem <- length(benefit) / k
      } else {
        stop(
          "Cannot infer `n_prem` for variable benefits supplied as a function.",
          call. = FALSE
        )
      }
    }
  } else {
    apv_benefits <- insurance_x(
      lt = lt,
      x = x,
      i = i,
      i_type = i_type,
      m = m,
      n = n,
      h = h,
      type = type,
      benefit = benefit,
      tidy = FALSE
    )

    if (is.null(n_prem)) {
      n_prem <- if (type == "whole") Inf else n
    }
  }

  # -------------------------------------------------------------------------
  # 1b. Premium-paying horizon validation
  # -------------------------------------------------------------------------

  if (type %in% c("term", "endowment", "variable_k") &&
      !is.infinite(n)) {
    if (is.infinite(n_prem)) {
      stop(
        "`n_prem` cannot be `Inf` for a finite-term insurance contract.",
        call. = FALSE
      )
    }

    premium_start_time <- if (premium_start == "issue") {
      0
    } else {
      h
    }

    premium_end_time <- premium_start_time + n_prem
    coverage_end_time <- h + n

    if (premium_end_time > coverage_end_time + 1e-10) {
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

  h_prem <- if (premium_start == "issue") {
    0
  } else {
    h
  }

  apv_premiums <- annuity_x(
    lt = lt,
    x = x,
    i = i,
    i_type = i_type,
    m = m,
    n = n_prem,
    h = h_prem,
    k = k,
    timing = timing,
    woolhouse = woolhouse,
    tidy = FALSE
  )

  if (!is.finite(apv_premiums) || apv_premiums <= 0) {
    stop("APV of premium annuity is nonpositive or not finite.", call. = FALSE)
  }

  P <- apv_benefits / apv_premiums

  if (!tidy) {
    return(P)
  }

  benefit_out <- if (is.function(benefit)) {
    NA_real_
  } else {
    suppressWarnings(as.numeric(benefit)[1])
  }

  tibble::tibble(
    x = x,
    age = x,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    type = type,
    insurance_type = type,
    benefit = benefit_out,
    n = n,
    term_years = n,
    h = h,
    deferral_years = h,
    k = k,
    payments_per_year = k,
    frac = frac,
    timing = timing,
    premium_timing = timing,
    premium_start = premium_start,
    n_prem = n_prem,
    premium_term_years = n_prem,
    woolhouse = woolhouse,
    P = P,
    premium = P,
    P_annual = k * P,
    premium_annual = k * P,
    apv_benefits = apv_benefits,
    a_premiums = apv_premiums,
    apv_premiums = apv_premiums
  )
}
