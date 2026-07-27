#' Net premium for single-life insurance by the equivalence principle
#'
#' Computes the net benefit premium of a single-life insurance contract from
#' the equation of equivalence:
#' \deqn{
#' \operatorname{APV}(\text{premiums})
#' =
#' \operatorname{APV}(\text{benefits}).
#' }
#'
#' For premiums payable \eqn{k} times per year, the function distinguishes
#' between the annualized premium \eqn{P^{(k)}} and the amount paid at each
#' installment, \eqn{P^{(k)} / k}.
#'
#' @param lt A life table containing at least columns \code{x} and \code{lx},
#'   or a single-life contract created with \code{\link{life_contract}}.
#' @param x Integer actuarial age at issue. Optional when \code{lt} is a
#'   single-life contract.
#' @param i Numeric scalar. Annual interest-rate input. Optional when
#'   \code{lt} is a single-life contract.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal interest-rate
#'   inputs.
#' @param type Insurance type. One of \code{"whole"}, \code{"term"},
#'   \code{"endowment"}, or \code{"variable_k"}.
#' @param benefit Benefit amount. For standard products, a single nonnegative
#'   numeric value. For \code{type = "variable_k"}, a numeric vector or a
#'   function of time may be supplied.
#' @param n Insurance term in years. Use \code{Inf} for whole-life insurance.
#' @param h Nonnegative integer deferment period in years.
#' @param k Positive integer. Number of premium payments per year.
#' @param frac Fractional-age assumption. One of \code{"UDD"}, \code{"CF"},
#'   \code{"CML"}, or \code{"Balducci"}.
#' @param timing Premium timing: \code{"due"} or \code{"immediate"}.
#' @param premium_start Start of premium payments: \code{"issue"} or
#'   \code{"deferred"}.
#' @param n_prem Premium-paying term in years. If finite, \code{n_prem * k}
#'   must be an integer.
#' @param woolhouse Woolhouse approximation order for the premium annuity:
#'   \code{"none"}, \code{"first"}, or \code{"second"}.
#' @param output Output level. \code{"value"} returns the annualized premium
#'   as a numeric scalar; \code{"summary"} returns a compact one-row executive
#'   summary; \code{"audit"} returns the components of the equivalence
#'   calculation in long format.
#' @param tidy Deprecated compatibility argument. \code{TRUE} maps to
#'   \code{output = "summary"} and \code{FALSE} maps to
#'   \code{output = "value"} when \code{output} is not supplied.
#' @param check Logical scalar. If \code{TRUE}, performs input validation.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{insurance_type}, \code{term_years}, \code{deferral_years},
#'   \code{payments_per_year}, \code{premium_timing}, and
#'   \code{premium_term_years}.
#'
#' @return
#' For \code{output = "value"}, a numeric scalar containing the annualized
#' premium \eqn{P^{(k)}}.
#'
#' For \code{output = "summary"}, a one-row tibble with six columns:
#' annualized premium, premium per payment, payment frequency, APV of
#' benefits, APV of the premium annuity, and equivalence residual.
#'
#' For \code{output = "audit"}, a long-format tibble containing the principal
#' inputs and calculated components.
#'
#' @details
#' Let \eqn{Z} be the present-value random variable of benefits and let
#' \eqn{Y^{(k)}} be the present-value random variable of a premium annuity
#' normalized to an annual payment rate of 1. The loss at issue is
#' \deqn{
#' L_0 = Z - P^{(k)}Y^{(k)}.
#' }
#'
#' The equivalence principle, \eqn{\operatorname{E}[L_0] = 0}, gives
#' \deqn{
#' P^{(k)}
#' =
#' \frac{\operatorname{APV}(\text{benefits})}
#'      {\operatorname{APV}(\text{premium annuity})}.
#' }
#'
#' Because \code{\link{annuity_x}} values a k-thly annuity with payments of
#' \eqn{1/k}, this quotient is the annualized premium. The actual installment
#' paid each fraction of the year is
#' \deqn{
#' P_{\text{per payment}} = \frac{P^{(k)}}{k}.
#' }
#'
#' A contract assembled with pipes may be valued directly:
#'
#' \preformatted{
#' life_contract(...) |>
#'   add_insurance(...) |>
#'   add_premium_schedule(...) |>
#'   premium_x(output = "summary")
#' }
#'
#' Explicit arguments supplied to \code{premium_x()} override values stored in
#' the contract components.
#'
#' @seealso \code{\link{life_contract}}, \code{\link{add_insurance}},
#'   \code{\link{add_premium_schedule}}, \code{\link{insurance_x}},
#'   \code{\link{insurance_variable_k}}, \code{\link{annuity_x}},
#'   \code{\link{premium_xy}}, \code{\link{premium_gross}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 40:90,
#'   lx = round(100000 * exp(-0.018 * (0:50)^1.35))
#' )
#' lt$lx[nrow(lt)] <- 0
#'
#' # Direct calculation: annualized premium
#' premium_x(
#'   lt = lt,
#'   x = 40,
#'   i = 0.05,
#'   type = "term",
#'   benefit = 100000,
#'   n = 20,
#'   k = 12,
#'   n_prem = 10
#' )
#'
#' # Executive result with annualized and monthly premiums
#' premium_x(
#'   lt = lt,
#'   x = 40,
#'   i = 0.05,
#'   type = "term",
#'   benefit = 100000,
#'   n = 20,
#'   k = 12,
#'   n_prem = 10,
#'   output = "summary"
#' )
#'
#' # Contract assembled with pipes
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
#'     n_prem = 10,
#'     timing = "due"
#'   ) |>
#'   premium_x(output = "summary")
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
    output = c("value", "summary", "audit"),
    tidy = NULL,
    check = TRUE,
    ...
) {
  lt_missing <- missing(lt)
  x_missing <- missing(x)
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
  woolhouse_missing <- missing(woolhouse)
  output_missing <- missing(output)

  dots <- list(...)

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
    "premium_term_years"
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
    if (!lt_missing) {
      stop(
        "Provide only one of `lt` or deprecated `mortality_table`.",
        call. = FALSE
      )
    }

    lt <- dots$mortality_table
    lt_missing <- FALSE
  }

  if (!is.null(dots$age)) {
    if (!x_missing) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }

    x <- dots$age
    x_missing <- FALSE
  }

  if (!is.null(dots$rate)) {
    if (!i_missing) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dots$rate
    i_missing <- FALSE
  }

  if (!is.null(dots$rate_type)) {
    if (!i_type_missing) {
      stop(
        "Provide only one of `i_type` or deprecated `rate_type`.",
        call. = FALSE
      )
    }

    i_type <- dots$rate_type
    i_type_missing <- FALSE
  }

  if (!is.null(dots$insurance_type)) {
    if (!type_missing) {
      stop(
        "Provide only one of `type` or deprecated `insurance_type`.",
        call. = FALSE
      )
    }

    type <- dots$insurance_type
    type_missing <- FALSE
  }

  if (!is.null(dots$term_years)) {
    if (!n_missing) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }

    n <- dots$term_years
    n_missing <- FALSE
  }

  if (!is.null(dots$deferral_years)) {
    if (!h_missing) {
      stop(
        "Provide only one of `h` or deprecated `deferral_years`.",
        call. = FALSE
      )
    }

    h <- dots$deferral_years
    h_missing <- FALSE
  }

  if (!is.null(dots$payments_per_year)) {
    if (!k_missing) {
      stop(
        "Provide only one of `k` or deprecated `payments_per_year`.",
        call. = FALSE
      )
    }

    k <- dots$payments_per_year
    k_missing <- FALSE
  }

  if (!is.null(dots$premium_timing)) {
    if (!timing_missing) {
      stop(
        "Provide only one of `timing` or deprecated `premium_timing`.",
        call. = FALSE
      )
    }

    timing <- dots$premium_timing
    timing_missing <- FALSE
  }

  if (!is.null(dots$premium_term_years)) {
    if (!n_prem_missing) {
      stop(
        "Provide only one of `n_prem` or deprecated `premium_term_years`.",
        call. = FALSE
      )
    }

    n_prem <- dots$premium_term_years
    n_prem_missing <- FALSE
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

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
  }

  contract <- NULL
  insurance_spec <- NULL
  premium_spec <- NULL

  if (!lt_missing &&
      exists(".as_life_contract", mode = "function") &&
      .as_life_contract(lt)) {
    contract <- lt

    if (!identical(contract$lives, "single")) {
      stop(
        "`premium_x()` currently supports only single-life contracts.",
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

    if (!is.null(insurance_spec)) {
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

  k_benefit <- if (!is.null(insurance_spec$k)) {
    insurance_spec$k
  } else {
    k
  }

  if (isTRUE(check)) {
    if (lt_missing) {
      stop("`lt` must be provided.", call. = FALSE)
    }

    if (x_missing) {
      stop("`x` must be provided.", call. = FALSE)
    }

    if (i_missing) {
      stop("`i` must be provided.", call. = FALSE)
    }

    if (!is.data.frame(lt)) {
      stop("`lt` must be a data frame or tibble.", call. = FALSE)
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

    if (!is.numeric(k_benefit) ||
        length(k_benefit) != 1L ||
        is.na(k_benefit) ||
        !is.finite(k_benefit) ||
        k_benefit < 1 ||
        abs(k_benefit - round(k_benefit)) > 1e-10) {
      stop("The benefit frequency must be a single positive integer.", call. = FALSE)
    }

    if (!is.numeric(n) ||
        length(n) != 1L ||
        is.na(n) ||
        n <= 0 ||
        (!is.infinite(n) &&
         (!is.finite(n) || abs(n - round(n)) > 1e-10))) {
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
        "`n_prem` must be NULL, Inf, or a positive value satisfying ",
        "`n_prem * k` integer.",
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
          "For standard products, `benefit` must be a single finite ",
          "nonnegative number.",
          call. = FALSE
        )
      }
    }
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))
  k_benefit <- as.integer(round(k_benefit))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (!is.null(n_prem) && !is.infinite(n_prem)) {
    n_prem <- round(n_prem * k) / k
  }

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
      k = k_benefit,
      frac = frac_benefit,
      tidy = FALSE,
      check = check
    )

    if (is.null(n_prem)) {
      if (!is.infinite(n)) {
        n_prem <- n
      } else if (!is.function(benefit)) {
        n_prem <- length(benefit) / k_benefit
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

  if (type %in% c("term", "endowment", "variable_k") &&
      !is.infinite(n)) {
    if (is.infinite(n_prem)) {
      stop(
        "`n_prem` cannot be `Inf` for a finite-term insurance contract.",
        call. = FALSE
      )
    }

    premium_start_time <- if (premium_start == "issue") 0 else h
    premium_end_time <- premium_start_time + n_prem
    coverage_end_time <- h + n

    if (premium_end_time > coverage_end_time + 1e-10) {
      stop(
        "Premium payments must not extend beyond the end of coverage. ",
        "Premiums end at time ", premium_end_time,
        ", while coverage ends at time ", coverage_end_time, ".",
        call. = FALSE
      )
    }
  }

  h_prem <- if (premium_start == "issue") 0 else h

  apv_premium_annuity <- annuity_x(
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
    frac = frac_premium,
    tidy = FALSE
  )

  if (!is.finite(apv_premium_annuity) || apv_premium_annuity <= 0) {
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

  insurance_term_value <- if (is.infinite(n)) Inf else as.numeric(n)
  premium_term_value <- if (is.infinite(n_prem)) Inf else as.numeric(n_prem)

  tibble::tibble(
    component = c(
      "apv_benefits",
      "apv_premium_annuity",
      "premium_annualized",
      "premium_per_payment",
      "payments_per_year",
      "insurance_term_years",
      "premium_term_years",
      "equivalence_residual"
    ),
    value = c(
      apv_benefits,
      apv_premium_annuity,
      premium_annualized,
      premium_per_payment,
      k,
      insurance_term_value,
      premium_term_value,
      equivalence_residual
    ),
    unit = c(
      "currency",
      "annuity factor",
      "currency per year",
      "currency per payment",
      "payments per year",
      "years",
      "years",
      "currency"
    )
  )
}
