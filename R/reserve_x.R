#' Benefit reserves for single-life insurance
#'
#' Computes terminal benefit reserves at selected integer policy durations by
#' the prospective or recursive method.
#'
#' The prospective method supports annual and true k-thly premiums payable in
#' advance. The premium input \code{P} is always interpreted as the annualized
#' premium \eqn{P^{(k)}}. The amount paid at each premium date is
#' \eqn{P^{(k)} / k}.
#'
#' @param lt A life table containing columns \code{x} and \code{lx}, or a
#'   single-life contract created with \code{\link{life_contract}}.
#' @param x Integer actuarial age at issue. Optional for a life contract.
#' @param i Numeric scalar. Annual interest-rate input. Optional for a life
#'   contract.
#' @param i_type Interest-rate type: \code{"effective"},
#'   \code{"nominal_interest"}, \code{"nominal_discount"}, or \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal interest rates.
#' @param type Insurance type: \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param n Insurance term in years. Use \code{Inf} for whole-life insurance.
#' @param h Nonnegative integer deferment period in years.
#' @param benefit Positive insurance benefit.
#' @param P Optional annualized premium \eqn{P^{(k)}}. If \code{NULL}, it is
#'   calculated with \code{\link{premium_x}} by the equivalence principle.
#' @param k Positive integer. Number of premium payments per year.
#' @param frac Fractional-age assumption for exact k-thly premium annuities.
#' @param timing Premium timing. Currently only \code{"due"} is supported,
#'   which gives reserves immediately before any premium payable at duration
#'   \code{t}.
#' @param premium_start Start of premium payments: \code{"issue"} or
#'   \code{"deferred"}.
#' @param n_prem Premium-paying term in years. A fractional term is permitted
#'   when \code{n_prem * k} is an integer.
#' @param woolhouse Woolhouse approximation order for the premium annuity:
#'   \code{"none"}, \code{"first"}, or \code{"second"}.
#' @param t Optional integer vector of policy durations. If \code{NULL}, all
#'   supported integer durations are returned.
#' @param method Reserve method: \code{"prospective"} or \code{"recursive"}.
#'   The recursive method currently supports annual premiums only
#'   (\code{k = 1}).
#' @param output Output level. \code{"summary"} returns a compact reserve
#'   schedule, \code{"value"} returns a named numeric vector, and
#'   \code{"audit"} returns the reserve components in long format.
#'   \code{"table"} is accepted as a deprecated alias for \code{"summary"}.
#' @param tidy Deprecated compatibility argument. \code{TRUE} maps to
#'   \code{output = "summary"} and \code{FALSE} to
#'   \code{output = "value"}.
#' @param check Logical scalar. If \code{TRUE}, validates inputs.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{insurance_type}, \code{term_years}, \code{premium},
#'   \code{premium_term_years}, and \code{durations}.
#'
#' @return
#' For \code{output = "summary"}, a tibble with at most six columns:
#' duration, attained age, reserve, annualized premium, premium per payment,
#' and method.
#'
#' For \code{output = "value"}, a named numeric vector.
#'
#' For \code{output = "audit"}, a long-format tibble containing future benefit
#' APV, future premium APV, reserve, and premium amounts at each duration.
#'
#' @details
#' At integer duration \eqn{t}, conditional on survival to age \eqn{x+t}, the
#' prospective reserve is
#' \deqn{
#' {}_tV_x =
#' \operatorname{APV}_t(\text{future benefits})
#' -
#' P^{(k)}
#' \operatorname{APV}_t(\text{future premium annuity}).
#' }
#'
#' The premium annuity is normalized to an annual payment rate of 1, so each
#' installment has amount \eqn{1/k}. Therefore, the multiplier in the reserve
#' formula is the annualized premium \eqn{P^{(k)}}, not the installment
#' \eqn{P^{(k)}/k}.
#'
#' Reserves are measured immediately before any premium payable at duration
#' \eqn{t}. This is the standard fully discrete terminal-reserve convention
#' for premiums payable in advance.
#'
#' For annual premiums, the recursive method uses
#' \deqn{
#' {}_{t+1}V_x =
#' \frac{
#' ({}_tV_x + P_t)(1+i) - b_{t+1}q_{x+t}
#' }{p_{x+t}}.
#' }
#'
#' @seealso \code{\link{premium_x}}, \code{\link{insurance_x}},
#'   \code{\link{annuity_x}}, \code{\link{life_contract}},
#'   \code{\link{add_insurance}}, \code{\link{add_premium_schedule}}
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
#' reserve_x(
#'   lt = lt,
#'   x = 40,
#'   i = 0.05,
#'   type = "term",
#'   n = 20,
#'   benefit = 100000,
#'   k = 12,
#'   n_prem = 10,
#'   t = c(0, 5, 10, 15, 20)
#' )
#'
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
#'     n_prem = 10
#'   ) |>
#'   reserve_x(
#'     t = c(0, 5, 10, 15, 20),
#'     output = "summary"
#'   )
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
    h = 0L,
    benefit = 1,
    P = NULL,
    k = 1L,
    frac = c("UDD", "CF", "CML", "Balducci"),
    timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    n_prem = NULL,
    woolhouse = c("none", "first", "second"),
    t = NULL,
    method = c("prospective", "recursive"),
    output = c("summary", "value", "audit", "table"),
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
  n_missing <- missing(n)
  h_missing <- missing(h)
  benefit_missing <- missing(benefit)
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
    "premium",
    "premium_term_years",
    "durations"
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

  if ("age" %in% names(dots)) {
    if (!x_missing) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }

    x <- dots[["age"]]
    x_missing <- FALSE
  }

  if ("rate" %in% names(dots)) {
    if (!i_missing) {
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

  if ("premium" %in% names(dots)) {
    if (!is.null(P)) {
      stop("Provide only one of `P` or deprecated `premium`.", call. = FALSE)
    }

    P <- dots[["premium"]]
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

  if ("durations" %in% names(dots)) {
    if (!is.null(t)) {
      stop(
        "Provide only one of `t` or deprecated `durations`.",
        call. = FALSE
      )
    }

    t <- dots[["durations"]]
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

  insurance_spec <- NULL
  premium_spec <- NULL

  if (!lt_missing &&
      exists(".as_life_contract", mode = "function") &&
      .as_life_contract(lt)) {
    contract <- lt

    if (!identical(contract$lives, "single")) {
      stop(
        "`reserve_x()` currently supports only single-life contracts.",
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

      if (n_missing) {
        n <- insurance_spec$n
        n_missing <- FALSE
      }

      if (h_missing) {
        h <- insurance_spec$h
        h_missing <- FALSE
      }

      if (benefit_missing) {
        benefit <- insurance_spec$benefit
        benefit_missing <- FALSE
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
  method <- match.arg(method)

  if (frac_missing) {
    frac <- if (!is.null(premium_spec$frac)) {
      premium_spec$frac
    } else if (!is.null(insurance_spec$frac)) {
      insurance_spec$frac
    } else {
      "UDD"
    }
  }

  frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))

  if (identical(frac, "CML")) {
    frac <- "CF"
  }

  if (!identical(timing, "due")) {
    stop(
      "`reserve_x()` currently supports premiums payable in advance only: ",
      "use `timing = 'due'`.",
      call. = FALSE
    )
  }

  if (identical(method, "recursive") && !identical(as.numeric(k), 1)) {
    stop(
      "The recursive method currently supports annual premiums only (`k = 1`). ",
      "Use `method = 'prospective'` for true k-thly premiums.",
      call. = FALSE
    )
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

    if (!is.numeric(n) ||
        length(n) != 1L ||
        is.na(n) ||
        n <= 0 ||
        (!is.infinite(n) &&
         (!is.finite(n) || abs(n - round(n)) > 1e-10))) {
      stop("`n` must be `Inf` or a single positive integer.", call. = FALSE)
    }

    if (type %in% c("term", "endowment") && is.infinite(n)) {
      stop("`n` must be finite for term and endowment insurance.", call. = FALSE)
    }

    if (!is.numeric(h) ||
        length(h) != 1L ||
        is.na(h) ||
        !is.finite(h) ||
        h < 0 ||
        abs(h - round(h)) > 1e-10) {
      stop("`h` must be a single nonnegative integer.", call. = FALSE)
    }

    if (!is.numeric(benefit) ||
        length(benefit) != 1L ||
        is.na(benefit) ||
        !is.finite(benefit) ||
        benefit <= 0) {
      stop("`benefit` must be a single positive finite number.", call. = FALSE)
    }

    if (!is.null(P) &&
        (!is.numeric(P) ||
         length(P) != 1L ||
         is.na(P) ||
         !is.finite(P))) {
      stop(
        "`P` must be `NULL` or a single finite annualized premium.",
        call. = FALSE
      )
    }

    if (!is.numeric(k) ||
        length(k) != 1L ||
        is.na(k) ||
        !is.finite(k) ||
        k < 1 ||
        abs(k - round(k)) > 1e-10) {
      stop("`k` must be a single positive integer.", call. = FALSE)
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
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (!is.null(n_prem) && !is.infinite(n_prem)) {
    n_prem <- round(n_prem * k) / k
  }

  lt <- lt[order(lt$x), , drop = FALSE]

  if (!is.numeric(lt$x) || !is.numeric(lt$lx)) {
    stop("Columns `x` and `lx` in `lt` must be numeric.", call. = FALSE)
  }

  if (any(is.na(lt$x)) ||
      any(!is.finite(lt$x)) ||
      any(abs(lt$x - round(lt$x)) > 1e-10)) {
    stop("Column `x` must contain finite integer ages.", call. = FALSE)
  }

  if (anyDuplicated(lt$x)) {
    stop("Life table ages in column `x` must be unique.", call. = FALSE)
  }

  if (any(is.na(lt$lx)) ||
      any(!is.finite(lt$lx)) ||
      any(lt$lx < 0)) {
    stop("Column `lx` must contain finite nonnegative values.", call. = FALSE)
  }

  ages <- as.integer(round(lt$x))
  lx_values <- as.numeric(lt$lx)

  if (!x %in% ages) {
    stop("Age `x` is not available in the life table.", call. = FALSE)
  }

  positive_ages <- ages[lx_values > 0]

  if (length(positive_ages) == 0L || x > max(positive_ages)) {
    stop("The life table has no positive exposure at age `x`.", call. = FALSE)
  }

  max_supported_age <- max(positive_ages)

  coverage_end <- if (identical(type, "whole")) {
    max_supported_age - x
  } else {
    h + n
  }

  if (!identical(type, "whole") && x + coverage_end > max(ages)) {
    stop(
      "The requested coverage horizon exceeds the available life table.",
      call. = FALSE
    )
  }

  if (coverage_end < 0) {
    stop("The contract horizon is not supported by the life table.", call. = FALSE)
  }

  if (is.null(n_prem)) {
    n_prem <- if (identical(type, "whole")) Inf else n
  }

  premium_start_time <- if (identical(premium_start, "issue")) 0 else h
  premium_end_time <- if (is.infinite(n_prem)) {
    Inf
  } else {
    premium_start_time + n_prem
  }

  if (!identical(type, "whole") &&
      premium_end_time > coverage_end + 1e-10) {
    stop(
      "Premium payments must not extend beyond the end of coverage.",
      call. = FALSE
    )
  }

  if (is.null(t)) {
    t_vec <- seq.int(0L, as.integer(coverage_end))
  } else {
    if (!is.numeric(t) ||
        length(t) < 1L ||
        any(is.na(t)) ||
        any(!is.finite(t)) ||
        any(abs(t - round(t)) > 1e-10)) {
      stop("`t` must be an integer vector of policy durations.", call. = FALSE)
    }

    t_vec <- sort(unique(as.integer(round(t))))

    if (any(t_vec < 0L) || any(t_vec > coverage_end)) {
      stop(
        "`t` must lie between 0 and the supported contract horizon.",
        call. = FALSE
      )
    }
  }

  future_benefit_apv <- function(tt) {
    current_age <- x + tt

    if (identical(type, "whole")) {
      remaining_h <- max(0L, h - tt)

      return(
        insurance_x(
          lt = lt,
          x = current_age,
          i = i,
          i_type = i_type,
          m = m,
          type = "whole",
          benefit = benefit,
          n = Inf,
          h = remaining_h,
          tidy = FALSE
        )
      )
    }

    if (tt < h) {
      return(
        insurance_x(
          lt = lt,
          x = current_age,
          i = i,
          i_type = i_type,
          m = m,
          type = type,
          benefit = benefit,
          n = n,
          h = h - tt,
          tidy = FALSE
        )
      )
    }

    elapsed_coverage <- tt - h
    remaining_coverage <- n - elapsed_coverage

    if (remaining_coverage < 0) {
      return(0)
    }

    if (remaining_coverage == 0) {
      return(if (identical(type, "endowment")) benefit else 0)
    }

    insurance_x(
      lt = lt,
      x = current_age,
      i = i,
      i_type = i_type,
      m = m,
      type = type,
      benefit = benefit,
      n = remaining_coverage,
      h = 0,
      tidy = FALSE
    )
  }

  future_premium_factor <- function(tt) {
    current_age <- x + tt

    if (tt < premium_start_time) {
      return(
        annuity_x(
          lt = lt,
          x = current_age,
          i = i,
          i_type = i_type,
          m = m,
          n = n_prem,
          h = premium_start_time - tt,
          k = k,
          timing = "due",
          woolhouse = woolhouse,
          frac = frac,
          tidy = FALSE
        )
      )
    }

    if (is.infinite(n_prem)) {
      return(
        annuity_x(
          lt = lt,
          x = current_age,
          i = i,
          i_type = i_type,
          m = m,
          n = Inf,
          h = 0,
          k = k,
          timing = "due",
          woolhouse = woolhouse,
          frac = frac,
          tidy = FALSE
        )
      )
    }

    remaining_premium_term <- premium_end_time - tt

    if (remaining_premium_term <= 0) {
      return(0)
    }

    annuity_x(
      lt = lt,
      x = current_age,
      i = i,
      i_type = i_type,
      m = m,
      n = remaining_premium_term,
      h = 0,
      k = k,
      timing = "due",
      woolhouse = woolhouse,
      frac = frac,
      tidy = FALSE
    )
  }

  premium_was_computed <- is.null(P)

  if (premium_was_computed) {
    P <- premium_x(
      lt = lt,
      x = x,
      i = i,
      i_type = i_type,
      m = m,
      type = type,
      benefit = benefit,
      n = n,
      h = h,
      k = k,
      frac = frac,
      timing = "due",
      premium_start = premium_start,
      n_prem = n_prem,
      woolhouse = woolhouse,
      output = "value",
      check = check
    )
  }

  apv_benefits_vec <- vapply(
    t_vec,
    future_benefit_apv,
    numeric(1L)
  )

  premium_factors_vec <- vapply(
    t_vec,
    future_premium_factor,
    numeric(1L)
  )

  apv_premiums_vec <- P * premium_factors_vec
  prospective_reserves <- apv_benefits_vec - apv_premiums_vec

  if (identical(method, "prospective")) {
    reserves <- prospective_reserves
  } else {
    i_effective <- standardize_interest(
      i_type = i_type,
      i = i,
      m = m
    )

    age_index <- stats::setNames(seq_along(ages), ages)

    get_lx <- function(age_value) {
      key <- as.character(age_value)

      if (key %in% names(age_index)) {
        return(lx_values[[age_index[[key]]]])
      }

      if (age_value == max(ages) + 1L) {
        return(0)
      }

      NA_real_
    }

    one_year_qx <- function(age_value) {
      l0 <- get_lx(age_value)
      l1 <- get_lx(age_value + 1L)

      if (is.na(l0) || is.na(l1) || l0 <= 0) {
        return(NA_real_)
      }

      1 - l1 / l0
    }

    full_t <- seq.int(0L, as.integer(coverage_end))
    full_reserves <- numeric(length(full_t))
    full_reserves[[1L]] <- future_benefit_apv(0) -
      P * future_premium_factor(0)

    if (length(full_t) > 1L) {
      for (idx in seq_len(length(full_t) - 1L)) {
        tt <- full_t[[idx]]
        current_age <- x + tt

        qx_t <- one_year_qx(current_age)

        if (is.na(qx_t)) {
          stop(
            "The life table does not support recursive reserve computation.",
            call. = FALSE
          )
        }

        px_t <- 1 - qx_t

        if (px_t <= 0) {
          stop(
            "The recursive reserve is undefined after survival probability ",
            "reaches zero.",
            call. = FALSE
          )
        }

        premium_due <- if (
          tt >= premium_start_time &&
          tt < premium_end_time
        ) {
          P
        } else {
          0
        }

        death_benefit <- if (identical(type, "whole")) {
          if (tt >= h) benefit else 0
        } else {
          if (tt >= h && tt < h + n) benefit else 0
        }

        full_reserves[[idx + 1L]] <-
          (
            (full_reserves[[idx]] + premium_due) *
              (1 + i_effective) -
              death_benefit * qx_t
          ) / px_t
      }
    }

    reserves <- full_reserves[match(t_vec, full_t)]
  }

  names(reserves) <- paste0("t=", t_vec)

  if (identical(output, "value")) {
    return(reserves)
  }

  premium_per_payment <- P / k

  if (identical(output, "summary")) {
    return(
      tibble::tibble(
        t = t_vec,
        age = x + t_vec,
        reserve = as.numeric(reserves),
        premium_annualized = P,
        premium_per_payment = premium_per_payment,
        method = method
      )
    )
  }

  tibble::tibble(
    t = rep(t_vec, each = 5L),
    age = rep(x + t_vec, each = 5L),
    component = rep(
      c(
        "apv_future_benefits",
        "apv_future_premiums",
        "reserve",
        "premium_annualized",
        "premium_per_payment"
      ),
      times = length(t_vec)
    ),
    value = as.numeric(
      unlist(
        Map(
          function(b, p, v) {
            c(
              b,
              p,
              v,
              P,
              premium_per_payment
            )
          },
          apv_benefits_vec,
          apv_premiums_vec,
          reserves
        ),
        use.names = FALSE
      )
    ),
    unit = rep(
      c(
        "currency",
        "currency",
        "currency",
        "currency per year",
        "currency per payment"
      ),
      times = length(t_vec)
    )
  )
}
