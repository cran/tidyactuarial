#' Benefit reserves for two-life insurance
#'
#' Computes prospective or recursive terminal reserves for joint-life and
#' last-survivor insurance. For k-thly premiums, `P` is the annualized premium
#' and each installment equals `P / k`.
#'
#' @param lt One life table, `list(lt_x, lt_y)`, or a two-life
#'   [life_contract()] object.
#' @param x,y Integer ages at issue. Optional for a life contract.
#' @param i Annual interest-rate input. Optional for a life contract.
#' @param i_type Interest-rate type.
#' @param m Conversion frequency for nominal rates.
#' @param type One of `"whole"`, `"term"`, `"endowment"`, or
#'   `"pure_endowment"`.
#' @param status `"joint"` or `"last"`.
#' @param n Insurance term in years after deferment.
#' @param h Nonnegative integer deferment in years.
#' @param benefit Positive benefit amount.
#' @param P Optional annualized premium. If `NULL`, [premium_xy()] is used.
#' @param n_prem Premium-paying term. Fractional values require `n_prem * k`
#'   to be an integer.
#' @param k Number of premium payments per year.
#' @param timing `"due"` or `"immediate"`.
#' @param premium_start `"issue"` or `"deferred"`.
#' @param frac Fractional-age assumption.
#' @param woolhouse `"none"`, `"first"`, or `"second"`.
#' @param t Integer policy durations.
#' @param method `"prospective"` or `"recursive"`. Recursion currently
#'   requires `h = 0`, `k = 1`, due premiums, and premiums starting at issue.
#' @param output `"summary"`, `"value"`, or `"audit"`. `"table"` is a
#'   compatibility alias for `"summary"`.
#' @param tidy Deprecated logical output selector.
#' @param check Logical input-check switch.
#' @param tol Numeric tolerance.
#' @param ... Deprecated argument aliases.
#'
#' @details
#' The prospective reserve is
#' \deqn{{}_tV =
#' \operatorname{APV}_t(\text{future benefits}) -
#' P^{(k)}\operatorname{APV}_t(\text{future premium annuity}).}
#' Because [annuity_xy()] assigns amount `1 / k` to each payment,
#' `P` is annualized and the actual installment is `P / k`.
#'
#' Reserves are measured immediately before a premium payable at duration
#' `t`. Pure endowments are valued as endowment APV minus term-insurance APV.
#'
#' @return A compact tibble, a named numeric vector, or a long audit tibble.
#'
#' @seealso [premium_xy()], [annuity_xy()], [insurance_xy()], [reserve_x()]
#' @family life-contingencies
#' @export
reserve_xy <- function(
    lt,
    x = NULL,
    y = NULL,
    i = NULL,
    i_type = "effective",
    m = 1L,
    type = c("whole", "term", "endowment", "pure_endowment"),
    status = c("joint", "last"),
    n = Inf,
    h = 0L,
    benefit = 1,
    P = NULL,
    n_prem = NULL,
    k = 1L,
    timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    frac = c("UDD", "CF", "CML", "Balducci"),
    woolhouse = c("none", "first", "second"),
    t = NULL,
    method = c("prospective", "recursive"),
    output = c("summary", "value", "audit", "table"),
    tidy = NULL,
    check = TRUE,
    tol = 1e-10,
    ...
) {
  lt_missing <- missing(lt)
  x_missing <- missing(x)
  y_missing <- missing(y)
  i_missing <- missing(i)
  i_type_missing <- missing(i_type)
  m_missing <- missing(m)
  type_missing <- missing(type)
  status_missing <- missing(status)
  n_missing <- missing(n)
  h_missing <- missing(h)
  benefit_missing <- missing(benefit)
  n_prem_missing <- missing(n_prem)
  k_missing <- missing(k)
  timing_missing <- missing(timing)
  premium_start_missing <- missing(premium_start)
  frac_missing <- missing(frac)
  woolhouse_missing <- missing(woolhouse)
  output_missing <- missing(output)

  dots <- list(...)
  allowed_old <- c(
    "mortality_table", "age_x", "age_y", "rate", "rate_type",
    "insurance_type", "cohort", "term_years", "deferment_years",
    "premium", "premium_term_years", "payments_per_year",
    "premium_timing", "at"
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
    if (!lt_missing) stop("Provide only one of `lt` or `mortality_table`.", call. = FALSE)
    lt <- dots$mortality_table
    lt_missing <- FALSE
  }
  if ("age_x" %in% names(dots)) {
    if (!x_missing && !is.null(x)) stop("Provide only one of `x` or `age_x`.", call. = FALSE)
    x <- dots$age_x
    x_missing <- FALSE
  }
  if ("age_y" %in% names(dots)) {
    if (!y_missing && !is.null(y)) stop("Provide only one of `y` or `age_y`.", call. = FALSE)
    y <- dots$age_y
    y_missing <- FALSE
  }
  if ("rate" %in% names(dots)) {
    if (!i_missing && !is.null(i)) stop("Provide only one of `i` or `rate`.", call. = FALSE)
    i <- dots$rate
    i_missing <- FALSE
  }
  if ("rate_type" %in% names(dots)) {
    if (!i_type_missing) stop("Provide only one of `i_type` or `rate_type`.", call. = FALSE)
    i_type <- dots$rate_type
    i_type_missing <- FALSE
  }
  if ("insurance_type" %in% names(dots)) {
    if (!type_missing) stop("Provide only one of `type` or `insurance_type`.", call. = FALSE)
    type <- dots$insurance_type
    type_missing <- FALSE
  }
  if ("cohort" %in% names(dots)) {
    if (!status_missing) stop("Provide only one of `status` or `cohort`.", call. = FALSE)
    old_cohort <- match.arg(dots$cohort, c("first", "last"))
    status <- if (old_cohort == "first") "joint" else "last"
    status_missing <- FALSE
  }
  if ("term_years" %in% names(dots)) {
    if (!n_missing) stop("Provide only one of `n` or `term_years`.", call. = FALSE)
    n <- dots$term_years
    n_missing <- FALSE
  }
  if ("deferment_years" %in% names(dots)) {
    if (!h_missing) stop("Provide only one of `h` or `deferment_years`.", call. = FALSE)
    h <- dots$deferment_years
    h_missing <- FALSE
  }
  if ("premium" %in% names(dots)) {
    if (!is.null(P)) stop("Provide only one of `P` or `premium`.", call. = FALSE)
    P <- dots$premium
  }
  if ("premium_term_years" %in% names(dots)) {
    if (!n_prem_missing) stop("Provide only one of `n_prem` or `premium_term_years`.", call. = FALSE)
    n_prem <- dots$premium_term_years
    n_prem_missing <- FALSE
  }
  if ("payments_per_year" %in% names(dots)) {
    if (!k_missing) stop("Provide only one of `k` or `payments_per_year`.", call. = FALSE)
    k <- dots$payments_per_year
    k_missing <- FALSE
  }
  if ("premium_timing" %in% names(dots)) {
    if (!timing_missing) stop("Provide only one of `timing` or `premium_timing`.", call. = FALSE)
    timing <- dots$premium_timing
    timing_missing <- FALSE
  }
  if ("at" %in% names(dots)) {
    if (!is.null(t)) stop("Provide only one of `t` or `at`.", call. = FALSE)
    t <- dots$at
  }

  if (!is.null(tidy)) {
    if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
      stop("`tidy` must be `NULL` or a logical scalar.", call. = FALSE)
    }
    if (!output_missing) stop("Provide only one of `output` or `tidy`.", call. = FALSE)
    output <- if (isTRUE(tidy)) "summary" else "value"
  }

  output <- match.arg(output)
  if (output == "table") output <- "summary"

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol < 0) {
    stop("`tol` must be a nonnegative finite scalar.", call. = FALSE)
  }

  contract_input <- NULL
  insurance_spec <- NULL
  premium_spec <- NULL

  if (!lt_missing &&
      exists(".as_life_contract", mode = "function") &&
      .as_life_contract(lt)) {
    contract_input <- lt

    if (!contract_input$lives %in% c("joint", "last_survivor")) {
      stop("`reserve_xy()` requires a two-life contract.", call. = FALSE)
    }

    lt <- contract_input$lt
    insurance_spec <- contract_input$insurance
    premium_spec <- contract_input$premium_schedule

    if (x_missing || is.null(x)) {
      x <- contract_input$x
      x_missing <- FALSE
    }
    if (y_missing || is.null(y)) {
      y <- contract_input$y
      y_missing <- FALSE
    }
    if (i_missing || is.null(i)) {
      i <- contract_input$i
      i_missing <- FALSE
    }
    if (i_type_missing || is.null(i_type)) {
      i_type <- contract_input$i_type
      i_type_missing <- FALSE
    }
    if (m_missing || is.null(m)) {
      m <- contract_input$m
      m_missing <- FALSE
    }
    if (status_missing) {
      status <- if (contract_input$lives == "joint") "joint" else "last"
      status_missing <- FALSE
    }

    if (!is.null(insurance_spec)) {
      if (insurance_spec$type == "variable_k") {
        stop("`reserve_xy()` does not support `type = 'variable_k'`.", call. = FALSE)
      }
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
      if (n_prem_missing) {
        n_prem <- premium_spec$n_prem
        n_prem_missing <- FALSE
      }
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
      if (woolhouse_missing) {
        woolhouse <- premium_spec$woolhouse
        woolhouse_missing <- FALSE
      }
    }
  }

  type <- match.arg(type)
  status <- match.arg(status)
  timing <- match.arg(timing)
  premium_start <- match.arg(premium_start)
  woolhouse <- match.arg(woolhouse)
  method <- match.arg(method)

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
    frac_benefit <- match.arg(frac)
    frac_premium <- frac_benefit
  }

  valid_frac <- c("UDD", "CF", "CML", "Balducci")
  frac_benefit <- match.arg(frac_benefit, valid_frac)
  frac_premium <- match.arg(frac_premium, valid_frac)
  if (frac_benefit == "CML") frac_benefit <- "CF"
  if (frac_premium == "CML") frac_premium <- "CF"

  scalar_integer <- function(z, lower = 0) {
    is.numeric(z) && length(z) == 1L && !is.na(z) &&
      is.finite(z) && z >= lower && abs(z - round(z)) <= tol
  }

  if (isTRUE(check)) {
    if (lt_missing) stop("`lt` must be provided.", call. = FALSE)
    if (x_missing || is.null(x) || !scalar_integer(x)) {
      stop("`x` must be a single integer age.", call. = FALSE)
    }
    if (y_missing || is.null(y) || !scalar_integer(y)) {
      stop("`y` must be a single integer age.", call. = FALSE)
    }
    if (i_missing || is.null(i) || !is.numeric(i) || length(i) != 1L ||
        is.na(i) || !is.finite(i)) {
      stop("`i` must be a single finite numeric value.", call. = FALSE)
    }
    valid_i_type <- c("effective", "nominal_interest", "nominal_discount", "force")
    if (!is.character(i_type) || length(i_type) != 1L ||
        is.na(i_type) || !i_type %in% valid_i_type) {
      stop("Invalid `i_type`.", call. = FALSE)
    }
    if (!scalar_integer(m, 1)) stop("`m` must be a positive integer.", call. = FALSE)
    if (!scalar_integer(h)) stop("`h` must be a nonnegative integer.", call. = FALSE)
    if (!scalar_integer(k, 1)) stop("`k` must be a positive integer.", call. = FALSE)
    if (!is.numeric(benefit) || length(benefit) != 1L ||
        is.na(benefit) || !is.finite(benefit) || benefit <= 0) {
      stop("`benefit` must be a positive finite scalar.", call. = FALSE)
    }
    if (!is.null(P) && (!is.numeric(P) || length(P) != 1L ||
        is.na(P) || !is.finite(P) || P < 0)) {
      stop("`P` must be `NULL` or a nonnegative annualized premium.", call. = FALSE)
    }
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0 ||
        (!is.infinite(n) && (!is.finite(n) || abs(n - round(n)) > tol))) {
      stop("`n` must be `Inf` or a positive integer.", call. = FALSE)
    }
    if (type %in% c("term", "endowment", "pure_endowment") && is.infinite(n)) {
      stop("`n` must be finite for this insurance type.", call. = FALSE)
    }
    if (!is.null(n_prem) &&
        (!is.numeric(n_prem) || length(n_prem) != 1L || is.na(n_prem) ||
         n_prem <= 0 ||
         (!is.infinite(n_prem) &&
          (!is.finite(n_prem) ||
           abs(n_prem * k - round(n_prem * k)) > tol)))) {
      stop("`n_prem * k` must be an integer.", call. = FALSE)
    }
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))
  if (!is.infinite(n)) n <- as.integer(round(n))
  if (!is.null(n_prem) && !is.infinite(n_prem)) {
    n_prem <- round(n_prem * k) / k
  }

  if (method == "recursive" &&
      (h != 0L || k != 1L || timing != "due" ||
       premium_start != "issue")) {
    stop(
      "The recursive method currently requires `h = 0`, `k = 1`, ",
      "`timing = 'due'`, and `premium_start = 'issue'`.",
      call. = FALSE
    )
  }

  validate_lifetable <- function(tab, label) {
    if (!is.data.frame(tab) || !all(c("x", "lx") %in% names(tab))) {
      stop("`", label, "` must be a life table with `x` and `lx`.", call. = FALSE)
    }
    if (!is.numeric(tab$x) || !is.numeric(tab$lx) ||
        anyNA(tab$x) || anyNA(tab$lx) ||
        any(!is.finite(tab$x)) || any(!is.finite(tab$lx)) ||
        any(abs(tab$x - round(tab$x)) > tol) ||
        any(tab$lx < 0) || anyDuplicated(tab$x)) {
      stop("Invalid life table in `", label, "`.", call. = FALSE)
    }
    tab <- tab[order(tab$x), , drop = FALSE]
    tab$x <- as.integer(round(tab$x))
    tab
  }

  if (is.data.frame(lt)) {
    table_x <- validate_lifetable(lt, "lt")
    table_y <- table_x
    table_use <- table_x
  } else if (is.list(lt) && length(lt) == 2L &&
             all(vapply(lt, is.data.frame, logical(1L)))) {
    table_x <- validate_lifetable(lt[[1L]], "lt[[1]]")
    table_y <- validate_lifetable(lt[[2L]], "lt[[2]]")
    table_use <- list(table_x, table_y)
  } else {
    stop("`lt` must be one life table or a list of two.", call. = FALSE)
  }

  if (!x %in% table_x$x) stop("`x` is absent from the first table.", call. = FALSE)
  if (!y %in% table_y$x) stop("`y` is absent from the second table.", call. = FALSE)

  horizon_x <- max(table_x$x) - x
  horizon_y <- max(table_y$x) - y
  status_horizon <- if (status == "joint") min(horizon_x, horizon_y) else max(horizon_x, horizon_y)

  if (h > status_horizon) stop("`h` exceeds the status horizon.", call. = FALSE)

  contract_horizon <- if (type == "whole") status_horizon else h + n
  if (contract_horizon > status_horizon) {
    stop("`h + n` exceeds the supported status horizon.", call. = FALSE)
  }
  contract_horizon <- as.integer(round(contract_horizon))

  if (is.null(n_prem)) n_prem <- if (type == "whole") Inf else n

  premium_start_time <- if (premium_start == "issue") 0 else h
  premium_end_time <- if (is.infinite(n_prem)) Inf else premium_start_time + n_prem

  if (type != "whole" && premium_end_time > contract_horizon + tol) {
    stop("Premium payments extend beyond coverage.", call. = FALSE)
  }

  if (is.null(t)) {
    t_vec <- seq.int(0L, contract_horizon)
  } else {
    if (!is.numeric(t) || length(t) < 1L || anyNA(t) ||
        any(!is.finite(t)) || any(abs(t - round(t)) > tol)) {
      stop("`t` must be an integer vector.", call. = FALSE)
    }
    t_vec <- sort(unique(as.integer(round(t))))
    if (any(t_vec < 0L) || any(t_vec > contract_horizon)) {
      stop("`t` lies outside the contract horizon.", call. = FALSE)
    }
  }

  insurance_value <- function(cx, cy, current_type, nn, hh) {
    args <- list(
      lt = table_use, x = cx, y = cy, i = i, i_type = i_type, m = m,
      status = status, n = nn, h = hh, benefit = benefit,
      frac = frac_benefit, tidy = FALSE
    )
    if (current_type == "pure_endowment") {
      args$type <- "endowment"
      endowment <- do.call(insurance_xy, args)
      args$type <- "term"
      value <- endowment - do.call(insurance_xy, args)
      return(if (value < 0 && abs(value) <= tol) 0 else value)
    }
    args$type <- current_type
    do.call(insurance_xy, args)
  }

  future_benefit_apv <- function(tt) {
    if (tt >= contract_horizon) {
      return(if (type %in% c("endowment", "pure_endowment") &&
                tt == contract_horizon) benefit else 0)
    }

    cx <- x + tt
    cy <- y + tt

    if (type == "whole") {
      return(insurance_value(cx, cy, "whole", Inf, max(0L, h - tt)))
    }

    if (tt < h) {
      nn <- n
      hh <- h - tt
    } else {
      nn <- h + n - tt
      hh <- 0L
    }

    if (nn <= 0) {
      return(if (type %in% c("endowment", "pure_endowment")) benefit else 0)
    }

    insurance_value(cx, cy, type, nn, hh)
  }

  future_premium_factor <- function(tt) {
    cx <- x + tt
    cy <- y + tt

    if (tt < premium_start_time) {
      return(annuity_xy(
        lt = table_use, x = cx, y = cy, i = i, i_type = i_type, m = m,
        status = status, n = n_prem, h = premium_start_time - tt,
        k = k, timing = timing, woolhouse = woolhouse,
        frac = frac_premium, tidy = FALSE
      ))
    }

    if (is.infinite(n_prem)) {
      return(annuity_xy(
        lt = table_use, x = cx, y = cy, i = i, i_type = i_type, m = m,
        status = status, n = Inf, h = 0, k = k, timing = timing,
        woolhouse = woolhouse, frac = frac_premium, tidy = FALSE
      ))
    }

    remaining <- premium_end_time - tt
    if (remaining <= 0) return(0)

    annuity_xy(
      lt = table_use, x = cx, y = cy, i = i, i_type = i_type, m = m,
      status = status, n = remaining, h = 0, k = k, timing = timing,
      woolhouse = woolhouse, frac = frac_premium, tidy = FALSE
    )
  }

  if (is.null(P)) {
    premium_args <- list(
      lt = if (!is.null(contract_input)) contract_input else table_use,
      x = x, y = y, i = i, i_type = i_type, m = m, type = type,
      status = status, n = n, h = h, benefit = benefit, n_prem = n_prem,
      k = k, timing = timing, premium_start = premium_start,
      woolhouse = woolhouse, output = "value", check = check, tol = tol
    )
    if (is.null(contract_input) || !frac_missing) {
      premium_args$frac <- frac_premium
    }
    P <- do.call(premium_xy, premium_args)
  }

  apv_benefits <- vapply(t_vec, future_benefit_apv, numeric(1L))
  premium_factors <- vapply(t_vec, future_premium_factor, numeric(1L))
  apv_premiums <- P * premium_factors
  prospective <- apv_benefits - apv_premiums

  if (method == "prospective") {
    reserves <- prospective
  } else {
    i_effective <- standardize_interest(i_type = i_type, i = i, m = m)

    one_year_status_p <- function(tt) {
      value <- t_pxy(
        lt = table_use, x = x + tt, y = y + tt, t = 1,
        frac = frac_benefit, status = status
      )
      if (is.na(value)) stop("Cannot compute status continuation.", call. = FALSE)
      value
    }

    full_t <- seq.int(0L, contract_horizon)
    full_reserves <- numeric(length(full_t))
    full_reserves[[1L]] <- future_benefit_apv(0) - P * future_premium_factor(0)

    if (length(full_t) > 1L) {
      for (idx in seq_len(length(full_t) - 1L)) {
        tt <- full_t[[idx]]
        p_status <- one_year_status_p(tt)

        if (p_status <= 0) {
          full_reserves[(idx + 1L):length(full_reserves)] <- 0
          break
        }

        premium_due <- if (tt >= premium_start_time &&
                           tt < premium_end_time) P else 0

        death_benefit <- if (type == "pure_endowment") {
          0
        } else if (type == "whole" || tt < n) {
          benefit
        } else {
          0
        }

        full_reserves[[idx + 1L]] <- (
          (full_reserves[[idx]] + premium_due) * (1 + i_effective) -
            death_benefit * (1 - p_status)
        ) / p_status
      }
    }

    if (type %in% c("endowment", "pure_endowment")) {
      full_reserves[[contract_horizon + 1L]] <- benefit
    }
    if (type == "term") full_reserves[[contract_horizon + 1L]] <- 0

    reserves <- full_reserves[match(t_vec, full_t)]
  }

  names(reserves) <- paste0("t=", t_vec)
  if (output == "value") return(reserves)

  premium_per_payment <- P / k

  if (output == "summary") {
    return(tibble::tibble(
      t = t_vec,
      age_x = x + t_vec,
      age_y = y + t_vec,
      reserve = as.numeric(reserves),
      premium_annualized = P,
      premium_per_payment = premium_per_payment
    ))
  }

  components <- c(
    "apv_future_benefits", "apv_future_premiums", "reserve",
    "premium_annualized", "premium_per_payment"
  )
  values <- cbind(
    apv_benefits, apv_premiums, as.numeric(reserves),
    rep(P, length(t_vec)), rep(premium_per_payment, length(t_vec))
  )

  tibble::tibble(
    t = rep(t_vec, each = length(components)),
    age_x = rep(x + t_vec, each = length(components)),
    age_y = rep(y + t_vec, each = length(components)),
    component = rep(components, times = length(t_vec)),
    value = as.vector(t(values)),
    unit = rep(
      c("currency", "currency", "currency",
        "currency per year", "currency per payment"),
      times = length(t_vec)
    )
  )
}
