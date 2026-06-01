#' Net premium for two-life insurance by the equivalence principle
#'
#' Computes the net benefit premium for a two-life insurance contract using the
#' equivalence principle and compact actuarial notation.
#'
#' The premium returned corresponds to one premium payment. For example, if
#' \code{k = 1}, it is an annual premium; if \code{k = 12}, it is a monthly
#' premium.
#'
#' The function separates:
#' \itemize{
#'   \item the insurance coverage period, controlled by \code{type},
#'   \code{n}, and \code{h};
#'   \item the premium-paying period, controlled by \code{n_prem},
#'   \code{premium_start}, \code{timing}, and \code{k}.
#' }
#'
#' @param lt Either a single life table used for both lives, a list of two life
#'   tables \code{list(lt_x, lt_y)}, or a \code{tidyact_life_contract} object
#'   created by \code{\link{life_contract}}. Each table must contain column
#'   \code{x} and at least one of \code{lx}, \code{px}, or \code{qx}.
#' @param x Integer actuarial age for the first life.
#' @param y Integer actuarial age for the second life.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param type Type of insurance: \code{"whole"}, \code{"term"},
#'   \code{"endowment"}, or \code{"pure_endowment"}.
#' @param benefit Benefit amount. Must be a single nonnegative number.
#' @param n Insurance term in years after deferment. Required as finite for
#'   \code{"term"}, \code{"endowment"}, and \code{"pure_endowment"}.
#' @param h Nonnegative integer deferment period in years.
#' @param k Positive integer. Number of premium payments per year.
#' @param frac Fractional-age assumption used for fractional premium payment
#'   times: \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#' @param timing Timing of premium payments: \code{"due"} for payments in
#'   advance or \code{"immediate"} for payments in arrears.
#' @param premium_start Start of premium payments: \code{"issue"} for time 0 or
#'   \code{"deferred"} for time \code{h}.
#' @param n_prem Optional premium-paying term in years, counted from
#'   \code{premium_start}. If \code{NULL}, defaults to the available status
#'   horizon for whole-life products and to \code{n} for finite products.
#' @param status Two-life status definition: \code{"joint"} for the joint-life
#'   status or \code{"last"} for the last-survivor status.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric premium. If
#'   \code{TRUE}, returns a one-row tibble with details.
#' @param check Logical. If \code{TRUE}, performs input validation.
#' @param tol Numeric tolerance for integer checks.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age_x}, \code{age_y}, \code{rate},
#'   \code{rate_type}, \code{insurance_type}, \code{term_years},
#'   \code{deferment_years}, \code{payments_per_year},
#'   \code{premium_timing}, \code{premium_term_years}, \code{cohort}, and
#'   \code{output}.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric net premium per payment.
#'
#' If \code{tidy = TRUE}, a one-row tibble with premium details.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table input, \code{x} and
#' \code{y} are the two actuarial ages, \code{i} is the interest-rate input,
#' \code{i_type} is the interest-rate type, \code{m} is the conversion
#' frequency for nominal rates, \code{n} is the insurance term, \code{h} is the
#' deferment period, and \code{k} is the premium payment frequency.
#'
#' The function assumes independent future lifetimes.
#'
#' Let \eqn{S(t)} denote the probability that the selected two-life status
#' survives \eqn{t} years from issue. For integer death benefits paid at the end
#' of the year of status failure, the benefit APV for a term insurance deferred
#' \eqn{h} years and lasting \eqn{n} years is
#' \deqn{
#' B \sum_{r=h}^{h+n-1} v^{r+1}\{S(r)-S(r+1)\}.
#' }
#'
#' For an endowment insurance, the pure endowment benefit
#' \eqn{B v^{h+n} S(h+n)} is added.
#'
#' Premiums are contingent on the selected two-life status being in force at the
#' premium payment time.
#'
#' @seealso \code{\link{premium_x}}, \code{\link{insurance_xy}},
#'   \code{\link{annuity_xy}}, \code{\link{reserve_xy}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:110,
#'   lx = seq(100000, 0, length.out = 51)
#' )
#'
#' premium_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   type = "term",
#'   n = 5,
#'   n_prem = 3,
#'   k = 12,
#'   status = "last",
#'   benefit = 100000,
#'   tidy = TRUE
#' )
#'
#' lt |>
#'  life_contract(lives = "joint", x = 60, y = 62, i = 0.05) |>
#'   premium_xy(
#'    type = "term",
#'    n = 5,
#'    n_prem = 3,
#'    k = 12,
#'    status = "joint",
#'    benefit = 100000)
#' @export
premium_xy <- function(
    lt,
    x = NULL,
    y = NULL,
    i = NULL,
    i_type = NULL,
    m = NULL,
    type = c("whole", "term", "endowment", "pure_endowment"),
    benefit = 1,
    n = Inf,
    h = 0L,
    k = 1L,
    frac = c("UDD", "CF", "CML", "Balducci"),
    timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    n_prem = NULL,
    status = c("joint", "last"),
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10,
    ...
) {
  dots <- list(...)
  status_missing <- missing(status)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "mortality_table",
    "age_x",
    "age_y",
    "rate",
    "rate_type",
    "insurance_type",
    "term_years",
    "deferment_years",
    "payments_per_year",
    "premium_timing",
    "premium_term_years",
    "cohort",
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

  if (!is.null(dots$age_x)) {
    if (!is.null(x)) {
      stop("Provide only one of `x` or deprecated `age_x`.", call. = FALSE)
    }

    x <- dots$age_x
  }

  if (!is.null(dots$age_y)) {
    if (!is.null(y)) {
      stop("Provide only one of `y` or deprecated `age_y`.", call. = FALSE)
    }

    y <- dots$age_y
  }

  if (!is.null(dots$rate)) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dots$rate
  }

  if (!is.null(dots$rate_type)) {
    if (!is.null(i_type)) {
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

  if (!is.null(dots$deferment_years)) {
    if (!identical(h, 0L) && !identical(h, 0)) {
      stop("Provide only one of `h` or deprecated `deferment_years`.", call. = FALSE)
    }

    h <- dots$deferment_years
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

  if (!is.null(dots$cohort)) {
    if (!status_missing) {
      stop("Provide only one of `status` or deprecated `cohort`.", call. = FALSE)
    }

    old_cohort <- match.arg(dots$cohort, c("first", "last"))
    status <- if (old_cohort == "first") "joint" else "last"
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
  status <- match.arg(status)

  if (frac == "CML") {
    frac <- "CF"
  }

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  # -------------------------------------------------------------------------
  # Resolve life_contract input
  # -------------------------------------------------------------------------

  if (!missing(lt) && inherits(lt, "tidyact_life_contract")) {
    contract <- lt

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop(
        "`premium_xy()` requires a two-life `life_contract()` object.",
        call. = FALSE
      )
    }

    lt <- contract$mortality_table

    if (is.null(x)) {
      x <- contract$x %||% contract$age_x
    }

    if (is.null(y)) {
      y <- contract$y %||% contract$age_y
    }

    if (is.null(i)) {
      i <- contract$i %||% contract$rate
    }

    if (is.null(i_type)) {
      i_type <- contract$i_type %||% contract$rate_type
    }

    if (is.null(m)) {
      m <- contract$m
    }

    if (contract$lives == "last_survivor" && status_missing) {
      status <- "last"
    }
  }

  if (is.null(i_type)) {
    i_type <- "effective"
  }

  if (is.null(m)) {
    m <- 1L
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      !is.finite(tol) ||
      tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  if (is.null(i) ||
      !is.numeric(i) ||
      length(i) != 1L ||
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

  if (is.null(x) ||
      !is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      abs(x - round(x)) > tol) {
    stop("`x` must be a single integer-valued age.", call. = FALSE)
  }

  if (is.null(y) ||
      !is.numeric(y) ||
      length(y) != 1L ||
      !is.finite(y) ||
      abs(y - round(y)) > tol) {
    stop("`y` must be a single integer-valued age.", call. = FALSE)
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      !is.finite(m) ||
      m <= 0 ||
      abs(m - round(m)) > tol) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(h) ||
      length(h) != 1L ||
      !is.finite(h) ||
      h < 0 ||
      abs(h - round(h)) > tol) {
    stop("`h` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(k) ||
      length(k) != 1L ||
      !is.finite(k) ||
      k <= 0 ||
      abs(k - round(k)) > tol) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      !is.finite(benefit) ||
      benefit < 0) {
    stop("`benefit` must be a single nonnegative finite number.", call. = FALSE)
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n < 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) || abs(n - round(n)) > tol))) {
    stop("`n` must be `Inf` or a single nonnegative integer.", call. = FALSE)
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (type %in% c("term", "endowment", "pure_endowment") &&
      is.infinite(n)) {
    stop(
      "`n` must be finite for term, endowment, and pure endowment insurance.",
      call. = FALSE
    )
  }

  if (!is.null(n_prem)) {
    if (!is.numeric(n_prem) ||
        length(n_prem) != 1L ||
        !is.finite(n_prem) ||
        n_prem < 0 ||
        abs(n_prem * k - round(n_prem * k)) > tol) {
      stop(
        "`n_prem` must be NULL or a single nonnegative value satisfying `n_prem * k` integer.",
        call. = FALSE
      )
    }

    n_prem <- round(n_prem * k) / k
  }

  if (type %in% c("term", "endowment", "pure_endowment") &&
      !is.null(n_prem) &&
      n_prem > n) {
    stop(
      "`n_prem` must not exceed `n` for finite two-life products.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Interest
  # -------------------------------------------------------------------------

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  v <- 1 / (1 + i_effective)

  # -------------------------------------------------------------------------
  # Life table preparation
  # -------------------------------------------------------------------------

  normalize_life_tables <- function(lt) {
    if (is.data.frame(lt)) {
      return(list(lt, lt))
    }

    if (!is.list(lt)) {
      stop(
        "`lt` must be a data.frame/tibble, a list of two life tables, ",
        "or a `tidyact_life_contract` object.",
        call. = FALSE
      )
    }

    if (length(lt) == 1L) {
      return(list(lt[[1L]], lt[[1L]]))
    }

    if (length(lt) != 2L) {
      stop("When `lt` is a list, its length must be 1 or 2.", call. = FALSE)
    }

    lt
  }

  validate_life_table <- function(tab, idx) {
    if (!is.data.frame(tab)) {
      stop(
        "Each element of `lt` must be a data.frame/tibble. ",
        "Problem at life ", idx, ".",
        call. = FALSE
      )
    }

    if (!("x" %in% names(tab))) {
      stop(
        "Each life table must contain column `x`. Problem at life ",
        idx, ".",
        call. = FALSE
      )
    }

    if (!any(c("lx", "px", "qx") %in% names(tab))) {
      stop(
        "Each life table must contain at least one of `lx`, `px`, or `qx`. ",
        "Problem at life ", idx, ".",
        call. = FALSE
      )
    }

    if (!is.numeric(tab$x) || any(!is.finite(tab$x))) {
      stop("Column `x` must be numeric and finite in every life table.", call. = FALSE)
    }

    invisible(TRUE)
  }

  prepare_life_table <- function(tab, idx) {
    if (isTRUE(check)) {
      validate_life_table(tab, idx)
    }

    tab <- tibble::as_tibble(tab)
    tab <- tab[order(tab$x), , drop = FALSE]

    if (anyDuplicated(tab$x)) {
      stop("Life table ages must be unique. Problem at life ", idx, ".", call. = FALSE)
    }

    if (any(abs(tab$x - round(tab$x)) > tol)) {
      stop("Life table ages must be integer-valued. Problem at life ", idx, ".", call. = FALSE)
    }

    tab$x <- as.integer(round(tab$x))

    if (!("px" %in% names(tab))) {
      if ("qx" %in% names(tab)) {
        tab$px <- 1 - as.numeric(tab$qx)
      } else if ("lx" %in% names(tab)) {
        lx <- as.numeric(tab$lx)
        px <- rep(NA_real_, length(lx))

        if (length(lx) >= 2L) {
          valid <- seq_len(length(lx) - 1L)
          px[valid] <- ifelse(lx[valid] > 0, lx[valid + 1L] / lx[valid], 0)
        }

        px[length(px)] <- 0
        tab$px <- px
      }
    }

    if (!("qx" %in% names(tab))) {
      tab$qx <- 1 - as.numeric(tab$px)
    }

    tab$px <- pmin(pmax(as.numeric(tab$px), 0), 1)
    tab$qx <- pmin(pmax(as.numeric(tab$qx), 0), 1)

    tab
  }

  table_list <- normalize_life_tables(lt)

  table_list <- list(
    prepare_life_table(table_list[[1L]], 1L),
    prepare_life_table(table_list[[2L]], 2L)
  )

  ages <- c(x, y)

  age_available <- vapply(seq_along(table_list), function(j) {
    ages[[j]] %in% table_list[[j]]$x
  }, logical(1L))

  if (!all(age_available)) {
    stop("Issue ages must be present in their corresponding life tables.", call. = FALSE)
  }

  max_age <- vapply(table_list, function(tab) max(tab$x, na.rm = TRUE), numeric(1L))
  horizon_life <- floor(max_age - ages)

  if (any(horizon_life < 0)) {
    stop("Life table horizon is insufficient for the requested issue ages.", call. = FALSE)
  }

  status_horizon <- if (status == "joint") {
    min(horizon_life)
  } else {
    max(horizon_life)
  }

  if (type %in% c("term", "endowment", "pure_endowment") &&
      (h + n) > status_horizon) {
    stop(
      "`n + h` exceeds the available two-life status horizon.",
      call. = FALSE
    )
  }

  if (type == "whole" && h > status_horizon) {
    stop(
      "`h` exceeds the available two-life status horizon.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Survival helpers
  # -------------------------------------------------------------------------

  single_life_survival <- function(tab, age0, tt, horizon) {
    if (tt <= 0) return(1)
    if (tt > horizon) return(0)

    n_int <- floor(tt)
    s <- tt - n_int

    if (abs(s) <= tol) s <- 0

    if (abs(s - 1) <= tol) {
      n_int <- n_int + 1L
      s <- 0
    }

    result <- 1

    if (n_int > 0L) {
      ages_int <- age0 + 0:(n_int - 1L)
      px_vec <- tab$px[match(ages_int, tab$x)]

      if (any(is.na(px_vec))) return(NA_real_)

      result <- prod(px_vec)
    }

    if (s > tol) {
      age_tail <- age0 + n_int
      idx_tail <- match(age_tail, tab$x)

      if (is.na(idx_tail)) return(NA_real_)

      q_tail <- tab$qx[[idx_tail]]
      p_tail <- tab$px[[idx_tail]]

      frac_surv <- if (frac == "UDD") {
        1 - s * q_tail
      } else if (frac == "CF") {
        p_tail^s
      } else {
        denominator <- 1 - (1 - s) * q_tail
        if (denominator <= 0) 0 else p_tail / denominator
      }

      result <- result * frac_surv
    }

    pmin(pmax(result, 0), 1)
  }

  life_survival <- function(j, tt) {
    single_life_survival(
      tab = table_list[[j]],
      age0 = ages[[j]],
      tt = tt,
      horizon = horizon_life[[j]]
    )
  }

  status_survival <- function(tt) {
    p <- vapply(1:2, life_survival, numeric(1L), tt = tt)

    if (anyNA(p)) {
      stop(
        "Could not compute two-life status survival at time ",
        tt,
        ".",
        call. = FALSE
      )
    }

    if (status == "joint") {
      prod(p)
    } else {
      1 - prod(1 - p)
    }
  }

  status_decrement_probability <- function(r) {
    pr <- status_survival(r)
    pr1 <- status_survival(r + 1)
    pmin(pmax(pr - pr1, 0), 1)
  }

  # -------------------------------------------------------------------------
  # APV of benefits
  # -------------------------------------------------------------------------

  term_insurance_apv <- function(start, term) {
    if (term <= 0L) return(0)

    years <- start:(start + term - 1L)

    sum(vapply(years, function(r) {
      benefit * (v^(r + 1)) * status_decrement_probability(r)
    }, numeric(1L)))
  }

  pure_endowment_apv <- function(time) {
    benefit * (v^time) * status_survival(time)
  }

  apv_benefits <- switch(
    type,
    whole = {
      if (h >= status_horizon) {
        0
      } else {
        term_insurance_apv(
          start = h,
          term = status_horizon - h
        )
      }
    },
    term = term_insurance_apv(
      start = h,
      term = n
    ),
    pure_endowment = pure_endowment_apv(
      time = h + n
    ),
    endowment = term_insurance_apv(
      start = h,
      term = n
    ) +
      pure_endowment_apv(
        time = h + n
      )
  )

  # -------------------------------------------------------------------------
  # Premium-paying period
  # -------------------------------------------------------------------------

  if (is.null(n_prem)) {
    n_prem <- if (type == "whole") {
      start_for_default <- if (premium_start == "issue") 0L else h
      max(0L, status_horizon - start_for_default)
    } else {
      n
    }
  }

  premium_start_time <- if (premium_start == "issue") 0L else h

  if (premium_start_time + n_prem > status_horizon) {
    stop(
      "Premium-paying period exceeds the available two-life status horizon.",
      call. = FALSE
    )
  }

  if (n_prem <= 0) {
    stop("`n_prem` must define a positive premium annuity.", call. = FALSE)
  }

  if (type %in% c("term", "endowment", "pure_endowment") &&
      n_prem > n) {
    stop(
      "`n_prem` must not exceed `n` for finite two-life products.",
      call. = FALSE
    )
  }

  n_payments_raw <- n_prem * k
  n_payments <- round(n_payments_raw)

  if (abs(n_payments_raw - n_payments) > tol) {
    stop("`n_prem * k` must be an integer.", call. = FALSE)
  }

  n_payments <- as.integer(n_payments)

  payment_times <- if (timing == "due") {
    premium_start_time + (0:(n_payments - 1L)) / k
  } else {
    premium_start_time + (1:n_payments) / k
  }

  a_premiums <- sum(vapply(payment_times, function(tt) {
    (v^tt) * status_survival(tt)
  }, numeric(1L)))

  if (!is.finite(a_premiums) || a_premiums <= 0) {
    stop("APV of premium annuity is nonpositive or not finite.", call. = FALSE)
  }

  P <- apv_benefits / a_premiums

  if (!tidy) {
    return(P)
  }

  tibble::tibble(
    x = x,
    y = y,
    age_x = x,
    age_y = y,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective,
    h = h,
    deferment_years = h,
    n = n,
    term_years = if (is.infinite(n)) NA_integer_ else n,
    type = type,
    insurance_type = type,
    status = status,
    cohort = if (status == "joint") "first" else "last",
    benefit = benefit,
    k = k,
    payments_per_year = k,
    frac = frac,
    timing = timing,
    premium_timing = timing,
    premium_start = premium_start,
    n_prem = n_prem,
    premium_term_years = n_prem,
    P = P,
    premium = P,
    P_annual = k * P,
    premium_annual = k * P,
    apv_benefits = apv_benefits,
    a_premiums = a_premiums,
    apv_premiums = a_premiums
  )
}
