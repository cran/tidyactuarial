#' Net premium for two-life insurance by the equivalence principle
#'
#' Computes the net benefit premium for a two-life insurance contract using the
#' equivalence principle.
#'
#' The premium returned corresponds to one premium payment. For example, if
#' `payments_per_year = 1`, it is an annual premium; if
#' `payments_per_year = 12`, it is a monthly premium.
#'
#' This function separates:
#' \itemize{
#'   \item the insurance coverage period, controlled by `insurance_type`,
#'   `term_years`, and `deferment_years`;
#'   \item the premium-paying period, controlled by `premium_term_years`,
#'   `premium_start`, `premium_timing`, and `payments_per_year`.
#' }
#'
#' @param mortality_table Either a single life table used for both lives, or a
#'   list of two life tables `list(table_x, table_y)`. Each table must contain
#'   column `x` and at least one of `lx`, `px`, or `qx`. A
#'   `tidyact_life_contract` object created by [life_contract()] is also
#'   accepted.
#' @param age_x Integer actuarial age for the first life.
#' @param age_y Integer actuarial age for the second life.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are `"effective"`, `"nominal_interest"`, `"nominal_discount"`, and `"force"`.
#' @param m Positive integer. Compounding frequency for nominal rates.
#' @param insurance_type Type of insurance: `"whole"`, `"term"`,
#'   `"endowment"`, or `"pure_endowment"`.
#' @param benefit Benefit amount. Must be a single nonnegative number.
#' @param term_years Optional insurance term in years after deferment. Required
#'   for `"term"`, `"endowment"`, and `"pure_endowment"`.
#' @param deferment_years Nonnegative integer deferral period in years.
#' @param payments_per_year Number of premium payments per year.
#' @param frac Fractional-age assumption used for fractional premium payment
#'   times: `"UDD"`, `"CF"`, `"CML"`, or `"Balducci"`.
#' @param premium_timing Timing of premium payments: `"due"` for payments in
#'   advance or `"immediate"` for payments in arrears.
#' @param premium_start Start of premium payments: `"issue"` for time 0 or
#'   `"deferred"` for time `deferment_years`.
#' @param premium_term_years Optional premium-paying term in years, counted from
#'   `premium_start`. If `NULL`, defaults to the available status horizon for
#'   whole-life products and to `term_years` for finite products.
#' @param cohort Status definition: `"first"` for joint-life status or `"last"`
#'   for last-survivor status.
#' @param output Character string. Use `"value"` for a numeric premium or
#'   `"table"` for a one-row tibble with details.
#' @param check Logical. If `TRUE`, performs input validation.
#' @param tol Numeric tolerance for integer checks.
#'
#' @return
#' If `output = "value"`, a numeric net premium per payment.
#' If `output = "table"`, a one-row tibble with premium details.
#'
#' @details
#' The function assumes independent future lifetimes.
#'
#' The supported two-life statuses are:
#' \itemize{
#'   \item `cohort = "first"`: joint-life status. The status survives while
#'   both lives survive.
#'   \item `cohort = "last"`: last-survivor status. The status survives while
#'   at least one life survives.
#' }
#'
#' Let `P(t)` denote the probability that the selected two-life status survives
#' `t` years from issue. For integer death benefits paid at the end of the year
#' of status failure, the benefit APV for a term insurance deferred `h` years
#' and lasting `n` years is
#' \deqn{
#' B \sum_{r=h}^{h+n-1} v^{r+1}\{P(r)-P(r+1)\}.
#' }
#'
#' For an endowment insurance, the pure endowment benefit
#' `B v^(h+n) P(h+n)` is added.
#'
#' Premiums are contingent on the selected two-life status being in force at
#' the premium payment time.
#'
#' @seealso [premium_x()], [insurance_xy()], [annuity_xy()], [reserve_xy()]
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
#'   mortality_table = lt,
#'   age_x = 60,
#'   age_y = 62,
#'   rate = 0.05,
#'   insurance_type = "term",
#'   term_years = 5,
#'   premium_term_years = 3,
#'   payments_per_year = 12,
#'   cohort = "last",
#'   benefit = 100000,
#'   output = "table"
#' )
#'
#' lt |>
#'   life_contract(lives = "joint", age_x = 60, age_y = 62, rate = 0.05) |>
#'   premium_xy(
#'     insurance_type = "term",
#'     term_years = 5,
#'     premium_term_years = 3,
#'     payments_per_year = 12,
#'     cohort = "first",
#'     benefit = 100000
#'   )
#'
#' @export
premium_xy <- function(
    mortality_table,
    age_x = NULL,
    age_y = NULL,
    rate = NULL,
    rate_type = NULL,
    m = NULL,
    insurance_type = c("whole", "term", "endowment", "pure_endowment"),
    benefit = 1,
    term_years = Inf,
    deferment_years = 0L,
    payments_per_year = 1L,
    frac = c("UDD", "CF", "CML", "Balducci"),
    premium_timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    premium_term_years = NULL,
    cohort = c("first", "last"),
    output = c("value", "table"),
    check = TRUE,
    tol = 1e-10
) {
  insurance_type <- match.arg(insurance_type)
  frac <- match.arg(frac)
  premium_timing <- match.arg(premium_timing)
  premium_start <- match.arg(premium_start)
  cohort <- match.arg(cohort)
  output <- match.arg(output)

  if (frac == "CML") {
    frac <- "CF"
  }

  # -------------------------------------------------------------------------
  # Resolve life_contract input
  # -------------------------------------------------------------------------

  if (inherits(mortality_table, "tidyact_life_contract")) {
    contract <- mortality_table

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop(
        "`premium_xy()` requires a two-life `life_contract()` object.",
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
  # Validation
  # -------------------------------------------------------------------------

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  if (is.null(rate) ||
      !is.numeric(rate) ||
      length(rate) != 1L ||
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
      !is.finite(age_x) ||
      abs(age_x - round(age_x)) > tol) {
    stop("`age_x` must be a single integer-valued age.", call. = FALSE)
  }

  if (is.null(age_y) ||
      !is.numeric(age_y) ||
      length(age_y) != 1L ||
      !is.finite(age_y) ||
      abs(age_y - round(age_y)) > tol) {
    stop("`age_y` must be a single integer-valued age.", call. = FALSE)
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      !is.finite(m) ||
      m <= 0 ||
      abs(m - round(m)) > tol) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(deferment_years) ||
      length(deferment_years) != 1L ||
      !is.finite(deferment_years) ||
      deferment_years < 0 ||
      abs(deferment_years - round(deferment_years)) > tol) {
    stop("`deferment_years` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(payments_per_year) ||
      length(payments_per_year) != 1L ||
      !is.finite(payments_per_year) ||
      payments_per_year <= 0 ||
      abs(payments_per_year - round(payments_per_year)) > tol) {
    stop("`payments_per_year` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      !is.finite(benefit) ||
      benefit < 0) {
    stop("`benefit` must be a single nonnegative finite number.", call. = FALSE)
  }

  if (!is.numeric(term_years) ||
      length(term_years) != 1L ||
      is.na(term_years) ||
      term_years < 0 ||
      (!is.infinite(term_years) &&
       (!is.finite(term_years) ||
        abs(term_years - round(term_years)) > tol))) {
    stop("`term_years` must be `Inf` or a single nonnegative integer.", call. = FALSE)
  }

  age_x <- as.integer(round(age_x))
  age_y <- as.integer(round(age_y))
  m <- as.integer(round(m))
  deferment_years <- as.integer(round(deferment_years))
  payments_per_year <- as.integer(round(payments_per_year))

  if (!is.infinite(term_years)) {
    term_years <- as.integer(round(term_years))
  }

  if (insurance_type %in% c("term", "endowment", "pure_endowment") &&
      is.infinite(term_years)) {
    stop(
      "`term_years` must be finite for term, endowment, and pure endowment insurance.",
      call. = FALSE
    )
  }

  if (!is.null(premium_term_years)) {
    if (!is.numeric(premium_term_years) ||
        length(premium_term_years) != 1L ||
        !is.finite(premium_term_years) ||
        premium_term_years < 0 ||
        abs(premium_term_years - round(premium_term_years)) > tol) {
      stop(
        "`premium_term_years` must be NULL or a single nonnegative integer.",
        call. = FALSE
      )
    }

    premium_term_years <- as.integer(round(premium_term_years))
  }

  if (insurance_type %in% c("term", "endowment", "pure_endowment") &&
      !is.null(premium_term_years) &&
      premium_term_years > term_years) {
    stop(
      "`premium_term_years` must not exceed `term_years` for finite two-life products.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Interest
  # -------------------------------------------------------------------------

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

  v <- 1 / (1 + i_effective)

  # -------------------------------------------------------------------------
  # Life table preparation
  # -------------------------------------------------------------------------

  normalize_life_tables <- function(mortality_table) {
    if (is.data.frame(mortality_table)) {
      return(list(mortality_table, mortality_table))
    }

    if (!is.list(mortality_table)) {
      stop(
        "`mortality_table` must be a data.frame/tibble, a list of two life tables, ",
        "or a `tidyact_life_contract` object.",
        call. = FALSE
      )
    }

    if (length(mortality_table) == 1L) {
      return(list(mortality_table[[1L]], mortality_table[[1L]]))
    }

    if (length(mortality_table) != 2L) {
      stop("When `mortality_table` is a list, its length must be 1 or 2.", call. = FALSE)
    }

    mortality_table
  }

  validate_life_table <- function(tab, idx) {
    if (!is.data.frame(tab)) {
      stop(
        "Each element of `mortality_table` must be a data.frame/tibble. ",
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

  table_list <- normalize_life_tables(mortality_table)

  table_list <- list(
    prepare_life_table(table_list[[1L]], 1L),
    prepare_life_table(table_list[[2L]], 2L)
  )

  ages <- c(age_x, age_y)

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

  status_horizon <- if (cohort == "first") {
    min(horizon_life)
  } else {
    max(horizon_life)
  }

  if (insurance_type %in% c("term", "endowment", "pure_endowment") &&
      (deferment_years + term_years) > status_horizon) {
    stop(
      "`term_years + deferment_years` exceeds the available two-life status horizon.",
      call. = FALSE
    )
  }

  if (insurance_type == "whole" && deferment_years > status_horizon) {
    stop(
      "`deferment_years` exceeds the available two-life status horizon.",
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

    if (cohort == "first") {
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
    insurance_type,
    whole = {
      if (deferment_years >= status_horizon) {
        0
      } else {
        term_insurance_apv(
          start = deferment_years,
          term = status_horizon - deferment_years
        )
      }
    },
    term = term_insurance_apv(
      start = deferment_years,
      term = term_years
    ),
    pure_endowment = pure_endowment_apv(
      time = deferment_years + term_years
    ),
    endowment = term_insurance_apv(
      start = deferment_years,
      term = term_years
    ) +
      pure_endowment_apv(
        time = deferment_years + term_years
      )
  )

  # -------------------------------------------------------------------------
  # Premium-paying period
  # -------------------------------------------------------------------------

  if (is.null(premium_term_years)) {
    premium_term_years <- if (insurance_type == "whole") {
      start_for_default <- if (premium_start == "issue") 0L else deferment_years
      max(0L, status_horizon - start_for_default)
    } else {
      term_years
    }
  }

  premium_start_time <- if (premium_start == "issue") 0L else deferment_years

  if (premium_start_time + premium_term_years > status_horizon) {
    stop(
      "Premium-paying period exceeds the available two-life status horizon.",
      call. = FALSE
    )
  }

  if (premium_term_years <= 0L) {
    stop("`premium_term_years` must define a positive premium annuity.", call. = FALSE)
  }

  if (insurance_type %in% c("term", "endowment", "pure_endowment") &&
      premium_term_years > term_years) {
    stop(
      "`premium_term_years` must not exceed `term_years` for finite two-life products.",
      call. = FALSE
    )
  }

  n_payments <- premium_term_years * payments_per_year

  payment_times <- if (premium_timing == "due") {
    premium_start_time + (0:(n_payments - 1L)) / payments_per_year
  } else {
    premium_start_time + (1:n_payments) / payments_per_year
  }

  apv_premiums <- sum(vapply(payment_times, function(tt) {
    (v^tt) * status_survival(tt)
  }, numeric(1L)))

  if (!is.finite(apv_premiums) || apv_premiums <= 0) {
    stop("APV of premium annuity is nonpositive or not finite.", call. = FALSE)
  }

  premium <- apv_benefits / apv_premiums

  if (output == "value") {
    return(premium)
  }

  tibble::tibble(
    age_x = age_x,
    age_y = age_y,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    deferment_years = deferment_years,
    term_years = if (is.infinite(term_years)) NA_integer_ else term_years,
    insurance_type = insurance_type,
    cohort = cohort,
    benefit = benefit,
    payments_per_year = payments_per_year,
    frac = frac,
    premium_timing = premium_timing,
    premium_start = premium_start,
    premium_term_years = premium_term_years,
    premium = premium,
    premium_annual = payments_per_year * premium,
    apv_benefits = apv_benefits,
    apv_premiums = apv_premiums
  )
}
