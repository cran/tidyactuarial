#' Actuarial present value of a two-life annuity
#'
#' Computes the actuarial present value of a discrete annuity contingent on two
#' independent lives. Supports joint-life, last-survivor, and state-based
#' reversionary-style payments.
#'
#' @param mortality_table Either a single life table used for both lives, a list
#'   of two life tables `list(table_x, table_y)`, or a `tidyact_life_contract`
#'   object created by [life_contract()]. Each life table must contain columns
#'   `x` and `lx`.
#' @param age_x Integer actuarial age for the first life.
#' @param age_y Integer actuarial age for the second life.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are `"effective"`, `"nominal_interest"`, `"nominal_discount"`, and `"force"`.
#' @param m Positive integer. Compounding frequency for nominal rates.
#' @param cohort Character string. Use `"first"` for joint-life payments while
#'   both lives are alive, or `"last"` for last-survivor payments while at least
#'   one life is alive. Ignored when `benefit` is supplied.
#' @param benefit Optional list with weights `both`, `x_only`, and `y_only`.
#' @param term_years Term in years. Use `Inf` for the maximum horizon allowed
#'   by the life tables.
#' @param deferment_years Integer deferment period in years.
#' @param payments_per_year Positive integer. Number of payments per year.
#' @param timing Payment timing. Use `"immediate"` or `"due"`.
#' @param woolhouse Woolhouse approximation for `payments_per_year > 1`:
#'   `"none"`, `"first"`, or `"second"`.
#' @param frac Fractional-age assumption for exact k-thly computation:
#'   `"UDD"`, `"CF"`, `"CML"`, or `"Balducci"`.
#' @param output Character string. Use `"value"` for a numeric APV or `"table"`
#'   for a one-row tibble.
#'
#' @details This function assumes independent future lifetimes.
#'
#' @return If `output = "value"`, a numeric scalar. If `output = "table"`,
#' a one-row tibble.
#'
#' @seealso [annuity_x()], [insurance_xy()], [premium_xy()], [t_pxy()]
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' annuity_xy(
#'   mortality_table = lt,
#'   age_x = 60,
#'   age_y = 62,
#'   rate = 0.05,
#'   cohort = "first",
#'   timing = "due"
#' )
#'
#' lt |>
#'   life_contract(lives = "joint", age_x = 60, age_y = 62, rate = 0.05) |>
#'   annuity_xy(cohort = "last", timing = "due")
#'
#' @export
annuity_xy <- function(
    mortality_table,
    age_x = NULL,
    age_y = NULL,
    rate = NULL,
    rate_type = NULL,
    m = NULL,
    cohort = c("first", "last"),
    benefit = NULL,
    term_years = Inf,
    deferment_years = 0L,
    payments_per_year = 1L,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    frac,
    output = c("value", "table")
) {
  cohort <- match.arg(cohort)
  timing <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)
  output <- match.arg(output)

  if (inherits(mortality_table, "tidyact_life_contract")) {
    contract <- mortality_table

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop("`annuity_xy()` requires a two-life `life_contract()` object.", call. = FALSE)
    }

    mortality_table <- contract$mortality_table
    if (is.null(age_x)) age_x <- contract$age_x
    if (is.null(age_y)) age_y <- contract$age_y
    if (is.null(rate)) rate <- contract$rate
    if (is.null(rate_type)) rate_type <- contract$rate_type
    if (is.null(m)) m <- contract$m

    if (contract$lives == "last_survivor") {
      cohort <- "last"
    }
  }

  if (is.null(rate_type)) rate_type <- "effective"
  if (is.null(m)) m <- 1L

  if (is.data.frame(mortality_table)) {
    if (!all(c("x", "lx") %in% names(mortality_table))) {
      stop("`mortality_table` must contain columns `x` and `lx`.", call. = FALSE)
    }
    table_x <- mortality_table
    table_y <- mortality_table
  } else if (
    is.list(mortality_table) &&
    length(mortality_table) == 2L &&
    all(vapply(mortality_table, is.data.frame, logical(1L)))
  ) {
    if (!all(c("x", "lx") %in% names(mortality_table[[1L]]))) {
      stop("The first life table must contain columns `x` and `lx`.", call. = FALSE)
    }
    if (!all(c("x", "lx") %in% names(mortality_table[[2L]]))) {
      stop("The second life table must contain columns `x` and `lx`.", call. = FALSE)
    }
    table_x <- mortality_table[[1L]]
    table_y <- mortality_table[[2L]]
  } else {
    stop(
      "`mortality_table` must be one life table, a list of two life tables, ",
      "or a `tidyact_life_contract` object.",
      call. = FALSE
    )
  }

  if (missing(frac)) {
    frac_x <- attr(table_x, "frac")
    frac_y <- attr(table_y, "frac")
    ok_x <- !is.null(frac_x) && frac_x %in% c("UDD", "CF", "Balducci")
    ok_y <- !is.null(frac_y) && frac_y %in% c("UDD", "CF", "Balducci")

    if (ok_x && ok_y) {
      if (!identical(frac_x, frac_y)) {
        stop("The two life tables carry different `frac` attributes. Supply `frac` explicitly.", call. = FALSE)
      }
      frac <- frac_x
    } else if (ok_x) {
      frac <- frac_x
    } else if (ok_y) {
      frac <- frac_y
    } else {
      frac <- "UDD"
    }
  } else {
    frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac == "CML") frac <- "CF"
  }

  if (is.null(rate) || !is.numeric(rate) || length(rate) != 1L ||
      is.na(rate) || !is.finite(rate)) {
    stop("`rate` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(rate_type) || length(rate_type) != 1L || is.na(rate_type)) {
    stop("`rate_type` must be a single character string.", call. = FALSE)
  }

  if (is.null(age_x) || !is.numeric(age_x) || length(age_x) != 1L ||
      is.na(age_x) || !is.finite(age_x) || abs(age_x - round(age_x)) > 1e-10) {
    stop("`age_x` must be a single integer age.", call. = FALSE)
  }

  if (is.null(age_y) || !is.numeric(age_y) || length(age_y) != 1L ||
      is.na(age_y) || !is.finite(age_y) || abs(age_y - round(age_y)) > 1e-10) {
    stop("`age_y` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      !is.finite(m) || m < 1 || abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(deferment_years) || length(deferment_years) != 1L ||
      is.na(deferment_years) || !is.finite(deferment_years) ||
      deferment_years < 0 || abs(deferment_years - round(deferment_years)) > 1e-10) {
    stop("`deferment_years` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(payments_per_year) || length(payments_per_year) != 1L ||
      is.na(payments_per_year) || !is.finite(payments_per_year) ||
      payments_per_year < 1 || abs(payments_per_year - round(payments_per_year)) > 1e-10) {
    stop("`payments_per_year` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(term_years) || length(term_years) != 1L ||
      is.na(term_years) || term_years < 0 ||
      (!is.infinite(term_years) &&
       (!is.finite(term_years) || abs(term_years - round(term_years)) > 1e-10))) {
    stop("`term_years` must be `Inf` or a single nonnegative integer.", call. = FALSE)
  }

  age_x <- as.integer(round(age_x))
  age_y <- as.integer(round(age_y))
  m <- as.integer(round(m))
  deferment_years <- as.integer(round(deferment_years))
  payments_per_year <- as.integer(round(payments_per_year))
  if (!is.infinite(term_years)) term_years <- as.integer(round(term_years))

  i_effective <- standardize_interest(type = rate_type, rate = rate, m = m)
  if (i_effective <= -1) {
    stop("The standardized annual effective interest rate must be greater than -1.", call. = FALSE)
  }

  v_fun <- function(tt) (1 + i_effective)^(-tt)

  omega_x <- max(table_x$x, na.rm = TRUE)
  omega_y <- max(table_y$x, na.rm = TRUE)

  horizon_x <- max(0L, omega_x - age_x)
  horizon_y <- max(0L, omega_y - age_y)

  if (is.null(benefit)) {
    horizon <- if (cohort == "first") min(horizon_x, horizon_y) else max(horizon_x, horizon_y)
  } else {
    if (!is.list(benefit) || !all(c("both", "x_only", "y_only") %in% names(benefit))) {
      stop("`benefit` must be a list with names `both`, `x_only`, and `y_only`.", call. = FALSE)
    }
    h_both <- if (!is.na(benefit$both) && benefit$both != 0) min(horizon_x, horizon_y) else 0L
    h_x_only <- if (!is.na(benefit$x_only) && benefit$x_only != 0) horizon_x else 0L
    h_y_only <- if (!is.na(benefit$y_only) && benefit$y_only != 0) horizon_y else 0L
    horizon <- max(h_both, h_x_only, h_y_only)
  }

  term_used <- if (is.infinite(term_years)) {
    max(0L, horizon - deferment_years)
  } else {
    term_years
  }

  if (term_used == 0L) {
    result <- 0
    if (output == "value") return(result)

    return(tibble::tibble(
      age_x = age_x,
      age_y = age_y,
      rate = rate,
      rate_type = rate_type,
      m = m,
      i_effective = i_effective,
      term_years = term_years,
      term_used = term_used,
      deferment_years = deferment_years,
      payments_per_year = payments_per_year,
      timing = timing,
      cohort = cohort,
      woolhouse = woolhouse,
      frac = frac,
      apv = result
    ))
  }

  if (is.null(benefit)) {
    benefit <- if (cohort == "first") {
      list(both = 1, x_only = 0, y_only = 0)
    } else {
      list(both = 1, x_only = 1, y_only = 1)
    }
  }

  expected_payment <- function(tt) {
    px_t <- t_px(table_x, x = age_x, t = tt, frac = frac, check = FALSE)
    py_t <- t_px(table_y, x = age_y, t = tt, frac = frac, check = FALSE)

    if (is.na(px_t) || is.na(py_t)) return(NA_real_)

    p_both <- px_t * py_t
    p_x_only <- px_t * (1 - py_t)
    p_y_only <- py_t * (1 - px_t)

    benefit$both * p_both +
      benefit$x_only * p_x_only +
      benefit$y_only * p_y_only
  }

  exact_apv <- function(nn, kk, tim) {
    N <- nn * kk
    if (N == 0L) return(0)

    j <- if (tim == "due") 0:(N - 1L) else 1:N
    u <- j / kk
    times <- deferment_years + u
    disc <- v_fun(times)
    ep <- vapply(times, expected_payment, numeric(1L))

    if (anyNA(ep)) return(NA_real_)

    sum((1 / kk) * disc * ep)
  }

  if (payments_per_year == 1L || woolhouse == "none") {
    result <- exact_apv(term_used, payments_per_year, timing)
  } else {
    annual_due <- exact_apv(term_used, 1L, "due")

    ep_start <- expected_payment(deferment_years)
    ep_end <- expected_payment(deferment_years + term_used)

    if (is.na(ep_start) || ep_start <= 0) ep_start <- 1
    if (is.na(ep_end)) ep_end <- 0

    nEx_status <- v_fun(term_used) * ep_end / ep_start

    adj1 <- (payments_per_year - 1) / (2 * payments_per_year) *
      (1 - nEx_status)

    if (woolhouse == "first") {
      due_k <- annual_due - adj1
    } else {
      delta <- log1p(i_effective)
      ep_m1 <- expected_payment(deferment_years + 1)

      if (is.na(ep_m1) || ep_start <= 0) {
        mu_start <- 0
      } else {
        p_status_1 <- ep_m1 / ep_start
        mu_start <- if (p_status_1 > 0) -log(p_status_1) else 0
      }

      ep_n1 <- expected_payment(deferment_years + term_used + 1)

      if (is.na(ep_n1) || is.na(ep_end) || ep_end <= 0) {
        mu_end <- 0
      } else {
        p_status_n1 <- ep_n1 / ep_end
        mu_end <- if (p_status_n1 > 0) -log(p_status_n1) else 0
      }

      adj2 <- (payments_per_year^2 - 1) /
        (12 * payments_per_year^2) *
        (delta + mu_start - nEx_status * (delta + mu_end))

      due_k <- annual_due - adj1 - adj2
    }

    result <- if (timing == "due") {
      due_k
    } else {
      due_k - (1 / payments_per_year) * (1 - nEx_status)
    }
  }

  if (output == "value") return(result)

  tibble::tibble(
    age_x = age_x,
    age_y = age_y,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    term_years = term_years,
    term_used = term_used,
    deferment_years = deferment_years,
    payments_per_year = payments_per_year,
    timing = timing,
    cohort = cohort,
    woolhouse = woolhouse,
    frac = frac,
    apv = result
  )
}
