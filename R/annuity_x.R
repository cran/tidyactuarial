#' Actuarial present value of a life annuity
#'
#' Computes the actuarial present value of a discrete life annuity using a life
#' table.
#'
#' The function supports:
#' \itemize{
#'   \item whole-life annuities,
#'   \item temporary annuities,
#'   \item integer deferral,
#'   \item annual or k-thly payments,
#'   \item exact fractional survival under UDD,
#'   \item first- and second-order Woolhouse approximations.
#' }
#'
#' @param mortality_table A life table as produced by \code{\link{lifetable}}.
#'   It must contain columns \code{x} and \code{lx}.
#' @param age Integer actuarial age.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Compounding frequency for nominal rates. Ignored
#'   for \code{rate_type = "effective"} and \code{rate_type = "force"}.
#' @param term_years Integer term in years. Use \code{Inf} for a whole-life
#'   annuity.
#' @param deferral_years Integer deferral period in years.
#' @param payments_per_year Positive integer. Number of payments per year.
#'   For example, use \code{12} for monthly payments.
#' @param timing Character string. Either \code{"immediate"} for payments at
#'   the end of each payment period or \code{"due"} for payments at the
#'   beginning of each payment period.
#' @param woolhouse Character string. For \code{payments_per_year > 1}, use
#'   \code{"none"} for exact UDD, \code{"first"} for the first-order Woolhouse
#'   approximation, or \code{"second"} for the second-order Woolhouse
#'   approximation.
#' @param output Character string. Use \code{"value"} to return a numeric APV
#'   or \code{"table"} to return a one-row tibble with intermediate quantities.
#'
#' @return
#' If \code{output = "value"}, a numeric scalar containing the actuarial present
#' value.
#'
#' If \code{output = "table"}, a one-row tibble with the main input values,
#' equivalent interest rate, deferral factor, pure endowment factor, and APV.
#'
#' @details
#' For annual annuities-due,
#' \deqn{\ddot{a}_{x:\overline{n}|} =
#' \sum_{j=0}^{n-1} v^j\,{}_jp_x.}
#'
#' For annual annuities-immediate,
#' \deqn{a_{x:\overline{n}|} =
#' \sum_{j=1}^{n} v^j\,{}_jp_x.}
#'
#' Deferral is handled through
#' \deqn{v^h\,{}_hp_x}
#' where \eqn{h} is \code{deferral_years}.
#'
#' For k-thly payments with \code{woolhouse = "none"}, fractional survival is
#' computed under UDD.
#'
#' @seealso \code{\link{insurance_x}}, \code{\link{premium_x}},
#'   \code{\link{reserve_x}}, \code{\link{t_px}}, \code{\link{t_Ex}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:65,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000)
#' )
#'
#' # Annual annuity-immediate
#' annuity_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   timing = "immediate"
#' )
#'
#' # Annual annuity-due
#' annuity_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   timing = "due"
#' )
#'
#' # Temporary annuity
#' annuity_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   term_years = 3,
#'   timing = "due"
#' )
#'
#' # Deferred annuity
#' annuity_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   deferral_years = 2,
#'   timing = "due"
#' )
#'
#' # Table output
#' annuity_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.06,
#'   term_years = 3,
#'   timing = "due",
#'   output = "table"
#' )
#'
#' @export
annuity_x <- function(
    mortality_table,
    age,
    rate,
    rate_type = "effective",
    m = 1L,
    term_years = Inf,
    deferral_years = 0L,
    payments_per_year = 1L,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    output = c("value", "table")
) {
  timing <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)
  output <- match.arg(output)

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (.as_life_contract(mortality_table)) {
    contract <- mortality_table

    if (!identical(contract$lives, "single")) {
      stop(
        "`annuity_x()` currently supports only single-life `life_contract()` objects.",
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

  if (!is.numeric(age) ||
      length(age) != 1L ||
      is.na(age) ||
      !is.finite(age) ||
      abs(age - round(age)) > 1e-10) {
    stop("`age` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(rate) ||
      length(rate) != 1L ||
      is.na(rate) ||
      !is.finite(rate)) {
    stop("`rate` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(rate_type) ||
      length(rate_type) != 1L ||
      is.na(rate_type)) {
    stop("`rate_type` must be a single character string.", call. = FALSE)
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      is.na(m) ||
      !is.finite(m) ||
      m < 1 ||
      abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(deferral_years) ||
      length(deferral_years) != 1L ||
      is.na(deferral_years) ||
      !is.finite(deferral_years) ||
      deferral_years < 0 ||
      abs(deferral_years - round(deferral_years)) > 1e-10) {
    stop("`deferral_years` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(payments_per_year) ||
      length(payments_per_year) != 1L ||
      is.na(payments_per_year) ||
      !is.finite(payments_per_year) ||
      payments_per_year < 1 ||
      abs(payments_per_year - round(payments_per_year)) > 1e-10) {
    stop("`payments_per_year` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(term_years) ||
      length(term_years) != 1L ||
      is.na(term_years) ||
      term_years < 0 ||
      (!is.infinite(term_years) &&
       (!is.finite(term_years) ||
        abs(term_years - round(term_years)) > 1e-10))) {
    stop(
      "`term_years` must be `Inf` or a single nonnegative integer.",
      call. = FALSE
    )
  }

  age <- as.integer(round(age))
  m <- as.integer(round(m))
  deferral_years <- as.integer(round(deferral_years))
  payments_per_year <- as.integer(round(payments_per_year))

  if (!is.infinite(term_years)) {
    term_years <- as.integer(round(term_years))
  }

  # -------------------------------------------------------------------------
  # Interest conversion
  # -------------------------------------------------------------------------

  i_effective <- standardize_interest(
    type = rate_type,
    rate = rate,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  v_pow <- function(tt) {
    (1 + i_effective)^(-tt)
  }

  # -------------------------------------------------------------------------
  # Life table preparation
  # -------------------------------------------------------------------------

  lt <- mortality_table[order(mortality_table$x), , drop = FALSE]

  if (!is.numeric(lt$x)) {
    stop("Column `x` in `mortality_table` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(lt$lx)) {
    stop("Column `lx` in `mortality_table` must be numeric.", call. = FALSE)
  }

  if (any(is.na(lt$x)) || any(!is.finite(lt$x))) {
    stop("Column `x` must contain finite non-missing values.", call. = FALSE)
  }

  if (any(abs(lt$x - round(lt$x)) > 1e-10)) {
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

  t_p_int <- function(current_age, tt) {
    if (tt == 0) {
      return(1)
    }

    l0 <- get_lx(current_age)
    l1 <- get_lx(current_age + tt)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    l1 / l0
  }

  t_p_udd <- function(current_age, u) {
    if (u < 0) {
      return(NA_real_)
    }

    if (u == 0) {
      return(1)
    }

    tt <- floor(u)
    s <- u - tt

    pt <- t_p_int(current_age, tt)

    if (is.na(pt)) {
      return(NA_real_)
    }

    if (s == 0) {
      return(pt)
    }

    y <- current_age + tt

    ly <- get_lx(y)
    ly1 <- get_lx(y + 1L)

    if (is.na(ly) || is.na(ly1) || ly <= 0) {
      return(NA_real_)
    }

    dy <- ly - ly1
    ps <- (ly - s * dy) / ly

    pt * ps
  }

  # -------------------------------------------------------------------------
  # Deferral
  # -------------------------------------------------------------------------

  deferment_factor <- v_pow(deferral_years) *
    t_p_int(age, deferral_years)

  if (is.na(deferment_factor)) {
    stop(
      "The deferral age `age + deferral_years` is outside the life table ",
      "or `lx(age)` is zero.",
      call. = FALSE
    )
  }

  start_age <- age + deferral_years

  # -------------------------------------------------------------------------
  # Term
  # -------------------------------------------------------------------------

  max_years <- max(0L, (omega + 1L) - start_age)

  if (is.infinite(term_years)) {
    term_used <- max_years
  } else {
    term_used <- term_years

    if (term_used > max_years) {
      stop(
        "`term_years` exceeds the horizon allowed by the life table. ",
        "The table must support ages up to `age + deferral_years + term_years`.",
        call. = FALSE
      )
    }
  }

  if (term_used == 0L) {
    result <- 0

    if (output == "value") {
      return(result)
    }

    return(tibble::tibble(
      age = age,
      rate = rate,
      rate_type = rate_type,
      m = m,
      i_effective = i_effective,
      term_years = term_years,
      term_used = term_used,
      deferral_years = deferral_years,
      start_age = start_age,
      payments_per_year = payments_per_year,
      timing = timing,
      woolhouse = woolhouse,
      deferment_factor = deferment_factor,
      pure_endowment_factor = 0,
      apv = result
    ))
  }

  pure_endowment_factor <- v_pow(term_used) *
    t_p_int(start_age, term_used)

  if (is.na(pure_endowment_factor)) {
    pure_endowment_factor <- 0
  }

  # -------------------------------------------------------------------------
  # Annual and k-thly computation helpers
  # -------------------------------------------------------------------------

  annual_exact <- function(current_age, nn, tim) {
    if (tim == "due") {
      times <- 0:(nn - 1L)
    } else {
      times <- 1:nn
    }

    survival <- vapply(
      times,
      function(tt) t_p_int(current_age, tt),
      numeric(1L)
    )

    if (anyNA(survival)) {
      stop(
        "The life table does not support the required ages for annual payments.",
        call. = FALSE
      )
    }

    sum(v_pow(times) * survival)
  }

  kthly_exact_udd <- function(current_age, nn, kk, tim) {
    if (tim == "due") {
      j <- 0:(kk * nn - 1L)
    } else {
      j <- 1:(kk * nn)
    }

    u <- j / kk

    survival <- vapply(
      u,
      function(uu) t_p_udd(current_age, uu),
      numeric(1L)
    )

    if (anyNA(survival)) {
      stop(
        "The life table does not support the required ages for UDD k-thly payments.",
        call. = FALSE
      )
    }

    sum((1 / kk) * v_pow(u) * survival)
  }

  # -------------------------------------------------------------------------
  # Main APV computation
  # -------------------------------------------------------------------------

  if (payments_per_year == 1L) {
    annuity_value_at_start <- annual_exact(
      current_age = start_age,
      nn = term_used,
      tim = timing
    )
  } else if (woolhouse == "none") {
    annuity_value_at_start <- kthly_exact_udd(
      current_age = start_age,
      nn = term_used,
      kk = payments_per_year,
      tim = timing
    )
  } else {
    annual_due <- annual_exact(
      current_age = start_age,
      nn = term_used,
      tim = "due"
    )

    adj1 <- (payments_per_year - 1) / (2 * payments_per_year) *
      (1 - pure_endowment_factor)

    if (woolhouse == "first") {
      due_k <- annual_due - adj1
    } else {
      delta <- log1p(i_effective)

      ly <- get_lx(start_age)
      ly1 <- get_lx(start_age + 1L)

      if (is.na(ly) || is.na(ly1) || ly <= 0) {
        stop("Cannot compute the force approximation at the starting age.", call. = FALSE)
      }

      p_y <- ly1 / ly
      mu_y <- if (!is.na(p_y) && p_y > 0) -log(p_y) else 0

      lyn <- get_lx(start_age + term_used)
      lyn1 <- get_lx(start_age + term_used + 1L)

      if (!is.na(lyn) && !is.na(lyn1) && lyn > 0) {
        p_yn <- lyn1 / lyn
        mu_yn <- if (p_yn > 0) -log(p_yn) else 0
      } else {
        mu_yn <- 0
      }

      adj2 <- (payments_per_year^2 - 1) /
        (12 * payments_per_year^2) *
        (
          delta + mu_y -
            pure_endowment_factor * (delta + mu_yn)
        )

      due_k <- annual_due - adj1 - adj2
    }

    annuity_value_at_start <- if (timing == "due") {
      due_k
    } else {
      due_k - (1 / payments_per_year) * (1 - pure_endowment_factor)
    }
  }

  result <- deferment_factor * annuity_value_at_start

  if (output == "value") {
    return(result)
  }

  tibble::tibble(
    age = age,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    term_years = term_years,
    term_used = term_used,
    deferral_years = deferral_years,
    start_age = start_age,
    payments_per_year = payments_per_year,
    timing = timing,
    woolhouse = woolhouse,
    deferment_factor = deferment_factor,
    pure_endowment_factor = pure_endowment_factor,
    annuity_value_at_start = annuity_value_at_start,
    apv = result
  )
}
