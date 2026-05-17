#' Net premium for two-life insurance by the equivalence principle
#'
#' @description
#' Computes the net benefit premium for a two-life insurance contract using the
#' equivalence principle:
#' \deqn{
#' P =
#' \frac{\mathrm{APV\ of\ benefits}}
#'      {\mathrm{APV\ of\ premium\ annuity}}.
#' }
#'
#' The premium returned corresponds to one premium payment. For example, if
#' \code{k = 1}, it is an annual premium; if \code{k = 12}, it is a monthly
#' premium.
#'
#' @details
#' This function is designed to be consistent with \code{\link{premium_x}}.
#' It separates:
#' \itemize{
#'   \item the insurance coverage period, controlled by \code{product},
#'   \code{n}, and \code{m};
#'   \item the premium-paying period, controlled by \code{n_prem},
#'   \code{prem_start}, \code{premium_timing}, and \code{k}.
#' }
#'
#' Therefore, a 5-year two-life term insurance may be paid with premiums during
#' only the first 3 years by using \code{n = 5} and \code{n_prem = 3}.
#'
#' The life-table input \code{lt} may be:
#' \itemize{
#'   \item a single life table, recycled for both lives; or
#'   \item a list of two life tables, one for each life.
#' }
#'
#' The supported two-life statuses are:
#' \itemize{
#'   \item \code{cohort = "first"}: joint-life status. The status survives while
#'   both lives survive.
#'   \item \code{cohort = "last"}: last-survivor status. The status survives
#'   while at least one life survives.
#' }
#'
#' Let \eqn{P(t)} denote the probability that the selected two-life status
#' survives \eqn{t} years from issue. For integer death benefits paid at the
#' end of the year of status failure, the benefit APV for a term insurance
#' deferred \eqn{m} years and lasting \eqn{n} years is
#' \deqn{
#' B \sum_{r=m}^{m+n-1} v^{r+1}\{P(r)-P(r+1)\}.
#' }
#'
#' For an endowment insurance, the pure endowment benefit
#' \eqn{Bv^{m+n}P(m+n)} is added.
#'
#' Premiums are contingent on the two-life status being in force at the premium
#' payment time. If \code{prem_start = "issue"}, premiums start at time 0. If
#' \code{prem_start = "deferred"}, premiums start at time \code{m}. The argument
#' \code{n_prem} is counted from the selected premium start time.
#'
#' @param lt A life table data frame with column \code{x} and at least one of
#'   \code{lx}, \code{px}, or \code{qx}; or a list of two such life tables.
#'   If a single life table is supplied, it is recycled for both lives.
#' @param x Integer actuarial age at issue for the first life.
#' @param y Integer actuarial age at issue for the second life.
#' @param i Effective annual interest rate. Must be greater than \code{-1}.
#' @param product Type of insurance: \code{"whole"}, \code{"term"},
#'   \code{"endowment"}, or \code{"pure_endowment"}.
#' @param type Deprecated alias for \code{product}. Included for compatibility.
#' @param benefit Benefit amount. Must be a single nonnegative number.
#' @param n Optional insurance term in years after deferment. Required for
#'   \code{"term"}, \code{"endowment"}, and \code{"pure_endowment"}.
#' @param m Nonnegative integer deferral period in years. Default is \code{0}.
#' @param k Number of premium payments per year. Default is \code{1}.
#' @param frac Fractional-age assumption used for fractional premium payment
#'   times: \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#'   The option \code{"CML"} is treated as constant force, equivalent to
#'   \code{"CF"}.
#' @param premium_timing Timing of premium payments:
#'   \code{"due"} for payments in advance or \code{"immediate"} for payments
#'   in arrears.
#' @param prem_start Start of premium payments:
#'   \code{"issue"} for time 0 or \code{"deferred"} for time \code{m}.
#' @param n_prem Optional premium-paying term in years, counted from
#'   \code{prem_start}. If \code{NULL}, defaults to the available status
#'   horizon for whole-life products and to \code{n} for temporary products.
#' @param cohort Status definition: \code{"first"} for joint-life status or
#'   \code{"last"} for last-survivor status.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble with details.
#' @param check Logical. If \code{TRUE}, performs input validation.
#' @param tol Numeric tolerance for integer checks.
#'
#' @return A numeric net premium per payment, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{premium_x}}, \code{\link{insurance_xy}},
#'   \code{\link{annuity_xy}}, \code{\link{reserve_xy}}, \code{\link{t_px}}
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:110,
#'   lx = seq(100000, 0, length.out = 51)
#' )
#'
#' # Joint-life whole life insurance, annual premiums
#' premium_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   product = "whole",
#'   cohort = "first",
#'   benefit = 100000
#' )
#'
#' # 5-year last-survivor term insurance paid over only 3 years
#' premium_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   product = "term",
#'   cohort = "last",
#'   n = 5,
#'   n_prem = 3,
#'   benefit = 100000
#' )
#'
#' # Monthly premiums
#' premium_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   product = "term",
#'   cohort = "first",
#'   n = 5,
#'   n_prem = 3,
#'   k = 12,
#'   benefit = 100000,
#'   tidy = TRUE
#' )
#'
#' @export
premium_xy <- function(
    lt,
    x,
    y,
    i,
    product = c("whole", "term", "endowment", "pure_endowment"),
    type = NULL,
    benefit = 1,
    n = NULL,
    m = 0,
    k = 1,
    frac = c("UDD", "CF", "CML", "Balducci"),
    premium_timing = c("due", "immediate"),
    prem_start = c("issue", "deferred"),
    n_prem = NULL,
    cohort = c("first", "last"),
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {
  if (!is.null(type)) {
    product <- type
  }

  product <- match.arg(product)
  frac <- match.arg(frac)
  premium_timing <- match.arg(premium_timing)
  prem_start <- match.arg(prem_start)
  cohort <- match.arg(cohort)

  if (missing(lt)) stop("`lt` is required.", call. = FALSE)
  if (missing(x)) stop("`x` is required.", call. = FALSE)
  if (missing(y)) stop("`y` is required.", call. = FALSE)
  if (missing(i)) stop("`i` is required.", call. = FALSE)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  if (!is.numeric(i) || length(i) != 1L || !is.finite(i) || i <= -1) {
    stop("`i` must be a single numeric value greater than -1.", call. = FALSE)
  }

  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      abs(x - round(x)) > tol) {
    stop("`x` must be a single integer-valued age.", call. = FALSE)
  }

  if (!is.numeric(y) || length(y) != 1L || !is.finite(y) ||
      abs(y - round(y)) > tol) {
    stop("`y` must be a single integer-valued age.", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
      k <= 0 || abs(k - round(k)) > tol) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || !is.finite(m) ||
      m < 0 || abs(m - round(m)) > tol) {
    stop("`m` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) || length(benefit) != 1L ||
      !is.finite(benefit) || benefit < 0) {
    stop("`benefit` must be a single nonnegative finite number.", call. = FALSE)
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  k <- as.integer(round(k))
  m <- as.integer(round(m))

  if (product %in% c("term", "endowment", "pure_endowment")) {
    if (is.null(n)) {
      stop("`n` must be provided for term, endowment, and pure_endowment insurance.", call. = FALSE)
    }

    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n < 0 || abs(n - round(n)) > tol) {
      stop("`n` must be a single nonnegative integer.", call. = FALSE)
    }

    n <- as.integer(round(n))
  }

  if (!is.null(n_prem)) {
    if (!is.numeric(n_prem) || length(n_prem) != 1L ||
        !is.finite(n_prem) || n_prem < 0 ||
        abs(n_prem - round(n_prem)) > tol) {
      stop("`n_prem` must be NULL or a single nonnegative integer.", call. = FALSE)
    }

    n_prem <- as.integer(round(n_prem))
  }

  normalize_life_tables <- function(lt) {
    if (is.data.frame(lt)) {
      return(list(lt, lt))
    }

    if (!is.list(lt)) {
      stop("`lt` must be a data.frame/tibble or a list of two life tables.", call. = FALSE)
    }

    if (length(lt) == 1L) {
      return(list(lt[[1]], lt[[1]]))
    }

    if (length(lt) != 2L) {
      stop("When `lt` is a list, its length must be 1 or 2.", call. = FALSE)
    }

    lt
  }

  validate_life_table <- function(tab, idx) {
    if (!is.data.frame(tab)) {
      stop("Each element of `lt` must be a data.frame/tibble. Problem at life ", idx, ".", call. = FALSE)
    }

    if (!("x" %in% names(tab))) {
      stop("Each life table must contain column `x`. Problem at life ", idx, ".", call. = FALSE)
    }

    if (!any(c("lx", "px", "qx") %in% names(tab))) {
      stop("Each life table must contain at least one of `lx`, `px`, or `qx`. Problem at life ", idx, ".", call. = FALSE)
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

    tab <- dplyr::as_tibble(tab)
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
        valid <- seq_len(length(lx) - 1L)
        px[valid] <- ifelse(lx[valid] > 0, lx[valid + 1L] / lx[valid], 0)
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

  lt_list <- normalize_life_tables(lt)
  lt_list <- list(
    prepare_life_table(lt_list[[1]], 1),
    prepare_life_table(lt_list[[2]], 2)
  )

  ages <- c(x, y)

  age_available <- vapply(seq_along(lt_list), function(j) {
    ages[j] %in% lt_list[[j]]$x
  }, logical(1))

  if (!all(age_available)) {
    stop("Issue ages must be present in their corresponding life tables.", call. = FALSE)
  }

  max_age <- vapply(lt_list, function(tab) max(tab$x, na.rm = TRUE), numeric(1))
  horizon_life <- floor(max_age - ages)

  if (any(horizon_life < 0)) {
    stop("Life table horizon is insufficient for the requested issue ages.", call. = FALSE)
  }

  status_horizon <- if (cohort == "first") {
    min(horizon_life)
  } else {
    max(horizon_life)
  }

  if (product %in% c("term", "endowment", "pure_endowment") &&
      (m + n) > status_horizon) {
    stop("Term plus deferment exceeds the available two-life status horizon.", call. = FALSE)
  }

  if (product == "whole" && m > status_horizon) {
    stop("Deferral exceeds the available two-life status horizon.", call. = FALSE)
  }

  v <- (1 + i)^(-1)

  single_life_survival <- function(tab, age0, tt, horizon) {
    if (tt <= 0) return(1)
    if (tt > horizon) return(0)

    n_int <- floor(tt)
    s <- tt - n_int

    if (abs(s) <= tol) {
      s <- 0
    }

    if (abs(s - 1) <= tol) {
      n_int <- n_int + 1L
      s <- 0
    }

    result <- 1

    if (n_int > 0L) {
      ages_int <- age0 + 0:(n_int - 1L)
      px_vec <- tab$px[match(ages_int, tab$x)]

      if (any(is.na(px_vec))) {
        return(NA_real_)
      }

      result <- prod(px_vec)
    }

    if (s > tol) {
      age_tail <- age0 + n_int
      idx_tail <- match(age_tail, tab$x)

      if (is.na(idx_tail)) {
        return(NA_real_)
      }

      q_tail <- tab$qx[idx_tail]
      p_tail <- tab$px[idx_tail]

      frac_surv <- if (frac == "UDD") {
        1 - s * q_tail
      } else if (frac %in% c("CF", "CML")) {
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
      tab = lt_list[[j]],
      age0 = ages[j],
      tt = tt,
      horizon = horizon_life[j]
    )
  }

  status_survival <- function(tt) {
    p <- vapply(1:2, life_survival, numeric(1), tt = tt)

    if (anyNA(p)) {
      stop("Could not compute two-life status survival at time ", tt, ".", call. = FALSE)
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

  term_insurance_apv <- function(start, term) {
    if (term <= 0L) return(0)

    years <- start:(start + term - 1L)

    sum(vapply(years, function(r) {
      benefit * (v^(r + 1)) * status_decrement_probability(r)
    }, numeric(1)))
  }

  pure_endowment_apv <- function(time) {
    benefit * (v^time) * status_survival(time)
  }

  if (product == "whole") {
    if (m >= status_horizon) {
      apv_benefits <- 0
    } else {
      apv_benefits <- term_insurance_apv(start = m, term = status_horizon - m)
    }
  }

  if (product == "term") {
    apv_benefits <- term_insurance_apv(start = m, term = n)
  }

  if (product == "pure_endowment") {
    apv_benefits <- pure_endowment_apv(time = m + n)
  }

  if (product == "endowment") {
    apv_benefits <- term_insurance_apv(start = m, term = n) +
      pure_endowment_apv(time = m + n)
  }

  if (is.null(n_prem)) {
    n_prem <- if (product == "whole") {
      start_for_default <- if (prem_start == "issue") 0L else m
      max(0L, status_horizon - start_for_default)
    } else {
      n
    }
  }

  premium_start_time <- if (prem_start == "issue") 0 else m

  if (premium_start_time + n_prem > status_horizon) {
    stop("Premium-paying period exceeds the available two-life status horizon.", call. = FALSE)
  }

  if (n_prem <= 0L) {
    stop("`n_prem` must define a positive premium annuity.", call. = FALSE)
  }

  n_payments <- n_prem * k

  payment_times <- if (premium_timing == "due") {
    premium_start_time + (0:(n_payments - 1L)) / k
  } else {
    premium_start_time + (1:n_payments) / k
  }

  apv_premiums <- sum(vapply(payment_times, function(tt) {
    (v^tt) * status_survival(tt)
  }, numeric(1)))

  if (!is.finite(apv_premiums) || apv_premiums <= 0) {
    stop("APV of premium annuity is nonpositive or not finite.", call. = FALSE)
  }

  premium <- apv_benefits / apv_premiums

  if (!isTRUE(tidy)) return(premium)

  tibble::tibble(
    x = x,
    y = y,
    m = m,
    n = if (is.null(n)) NA_integer_ else n,
    product = product,
    cohort = cohort,
    benefit = benefit,
    k = k,
    frac = frac,
    premium_timing = premium_timing,
    prem_start = prem_start,
    n_prem = n_prem,
    premium = premium,
    premium_annual = k * premium,
    apv_benefits = apv_benefits,
    apv_premiums = apv_premiums
  )
}
