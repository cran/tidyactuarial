#' Net premium for two-life insurance by the equivalence principle
#'
#' Computes the net benefit premium of a two-life insurance contract
#' (independent lives) using the equivalence principle
#' (Finan, Sections 41 and 58--59):
#' \deqn{P = \frac{\text{APV of benefits}}{\text{APV of premium annuity}}.}
#'
#' The premium returned corresponds to one premium payment
#' (annual if \code{k = 1}, monthly if \code{k = 12}, etc.).
#'
#' @param lt A life table data frame with columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for life 1 at issue.
#' @param y Integer actuarial age for life 2 at issue.
#' @param i Annual effective interest rate (must be \code{> -1}).
#' @param type Insurance type: \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param cohort Status cohort (benefit trigger): \code{"first"}
#'   (joint-life) or \code{"last"} (last-survivor).
#' @param benefit Benefit amount (single nonnegative numeric value).
#' @param n Optional insurance term in years. Required for \code{"term"}
#'   and \code{"endowment"}.
#' @param m Nonnegative integer deferral period in years (default \code{0}).
#' @param k Number of premium payments per year (default \code{1}).
#' @param premium_timing Timing of premium payments: \code{"due"}
#'   (in advance, default) or \code{"immediate"} (in arrears).
#' @param prem_start Start of premium payments: \code{"issue"} (start at
#'   time 0) or \code{"deferred"} (start at time \code{m}).
#' @param n_prem Optional premium-paying term in years, counted from
#'   \code{prem_start}. If \code{NULL}, defaults to whole-life for
#'   \code{"whole"} or insurance term \code{n} for temporary products.
#' @param woolhouse Woolhouse order for the premium annuity when
#'   \code{k > 1}: \code{"none"} (exact k-thly, default), \code{"first"},
#'   or \code{"second"}. Passed to \code{\link{annuity_xy}}.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#' @param check Logical. If \code{TRUE}, performs basic input validation.
#'
#' @details
#' The APV of benefits is computed via \code{\link{insurance_xy}},
#' valued at issue (time 0, age \code{x}).
#'
#' The APV of the premium annuity is computed via
#' \code{\link{annuity_xy}}, also valued at time 0. Both APVs are at the
#' same valuation point, ensuring the equivalence principle is applied
#' correctly.
#'
#' **Premium payment start**:
#' \itemize{
#'   \item \code{prem_start = "issue"}: premiums start at time 0 (ages
#'     \code{x}, \code{y}). The premium annuity is
#'     \eqn{\ddot{a}_{xy:\overline{n_p}|}}.
#'   \item \code{prem_start = "deferred"}: premiums start at time \code{m}
#'     (ages \code{x+m}, \code{y+m}). The premium annuity, valued at
#'     time 0, is a deferred annuity with \code{m} years of deferral:
#'     \eqn{v^m \cdot {}_mp_{xy} \cdot
#'     \ddot{a}_{x+m,y+m:\overline{n_p}|}}.
#' }
#'
#' The premium annuity uses the same two-life status as the benefit.
#'
#' @return A numeric net premium per payment, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{insurance_xy}} for the APV of two-life benefits,
#'   \code{\link{annuity_xy}} for the APV of the premium annuity,
#'   \code{\link{premium_x}} for single-life premiums,
#'   \code{\link{premium_gross}} for gross premiums with expenses.
#'
#' @examples
#' lt <- data.frame(x = 60:110, lx = seq(100000, 0, length.out = 51))
#'
#' # Joint-life whole insurance, annual premium from issue
#' premium_xy(lt, x = 60, y = 62, i = 0.05,
#'            type = "whole", cohort = "first",
#'            benefit = 100000)
#'
#' # 4-year term, last-survivor
#' premium_xy(lt, x = 60, y = 62, i = 0.05,
#'            type = "term", cohort = "last",
#'            n = 4, benefit = 100000)
#'
#' @export
premium_xy <- function(
    lt, x, y, i,
    type = c("whole", "term", "endowment"),
    cohort = c("first", "last"),
    benefit,
    n = NULL,
    m = 0,
    k = 1,
    premium_timing = c("due", "immediate"),
    prem_start = c("issue", "deferred"),
    n_prem = NULL,
    woolhouse = c("none", "first", "second"),
    tidy = FALSE,
    check = TRUE
) {
  type           <- match.arg(type)
  cohort         <- match.arg(cohort)
  premium_timing <- match.arg(premium_timing)
  prem_start     <- match.arg(prem_start)
  woolhouse      <- match.arg(woolhouse)

  if (isTRUE(check)) {
    if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
    if (!all(c("x", "lx") %in% names(lt))) {
      stop("Life table must contain columns 'x' and 'lx'.")
    }
    if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
        abs(x - round(x)) > 1e-10) {
      stop("'x' must be a single integer age.")
    }
    if (!is.numeric(y) || length(y) != 1L || is.na(y) ||
        abs(y - round(y)) > 1e-10) {
      stop("'y' must be a single integer age.")
    }
    if (!is.numeric(i) || length(i) != 1L || is.na(i) ||
        i <= -1) {
      stop("'i' must be a single numeric value > -1.")
    }
    if (!is.numeric(k) || length(k) != 1L || is.na(k) ||
        k <= 0 || abs(k - round(k)) > 1e-10) {
      stop("'k' must be a single positive integer.")
    }
    if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
        m < 0 || abs(m - round(m)) > 1e-10) {
      stop("'m' must be a single nonnegative integer.")
    }
    if (type %in% c("term", "endowment") && is.null(n)) {
      stop("'n' required for term or endowment insurance.")
    }
    if (missing(benefit)) stop("'benefit' must be provided.")
    if (!is.numeric(benefit) || length(benefit) != 1L ||
        is.na(benefit) || benefit < 0) {
      stop("'benefit' must be a single nonnegative number.")
    }
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))
  k <- as.integer(round(k))
  if (!is.null(n)) n <- as.integer(round(n))
  if (!is.null(n_prem)) n_prem <- as.integer(round(n_prem))

  # --------------------------------------------------
  # 1) APV of benefits (valued at time 0)
  # --------------------------------------------------
  apv_unit <- insurance_xy(
    lt = lt, x = x, y = y, i = i,
    n = n, m = m, type = type, cohort = cohort
  )
  apv_benefits <- benefit * apv_unit

  # --------------------------------------------------
  # 2) Default premium-paying term
  # --------------------------------------------------
  if (is.null(n_prem)) {
    omega <- max(lt$x, na.rm = TRUE)
    if (type == "whole") {
      x_start <- if (prem_start == "issue") x else (x + m)
      y_start <- if (prem_start == "issue") y else (y + m)
      # Horizon depends on status
      if (cohort == "first") {
        n_prem <- max(0L, omega - max(x_start, y_start))
      } else {
        n_prem <- max(0L, omega - min(x_start, y_start))
      }
    } else {
      n_prem <- n
    }
  }

  # --------------------------------------------------
  # 3) APV of premium annuity (valued at time 0)
  #    Delegate to annuity_xy for consistency & Woolhouse
  # --------------------------------------------------
  # prem_start = "issue"  -> annuity at (x,y) with m=0
  # prem_start = "deferred" -> annuity at (x,y) with m=m
  #   (annuity_xy includes deferral factor automatically)
  m_prem <- if (prem_start == "issue") 0L else m

  apv_prem <- annuity_xy(
    lt = lt, x = x, y = y, i = i,
    cohort = cohort, n = n_prem, m = m_prem, k = k,
    timing = premium_timing, woolhouse = woolhouse
  )

  if (!is.finite(apv_prem) || apv_prem <= 0) {
    stop("APV of premium annuity is nonpositive or not finite.")
  }

  # P per payment = APV(benefits) / APV(premium annuity)
  premium <- apv_benefits / apv_prem

  if (!isTRUE(tidy)) return(premium)

  tibble::tibble(
    x              = x,
    y              = y,
    m              = m,
    n              = n,
    type           = type,
    cohort         = cohort,
    benefit        = benefit,
    k              = k,
    premium_timing = premium_timing,
    prem_start     = prem_start,
    n_prem         = n_prem,
    woolhouse      = woolhouse,
    premium        = premium,
    premium_annual = k * premium,
    apv_benefits   = apv_benefits,
    apv_premiums   = apv_prem
  )
}
