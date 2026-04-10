#' Net premium for life insurance by the equivalence principle
#'
#' Computes the net benefit premium of a life insurance contract using the
#' equivalence principle (Finan, Section 41):
#' \deqn{P = \frac{\text{APV of benefits}}{\text{APV of premium annuity}}.}
#'
#' The premium returned corresponds to one premium payment
#' (annual if \code{k = 1}, monthly if \code{k = 12}, etc.).
#'
#' @param lt A life table data frame containing at least columns \code{x}
#'   and \code{lx}.
#' @param x Integer actuarial age at issue.
#' @param i Effective annual interest rate (must be \code{> -1}).
#' @param product Type of insurance: \code{"whole"}, \code{"term"},
#'   \code{"endowment"}, or \code{"variable_k"}.
#' @param benefit Benefit amount. For standard products, a single numeric
#'   value. For \code{"variable_k"}, a numeric vector or a function of time.
#' @param n Optional insurance term in years. Required for \code{"term"}
#'   and \code{"endowment"}.
#' @param m Nonnegative integer deferral period in years (default \code{0}).
#' @param k Number of premium payments per year (default \code{1}).
#' @param frac Fractional-age assumption for survival probabilities:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#' @param premium_timing Timing of premium payments:
#'   \code{"due"} (in advance, default) or \code{"immediate"} (in arrears).
#' @param prem_start Start of premium payments:
#'   \code{"issue"} (start at issue, time 0) or \code{"deferred"} (start at
#'   time \code{m}).
#' @param n_prem Optional premium-paying term in years, counted from
#'   \code{prem_start}. If \code{NULL}, defaults to whole-life to the end of
#'   the table (for \code{"whole"}) or to the insurance term \code{n}
#'   (for temporary products).
#' @param woolhouse Woolhouse order for the premium annuity when \code{k > 1}:
#'   \code{"none"} (exact UDD, default), \code{"first"}, or \code{"second"}.
#'   Passed to \code{\link{annuity_x}}.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble with details.
#' @param check Logical. If \code{TRUE}, performs basic input validation.
#'
#' @details
#' The benefit premium is the level payment that satisfies the equivalence
#' principle: the APV of premiums equals the APV of benefits at issue.
#'
#' For standard products, the APV of benefits is (Finan, Sections 27 and 41):
#' \itemize{
#'   \item Whole life: \eqn{A_x} (Finan, Sec. 41.1)
#'   \item n-year term: \eqn{A^1_{x:\overline{n}|}} (Finan, Sec. 41.2)
#'   \item n-year endowment: \eqn{A_{x:\overline{n}|}} (Finan, Sec. 41.4)
#' }
#'
#' The APV of the premium annuity is computed via \code{\link{annuity_x}},
#' supporting \code{k}-thly payments and Woolhouse approximations.
#'
#' **Premium payment start** (Finan, Sec. 41.5):
#' \itemize{
#'   \item \code{prem_start = "issue"}: premiums start at age \eqn{x}
#'     (time 0). The premium annuity is \eqn{\ddot{a}_{x:\overline{n_p}|}}.
#'   \item \code{prem_start = "deferred"}: premiums start at age \eqn{x+m}
#'     (after the deferral period). The premium annuity, valued at time 0,
#'     is \eqn{v^m \cdot {}_mp_x \cdot \ddot{a}_{x+m:\overline{n_p}|}}
#'     which equals a deferred annuity with \code{m} years of deferral.
#' }
#'
#' @return A numeric net premium per payment, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{insurance_x}} for the APV of benefits,
#'   \code{\link{annuity_x}} for the APV of the premium annuity,
#'   \code{\link{premium_xy}} for two-life premiums,
#'   \code{\link{premium_gross}} for gross premiums with expenses.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Whole life insurance, annual premium (Finan, Sec. 41.1):
#' # P(A_x) = A_x / \ddot{a}_x
#' premium_x(lt, x = 60, i = 0.05, product = "whole", benefit = 100000)
#'
#' # Verify manually: P = A / \ddot{a}
#' A  <- insurance_x(lt, x = 60, i = 0.05, type = "whole")
#' ad <- annuity_x(lt, x = 60, i = 0.05, timing = "due")
#' 100000 * A / ad
#'
#' # 5-year term insurance, annual premium (Finan, Sec. 41.2):
#' # P(A^1_{x:n}) = A^1_{x:n} / \ddot{a}_{x:n}
#' premium_x(lt, x = 60, i = 0.05, product = "term", n = 5, benefit = 100000)
#'
#' # 5-year endowment insurance (Finan, Sec. 41.4):
#' # P(A_{x:n}) = A_{x:n} / \ddot{a}_{x:n}
#' premium_x(lt, x = 60, i = 0.05, product = "endowment", n = 5,
#'           benefit = 100000)
#'
#' # Deferred whole life: premiums start after deferral (Finan, Sec. 41.5)
#' premium_x(lt, x = 60, i = 0.05, product = "whole", benefit = 100000,
#'           m = 2, prem_start = "deferred")
#'
#' # t-payment whole life: 3 years of premiums for whole life coverage
#' premium_x(lt, x = 60, i = 0.05, product = "whole", benefit = 100000,
#'           n_prem = 3)
#'
#' # Monthly premiums with Woolhouse
#' premium_x(lt, x = 60, i = 0.05, product = "whole", benefit = 100000,
#'           k = 12, woolhouse = "first")
#'
#' # Tidy output with all details
#' premium_x(lt, x = 60, i = 0.05, product = "term", n = 5,
#'           benefit = 100000, tidy = TRUE)
#'
#' @export
premium_x <- function(
    lt,
    x,
    i,
    product = c("whole", "term", "endowment", "variable_k"),
    benefit,
    n = NULL,
    m = 0,
    k = 1,
    frac = c("UDD", "CF", "CML", "Balducci"),
    premium_timing = c("due", "immediate"),
    prem_start = c("issue", "deferred"),
    n_prem = NULL,
    woolhouse = c("none", "first", "second"),
    tidy = FALSE,
    check = TRUE
) {
  product        <- match.arg(product)
  frac           <- match.arg(frac)
  premium_timing <- match.arg(premium_timing)
  prem_start     <- match.arg(prem_start)
  woolhouse      <- match.arg(woolhouse)

  if (isTRUE(check)) {
    if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
    if (!all(c("x", "lx") %in% names(lt))) {
      stop("Life table must contain columns 'x' and 'lx'.")
    }
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || abs(x - round(x)) > 1e-10) {
      stop("'x' must be a single integer age.")
    }
    if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) {
      stop("'i' must be a single numeric value > -1.")
    }
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0 || abs(k - round(k)) > 1e-10) {
      stop("'k' must be a single positive integer.")
    }
    if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 0 || abs(m - round(m)) > 1e-10) {
      stop("'m' must be a single nonnegative integer.")
    }
    if (product %in% c("term", "endowment") && is.null(n)) {
      stop("'n' must be provided for term or endowment insurance.")
    }
    if (missing(benefit)) stop("'benefit' must be provided.")
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  k <- as.integer(round(k))
  if (!is.null(n)) n <- as.integer(round(n))
  if (!is.null(n_prem)) n_prem <- as.integer(round(n_prem))

  # --------------------------------------------------
  # 1) APV of benefits (valued at time 0, age x)
  # --------------------------------------------------
  if (product == "variable_k") {

    apv_benefits <- insurance_variable_k(
      lt = lt, x = x, i = i,
      benefit = benefit,
      n = n, m = m, k = k, frac = frac
    )

    if (is.null(n_prem)) {
      if (!is.null(n)) {
        n_prem <- n
      } else if (!is.function(benefit)) {
        n_prem <- as.integer(length(benefit) / k)
      } else {
        stop("Cannot infer 'n_prem' for variable benefits supplied as function.")
      }
    }

  } else {

    if (!is.numeric(benefit) || length(benefit) != 1L || is.na(benefit) || benefit < 0) {
      stop("For standard products, 'benefit' must be a single nonnegative number.")
    }

    # insurance_x returns APV per unit benefit, valued at time 0
    apv_unit <- insurance_x(
      lt = lt, x = x, i = i,
      n = n, m = m,
      type = product
    )

    apv_benefits <- benefit * apv_unit

    # Default premium-paying term
    if (is.null(n_prem)) {
      if (product == "whole") {
        x_start <- if (prem_start == "issue") x else (x + m)
        n_prem <- max(lt$x, na.rm = TRUE) - x_start
      } else {
        n_prem <- n
      }
    }
  }

  # --------------------------------------------------
  # 2) APV of premium annuity (valued at time 0, age x)
  #    Using annuity_x for consistency and Woolhouse support
  # --------------------------------------------------
  # Both APVs must be at the same valuation point (time 0, age x).
  #
  # prem_start = "issue":
  #   Premiums start at age x -> annuity_x(x, n=n_prem, m=0, ...)
  #
  # prem_start = "deferred":
  #   Premiums start at age x+m -> annuity_x(x, n=n_prem, m=m, ...)
  #   This includes the deferral factor v^m * m_p_x automatically.

  m_prem <- if (prem_start == "issue") 0L else m

  apv_prem <- annuity_x(
    lt = lt, x = x, i = i,
    n = n_prem, m = m_prem, k = k,
    timing = premium_timing,
    woolhouse = woolhouse
  )

  if (!is.finite(apv_prem) || apv_prem <= 0) {
    stop("APV of premium annuity is nonpositive or not finite.")
  }

  # P per payment = APV(benefits) / APV(premium annuity)
  premium <- apv_benefits / apv_prem

  if (!isTRUE(tidy)) return(premium)

  tibble::tibble(
    x              = x,
    m              = m,
    n              = n,
    product        = product,
    benefit        = if (is.function(benefit)) NA_real_ else suppressWarnings(as.numeric(benefit)[1]),
    k              = k,
    frac           = frac,
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
