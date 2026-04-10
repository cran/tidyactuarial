#' Actuarial present value of a multi-life annuity (n lives)
#'
#' Computes the APV of a discrete annuity contingent on multiple independent lives.
#' Supports status-based annuities (joint-life / last-survivor) and a
#' joint-and-survivor style annuity ("reversionary") that pays 1 while all
#' lives are alive and then pays a fraction \eqn{\alpha} while at least one
#' life remains alive.
#'
#' @param lt Life table with column \code{x} and at least one of \code{lx},
#'   \code{px}, \code{qx}.
#' @param ages Integer vector of actuarial ages for the lives at issue.
#' @param annuity Type of annuity logic: \code{"cohort"} uses the status
#'   defined in \code{cohort}; \code{"reversionary"} uses the \eqn{\alpha}
#'   fractional reduction.
#' @param cohort Survival status: \code{"first"} (joint-life, pays while all
#'   are alive) or \code{"last"} (last-survivor, pays while at least one
#'   remains alive). Used only when \code{annuity = "cohort"}.
#' @param alpha Reversionary fraction (typically \eqn{0 \le \alpha \le 1}).
#'   Used only when \code{annuity = "reversionary"}.
#'   Note: \code{alpha = 0} matches joint-life; \code{alpha = 1}
#'   matches last-survivor.
#' @param n Integer term in years after deferment. If \code{NULL}, runs to
#'   the end of the table.
#' @param m Integer deferment in years.
#' @param k Integer payments per year. If \code{k > 1}, Woolhouse approximations
#'   may be applied.
#' @param timing \code{"immediate"} or \code{"due"}.
#' @param woolhouse \code{"none"}, \code{"first"}, or \code{"second"}.
#' @param i Annual effective interest rate.
#'
#' @details
#' Under the assumption of independent future lifetimes (Finan, Section 51),
#' the survival probability for the status is calculated as:
#' \itemize{
#'   \item \strong{Joint-life (first-death):}
#'     \eqn{{}_t p_{x_1 x_2 \dots x_n} = \prod_{j=1}^n {}_t p_{x_j}}
#'   \item \strong{Last-survivor:}
#'     \eqn{{}_t p_{\overline{x_1 x_2 \dots x_n}} = 1 - \prod_{j=1}^n (1 - {}_t p_{x_j})}
#' }
#'
#' For \code{annuity = "reversionary"}, the APV is a weighted combination
#' of the two statuses (Finan, Section 53.3):
#' \deqn{APV = APV(\text{joint-life}) + \alpha [APV(\text{last-survivor}) - APV(\text{joint-life})]}
#'
#' @return A single numeric value representing the APV.
#'
#' @seealso \code{\link{annuity_x}} for single-life annuities,
#'   \code{\link{t_px}} for survival probabilities.
#' @export
annuity_multi <- function(
    lt, ages,
    annuity = c("cohort", "reversionary"),
    cohort = c("first", "last"),
    alpha = NULL,
    n = NULL,
    m = 0L,
    k = 1L,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    i
) {
  annuity <- match.arg(annuity)
  cohort <- match.arg(cohort)
  timing <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)

  if (missing(i) || !is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) {
    stop("'i' must be a single numeric rate greater than -1.")
  }
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!("x" %in% names(lt))) stop("Life table must contain column 'x'.")
  if (!("lx" %in% names(lt)) && !("px" %in% names(lt)) && !("qx" %in% names(lt))) {
    stop("Life table must contain 'lx', 'px', or 'qx'.")
  }

  if (!is.numeric(ages) || length(ages) < 1) stop("'ages' must be a non-empty numeric vector.")
  ages <- as.integer(round(ages))

  m <- as.integer(round(m))
  k <- as.integer(round(k))
  if (m < 0) stop("'m' must be nonnegative.")
  if (k < 1) stop("'k' must be >= 1.")

  if (annuity == "reversionary") {
    if (is.null(alpha) || !is.numeric(alpha) || length(alpha) != 1L || is.na(alpha)) {
      stop("'alpha' must be a single numeric value.")
    }
  }

  lt <- lt[order(lt$x), ]
  omega <- max(lt$x, na.rm = TRUE)

  # term after deferment
  if (is.null(n)) {
    n <- max(0L, omega - max(ages) - m)
  } else {
    n <- as.integer(round(n))
    if (n < 0) stop("'n' must be nonnegative or NULL.")
  }
  if (n == 0L) return(0)

  v <- function(t) (1 + i)^(-t)

  # integer-year survival for one life: {}_t p_age
  t_px1 <- function(age, t) {
    if (t == 0) return(1)

    if ("lx" %in% names(lt)) {
      l0 <- lt$lx[match(age, lt$x)]
      l1 <- lt$lx[match(age + t, lt$x)]
      if (is.na(l0) || is.na(l1) || l0 <= 0) return(NA_real_)
      return(l1 / l0)
    }

    if ("px" %in% names(lt)) {
      pxv <- lt$px[match(age:(age + t - 1), lt$x)]
      if (anyNA(pxv)) return(NA_real_)
      return(prod(pxv))
    }

    qxv <- lt$qx[match(age:(age + t - 1), lt$x)]
    if (anyNA(qxv)) return(NA_real_)
    prod(1 - qxv)
  }

  # Prob(all alive at time u) and Prob(any alive at time u), from ISSUE ages
  P_all_u <- function(u) {
    p_vec <- vapply(ages, function(a) t_px1(a, u), numeric(1))
    if (anyNA(p_vec)) return(NA_real_)
    prod(p_vec)
  }
  P_any_u <- function(u) {
    p_vec <- vapply(ages, function(a) t_px1(a, u), numeric(1))
    if (anyNA(p_vec)) return(NA_real_)
    1 - prod(1 - p_vec)
  }

  # Expected payment at absolute time u (years since issue)
  E_pay_u <- function(u) {
    pa <- P_all_u(u)
    pn <- P_any_u(u)
    if (is.na(pa) || is.na(pn)) return(NA_real_)

    if (annuity == "cohort") {
      if (cohort == "first") return(pa)
      return(pn)
    }

    # reversionary: 1 while all alive; alpha while partially alive (at least one alive)
    pa + alpha * (pn - pa)
  }

  # Annual APV (k=1 base), payments at u = m + t, t integer times depending on timing
  annual_apv <- function(n, timing) {
    tt <- if (timing == "due") 0:(n - 1L) else 1:n
    u <- m + tt
    ep <- vapply(u, E_pay_u, numeric(1))
    if (anyNA(ep)) stop("Life table does not support required ages for this annuity.")
    sum(v(u) * ep)
  }

  # --- exact annual or exact k-thly not implemented here (you are using Woolhouse/annual base) ---
  if (k == 1L || woolhouse == "none") {
    return(annual_apv(n, timing))
  }

  # Woolhouse: adjust annual annuity-due value; scale by PV of expected first (due) payment at time m
  adue <- annual_apv(n, "due")

  # PV of expected payment at the first due date (time m)
  pv_first <- v(m) * E_pay_u(m)

  adj1 <- (k - 1) / (2 * k)

  if (woolhouse == "first") {
    adue_k <- adue - pv_first * adj1
  } else {
    delta <- log(1 + i)

    # Approximate mu of the STATUS around time m using one-year survival ratio:
    # p_status1 \approx status_prob(m+1) / status_prob(m)
    p_m  <- E_pay_u(m)   # expected payment at m (proportional to being "in force")
    p_m1 <- E_pay_u(m + 1L)

    if (!is.finite(p_m) || p_m <= 0 || !is.finite(p_m1) || p_m1 < 0) {
      stop("Cannot compute second-order Woolhouse adjustment (status probability at m is nonpositive).")
    }

    p_status1 <- p_m1 / p_m
    if (p_status1 <= 0) stop("Cannot compute mu: derived status one-year survival <= 0.")

    mu <- -log(p_status1)

    adj2 <- (k^2 - 1) / (12 * k^2) * (delta + mu)

    adue_k <- adue - pv_first * (adj1 + adj2)
  }

  if (timing == "due") return(adue_k)

  # Convert due -> immediate (approx): remove 1/k of the first due payment PV
  adue_k - pv_first * (1 / k)
}
