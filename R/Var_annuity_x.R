#' Variance of the actuarial present value of a life annuity
#'
#' Computes the variance of the present value random variable of a discrete
#' life annuity with \eqn{k}-thly payments, using the exact second moment
#' under the UDD (Uniform Distribution of Deaths) assumption.
#'
#' @param lt A life table data frame containing columns \code{x} and
#'   \code{lx}. Requires \code{lx} to apply exact UDD at fractional ages.
#' @param x Integer actuarial age.
#' @param i Effective annual interest rate (must be \code{> -1}).
#' @param n Integer term in years. If \code{NULL}, whole life to end of table.
#' @param m Integer deferment in years (default \code{0}).
#' @param k Integer payments per year (default \code{1}). Example:
#'   \code{k = 12} for monthly.
#' @param timing Payment timing: \code{"immediate"} (arrears) or
#'   \code{"due"} (advance).
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble with details.
#'
#' @details
#' Let payments be of amount \eqn{1/k} at times \eqn{m + u_j}, where
#' \eqn{u_j = j/k} with \eqn{j = 0, \ldots, kn-1} (due) or
#' \eqn{j = 1, \ldots, kn} (immediate).
#'
#' Define \eqn{Z = \sum_j \frac{1}{k} v^{m+u_j} I(T_{x+m} > u_j)}.
#' Then:
#' \deqn{\mathrm{Var}(Z) = E[Z^2] - (E[Z])^2}
#' with
#' \deqn{E[Z^2] = \sum_j \sum_\ell \frac{1}{k^2}
#'   v^{2m + u_j + u_\ell} \cdot {}_{\max(u_j, u_\ell)}p_{x+m}.}
#'
#' This uses the fact that payments \eqn{j} and \eqn{\ell} are both made
#' if and only if the annuitant survives to \eqn{\max(u_j, u_\ell)}.
#' The double sum is computed via a vectorized matrix product
#' (\code{outer}) for efficiency.
#'
#' For the special case of annual payments (\code{k = 1}) without deferral,
#' the variance satisfies the classical identity (Finan, Section 37.1,
#' Example 37.5):
#' \deqn{\mathrm{Var}(\ddot{Y}_{x:\overline{n}|}) =
#'   \frac{{}^2A_{x:\overline{n}|} - (A_{x:\overline{n}|})^2}{d^2}}
#' where \eqn{d = i/(1+i)} and \eqn{{}^2A} is the second moment of the
#' endowment insurance (Finan, Section 27). This identity can be used to
#' cross-validate results.
#'
#' Fractional survival probabilities are exact under UDD
#' (Finan, Section 24.1).
#'
#' @return A numeric variance, or a one-row tibble if \code{tidy = TRUE}
#'   with columns \code{x}, \code{n}, \code{m}, \code{k}, \code{timing},
#'   \code{i}, \code{EZ}, \code{EZ2}, \code{variance}.
#'
#' @seealso \code{\link{annuity_x}} for the expected APV (first moment),
#'   \code{\link{Var_insurance_x}} for the variance of life insurance APV,
#'   \code{\link{insurance_x}} for life insurance APVs.
#'
#' @examples
#' lt <- data.frame(x = 60:110, lx = seq(100000, 0, length.out = 51))
#'
#' # Annual variance, annuity-due
#' Var_annuity_x(lt, x = 60, i = 0.06, timing = "due")
#'
#' # Monthly (k=12) exact under UDD for a 10-year term
#' Var_annuity_x(lt, x = 60, i = 0.06, n = 10, k = 12, timing = "immediate")
#'
#' # Cross-validate with Finan identity (Sec. 37.1):
#' # Var(\ddot{Y}_x) = (2A_x - A_x^2) / d^2
#' v_out <- Var_annuity_x(lt, x = 60, i = 0.06, timing = "due", tidy = TRUE)
#' v_out
#'
#' # 5-year temporary with tidy output
#' Var_annuity_x(lt, x = 60, i = 0.06, n = 5, k = 4,
#'               timing = "due", tidy = TRUE)
#'
#' # Deferred annuity variance
#' Var_annuity_x(lt, x = 60, i = 0.06, m = 5, timing = "due")
#'
#' @export
Var_annuity_x <- function(
    lt,
    x,
    i,
    n = NULL,
    m = 0L,
    k = 1L,
    timing = c("immediate", "due"),
    tidy = FALSE
) {
  timing <- match.arg(timing)

  # --- checks ---
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) stop("`lt` must contain columns `x` and `lx`.")
  if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) stop("'i' must be > -1.")
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || abs(x - round(x)) > 1e-10) stop("'x' must be an integer.")
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 0 || abs(m - round(m)) > 1e-10) stop("'m' must be a nonnegative integer.")
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 || abs(k - round(k)) > 1e-10) stop("'k' must be a positive integer.")
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || abs(n - round(n)) > 1e-10) {
      stop("'n' must be a nonnegative integer or NULL.")
    }
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  k <- as.integer(round(k))
  if (!is.null(n)) n <- as.integer(round(n))

  lt <- lt[order(lt$x), ]
  if (anyDuplicated(lt$x)) stop("Life table ages `x` must be unique.")
  if (any(is.na(lt$lx)) || any(lt$lx < 0)) stop("`lx` must be nonnegative and not NA.")

  ages <- as.integer(lt$x)
  lx   <- as.numeric(lt$lx)
  omega <- max(ages)
  y <- x + m

  # Default n (whole life)
  max_years <- max(0L, (omega + 1L) - y)
  if (is.null(n)) {
    n <- max_years
  } else {
    if (n > max_years) stop("`n` exceeds the horizon allowed by the life table.")
  }

  if (n == 0L) {
    if (!isTRUE(tidy)) return(0)
    return(tibble::tibble(x = x, n = n, m = m, k = k, timing = timing,
                          i = i, EZ = 0, EZ2 = 0, variance = 0))
  }

  # --- helpers ---
  get_lx <- function(age) {
    idx <- match(age, ages)
    if (!is.na(idx)) return(lx[idx])
    if (age == (omega + 1L)) return(0)
    NA_real_
  }

  t_p_int <- function(age, tt) {
    if (tt == 0) return(1)
    l0 <- get_lx(age)
    l1 <- get_lx(age + tt)
    if (is.na(l0) || is.na(l1) || l0 <= 0) return(NA_real_)
    l1 / l0
  }

  t_p_udd <- function(age, u) {
    if (u < 0) return(NA_real_)
    if (u == 0) return(1)
    tt <- floor(u)
    s  <- u - tt
    pt <- t_p_int(age, tt)
    if (is.na(pt)) return(NA_real_)
    if (s == 0) return(pt)
    yy  <- age + tt
    ly  <- get_lx(yy)
    ly1 <- get_lx(yy + 1L)
    if (is.na(ly) || is.na(ly1) || ly <= 0) return(NA_real_)
    dy  <- ly - ly1
    ps  <- (ly - s * dy) / ly
    pt * ps
  }

  v_fun <- function(tt) (1 + i)^(-tt)

  # --- payment times ---
  N <- n * k
  j <- if (timing == "due") 0:(N - 1L) else 1:N
  u <- j / k

  # --- survival at each payment time from age y ---
  p_u <- vapply(u, function(uu) t_p_udd(y, uu), numeric(1))
  if (anyNA(p_u)) stop("Life table does not support required ages for UDD k-thly valuation.")

  # --- E[Z]: first moment ---
  v_vec <- v_fun(m + u)
  EZ <- sum((1 / k) * v_vec * p_u)

  # --- E[Z^2]: second moment via vectorized outer product ---
  # Precompute survival at all steps r/k for r = 0..max(j)
  r_max <- max(j)
  p_r <- vapply(0:r_max, function(r) t_p_udd(y, r / k), numeric(1))
  if (anyNA(p_r)) stop("Life table does not support required ages for UDD second moment.")

  # P_mat[a,b] = survival to later of two payment times
  # V_mat[a,b] = product of two discount factors
  P_mat <- p_r[outer(j, j, pmax) + 1L]
  V_mat <- outer(v_vec, v_vec)
  EZ2   <- sum(V_mat * P_mat) / k^2

  var_val <- EZ2 - EZ^2
  # Numerical guard for floating-point noise
  if (var_val < 0 && var_val > -1e-12) var_val <- 0

  if (!isTRUE(tidy)) return(var_val)

  tibble::tibble(
    x = x, n = n, m = m, k = k, timing = timing, i = i,
    EZ = EZ, EZ2 = EZ2, variance = var_val
  )
}
