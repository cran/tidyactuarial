#' Variance of the actuarial present value of a life insurance
#'
#' Computes the variance of the present value random variable of a discrete
#' life insurance under the exact UDD assumption on a \eqn{1/k}-year grid.
#'
#' @param lt Life table data frame containing columns \code{x} and \code{lx}.
#'   Requires \code{lx} to apply exact UDD at fractional ages.
#' @param x Integer actuarial age at issue.
#' @param i Effective annual interest rate (must be \code{> -1}).
#' @param product Insurance type: \code{"whole"}, \code{"term"},
#'   \code{"endowment"}, or \code{"pure_endowment"}.
#' @param benefit Numeric scalar benefit amount (default \code{1}).
#' @param n Integer term in years after deferment. Required for
#'   \code{"term"}, \code{"endowment"}, and \code{"pure_endowment"}.
#' @param m Integer deferment in years (default \code{0}).
#' @param k Integer grid frequency per year (default \code{1}). Example:
#'   \code{k = 12} for monthly benefit timing.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble with details.
#'
#' @details
#' The benefit is paid at the end of the \eqn{1/k}-year subperiod in which
#' death occurs (or at time \eqn{m + n} for the pure endowment).
#'
#' Under UDD (Finan, Section 24.1), fractional survival is computed by
#' linear interpolation of \eqn{\ell_{x+s}} within each year. The
#' death-in-subinterval probability is:
#' \deqn{\Pr(T \in (u_{r-1}, u_r]) = S(u_{r-1}) - S(u_r)}
#' where \eqn{S(u) = {}_up_{x+m}}.
#'
#' **Whole life insurance** (Finan, Section 27):
#' \deqn{\mathrm{Var}(Z_x) = {}^2A_x - (A_x)^2}
#' where \eqn{{}^2A_x = \sum_{k=0}^{\infty} v^{2(k+1)} \cdot {}_kp_x
#' \cdot q_{x+k}}.
#'
#' **Term insurance** (Finan, Section 27):
#' \deqn{\mathrm{Var}(Z^1_{x:\overline{n}|}) =
#'   {}^2A^1_{x:\overline{n}|} - (A^1_{x:\overline{n}|})^2.}
#'
#' **Pure endowment** (Finan, Section 26.3.1):
#' \deqn{\mathrm{Var}(\bar{Z}^{\phantom{1}}_{\phantom{1}x:\overline{n}|})
#'   = v^{2n} \cdot {}_np_x \cdot {}_nq_x.}
#'
#' **Endowment insurance** (Finan, Example 26.15 and Section 26.3.2):
#' Since the term component and the pure endowment are mutually exclusive
#' (\eqn{Z^1_{x:\overline{n}|} \cdot
#' \bar{Z}^{\phantom{1}}_{\phantom{1}x:\overline{n}|} = 0}), the
#' covariance is \eqn{-A^1_{x:\overline{n}|} \cdot A^{\phantom{1}}_{\phantom{1}x:\overline{n}|}}
#' and therefore:
#' \deqn{\mathrm{Var}(Z_{x:\overline{n}|}) =
#'   {}^2A^1_{x:\overline{n}|} + {}^2A^{\phantom{1}}_{\phantom{1}x:\overline{n}|}
#'   - (A^1_{x:\overline{n}|} + A^{\phantom{1}}_{\phantom{1}x:\overline{n}|})^2.}
#'
#' For the k-thly case, the sums run over all \eqn{1/k}-year subperiods,
#' giving the exact UDD result.
#'
#' @return Numeric variance, or a one-row tibble if \code{tidy = TRUE}
#'   with columns \code{x}, \code{m}, \code{n}, \code{k}, \code{product},
#'   \code{i}, \code{benefit}, \code{EZ}, \code{EZ2}, \code{variance}.
#'
#' @seealso \code{\link{insurance_x}} for the expected APV (first moment),
#'   \code{\link{Var_annuity_x}} for the variance of life annuity APV,
#'   \code{\link{annuity_x}} for life annuity APVs.
#'
#' @examples
#' lt <- data.frame(x = 60:110, lx = seq(100000, 0, length.out = 51))
#'
#' # Whole life, annual
#' Var_insurance_x(lt, x = 60, i = 0.06, product = "whole")
#'
#' # 10-year term, monthly (k=12)
#' Var_insurance_x(lt, x = 60, i = 0.06, product = "term", n = 10, k = 12)
#'
#' # Finan Section 27 style: Var(Z) = 2A - A^2
#' # Verify by comparing tidy output EZ2 - EZ^2
#' Var_insurance_x(lt, x = 60, i = 0.06, product = "whole", tidy = TRUE)
#'
#' # 10-year term with 5-year deferral, monthly
#' Var_insurance_x(lt, x = 60, i = 0.06, product = "term",
#'                  n = 10, m = 5, k = 12)
#'
#' # Pure endowment variance (Finan, Sec. 26.3.1):
#' # Var = v^{2n} * n_p_x * n_q_x
#' Var_insurance_x(lt, x = 60, i = 0.06, product = "pure_endowment",
#'                  n = 10, tidy = TRUE)
#'
#' # Endowment = term + pure endowment (Finan, Sec. 26.3.2)
#' Var_insurance_x(lt, x = 60, i = 0.06, product = "endowment",
#'                  n = 10, tidy = TRUE)
#'
#' # Benefit of $100,000
#' Var_insurance_x(lt, x = 60, i = 0.06, product = "whole",
#'                  benefit = 100000)
#'
#' @export
Var_insurance_x <- function(
    lt,
    x,
    i,
    product = c("whole", "term", "endowment", "pure_endowment"),
    benefit = 1,
    n = NULL,
    m = 0L,
    k = 1L,
    tidy = FALSE
) {
  product <- match.arg(product)

  # --- checks ---
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) stop("`lt` must contain columns `x` and `lx`.")
  if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) stop("'i' must be > -1.")
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || abs(x - round(x)) > 1e-10) stop("'x' must be an integer.")
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 0 || abs(m - round(m)) > 1e-10) stop("'m' must be a nonnegative integer.")
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 || abs(k - round(k)) > 1e-10) stop("'k' must be a positive integer.")
  if (!is.numeric(benefit) || length(benefit) != 1L || is.na(benefit) || benefit < 0) {
    stop("'benefit' must be a single nonnegative number.")
  }

  if (product %in% c("term", "endowment", "pure_endowment")) {
    if (is.null(n)) stop("'n' must be provided for term/endowment/pure_endowment.")
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || abs(n - round(n)) > 1e-10) {
      stop("'n' must be a single nonnegative integer.")
    }
    n <- as.integer(round(n))
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  k <- as.integer(round(k))

  lt <- lt[order(lt$x), ]
  if (anyDuplicated(lt$x)) stop("Life table ages `x` must be unique.")
  if (any(is.na(lt$lx)) || any(lt$lx < 0)) stop("`lx` must be nonnegative and not NA.")

  ages <- as.integer(lt$x)
  lx   <- as.numeric(lt$lx)
  omega <- max(ages)
  y <- x + m

  # horizon (assume lx(omega+1) = 0)
  max_years <- max(0L, (omega + 1L) - y)
  if (product == "whole") {
    n_use <- max_years
  } else {
    n_use <- min(n, max_years)
  }
  if (n_use < 0L) n_use <- 0L

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

  # --- Pure endowment moments (Finan, Sec. 26.3.1) ---
  # Var = v^{2n} * n_p_y * n_q_y
  pure_endow_moments <- function(n_years) {
    S <- t_p_int(y, n_years)
    if (is.na(S)) stop("Cannot compute survival to m+n (check life table horizon).")
    tpay <- m + n_years
    EZ  <- benefit * v_fun(tpay) * S
    EZ2 <- (benefit^2) * v_fun(2 * tpay) * S
    c(EZ = EZ, EZ2 = EZ2)
  }

  # --- Term/Whole death-benefit moments (vectorized, Finan Sec. 27) ---
  death_moments <- function(n_years) {
    N <- as.integer(n_years * k)
    if (N <= 0L) return(c(EZ = 0, EZ2 = 0))

    # Grid: u_0, u_1, ..., u_N where u_r = r/k
    u0 <- (0:N) / k
    Su <- vapply(u0, function(uu) t_p_udd(y, uu), numeric(1))
    if (anyNA(Su)) stop("Cannot compute UDD survival on required grid (check life table horizon).")

    # Death in subinterval (u_{r-1}, u_r]: dq_r = S(u_{r-1}) - S(u_r)
    dq <- Su[1:N] - Su[2:(N + 1L)]

    # Payment at end of subperiod: time = m + u_r for r = 1..N
    ur    <- (1:N) / k
    disc1 <- v_fun(m + ur)
    disc2 <- v_fun(2 * (m + ur))

    EZ  <- benefit * sum(disc1 * dq)
    EZ2 <- (benefit^2) * sum(disc2 * dq)

    c(EZ = EZ, EZ2 = EZ2)
  }

  # --- combine by product ---
  if (product == "pure_endowment") {
    moms <- pure_endow_moments(n_use)
  } else if (product == "term" || product == "whole") {
    moms <- death_moments(n_use)
  } else {
    # endowment = term + pure endowment
    # Cov = 0 by mutual exclusivity (Finan, Example 26.15)
    # so E[Z] = E[Z_term] + E[Z_pe], E[Z^2] = E[Z^2_term] + E[Z^2_pe]
    moms_term <- death_moments(n_use)
    moms_pe   <- pure_endow_moments(n_use)
    moms <- moms_term + moms_pe
  }

  var_val <- moms["EZ2"] - moms["EZ"]^2
  if (var_val < 0 && var_val > -1e-12) var_val <- 0

  if (!isTRUE(tidy)) return(as.numeric(var_val))

  tibble::tibble(
    x        = x,
    m        = m,
    n        = if (product == "whole") NA_integer_ else n_use,
    k        = k,
    product  = product,
    i        = i,
    benefit  = benefit,
    EZ       = as.numeric(moms["EZ"]),
    EZ2      = as.numeric(moms["EZ2"]),
    variance = as.numeric(var_val)
  )
}
