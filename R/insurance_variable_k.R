#' Actuarial present value of a life insurance with variable k-thly benefits
#'
#' Computes the actuarial present value of a life insurance where the death
#' benefit may vary by subperiod (\eqn{k} payments per year) and is payable
#' at the end of the subperiod of death.
#'
#' @param lt A life table object or data frame containing at least \code{x}
#'   and \code{lx}.
#' @param x Integer actuarial age at issue.
#' @param i Effective annual interest rate (must be \code{> -1}).
#' @param benefit Numeric vector of benefits by subperiod, or a function of
#'   time returning the benefit at time \eqn{t}.
#' @param n Optional term in years. If \code{NULL}, the term is inferred
#'   from the length of \code{benefit} (when numeric).
#' @param m Nonnegative integer deferral period in years (default \code{0}).
#' @param k Number of subperiods per year (default \code{12}).
#' @param frac Fractional-age assumption used in survival probabilities:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#'   If not specified and \code{lt} carries a \code{frac} attribute, that
#'   value is used.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#'
#' @details
#' Let \eqn{k} be the number of subperiods per year and \eqn{N = nk} the
#' total number of subperiods. With deferral \eqn{m}, the actuarial present
#' value at age \eqn{x} is (generalizing Finan, Section 29):
#' \deqn{
#'   \text{APV} = \sum_{j=1}^{N} v^{m + t_j} \cdot b_j \cdot
#'   \left({{}_{m + t_{j-1}}}p_x - {{}_{m + t_j}}p_x\right)
#' }
#' where \eqn{t_j = j/k} and \eqn{b_j} is the benefit payable if death
#' occurs in the interval \eqn{(m + t_{j-1},\, m + t_j]}.
#'
#' This is the general form that encompasses:
#' \itemize{
#'   \item \strong{Level insurance}: constant \eqn{b_j = 1} reduces to
#'     \eqn{A^{(k)1}_{x:\overline{n}|}} (Finan, Section 31).
#'   \item \strong{Increasing insurance}: \eqn{b_j = \lceil j/k \rceil}
#'     gives \eqn{(IA)^{(k)}_{x:\overline{n}|}} (Finan, Section 29.3).
#'   \item \strong{Decreasing insurance}: \eqn{b_j = n - \lfloor(j-1)/k\rfloor}
#'     gives \eqn{(DA)^{(k)}_{x:\overline{n}|}} (Finan, Section 29.3).
#'   \item \strong{Credit insurance}: \eqn{b_j = B(t_j)} where \eqn{B(t)}
#'     is the outstanding loan balance.
#' }
#'
#' The benefit may be supplied either as a numeric vector indexed by
#' subperiod or as a function of time.
#'
#' Fractional survival probabilities are computed via \code{\link{t_px}}
#' under the selected assumption (Finan, Section 24).
#'
#' @return A numeric actuarial present value, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{insurance_x}} for level-benefit life insurance,
#'   \code{\link{annuity_x}} for life annuity APVs,
#'   \code{\link{t_px}} for survival probabilities,
#'   Monte Carlo simulation tools may be used for empirical variance estimation.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Monthly insurance with increasing benefits (Finan, Sec. 29.3 style)
#' insurance_variable_k(
#'   lt, x = 60, i = 0.05,
#'   benefit = seq(100, 1200, length.out = 12),
#'   n = 1, k = 12
#' )
#'
#' # Credit-style insurance with declining outstanding balance
#' balance <- function(t) 2000 * exp(-0.3 * t)
#' insurance_variable_k(
#'   lt, x = 60, i = 0.05,
#'   benefit = balance,
#'   n = 1, k = 12
#' )
#'
#' # Level benefit = 1 should match a term insurance
#' # (approximately, since insurance_x uses annual and this uses k-thly)
#' insurance_variable_k(lt, x = 60, i = 0.05, benefit = 1, n = 5, k = 1)
#' insurance_x(mortality_table = lt, age = 60, rate = 0.05, term_years = 5, insurance_type = "term")
#'
#' # 2-year deferred, 3-year term with monthly varying benefits
#' insurance_variable_k(
#'   lt, x = 60, i = 0.05,
#'   benefit = rep(1000, 36),
#'   n = 3, m = 2, k = 12
#' )
#'
#' # Tidy output
#' insurance_variable_k(
#'   lt, x = 60, i = 0.05,
#'   benefit = rep(1000, 12), n = 1, k = 12, tidy = TRUE
#' )
#'
#' @export
insurance_variable_k <- function(
    lt,
    x,
    i,
    benefit,
    n = NULL,
    m = 0,
    k = 12,
    frac,
    tidy = FALSE,
    check = TRUE
) {

  # --- inherit frac from lifetable attribute if not supplied ---
  if (missing(frac)) {
    lt_frac <- attr(lt, "frac")
    if (!is.null(lt_frac) && lt_frac %in% c("UDD", "CF", "Balducci")) {
      frac <- lt_frac
    } else {
      frac <- "UDD"
    }
  } else {
    frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac == "CML") frac <- "CF"
  }

  if (isTRUE(check)) {
    if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
    if (!all(c("x", "lx") %in% names(lt))) {
      stop("Life table must contain columns 'x' and 'lx'.")
    }
    if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) {
      stop("'i' must be a single numeric value > -1.")
    }
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || abs(x - round(x)) > 1e-10) {
      stop("'x' must be a single integer age.")
    }
    if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 0 || abs(m - round(m)) > 1e-10) {
      stop("'m' must be a single nonnegative integer.")
    }
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0 || abs(k - round(k)) > 1e-10) {
      stop("'k' must be a single positive integer.")
    }
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  k <- as.integer(round(k))

  # --- determine N (total subperiods) ---
  if (is.null(n)) {
    if (is.function(benefit)) stop("Provide 'n' when benefit is a function.")
    N <- length(benefit)
    n <- N / k
  } else {
    n <- as.integer(round(n))
    N <- n * k
  }

  # --- build benefit vector ---
  if (is.function(benefit)) {
    times <- (1:N) / k
    bvec  <- benefit(times)
  } else {
    if (length(benefit) == 1L) benefit <- rep(benefit, N)
    if (length(benefit) != N) stop("'benefit' length must equal n * k = ", N, ".")
    bvec <- benefit
  }

  # --- compute APV (vectorized) ---
  # Payment times relative to age x: m + j/k
  # Survival from age x to m + t_{j-1} and m + t_j
  t0_vec <- m + (0:(N - 1L)) / k   # start of each subperiod from age x
  t1_vec <- m + (1:N) / k           # end of each subperiod from age x

  # Survival probabilities from age x
  S_t0 <- t_px(lt, x = x, t = t0_vec, frac = frac, check = FALSE)
  S_t1 <- t_px(lt, x = x, t = t1_vec, frac = frac, check = FALSE)

  # Death probability in each subperiod: P(death in (t0, t1]) = S(t0) - S(t1)
  dq <- S_t0 - S_t1

  # Discount to time 0 at age x: v^(m + j/k)
  disc <- (1 + i)^(-t1_vec)

  # APV = sum b_j * v^{m+t_j} * [S(m+t_{j-1}) - S(m+t_j)]
  apv <- sum(bvec * disc * dq)

  if (!isTRUE(tidy)) return(apv)

  tibble::tibble(
    x = x, m = m, n = n, k = k, i = i,
    frac = frac, apv = apv
  )
}
