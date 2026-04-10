#' Actuarial present value of a life annuity
#'
#' Computes the APV of a discrete life annuity at actuarial age \eqn{x} using
#' a life table. Supports term \eqn{n} (temporary), integer deferral \eqn{m},
#' \eqn{k}-thly payments (exact under UDD), and Woolhouse approximations up
#' to second order.
#'
#' @param lt A lifetable object as produced by \code{\link{lifetable}}.
#'   Must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age.
#' @param i Effective annual interest rate (must satisfy \code{i > -1}).
#' @param n Integer term in years. If \code{NULL} (default), whole life to
#'   end of table.
#' @param m Integer deferral in years (default \code{0}).
#' @param k Integer payments per year (default \code{1}). Example:
#'   \code{k = 12} for monthly.
#' @param timing \code{"immediate"} (payments at end of period) or
#'   \code{"due"} (beginning). Default \code{"immediate"}.
#' @param woolhouse For \code{k > 1}: \code{"none"} (exact UDD),
#'   \code{"first"} (2-term Woolhouse), or \code{"second"} (3-term Woolhouse).
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#'
#' @details
#' **Annual annuity-due** (Finan, Section 37.2, Example 37.9):
#' \deqn{\ddot{a}_{x:\overline{n}|} = \sum_{j=0}^{n-1} v^j \times {}_j p_x}
#'
#' **Annual annuity-immediate** (Finan, Section 37.5):
#' \deqn{a_{x:\overline{n}|} = \sum_{j=1}^{n} v^j \times {}_j p_x}
#'
#' **Deferral** (Finan, Section 37.3):
#' \deqn{{}_{m|}\ddot{a}_x = v^m \times {}_m p_x \times \ddot{a}_{x+m}}
#'
#' **k-thly exact under UDD** (Finan, Section 38): each \eqn{1/k}-year
#' fractional survival is computed under UDD, giving an exact result:
#' \deqn{\ddot{a}^{(k)}_{x:\overline{n}|} = \frac{1}{k}\sum_{j=0}^{kn-1}
#'   v^{j/k} \times {}_{j/k} p_x}
#'
#' **Woolhouse approximations** (Finan, Problems 38.9-38.10):
#' \itemize{
#'   \item 2-term (first order):
#'     \eqn{\ddot{a}^{(k)}_{x:\overline{n}|} \approx \ddot{a}_{x:\overline{n}|} - \frac{k-1}{2k}(1 - {}_n E_x)}
#'   \item 3-term (second order):
#'     \eqn{\ddot{a}^{(k)}_{x:\overline{n}|} \approx \ddot{a}_{x:\overline{n}|} - \frac{k-1}{2k}(1 - {}_n E_x) - \frac{k^2-1}{12k^2}[\delta + \mu_x - {}_n E_x(\delta + \mu_{x+n})]}
#' }
#'
#' Conversion from due to immediate for k-thly temporaries
#' (Finan, Problem 38.6):
#' \deqn{a^{(k)}_{x:\overline{n}|} = \ddot{a}^{(k)}_{x:\overline{n}|} - \frac{1}{k}(1 - {}_n E_x)}
#'
#' @return A single numeric APV value, or a one-row tibble if
#'   \code{tidy = TRUE}.
#' @export
annuity_x <- function(
    lt, x, i,
    n = NULL,
    m = 0L,
    k = 1L,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    tidy = FALSE
) {
  timing    <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)

  # --- checks ---
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) stop("Life table must contain columns 'x' and 'lx'.")
  if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) stop("'i' must be > -1.")
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || abs(x - round(x)) > 1e-10) stop("'x' must be an integer.")
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 0 || abs(m - round(m)) > 1e-10) stop("'m' must be a nonnegative integer.")
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 || abs(k - round(k)) > 1e-10) stop("'k' must be a positive integer.")

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  k <- as.integer(round(k))

  # sort and basic integrity
  lt <- lt[order(lt$x), ]
  if (anyDuplicated(lt$x)) stop("Life table ages 'x' must be unique.")
  if (any(is.na(lt$lx)) || any(lt$lx < 0)) stop("'lx' must be nonnegative and not NA.")

  ages <- as.integer(lt$x)
  lx   <- as.numeric(lt$lx)
  omega <- max(ages)

  # helper: get lx at integer age; assume lx(omega+1) = 0
  get_lx <- function(age) {
    idx <- match(age, ages)
    if (!is.na(idx)) return(lx[idx])
    if (age == omega + 1L) return(0)
    NA_real_
  }

  # integer-year survival
  t_p_int <- function(age, tt) {
    if (tt == 0) return(1)
    l0 <- get_lx(age)
    l1 <- get_lx(age + tt)
    if (is.na(l0) || is.na(l1) || l0 <= 0) return(NA_real_)
    l1 / l0
  }

  # UDD fractional survival
  t_p_udd <- function(age, u) {
    if (u < 0) return(NA_real_)
    if (u == 0) return(1)
    tt <- floor(u)
    s  <- u - tt
    pt <- t_p_int(age, tt)
    if (is.na(pt)) return(NA_real_)
    if (s == 0) return(pt)
    y   <- age + tt
    ly  <- get_lx(y)
    ly1 <- get_lx(y + 1L)
    if (is.na(ly) || is.na(ly1) || ly <= 0) return(NA_real_)
    dy <- ly - ly1
    ps <- (ly - s * dy) / ly
    pt * ps
  }

  # discount
  v_pow <- function(tt) (1 + i)^(-tt)

  # --- deferral factor: v^m * m_p_x ---
  defer <- v_pow(m) * t_p_int(x, m)
  if (is.na(defer)) stop("Deferral age x+m is outside the life table (or lx(x)=0).")

  y <- x + m  # annuity starts at age y

  # --- determine n (whole life default) ---
  max_years <- max(0L, (omega + 1L) - y)
  if (is.null(n)) {
    n <- max_years
  } else {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || abs(n - round(n)) > 1e-10) {
      stop("'n' must be a nonnegative integer or NULL.")
    }
    n <- as.integer(round(n))
    if (n > max_years) stop("'n' exceeds the horizon allowed by the life table (need lx up to x+m+n).")
  }

  if (n == 0) {
    result <- 0
    if (isTRUE(tidy)) {
      return(tibble::tibble(x = x, i = i, n = n, m = m, k = k,
                            timing = timing, woolhouse = woolhouse, apv = 0))
    }
    return(0)
  }

  # --- nEx = v^n * n_p_y (pure endowment at y for term n) ---
  nEx <- v_pow(n) * t_p_int(y, n)
  if (is.na(nEx)) nEx <- 0

  # --- annual exact ---
  annual_exact <- function(age, nn, tim) {
    if (tim == "due") {
      tt <- 0:(nn - 1L)
    } else {
      tt <- 1:nn
    }
    sp <- vapply(tt, function(t) t_p_int(age, t), numeric(1))
    if (anyNA(sp)) stop("Life table does not support required ages for annual payments.")
    sum(v_pow(tt) * sp)
  }

  # --- k-thly exact under UDD ---
  kthly_exact_udd <- function(age, nn, kk, tim) {
    if (tim == "due") {
      j <- 0:(kk * nn - 1L)
    } else {
      j <- 1:(kk * nn)
    }
    u  <- j / kk
    sp <- vapply(u, function(uu) t_p_udd(age, uu), numeric(1))
    if (anyNA(sp)) stop("Life table does not support required ages for UDD k-thly payments.")
    sum((1 / kk) * v_pow(u) * sp)
  }

  # --- main computation ---
  if (k == 1L) {
    result <- defer * annual_exact(y, n, timing)
  } else if (woolhouse == "none") {
    # exact under UDD
    result <- defer * kthly_exact_udd(y, n, k, timing)
  } else {
    # --- Woolhouse approximation (Finan, Problems 38.9-38.10) ---
    # Always start from annual annuity-due
    adue <- annual_exact(y, n, "due")

    # 2-term Woolhouse (Finan, Problem 38.9a):
    # \ddot{a}^(m)_{y:n} \approx \ddot{a}_{y:n} - (m-1)/(2m) * (1 - nEx)
    adj1 <- (k - 1) / (2 * k) * (1 - nEx)

    if (woolhouse == "first") {
      adue_k <- adue - adj1
    } else {
      # 3-term Woolhouse (Finan, Problem 38.10a):
      # \ddot{a}^(m)_{y:n} \approx \ddot{a}_{y:n} - adj1 - (m^2-1)/(12m^2) * [\delta+\mu(y) - nEx*(\delta+\mu(y+n))]
      delta <- log(1 + i)

      # \mu(y) approximated from p_y
      ly  <- get_lx(y)
      ly1 <- get_lx(y + 1L)
      if (is.na(ly) || is.na(ly1) || ly <= 0) stop("Cannot compute mu at age y from lx.")
      p_y <- ly1 / ly
      mu_y <- if (!is.na(p_y) && p_y > 0) -log(p_y) else 0

      # \mu(y+n) approximated from p_{y+n} (for temporary correction)
      lyn  <- get_lx(y + n)
      lyn1 <- get_lx(y + n + 1L)
      if (!is.na(lyn) && !is.na(lyn1) && lyn > 0) {
        p_yn <- lyn1 / lyn
        mu_yn <- if (p_yn > 0) -log(p_yn) else 0
      } else {
        mu_yn <- 0
      }

      adj2 <- (k^2 - 1) / (12 * k^2) * (delta + mu_y - nEx * (delta + mu_yn))
      adue_k <- adue - adj1 - adj2
    }

    # Convert due -> immediate for k-thly (Finan, Problem 38.6):
    # a^(m)_{y:n} = \ddot{a}^(m)_{y:n} - (1/m)(1 - nEx)
    a_k <- if (timing == "due") adue_k else (adue_k - (1 / k) * (1 - nEx))

    result <- defer * a_k
  }

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x, i = i, n = as.integer(n), m = m, k = k,
      timing = timing, woolhouse = woolhouse, apv = result
    ))
  }
  result
}
