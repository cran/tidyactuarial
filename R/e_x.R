#' Expected future lifetime from an annual life table
#'
#' Computes the curtate or complete expected future lifetime at integer age
#' \eqn{x}, optionally restricted to a temporary horizon of \eqn{t} years.
#'
#' @param lt A life table object as produced by \code{\link{lifetable}}
#'   (must contain columns \code{x} and \code{lx}).
#' @param x Integer age(s).
#' @param t Optional nonnegative numeric duration(s). If \code{NULL} (default),
#'   the whole-life expectancy is computed (i.e., horizon extends to
#'   \eqn{\omega - x}). If a numeric value is provided, the \eqn{t}-year
#'   temporary life expectancy is returned.
#' @param type Character: \code{"curtate"} (default) or \code{"complete"}.
#' @param frac Fractional-age assumption for \code{type = "complete"}:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"} (alias of CF), or
#'   \code{"Balducci"}. If not specified and \code{lt} carries a \code{frac}
#'   attribute (set by \code{\link{lifetable}}), that value is used.
#' @param tidy Logical. If \code{TRUE}, returns a tibble.
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#' @param tol Numeric tolerance for integer checks.
#'
#' @details
#' **Curtate life expectancy** (Finan, Section 23.7):
#' \deqn{e_x = \sum_{k=1}^{\omega - x} {}_k p_x = \frac{1}{\ell_x}
#'   \sum_{k=1}^{\omega - x} \ell_{x+k}.}
#'
#' The \eqn{t}-year temporary curtate expectancy is (Finan, Sec. 23.7):
#' \deqn{e_{x:\overline{t}|} = \sum_{k=1}^{t} {}_k p_x.}
#'
#' **Complete life expectancy** (Finan, Section 23.3):
#' \deqn{\breve{e}_x = \int_0^{\omega - x} {}_t p_x \, dt =
#'   \frac{T_x}{\ell_x}.}
#'
#' The integral is decomposed year-by-year. Within each year, the
#' within-year survival integral \eqn{\int_0^s {}_u p_y \, du} is evaluated
#' in closed form under the selected fractional-age assumption
#' (Finan, Section 24):
#' \itemize{
#'   \item UDD (Sec. 24.1): \eqn{\int_0^s {}_u p_y \, du = s - \frac{1}{2} s^2 q_y}
#'   \item CF (Sec. 24.2): \eqn{\int_0^s {}_u p_y \, du = (1 - p_y^s) / (-\ln p_y)}
#'   \item Balducci (Sec. 24.3): \eqn{\int_0^s {}_u p_y \, du = \frac{p_y}{q_y} \ln \left( \frac{p_y + q_y s}{p_y} \right)}
#' }
#'
#' Under UDD, the complete expectancy satisfies the well-known approximation
#' (Finan, Example 20.24):
#' \deqn{\breve{e}_x \approx e_x + \frac{1}{2}.}
#'
#' @return A numeric vector of expected future lifetimes, or a tibble if
#'   \code{tidy = TRUE} with columns \code{x}, \code{t}, \code{type},
#'   \code{frac}, \code{ex}.
#' @export
e_x <- function(
    lt,
    x,
    t = NULL,
    type = c("curtate", "complete"),
    frac,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {
  type <- match.arg(type)

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

  if (missing(lt)) stop("`lt` is required.")
  if (!("x" %in% names(lt)) || !("lx" %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.")
  }
  if (missing(x)) stop("`x` is required.")

  x <- as.numeric(x)
  if (check) {
    if (any(!is.finite(x))) stop("`x` must be finite.")
    if (any(abs(x - round(x)) > tol)) stop("`x` must be integer ages.")
  }
  x <- as.integer(round(x))

  # --- Duration vector handling ---
  # NULL = whole-life (sentinel NA internally); numeric = temporary
  if (is.null(t)) {
    t <- rep(NA_real_, length(x))
  } else {
    t <- as.numeric(t)
    if (check) {
      if (any(!is.finite(t))) stop("`t` must be finite when provided.")
      if (any(t < 0)) stop("`t` must be nonnegative.")
    }
  }

  # --- recycle x, t ---
  nx <- length(x); nt <- length(t)
  n <- max(nx, nt)
  if (!((nx == 1L || nx == n) && (nt == 1L || nt == n))) {
    stop("`x` and `t` must have the same length or one of them must have length 1.")
  }
  if (nx == 1L && n > 1L) x <- rep(x, n)
  if (nt == 1L && n > 1L) t <- rep(t, n)

  ages <- as.integer(lt$x)
  lx   <- as.numeric(lt$lx)
  omega <- max(ages)

  qx_col <- if ("qx" %in% names(lt)) as.numeric(lt$qx) else rep(NA_real_, length(ages))
  px_col <- if ("px" %in% names(lt)) as.numeric(lt$px) else rep(NA_real_, length(ages))

  # --- Within-year survival integral: int_0^s {}_u p_y du ---
  # Closed-form under each fractional-age assumption (Finan, Section 24)
  int_ps <- function(py, qy, s, frac_arg) {
    if (is.na(s) || s <= 0) return(0)
    if (is.na(py) || is.na(qy)) return(NA_real_)
    if (s >= 1) s <- 1

    if (frac_arg == "UDD") {
      # Sec. 24.1: int_0^s (1 - u*q) du = s - 0.5*s^2*q
      return(s - 0.5 * s * s * qy)
    }

    if (frac_arg == "CF") {
      # Sec. 24.2: int_0^s p^u du = (1 - p^s) / (-ln p)
      if (py >= 1) return(s)
      if (py <= 0) return(0)
      a <- -log(py)
      return((1 - py^s) / a)
    }

    # Balducci (Sec. 24.3): int_0^s p/(1-(1-u)q) du = (p/q)*ln((p+qs)/p)
    if (qy <= 0) return(s)
    (py / qy) * log((py + qy * s) / py)
  }

  ex <- rep(NA_real_, n)

  for (m in seq_len(n)) {
    xm <- x[m]
    tm <- t[m]

    i0 <- match(xm, ages)
    if (is.na(i0) || lx[i0] <= 0) {
      ex[m] <- NA_real_
      next
    }

    max_k <- omega - xm
    if (max_k <= 0) {
      ex[m] <- 0
      next
    }

    # Determine horizon: NA = whole-life, numeric = temporary
    if (is.na(tm)) {
      tm_use <- max_k
    } else {
      tm_use <- min(tm, max_k)
    }

    n_int <- floor(tm_use)
    s <- tm_use - n_int

    # Build k_p_x = l_{x+k}/l_x for k = 0..n_int
    idx_k <- match(xm + (0:n_int), ages)
    if (any(is.na(idx_k))) {
      ex[m] <- NA_real_
      next
    }
    p_k <- lx[idx_k] / lx[i0]

    # --- Curtate: e_{x:n} = sum_{k=1}^{n} k_p_x (Finan, Sec. 23.7) ---
    if (type == "curtate") {
      ex[m] <- if (n_int >= 1) sum(p_k[-1]) else 0
      next
    }

    # --- Complete: \breve{e}_{x:n} = int_0^n t_p_x dt (Finan, Sec. 23.3) ---
    if (n_int == 0) {
      if (s == 0) {
        ex[m] <- 0
        next
      }
      y <- xm
      iy <- match(y, ages)
      qy <- qx_col[iy]; py <- px_col[iy]
      if (is.na(qy) || is.na(py)) {
        iy1 <- match(y + 1L, ages)
        if (!is.na(iy1) && lx[iy] > 0) {
          qy <- (lx[iy] - lx[iy1]) / lx[iy]
          py <- 1 - qy
        }
      }
      ex[m] <- int_ps(py, qy, s, frac)
      next
    }

    # Full years: sum_{k=0}^{n-1} k_p_x * int_0^1 u_p_{x+k} du
    ys <- xm + (0:(n_int - 1))
    iy <- match(ys, ages)

    qy <- qx_col[iy]
    py <- px_col[iy]

    need <- is.na(qy) | is.na(py)
    if (any(need)) {
      iy1 <- match(ys + 1L, ages)
      okd <- need & !is.na(iy1) & (lx[iy] > 0)
      qy[okd] <- (lx[iy[okd]] - lx[iy1[okd]]) / lx[iy[okd]]
      py[okd] <- 1 - qy[okd]
    }

    I1 <- vapply(seq_along(py), function(j) int_ps(py[j], qy[j], 1, frac), numeric(1))
    part1 <- sum(p_k[seq_along(I1)] * I1)

    # Partial last year at age x+n
    if (s == 0) {
      part2 <- 0
    } else {
      yN <- xm + n_int
      iN <- match(yN, ages)
      qN <- qx_col[iN]; pN <- px_col[iN]
      if (is.na(qN) || is.na(pN)) {
        iN1 <- match(yN + 1L, ages)
        if (!is.na(iN1) && lx[iN] > 0) {
          qN <- (lx[iN] - lx[iN1]) / lx[iN]
          pN <- 1 - qN
        }
      }
      part2 <- p_k[length(p_k)] * int_ps(pN, qN, s, frac)
    }

    ex[m] <- part1 + part2
  }

  if (isTRUE(tidy)) {
    t_out <- t
    t_out[is.na(t_out)] <- NA_real_
    return(tibble::tibble(x = x, t = t_out, type = type, frac = frac, ex = ex))
  }
  ex
}
