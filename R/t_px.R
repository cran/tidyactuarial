#' t-year survival probability from a life table
#'
#' Computes the t-year survival probability
#' \deqn{{}_t p_x = P[T(x) > t]}
#' using an annual life table, allowing for fractional ages under standard
#' actuarial assumptions.
#'
#' @param lt A lifetable object as produced by \code{\link{lifetable}}.
#'   Must contain columns \code{x} and \code{lx}. Columns \code{qx} or
#'   \code{px} are used if present.
#' @param x Integer age(s) at which survival starts.
#' @param t Nonnegative numeric duration(s) in years (can be fractional).
#' @param frac Fractional-age assumption:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"} (alias of CF), or
#'   \code{"Balducci"}. If not specified and \code{lt} carries a \code{frac}
#'   attribute (set by \code{\link{lifetable}}), that value is used.
#' @param tidy Logical. If \code{TRUE}, returns a tibble with columns
#'   \code{x}, \code{t}, \code{frac}, \code{tpx}.
#' @param check Logical. If \code{TRUE}, performs validity checks.
#' @param tol Numeric tolerance for integer checks on \code{x}.
#'
#' @details
#' The integer-year survival is obtained directly from the life table
#' (Finan, Section 22):
#' \deqn{{}_n p_x = \frac{\ell_{x+n}}{\ell_x}}
#'
#' For non-integer durations, let \eqn{t = n + s} with \eqn{n = \lfloor t \rfloor}
#' and \eqn{s \in [0,1)}. Then (Finan, Section 24):
#' \deqn{{}_t p_x = {}_n p_x \times {}_s p_{x+n}}
#'
#' The fractional-year factor \eqn{{}_s p_y} depends on the assumption:
#' \itemize{
#'   \item UDD (Finan, Sec. 24.1): \eqn{{}_s p_y = 1 - s \times q_y}
#'   \item CF (Finan, Sec. 24.2): \eqn{{}_s p_y = (p_y)^s}
#'   \item Balducci (Finan, Sec. 24.3):
#'     \eqn{{}_s p_y = \frac{p_y}{1 - (1 - s) \times q_y}}
#' }
#'
#' If \eqn{x + t > \omega} (the terminal age), the function returns 0
#' since no survival is possible beyond the table's limiting age.
#'
#' @return Numeric vector of \eqn{{}_t p_x}, or a tibble if
#'   \code{tidy = TRUE}.
#' @export
t_px <- function(
    lt,
    x,
    t,
    frac,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
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

  if (missing(lt)) stop("`lt` is required.")
  if (!("x" %in% names(lt)) || !("lx" %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.")
  }
  if (missing(x) || missing(t)) stop("`x` and `t` are required.")

  x <- as.numeric(x)
  t <- as.numeric(t)

  if (check) {
    if (any(!is.finite(x)) || any(!is.finite(t))) stop("`x` and `t` must be finite.")
    if (any(abs(x - round(x)) > tol)) stop("`x` must be integer ages.")
    if (any(t < 0)) stop("`t` must be nonnegative.")
  }

  x <- as.integer(round(x))

  # --- recycle x, t ---
  nx <- length(x); nt <- length(t)
  n <- max(nx, nt)
  if (!((nx == 1L || nx == n) && (nt == 1L || nt == n))) {
    stop("`x` and `t` must have the same length or one of them must have length 1.")
  }
  if (nx == 1L && n > 1L) x <- rep(x, n)
  if (nt == 1L && n > 1L) t <- rep(t, n)

  # --- extract table vectors ---
  ages <- as.integer(lt$x)
  lx   <- as.numeric(lt$lx)
  omega <- max(ages)

  qx_col <- if ("qx" %in% names(lt)) as.numeric(lt$qx) else NULL
  px_col <- if ("px" %in% names(lt)) as.numeric(lt$px) else NULL

  out <- rep(NA_real_, n)

  # --- fast path: t = 0 always returns 1 ---
  zero_t <- (t == 0)
  out[zero_t] <- 1

  # --- fast path: x + t beyond omega returns 0 ---
  beyond <- (!zero_t) & (x + t > omega + tol)
  out[beyond] <- 0

  # --- remaining cases ---
  todo <- which(is.na(out))
  if (length(todo) == 0L) {
    if (isTRUE(tidy)) {
      return(tibble::tibble(x = x, t = t, frac = frac, tpx = out))
    }
    return(out)
  }

  # Split t = n_int + s
  n_int <- floor(t[todo])
  s     <- t[todo] - n_int

  x_todo <- x[todo]

  i0 <- match(x_todo, ages)
  i1 <- match(x_todo + as.integer(n_int), ages)

  ok <- !is.na(i0) & !is.na(i1) & (lx[i0] > 0)

  # Integer survival part: n_p_x = l_{x+n} / l_x
  p_int <- rep(NA_real_, length(todo))
  p_int[ok] <- lx[i1[ok]] / lx[i0[ok]]

  # --- fractional part: s_p_{x+n} ---
  # Get q_y and p_y at y = x + n_int
  y  <- x_todo + as.integer(n_int)
  iy <- match(y, ages)

  qy <- py <- rep(NA_real_, length(todo))

  if (!is.null(qx_col)) qy[ok] <- qx_col[iy[ok]]
  if (!is.null(px_col)) py[ok] <- px_col[iy[ok]]

  # Derive qy from lx if not available
  need <- ok & is.na(qy)
  if (any(need)) {
    iy1 <- match(y + 1L, ages)
    good <- need & !is.na(iy1)
    qy[good] <- (lx[iy[good]] - lx[iy1[good]]) / lx[iy[good]]
  }
  py[is.na(py)] <- 1 - qy[is.na(py)]

  # Compute fractional survival (Finan, Section 24)
  ps <- rep(NA_real_, length(todo))
  ok2 <- ok & !is.na(py) & !is.na(qy)

  if (any(ok2)) {
    s_ok  <- s[ok2]
    py_ok <- py[ok2]
    qy_ok <- qy[ok2]

    if (frac == "UDD") {
      # Sec. 24.1: s_p_y = 1 - s * q_y
      ps[ok2] <- 1 - s_ok * qy_ok
    } else if (frac == "CF") {
      # Sec. 24.2: s_p_y = p_y^s
      ps[ok2] <- py_ok^s_ok
    } else if (frac == "Balducci") {
      # Sec. 24.3: s_p_y = p_y / (1 - (1-s) * q_y)
      denom <- 1 - (1 - s_ok) * qy_ok
      ps[ok2] <- ifelse(denom > 0, py_ok / denom, 0)
    }
  }

  # Combine integer and fractional parts
  result <- p_int * ps
  # Clamp to [0, 1]
  result <- pmin(pmax(result, 0), 1)

  out[todo] <- result

  if (isTRUE(tidy)) {
    return(tibble::tibble(x = x, t = t, frac = frac, tpx = out))
  }
  out
}
