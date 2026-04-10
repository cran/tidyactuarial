#' Expected future lifetime for two independent lives
#'
#' Computes the expected future lifetime for two independent lives aged
#' \code{x} and \code{y}, for either joint-life (first death) or
#' last-survivor (second death).
#'
#' @param lt A life table data frame with columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for life 1.
#' @param y Integer actuarial age for life 2.
#' @param t Optional nonnegative numeric duration(s). If \code{NULL}, uses
#'   the maximum horizon allowed by the table.
#' @param type Character: \code{"curtate"} or \code{"complete"}.
#' @param frac Fractional-age assumption for \code{type = "complete"},
#'   passed to \code{\link{t_pxy}}: \code{"UDD"}, \code{"CF"},
#'   \code{"CML"}, or \code{"Balducci"}. If not specified and \code{lt}
#'   carries a \code{frac} attribute, that value is used.
#' @param cohort Two-life cohort: \code{"first"} (joint-life) or
#'   \code{"last"} (last survivor).
#' @param tidy Logical. If \code{TRUE}, returns a tibble.
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#' @param tol Numeric tolerance for integer checks.
#'
#' @details
#' **Curtate expectation** (Finan, Section 56.4 / Section 57):
#' \deqn{e_{xy} = \sum_{k=1}^{\infty} {}_kp_{xy}, \quad
#'   e_{\overline{xy}} = \sum_{k=1}^{\infty} {}_kp_{\overline{xy}}.}
#'
#' **Complete expectation** (Finan, Section 56.4):
#' \deqn{\mathring{e}_{xy} = \int_0^{\infty} {}_tp_{xy} \, dt.}
#'
#' The integral is decomposed year-by-year. Within each year, the survival
#' integral for the two-life status is computed numerically via composite
#' trapezoid (80-point grid), since closed-form expressions for joint/last
#' survivor under fractional-age assumptions are complex.
#'
#' **Key identity** (Finan, Example 57.4):
#' \deqn{\mathring{e}_{\overline{xy}} = \mathring{e}_x + \mathring{e}_y
#'   - \mathring{e}_{xy}.}
#'
#' This can be used to cross-validate results.
#'
#' @return Numeric vector, or tibble if \code{tidy = TRUE}.
#'
#' @seealso \code{\link{e_x}} for single-life expectancy,
#'   \code{\link{t_pxy}} for two-life survival probabilities,
#'   \code{\link{annuity_xy}} for two-life annuity APVs.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Curtate joint-life expectancy (Finan, Sec. 56.4)
#' e_xy(lt, x = 60, y = 62, type = "curtate", cohort = "first")
#'
#' # Curtate last-survivor expectancy (Finan, Sec. 57)
#' e_xy(lt, x = 60, y = 62, type = "curtate", cohort = "last")
#'
#' # Verify identity (Finan, Example 57.4):
#' # e_{xy-bar} = e_x + e_y - e_xy
#' e_joint <- e_xy(lt, x = 60, y = 62, type = "curtate", cohort = "first")
#' e_last  <- e_xy(lt, x = 60, y = 62, type = "curtate", cohort = "last")
#' e_x_val <- e_x(lt, x = 60, type = "curtate")
#' e_y_val <- e_x(lt, x = 62, type = "curtate")
#' c(last_surv = e_last, sum_minus_joint = e_x_val + e_y_val - e_joint)
#'
#' # Complete joint-life expectancy under UDD
#' e_xy(lt, x = 60, y = 62, type = "complete", frac = "UDD", cohort = "first")
#'
#' # Temporary: 3-year curtate joint-life
#' e_xy(lt, x = 60, y = 62, t = 3, type = "curtate", cohort = "first")
#'
#' # Tidy output
#' e_xy(lt, x = 60, y = 62, type = "curtate", cohort = "first", tidy = TRUE)
#'
#' @export
e_xy <- function(
    lt,
    x, y,
    t = NULL,
    type = c("curtate", "complete"),
    frac,
    cohort = c("first", "last"),
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {
  type   <- match.arg(type)
  cohort <- match.arg(cohort)

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
  if (!is.data.frame(lt)) stop("`lt` must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) stop("`lt` must contain columns `x` and `lx`.")
  if (missing(x) || missing(y)) stop("`x` and `y` are required.")

  x <- as.numeric(x); y <- as.numeric(y)

  if (check) {
    if (any(!is.finite(x)) || any(!is.finite(y))) stop("`x` and `y` must be finite.")
    if (any(abs(x - round(x)) > tol) || any(abs(y - round(y)) > tol)) {
      stop("`x` and `y` must be integer ages.")
    }
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))

  # --- duration handling ---
  if (is.null(t)) {
    t <- rep(NA_real_, max(length(x), length(y)))
  } else {
    t <- as.numeric(t)
    if (check) {
      if (any(!is.finite(t))) stop("`t` must be finite when provided.")
      if (any(t < 0)) stop("`t` must be nonnegative.")
    }
  }

  # --- recycle x, y, t ---
  nx <- length(x); ny <- length(y); nt <- length(t)
  nout <- max(nx, ny, nt)
  ok <- (nx == 1L || nx == nout) && (ny == 1L || ny == nout) && (nt == 1L || nt == nout)
  if (!ok) stop("`x`, `y`, and `t` must have the same length or be length 1.")

  if (nx == 1L && nout > 1L) x <- rep(x, nout)
  if (ny == 1L && nout > 1L) y <- rep(y, nout)
  if (nt == 1L && nout > 1L) t <- rep(t, nout)

  ages  <- as.integer(lt$x)
  omega <- max(ages, na.rm = TRUE)

  status <- if (cohort == "first") "joint" else "last"

  # --- within-year survival integral via composite trapezoid (80 points) ---
  # Closed-form for joint/last under each frac is complex;
  # numerical integration is robust and accurate for annual life tables.
  int_surv <- function(x0, y0, s) {
    if (is.na(s) || s <= 0) return(0)
    if (s >= 1) s <- 1
    grid_n <- 80L
    u  <- seq(0, s, length.out = grid_n)
    Su <- vapply(u, function(uu) {
      t_pxy(lt, x0, y0, t = uu, frac = frac, status = status)
    }, numeric(1))
    if (anyNA(Su)) return(NA_real_)
    sum((head(Su, -1) + tail(Su, -1)) / 2) * (s / (grid_n - 1))
  }

  ex <- rep(NA_real_, nout)

  for (j in seq_len(nout)) {
    xj <- x[j]; yj <- y[j]; tj <- t[j]

    # --- maximum horizon depends on status ---
    # Joint life: limited by the older life (both must be alive)
    # Last survivor: limited by the younger life (at least one alive)
    if (cohort == "first") {
      max_k <- omega - max(xj, yj)
    } else {
      max_k <- omega - min(xj, yj)
    }
    if (max_k <= 0) {
      ex[j] <- 0
      next
    }

    if (is.na(tj)) {
      tj_use <- max_k
    } else {
      tj_use <- min(tj, max_k)
    }

    n_int <- floor(tj_use)
    s <- tj_use - n_int

    # --- curtate: e = sum_{k=1}^{n} k_p (Finan, Sec. 56.4) ---
    if (type == "curtate") {
      if (n_int <= 0) {
        ex[j] <- 0
      } else {
        pk <- vapply(1:n_int, function(k) {
          t_pxy(lt, xj, yj, t = k, frac = frac, status = status)
        }, numeric(1))
        ex[j] <- if (anyNA(pk)) NA_real_ else sum(pk)
      }
      next
    }

    # --- complete: integral decomposed year-by-year ---
    if (n_int == 0) {
      ex[j] <- int_surv(xj, yj, s)
      next
    }

    # Full years: sum_{k=0}^{n-1} k_p * int_0^1 u_p_{x+k, y+k} du
    p0k <- vapply(0:(n_int - 1), function(k) {
      t_pxy(lt, xj, yj, t = k, frac = frac, status = status)
    }, numeric(1))
    I1 <- vapply(0:(n_int - 1), function(k) {
      int_surv(xj + k, yj + k, 1)
    }, numeric(1))
    if (anyNA(p0k) || anyNA(I1)) {
      ex[j] <- NA_real_
      next
    }
    part1 <- sum(p0k * I1)

    # Last partial year
    part2 <- if (s == 0) 0 else {
      pn <- t_pxy(lt, xj, yj, t = n_int, frac = frac, status = status)
      In <- int_surv(xj + n_int, yj + n_int, s)
      if (is.na(pn) || is.na(In)) NA_real_ else pn * In
    }

    ex[j] <- part1 + part2
  }

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x, y = y, t = t, type = type, frac = frac,
      cohort = cohort, ex = ex
    ))
  }
  ex
}
