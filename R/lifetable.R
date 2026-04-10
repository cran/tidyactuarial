#' Build an annual life table (tidy tibble) from lx, qx, px, or mx
#'
#' Creates an annual life table with **integer, consecutive ages** and returns a
#' **tibble** (tidyverse-friendly) with class `"lifetable"`.
#'
#' The table can be built from exactly one of:
#' \itemize{
#'   \item \code{lx} (survivors), or
#'   \item \code{qx} (one-year death probabilities), or
#'   \item \code{px} (one-year survival probabilities), or
#'   \item \code{mx} (central death rates),
#' }
#' and the function will compute the remaining columns consistently: \code{dx},
#' \code{qx}, \code{px}, and \code{mx}.
#'
#' When multiple inputs are provided, priority is:
#' \code{lx} > \code{qx} > \code{px} > \code{mx}.
#' If \code{lx} is provided together with \code{qx}, cross-consistency is
#' validated (both must agree via \eqn{q_x = (\ell_x - \ell_{x+1}) / \ell_x}).
#'
#' By default, the table is actuarially closed at \code{omega}:
#' \deqn{\ell_{\omega+1} = 0 \Rightarrow d_{\omega} = \ell_{\omega}
#'       \Rightarrow q_{\omega} = 1 \Rightarrow p_{\omega} = 0.}
#'
#' @details
#' The life table follows the standard actuarial construction described in
#' Finan, Sections 22--24 (Exam MLC preparation).
#'
#' The basic identities are (Finan, Section 22):
#' \deqn{\ell_x = \ell_0 \cdot s(x), \quad d_x = \ell_x - \ell_{x+1},
#'       \quad q_x = d_x / \ell_x, \quad p_x = \ell_{x+1} / \ell_x.}
#'
#' The central death rate \code{mx} is computed via the discrete approximation
#' (Finan, Section 23.9):
#' \deqn{m_x = \frac{q_x}{1 - a_x \cdot q_x}}
#' which under UDD (\code{ax = 0.5}) reduces to the classical formula
#' \eqn{m_x = q_x / (1 - 0.5 \, q_x)} (Finan, Section 24.1). This arises
#' because under UDD, \eqn{L_x = \ell_x - \tfrac{1}{2} d_x}, and therefore
#' \eqn{m_x = d_x / L_x}.
#'
#' At the terminal age \eqn{\omega} with \code{close = TRUE}, closure forces
#' \eqn{q_\omega = 1}, \eqn{p_\omega = 0}, and \eqn{d_\omega = \ell_\omega}.
#' The corresponding \eqn{m_\omega} equals \eqn{1/(1 - a_x)}, which is 2
#' under UDD (\code{ax = 0.5}). If \code{ax = 1}, \eqn{m_\omega = \infty}.
#'
#' @param x Numeric vector of ages. Must be **integer** and **consecutive**
#'   (annual table), e.g. \code{0:110}.
#' @param lx Optional numeric vector of survivors \eqn{\ell_x}. Must be
#'   nonnegative and nonincreasing.
#' @param qx Optional numeric vector of one-year death probabilities \eqn{q_x}.
#'   \code{NA} values are allowed (useful at the last age when \code{close=TRUE}),
#'   but \code{Inf}/\code{NaN} are not allowed.
#' @param px Optional numeric vector of one-year survival probabilities \eqn{p_x}.
#'   \code{NA} values are allowed, but \code{Inf}/\code{NaN} are not allowed.
#'   If provided, \code{qx = 1 - px}.
#' @param mx Optional numeric vector of central death rates \eqn{m_x}.
#'   \code{NA} values are allowed, but \code{Inf}/\code{NaN} are not allowed.
#'   If provided, converted to \code{qx} using \code{ax}.
#' @param radix Optional positive scalar. Required if building \code{lx} from
#'   (\code{qx}/\code{px}/\code{mx}) and \code{lx} is not provided.
#' @param omega Optional integer limiting age. If \code{omega < max(x)}, the
#'   table is truncated to \code{omega}. If \code{omega > max(x)}, an error is
#'   raised (the function will not invent missing ages).
#' @param close Logical. If \code{TRUE} (default), closes the table at
#'   \code{omega} (forces terminal conditions).
#' @param ax Scalar in \code{[0,1]}. Average fraction of the year lived by
#'   those who die in the interval \eqn{[x, x+1)}. Under UDD (Finan, Sec. 24.1),
#'   \code{ax = 0.5}. Under constant force, \eqn{a_x = 1/\mu - 1/(\exp(\mu)-1)}.
#'   At the terminal age with \code{close = TRUE}, \code{mx} equals
#'   \eqn{1/(1 - a_x)}, which is 2 for \code{ax = 0.5}. Default is \code{0.5}.
#' @param type Character. \code{"ultimate"} or \code{"select"} (metadata).
#'   Stored as an attribute and used by downstream functions.
#' @param frac Character. \code{"UDD"}, \code{"CF"}, or \code{"Balducci"}
#'   (metadata). Stored as an attribute and used by fractional-age functions
#'   such as \code{\link{t_px}}.
#' @param check Logical. If \code{TRUE} (default), performs strict validity
#'   and consistency checks.
#' @param tol Numeric tolerance for integer checks and consistency checks.
#'
#' @return A \code{tibble} with class
#'   \code{c("lifetable","tbl_df","tbl","data.frame")} and columns:
#' \itemize{
#'   \item \code{x}: integer ages
#'   \item \code{lx}: survivors at exact age x
#'   \item \code{dx}: deaths in \eqn{[x, x+1)}
#'   \item \code{qx}: probability of death in \eqn{[x, x+1)}
#'   \item \code{px}: probability of survival to \eqn{x+1}
#'   \item \code{mx}: central death rate (derived using \code{ax}). At the
#'     terminal age with \code{close = TRUE}, \code{mx} equals
#'     \eqn{1/(1 - a_x)} and may be \code{Inf} if \code{ax = 1}.
#' }
#' Attributes include: \code{radix}, \code{omega}, \code{type}, \code{frac},
#' \code{closed}, \code{ax}.
#'
#' @seealso
#' \code{\link{km_lifetable}} for Kaplan--Meier construction,
#' \code{\link{t_px}} and \code{\link{t_qx}} for survival and death
#' probabilities (including fractional ages),
#' \code{\link{e_x}} for curtate and complete life expectancy,
#' \code{\link{annuity_x}} and \code{\link{insurance_x}} for life
#' contingency valuations that consume a life table.
#'
#' @examples
#' # Example 1: build from lx (Finan, Section 22 style)
#' x  <- 0:5
#' lx <- c(100000, 99500, 99000, 98200, 97000, 95000)
#' lt1 <- lifetable(x = x, lx = lx, omega = 5, close = TRUE)
#' lt1
#'
#' # Example 2: build from qx (radix required)
#' qx <- c(0.005, 0.005, 0.008, 0.012, 0.020, 1)
#' lt2 <- lifetable(x = x, qx = qx, radix = 100000, omega = 5, close = TRUE)
#' lt2
#'
#' # Example 3: build from px
#' px <- 1 - c(0.005, 0.005, 0.008, 0.012, 0.020, 1)
#' lt3 <- lifetable(x = x, px = px, radix = 100000, omega = 5, close = TRUE)
#' lt3
#'
#' # Example 4: build from mx
#' mx <- c(0.005, 0.006, 0.008, 0.012, 0.020, 0.030)
#' lt4 <- lifetable(x = x, mx = mx, radix = 100000, omega = 5, close = TRUE, ax = 0.5)
#' lt4
#'
#' # Example 5: truncate to a smaller omega
#' lt5 <- lifetable(x = 0:10, lx = 100000 * exp(-0.01 * (0:10)), omega = 7, close = TRUE)
#' lt5
#'
#' # Example 6: Finan Example 22.1 - exponential survival s(x) = exp(-0.005x)
#' lt_exp <- lifetable(
#'   x  = 0:7,
#'   lx = 1000 * exp(-0.005 * (0:7)),
#'   close = TRUE
#' )
#' lt_exp
#'
#' # Example 7: verify survival identity (Finan, Section 22)
#' # 2_p_2 = l_4 / l_2 = 97000 / 99000
#' lt1$lx[lt1$x == 4] / lt1$lx[lt1$x == 2]
#'
#' # Example 8: without closure - qx at omega is not forced to 1
#' lt_open <- lifetable(x = 0:3, lx = c(1000, 900, 750, 500), close = FALSE)
#' lt_open$qx  # last element is NA
#'
#' # Example 9: access table metadata
#' attr(lt1, "omega")   # 5
#' attr(lt1, "closed")  # TRUE
#' attr(lt1, "frac")    # "UDD"
#' attr(lt1, "ax")      # 0.5
#'
#' @export
lifetable <- function(
    x,
    lx = NULL,
    qx = NULL,
    px = NULL,
    mx = NULL,
    radix = NULL,
    omega = NULL,
    close = TRUE,
    ax = 0.5,
    type = c("ultimate", "select"),
    frac = c("UDD", "CF", "Balducci"),
    check = TRUE,
    tol = 1e-10
) {
  type <- match.arg(type)
  frac <- match.arg(frac)

  # --- Basic checks on x ---
  if (missing(x)) stop("`x` is required.")
  if (!is.numeric(x)) stop("`x` must be numeric.")
  x <- as.numeric(x)

  if (check) {
    if (any(!is.finite(x))) stop("`x` must be finite.")
    if (any(abs(x - round(x)) > tol)) stop("`x` must contain integer ages (annual table).")
    # Bug fix: use tolerance for consecutive-age check (floating-point safe)
    if (any(abs(diff(x) - 1) > tol)) stop("`x` must be consecutive annual ages (diff(x) == 1).")
  }

  # --- omega handling ---
  max_x <- max(x)
  if (is.null(omega)) {
    omega <- max_x
  } else {
    if (!is.numeric(omega) || length(omega) != 1L || !is.finite(omega)) stop("`omega` must be a single finite number.")
    if (abs(omega - round(omega)) > tol) stop("`omega` must be an integer age.")
    omega <- as.numeric(omega)

    if (omega > max_x + tol) {
      stop("`omega` cannot exceed max(x). Provide ages up to omega (or extend the table with a separate function).")
    }
    if (omega < max_x - tol) {
      keep <- x <= omega
      x <- x[keep]
      if (!is.null(lx)) lx <- lx[keep]
      if (!is.null(qx)) qx <- qx[keep]
      if (!is.null(px)) px <- px[keep]
      if (!is.null(mx)) mx <- mx[keep]
      max_x <- max(x)
    }
  }

  n <- length(x)
  if (n < 2) stop("Need at least two ages after applying omega.")

  # --- ax checks ---
  ax <- as.numeric(ax)
  if (check) {
    if (!is.finite(ax) || length(ax) != 1L) stop("`ax` must be a single finite number.")
    if (ax < 0 - tol || ax > 1 + tol) stop("`ax` must be in [0,1].")
  }

  # --- Normalize vectors / lengths ---
  norm_vec <- function(z, name) {
    if (is.null(z)) return(NULL)
    if (!is.numeric(z)) stop("`", name, "` must be numeric.")
    z <- as.numeric(z)
    if (length(z) != n) stop("`", name, "` must have the same length as `x`.")

    if (check) {
      # Allow NA, but reject Inf/NaN
      bad <- !is.finite(z) & !is.na(z)
      if (any(bad)) stop("`", name, "` must not contain Inf/NaN.")
    }
    z
  }

  lx <- norm_vec(lx, "lx")
  qx <- norm_vec(qx, "qx")
  px <- norm_vec(px, "px")
  mx <- norm_vec(mx, "mx")

  # --- Determine radix if needed ---
  if (is.null(radix)) {
    radix <- if (!is.null(lx)) lx[1] else NA_real_
  }
  if (check && !is.na(radix)) {
    if (!is.numeric(radix) || length(radix) != 1L || !is.finite(radix)) stop("`radix` must be a single finite number.")
    if (radix <= 0) stop("`radix` must be positive.")
  }

  # --- Build qx from px or mx if needed ---
  # Priority: lx > qx > px > mx
  if (is.null(qx) && !is.null(px)) {
    qx <- 1 - px
  }

  # mx -> qx: q = m / (1 + ax*m)   (Finan, Section 23.9 inverse)
  if (is.null(qx) && !is.null(mx)) {
    if (check && any(mx < -tol, na.rm = TRUE)) stop("`mx` must be nonnegative.")
    qx <- mx / (1 + ax * mx)
  }

  if (!is.null(qx) && is.null(px)) {
    px <- 1 - qx
  }

  # --- If lx is missing, build it from qx (needs radix) ---
  if (is.null(lx)) {
    if (is.null(qx) && is.null(px)) stop("Provide either `lx` or one of `qx`, `px`, `mx`.")
    if (is.na(radix)) stop("`radix` is required when building `lx` from (`qx`/`px`/`mx`).")

    # Vectorized construction: lx[k] = radix * prod(px[1], ..., px[k-1])
    lx <- c(radix, radix * cumprod(px[-n]))
  }

  if (check) {
    if (any(lx < -tol, na.rm = TRUE)) stop("`lx` must be nonnegative.")
    if (any(diff(lx) > tol, na.rm = TRUE)) stop("`lx` must be nonincreasing.")
  }

  # --- Derive dx from lx ---
  # Actuarial convention: l_{omega+1} = 0 (Finan, Section 22)
  lx_next <- c(lx[-1], 0)
  dx <- lx - lx_next

  # --- Derive / validate qx from lx for ages < omega ---
  if (is.null(qx)) {
    qx <- dx / lx
    px <- 1 - qx
  } else if (check) {
    idx <- which(seq_len(n) < n & lx > 0)
    qx_lx <- dx / lx
    if (any(abs(qx[idx] - qx_lx[idx]) > 10 * tol, na.rm = TRUE)) {
      stop("Inconsistent inputs: provided `qx` does not match `lx` via qx = (lx - l_{x+1})/lx.")
    }
  }

  # --- Actuarial closure at omega (default TRUE) ---
  if (check && abs(max(x) - omega) > tol) stop("Internal error: max(x) must equal omega after processing.")
  j <- which(x == omega)

  if (close) {
    dx[j] <- lx[j]
    qx[j] <- 1
    px[j] <- 0
  } else {
    if (is.na(qx[j])) {
      dx[j] <- NA_real_
      qx[j] <- NA_real_
      px[j] <- NA_real_
    } else {
      dx[j] <- lx[j] * qx[j]
      px[j] <- 1 - qx[j]
    }
  }

  # --- Compute mx from qx using ax: m = q / (1 - ax*q) ---
  # Under UDD (ax=0.5): mx = qx / (1 - 0.5*qx) (Finan, Section 24.1)
  # At omega with close=TRUE: qx=1, mx = 1/(1-ax). If ax=1, mx=Inf.
  mx_out <- ifelse(qx == 1 & abs(ax - 1) < tol, Inf, qx / (1 - ax * qx))

  out <- tibble::tibble(
    x  = as.integer(round(x)),
    lx = as.numeric(lx),
    dx = as.numeric(dx),
    qx = as.numeric(qx),
    px = as.numeric(px),
    mx = as.numeric(mx_out)
  )

  class(out) <- c("lifetable", class(out))

  attr(out, "radix")  <- out$lx[1]
  attr(out, "omega")  <- as.integer(round(omega))
  attr(out, "type")   <- type
  attr(out, "frac")   <- frac
  attr(out, "closed") <- isTRUE(close)
  attr(out, "ax")     <- ax

  if (check) {
    if (any(out$qx < -tol | out$qx > 1 + tol, na.rm = TRUE)) stop("`qx` must be in [0,1].")
    if (any(out$px < -tol | out$px > 1 + tol, na.rm = TRUE)) stop("`px` must be in [0,1].")
  }

  out
}
