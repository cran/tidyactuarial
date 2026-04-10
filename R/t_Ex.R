#' Pure endowment (discounted survival): \eqn{{}_tE_x}
#'
#' Computes the actuarial present value of a pure endowment, i.e., the
#' expected present value of a payment of 1 made at time \eqn{t} if and only
#' if a life aged \eqn{x} survives to age \eqn{x + t}:
#' \deqn{{}_tE_x = v^t \cdot {}_tp_x = (1 + i)^{-t} \cdot
#'   \frac{\ell_{x+t}}{\ell_x}.}
#'
#' @param lt A lifetable object as produced by \code{\link{lifetable}}.
#'   Must contain columns \code{x} and \code{lx}.
#' @param x Integer age(s) at which the endowment starts.
#' @param t Nonnegative numeric duration(s) in years (can be fractional).
#' @param i Annual effective interest rate(s). Must satisfy \code{i > -1}.
#' @param frac Fractional-age assumption passed to \code{\link{t_px}}:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"} (alias of CF), or
#'   \code{"Balducci"}. If not specified and \code{lt} carries a \code{frac}
#'   attribute, that value is used.
#' @param tidy Logical. If \code{TRUE}, returns a tibble with columns
#'   \code{x}, \code{t}, \code{i}, \code{frac}, \code{nEx}.
#' @param check Logical. If \code{TRUE}, performs validity checks.
#' @param tol Numeric tolerance for integer checks on \code{x}.
#'
#' @details
#' The pure endowment is a fundamental building block in life contingency
#' mathematics (Finan, Section 26.3.1). It serves as the actuarial discount
#' factor, combining financial discounting with mortality:
#' \deqn{{}_nE_x = v^n \cdot {}_np_x.}
#'
#' Key identities involving \eqn{{}_nE_x}:
#' \itemize{
#'   \item Actuarial accumulated value:
#'     \eqn{\ddot{s}_{x:\overline{n}|} = \ddot{a}_{x:\overline{n}|} /
#'     {}_nE_x} (Finan, Sec. 34).
#'   \item Endowment insurance decomposition:
#'     \eqn{A_{x:\overline{n}|} = A^1_{x:\overline{n}|} + {}_nE_x}
#'     (Finan, Sec. 26.3.2).
#'   \item Deferred annuity:
#'     \eqn{{}_{n|}\ddot{a}_x = {}_nE_x \cdot \ddot{a}_{x+n}}
#'     (Finan, Sec. 35).
#' }
#'
#' For a constant force of mortality \eqn{\mu} and force of interest
#' \eqn{\delta}:
#' \deqn{{}_nE_x = e^{-n(\mu + \delta)}}
#' (Finan, Example 26.14).
#'
#' The variance of the pure endowment random variable is
#' (Finan, Sec. 26.3.1):
#' \deqn{\mathrm{Var}(\bar{Z}^{\phantom{1}}_{\phantom{1}x:\overline{n}|})
#'   = v^{2n} \cdot {}_np_x \cdot {}_nq_x.}
#'
#' @return Numeric vector of \eqn{{}_tE_x}, or a tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{t_px}} for survival probabilities (without
#'   discounting), \code{\link{insurance_x}} for term, whole life, and
#'   endowment insurance APVs, \code{\link{annuity_x}} for life annuity
#'   APVs that use \eqn{{}_nE_x} internally.
#'
#' @examples
#' x  <- 0:5
#' lx <- c(100000, 99500, 99000, 98200, 97000, 95000)
#' lt <- lifetable(x = x, lx = lx, omega = 5, close = TRUE)
#'
#' # Basic pure endowment: 3_E_0 = v^3 * 3_p_0
#' t_Ex(lt, x = 0, t = 3, i = 0.06)
#' # Verify manually:
#' (1.06)^(-3) * t_px(lt, x = 0, t = 3)
#'
#' # Finan Example 26.14 style: constant force mu, delta
#' # For exponential survival with mu = 0.05, delta = 0.10:
#' # 10_E_30 = exp(-10*(0.05 + 0.10)) = exp(-1.5)
#' lt_exp <- lifetable(x = 0:50, lx = 100000 * exp(-0.05 * (0:50)))
#' t_Ex(lt_exp, x = 30, t = 10, i = exp(0.10) - 1)
#' exp(-1.5)  # theoretical value
#'
#' # Finan Problem 26.16: 5-year pure endowment for (30), i = 6%
#' lt_ilt <- lifetable(
#'   x  = 30:35,
#'   lx = c(9501381, 9486854, 9471591, 9455522, 9438571, 9420657)
#' )
#' t_Ex(lt_ilt, x = 30, t = 5, i = 0.06)
#' # = v^5 * l_35 / l_30 = (1.06)^(-5) * 9420657/9501381
#'
#' # Vectorized: multiple ages at once
#' t_Ex(lt, x = c(0, 1, 2), t = 3, i = 0.05, tidy = TRUE)
#'
#' # Use in a tidy pipeline
#' if (requireNamespace("dplyr", quietly = TRUE)) {
#'   tibble::tibble(age = 0:4, term = c(5, 4, 3, 2, 1)) |>
#'     dplyr::mutate(pure_endow = t_Ex(lt, x = age, t = term, i = 0.06))
#' }
#'
#' @export
t_Ex <- function(
    lt,
    x,
    t,
    i,
    frac,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {

  # --- checks ---
  if (missing(lt)) stop("`lt` is required.")
  if (!("x" %in% names(lt)) || !("lx" %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.")
  }
  if (missing(x) || missing(t) || missing(i)) {
    stop("`x`, `t`, and `i` are required.")
  }

  x <- as.numeric(x)
  t <- as.numeric(t)
  i <- as.numeric(i)

  if (check) {
    if (any(!is.finite(x))) stop("`x` must be finite.")
    if (any(abs(x - round(x)) > tol)) stop("`x` must be integer ages.")
    if (any(!is.finite(t))) stop("`t` must be finite.")
    if (any(t < 0)) stop("`t` must be nonnegative.")
    if (any(!is.finite(i))) stop("`i` must be finite.")
    if (any(i <= -1)) stop("`i` must be greater than -1.")
  }

  x <- as.integer(round(x))

  # --- recycle x, t, i to common length ---
  lengths <- c(length(x), length(t), length(i))
  n <- max(lengths)
  if (!all(lengths %in% c(1L, n))) {
    stop("`x`, `t`, and `i` must have the same length or length 1.")
  }
  if (length(x) == 1L && n > 1L) x <- rep(x, n)
  if (length(t) == 1L && n > 1L) t <- rep(t, n)
  if (length(i) == 1L && n > 1L) i <- rep(i, n)

  # --- compute survival probability via t_px ---
  if (missing(frac)) {
    tpx <- t_px(lt = lt, x = x, t = t, tidy = FALSE, check = FALSE, tol = tol)
    # recover frac for tidy output
    lt_frac <- attr(lt, "frac")
    frac_used <- if (!is.null(lt_frac) && lt_frac %in% c("UDD", "CF", "Balducci")) {
      lt_frac
    } else {
      "UDD"
    }
  } else {
    frac_used <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac_used == "CML") frac_used <- "CF"
    tpx <- t_px(lt = lt, x = x, t = t, frac = frac_used,
                 tidy = FALSE, check = FALSE, tol = tol)
  }

  # --- pure endowment: nEx = v^t * t_px (Finan, Section 26.3.1) ---
  v <- (1 + i)^(-t)
  nEx <- v * tpx

  # Clamp to [0, 1] (it's a probability-weighted discount)
  nEx <- pmin(pmax(nEx, 0), 1)

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x    = x,
      t    = t,
      i    = i,
      frac = frac_used,
      nEx  = nEx
    ))
  }
  nEx
}
