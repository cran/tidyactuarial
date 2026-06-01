#' Pure endowment (discounted survival): \eqn{{}_tE_x}
#'
#' Computes the actuarial present value of a pure endowment, i.e., the expected
#' present value of a payment of 1 made at time \eqn{t} if and only if a life
#' aged \eqn{x} survives to age \eqn{x + t}:
#' \deqn{{}_tE_x = v^t \cdot {}_tp_x.}
#'
#' @param lt A lifetable object as produced by \code{\link{lifetable}}.
#'   Must contain columns \code{x} and \code{lx}.
#' @param x Integer age(s) at which the endowment starts.
#' @param t Nonnegative numeric duration(s) in years. Fractional durations are
#'   allowed and are handled through \code{\link{t_px}}.
#' @param i Annual interest-rate input(s).
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param frac Fractional-age assumption passed to \code{\link{t_px}}:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"} (alias of CF), or
#'   \code{"Balducci"}. If not specified and \code{lt} carries a \code{frac}
#'   attribute, that value is used.
#' @param tidy Logical. If \code{TRUE}, returns a tibble with columns
#'   \code{x}, \code{t}, \code{i}, \code{i_type}, \code{m}, \code{i_effective},
#'   \code{frac}, and \code{nEx}.
#' @param check Logical. If \code{TRUE}, performs validity checks.
#' @param tol Numeric tolerance for integer checks on \code{x}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age,
#' \code{t} is the duration, \code{i} is the interest-rate input,
#' \code{i_type} is the interest-rate type, and \code{m} is the conversion
#' frequency for nominal rates.
#'
#' The pure endowment is a fundamental building block in life contingency
#' mathematics. It serves as the actuarial discount factor, combining financial
#' discounting with mortality:
#' \deqn{{}_tE_x = v^t \cdot {}_tp_x.}
#'
#' The interest-rate input is converted to an annual effective rate before
#' applying the discount factor:
#' \deqn{v = (1+i_e)^{-1}.}
#'
#' Key identities involving \eqn{{}_nE_x}:
#' \itemize{
#'   \item Actuarial accumulated value:
#'     \eqn{\ddot{s}_{x:\overline{n}|} = \ddot{a}_{x:\overline{n}|} /
#'     {}_nE_x}.
#'   \item Endowment insurance decomposition:
#'     \eqn{A_{x:\overline{n}|} = A^1_{x:\overline{n}|} + {}_nE_x}.
#'   \item Deferred annuity:
#'     \eqn{{}_{n|}\ddot{a}_x = {}_nE_x \cdot \ddot{a}_{x+n}}.
#' }
#'
#' For a constant force of mortality \eqn{\mu} and force of interest
#' \eqn{\delta}:
#' \deqn{{}_nE_x = e^{-n(\mu + \delta)}.}
#'
#' The variance of the pure endowment random variable is:
#' \deqn{\mathrm{Var}(Z) = v^{2n} \cdot {}_np_x \cdot {}_nq_x.}
#'
#' @return Numeric vector of \eqn{{}_tE_x}, or a tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{t_px}} for survival probabilities without discounting,
#'   \code{\link{insurance_x}} for insurance APVs, \code{\link{annuity_x}} for
#'   life annuity APVs.
#'
#' @examples
#' x  <- 0:5
#' lx <- c(100000, 99500, 99000, 98200, 97000, 95000)
#' lt <- lifetable(x = x, lx = lx, omega = 5, close = TRUE)
#'
#' # Basic pure endowment: 3_E_0 = v^3 * 3_p_0
#' t_Ex(lt, x = 0, t = 3, i = 0.06)
#'
#' # Verify manually:
#' (1.06)^(-3) * t_px(lt, x = 0, t = 3)
#'
#' # Constant force of interest
#' lt_exp <- lifetable(x = 0:50, lx = 100000 * exp(-0.05 * (0:50)))
#' t_Ex(
#'   lt_exp,
#'   x = 30,
#'   t = 10,
#'   i = 0.10,
#'   i_type = "force"
#' )
#' exp(-1.5)
#'
#' # Nominal annual interest rate convertible monthly
#' t_Ex(
#'   lt,
#'   x = 0,
#'   t = 3,
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12
#' )
#'
#' # Vectorized: multiple ages at once
#' t_Ex(lt, x = c(0, 1, 2), t = 3, i = 0.05, tidy = TRUE)
#'
#' # Use in a tidy pipeline
#' if (requireNamespace("dplyr", quietly = TRUE)) {
#'   tibble::tibble(x = 0:4, t = c(5, 4, 3, 2, 1)) |>
#'     dplyr::mutate(pure_endow = t_Ex(lt, x = x, t = t, i = 0.06))
#' }
#'
#' @export
t_Ex <- function(
    lt,
    x,
    t,
    i,
    i_type = "effective",
    m = 1L,
    frac,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {

  # --- checks ---
  if (missing(lt)) stop("`lt` is required.", call. = FALSE)

  if (!("x" %in% names(lt)) || !("lx" %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
  }

  if (missing(x) || missing(t) || missing(i)) {
    stop("`x`, `t`, and `i` are required.", call. = FALSE)
  }

  x <- as.numeric(x)
  t <- as.numeric(t)
  i <- as.numeric(i)
  m <- as.numeric(m)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  if (check) {
    if (any(!is.finite(x))) stop("`x` must be finite.", call. = FALSE)
    if (any(abs(x - round(x)) > tol)) {
      stop("`x` must be integer ages.", call. = FALSE)
    }

    if (any(!is.finite(t))) stop("`t` must be finite.", call. = FALSE)
    if (any(t < 0)) stop("`t` must be nonnegative.", call. = FALSE)

    if (any(!is.finite(i))) stop("`i` must be finite.", call. = FALSE)

    if (!is.character(i_type)) {
      stop("`i_type` must be a character vector.", call. = FALSE)
    }

    valid_i_type <- c(
      "effective",
      "nominal_interest",
      "nominal_discount",
      "force"
    )

    bad_i_type <- !is.na(i_type) & !(i_type %in% valid_i_type)
    if (any(bad_i_type)) {
      stop(
        "`i_type` must contain only: ",
        paste(sprintf("'%s'", valid_i_type), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (any(!is.finite(m))) stop("`m` must be finite.", call. = FALSE)
    if (any(abs(m - round(m)) > tol)) {
      stop("`m` must contain integer values.", call. = FALSE)
    }
    if (any(m < 1)) stop("`m` must contain positive integers.", call. = FALSE)
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))

  # --- recycle x, t, i, i_type, m to common length ---
  lengths <- c(length(x), length(t), length(i), length(i_type), length(m))
  n <- max(lengths)

  if (!all(lengths %in% c(1L, n))) {
    stop(
      "`x`, `t`, `i`, `i_type`, and `m` must have the same length or length 1.",
      call. = FALSE
    )
  }

  if (length(x) == 1L && n > 1L) x <- rep(x, n)
  if (length(t) == 1L && n > 1L) t <- rep(t, n)
  if (length(i) == 1L && n > 1L) i <- rep(i, n)
  if (length(i_type) == 1L && n > 1L) i_type <- rep(i_type, n)
  if (length(m) == 1L && n > 1L) m <- rep(m, n)

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (any(is.na(i_effective)) ||
      any(!is.finite(i_effective)) ||
      any(i_effective <= -1)) {
    stop(
      "The standardized annual effective interest rates must be greater than -1.",
      call. = FALSE
    )
  }

  # --- compute survival probability via t_px ---
  if (missing(frac)) {
    tpx <- t_px(
      lt = lt,
      x = x,
      t = t,
      tidy = FALSE,
      check = FALSE,
      tol = tol
    )

    lt_frac <- attr(lt, "frac")

    frac_used <- if (!is.null(lt_frac) && lt_frac %in% c("UDD", "CF", "Balducci")) {
      lt_frac
    } else {
      "UDD"
    }
  } else {
    frac_used <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac_used == "CML") frac_used <- "CF"

    tpx <- t_px(
      lt = lt,
      x = x,
      t = t,
      frac = frac_used,
      tidy = FALSE,
      check = FALSE,
      tol = tol
    )
  }

  # --- pure endowment: tEx = v^t * t_px ---
  v <- (1 + i_effective)^(-t)
  nEx <- v * tpx

  # Preserve the previous behavior: clamp to [0, 1].
  nEx <- pmin(pmax(nEx, 0), 1)

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x,
      t = t,
      i = i,
      i_type = i_type,
      m = m,
      i_effective = i_effective,
      frac = frac_used,
      nEx = nEx
    ))
  }

  nEx
}
