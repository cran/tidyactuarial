#' Duration-based immunization with multiple assets
#'
#' Computes asset weights that duration-immunize a stream of liabilities using
#' two or more assets, using compact actuarial notation.
#'
#' The method enforces:
#' \enumerate{
#'   \item present value of assets = present value of liabilities;
#'   \item Macaulay duration of assets = Macaulay duration of liabilities.
#' }
#'
#' For exactly two assets, a closed-form solution is used. For three or more
#' assets, a minimum-norm solution is computed by linear algebra.
#'
#' @param L Numeric vector with liability payments.
#' @param t Numeric vector of the same length as \code{L}, giving the times at
#'   which each liability payment occurs.
#' @param P Numeric vector with present values or prices of the immunizing
#'   assets, evaluated on the same yield basis.
#' @param D Numeric vector with the Macaulay duration of each asset, expressed
#'   in the same time units as \code{t}.
#' @param i Numeric scalar. Interest-rate input used to discount the
#'   liabilities.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#'
#' @return A tibble with:
#' \describe{
#'   \item{w}{Numeric vector of asset weights or units.}
#'   \item{PV_L}{Present value of the liabilities.}
#'   \item{D_L}{Macaulay duration of the liabilities.}
#'   \item{PV_A}{Present value of the immunized asset portfolio.}
#'   \item{D_A}{Macaulay duration of the asset portfolio.}
#'   \item{n_assets}{Number of assets used.}
#' }
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{L} denotes liabilities, \code{t} denotes payment
#' times, \code{P} denotes asset prices or present values, \code{D} denotes
#' asset durations, \code{i} denotes the interest rate, \code{i_type} denotes
#' the interest-rate type, and \code{m} denotes the conversion frequency for
#' nominal rates.
#'
#' Let \eqn{PV_L} and \eqn{D_L} be the present value and Macaulay duration of
#' the liability stream at yield \eqn{i}. Let \eqn{P_j} and \eqn{D_j} be the
#' price and duration of asset \eqn{j}. The weights \eqn{w_j} are chosen so
#' that:
#'
#' \deqn{\sum_j w_j P_j = PV_L}
#'
#' and
#'
#' \deqn{
#'   \frac{\sum_j w_j P_j D_j}{\sum_j w_j P_j} = D_L.
#' }
#'
#' For two assets, the closed-form solution is:
#'
#' \deqn{
#' w_1 =
#' \frac{PV_L(D_L - D_2)}{P_1(D_1 - D_2)}, \qquad
#' w_2 = \frac{PV_L - w_1 P_1}{P_2}.
#' }
#'
#' For three or more assets, the minimum-norm solution of the linear system
#' \eqn{Aw = b} is computed, where \eqn{A} is a \eqn{2 \times r} matrix with
#' rows
#'
#' \deqn{(P_1,\ldots,P_r)}
#'
#' and
#'
#' \deqn{(P_1D_1,\ldots,P_rD_r),}
#'
#' and
#'
#' \deqn{b = (PV_L, PV_LD_L)^T.}
#'
#' @seealso \code{\link{immunize_duration_convexity}},
#'   \code{\link{bond_duration}}, \code{\link{bond_convexity}}
#'
#' @family immunization
#'
#' @examples
#' # Two-asset immunization
#' immunize_duration(
#'   L = c(5000, 8000),
#'   t = c(3, 7),
#'   P = c(100, 200),
#'   D = c(3, 7),
#'   i = 0.05
#' )
#'
#' # Three-asset immunization: minimum-norm solution
#' immunize_duration(
#'   L = c(5000, 8000),
#'   t = c(3, 7),
#'   P = c(100, 150, 200),
#'   D = c(2, 5, 8),
#'   i = 0.05
#' )
#'
#' # Nominal annual interest rate convertible monthly
#' immunize_duration(
#'   L = c(5000, 8000),
#'   t = c(3, 7),
#'   P = c(100, 200),
#'   D = c(3, 7),
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12
#' )
#'
#' @export
immunize_duration <- function(
    L,
    t,
    P,
    D,
    i,
    i_type = "effective",
    m = 1L
) {
  # --- Basic checks ---
  if (missing(L) || !is.numeric(L)) {
    stop("`L` must be a numeric vector of liabilities.", call. = FALSE)
  }

  if (missing(t) || !is.numeric(t)) {
    stop("`t` must be a numeric vector of liability times.", call. = FALSE)
  }

  if (length(L) != length(t)) {
    stop("`L` and `t` must have the same length.", call. = FALSE)
  }

  if (length(L) == 0L) {
    stop("`L` and `t` must have positive length.", call. = FALSE)
  }

  if (any(is.na(L)) || any(!is.finite(L))) {
    stop("`L` must contain only finite numeric values.", call. = FALSE)
  }

  if (any(is.na(t)) || any(!is.finite(t)) || any(t < 0)) {
    stop("`t` must contain only finite values greater than or equal to 0.", call. = FALSE)
  }

  if (missing(P) || !is.numeric(P)) {
    stop("`P` must be a numeric vector of asset prices or present values.", call. = FALSE)
  }

  if (missing(D) || !is.numeric(D)) {
    stop("`D` must be a numeric vector of asset durations.", call. = FALSE)
  }

  if (length(P) != length(D)) {
    stop("`P` and `D` must have the same length.", call. = FALSE)
  }

  n_assets <- length(P)

  if (n_assets < 2L) {
    stop("At least two assets are required for duration immunization.", call. = FALSE)
  }

  if (any(is.na(P)) || any(!is.finite(P)) || any(P <= 0)) {
    stop("`P` must contain only finite positive values.", call. = FALSE)
  }

  if (any(is.na(D)) || any(!is.finite(D))) {
    stop("`D` must contain only finite numeric values.", call. = FALSE)
  }

  if (missing(i) || !is.numeric(i) || length(i) != 1L ||
      is.na(i) || !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
    stop("`i_type` must be a single character string.", call. = FALSE)
  }

  valid_i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  if (!i_type %in% valid_i_type) {
    stop(
      "`i_type` must be one of: ",
      paste(sprintf("'%s'", valid_i_type), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      !is.finite(m) || m < 1 || abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  m <- as.integer(round(m))

  i_effective <- standardize_interest(
    type = i_type,
    rate = i,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  # --- Present value and duration of liabilities ---
  v <- 1 / (1 + i_effective)
  pv_L <- sum(L * v^t)

  if (abs(pv_L) < .Machine$double.eps * 100) {
    stop("Present value of liabilities is zero; immunization is not meaningful.", call. = FALSE)
  }

  dur_L <- sum(t * L * v^t) / pv_L

  # --- Build linear system A w = b ---
  # Row 1: PV constraint
  # Row 2: PV * duration constraint
  A <- rbind(
    P,
    P * D
  )
  b <- c(pv_L, pv_L * dur_L)

  # --- Solve for weights ---
  if (n_assets == 2L) {
    # Closed-form solution for two assets
    P1 <- P[1L]
    P2 <- P[2L]
    D1 <- D[1L]
    D2 <- D[2L]

    den <- P1 * (D1 - D2)
    if (abs(den) < .Machine$double.eps * 100) {
      stop("Degenerate case: assets have equal duration; cannot solve uniquely.", call. = FALSE)
    }

    w1 <- pv_L * (dur_L - D2) / den
    w2 <- (pv_L - w1 * P1) / P2
    w <- c(w1, w2)
  } else {
    # Minimum-norm solution for n_assets >= 3:
    # w = t(A) %*% solve(A %*% t(A), b)
    At <- t(A)
    M <- A %*% At
    M_solve <- tryCatch(
      solve(M, b),
      error = function(e) {
        stop("System matrix is singular; cannot compute minimum-norm solution.", call. = FALSE)
      }
    )
    w <- as.numeric(At %*% M_solve)
  }

  # --- Check resulting portfolio ---
  pv_A <- sum(w * P)
  dur_A <- sum(w * P * D) / pv_A

  tibble::tibble(
    w = w,
    PV_L = pv_L,
    D_L = dur_L,
    PV_A = pv_A,
    D_A = dur_A,
    n_assets = n_assets
  )
}
