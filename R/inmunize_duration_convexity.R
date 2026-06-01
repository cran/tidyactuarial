#' Duration and convexity immunization with multiple assets
#'
#' Computes asset weights that immunize a stream of liabilities using three or
#' more assets, enforcing:
#' \enumerate{
#'   \item present value of assets = present value of liabilities;
#'   \item Macaulay duration of assets = Macaulay duration of liabilities;
#'   \item convexity of assets = convexity of liabilities.
#' }
#'
#' For exactly three assets, the system is solved directly. For four or more
#' assets, a minimum-norm solution is computed by linear algebra.
#'
#' @param L Numeric vector with liability payments.
#' @param t Numeric vector of the same length as \code{L}, giving the times at
#'   which each liability payment occurs.
#' @param P Numeric vector with present values or prices of the immunizing
#'   assets, evaluated on the same yield basis.
#' @param D Numeric vector with the Macaulay duration of each asset, expressed
#'   in the same time units as \code{t}.
#' @param C Numeric vector with the discrete convexity of each asset, expressed
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
#'   \item{C_L}{Discrete convexity of the liabilities.}
#'   \item{PV_A}{Present value of the asset portfolio.}
#'   \item{D_A}{Macaulay duration of the asset portfolio.}
#'   \item{C_A}{Discrete convexity of the asset portfolio.}
#'   \item{n_assets}{Number of assets used.}
#' }
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{L} denotes liabilities, \code{t} denotes payment
#' times, \code{P} denotes asset prices or present values, \code{D} denotes
#' asset durations, \code{C} denotes asset convexities, \code{i} denotes the
#' interest rate, \code{i_type} denotes the interest-rate type, and \code{m}
#' denotes the conversion frequency for nominal rates.
#'
#' Let \eqn{PV_L}, \eqn{D_L}, and \eqn{C_L} be the present value, Macaulay
#' duration, and discrete convexity of the liability stream at yield \eqn{i}.
#' The discrete convexity of the liabilities is computed as:
#' \deqn{
#' C_L =
#' \frac{\sum_t L_t\,t(t+1)\,v^{t+2}}{PV_L},
#' }
#' where \eqn{v = 1/(1+i)} after converting \code{i} to the equivalent
#' effective rate.
#'
#' The weights \eqn{w_j} satisfy the \eqn{3 \times r} system \eqn{Aw = b},
#' where the rows of \eqn{A} are
#' \deqn{(P_1,\ldots,P_r),}
#' \deqn{(P_1D_1,\ldots,P_rD_r),}
#' and
#' \deqn{(P_1C_1,\ldots,P_rC_r),}
#' and
#' \deqn{
#' b = (PV_L,\; PV_LD_L,\; PV_LC_L)^T.
#' }
#'
#' For three assets, the system is square and solved directly. For four or more
#' assets, the minimum-norm solution
#' \deqn{
#' w = A^T(AA^T)^{-1}b
#' }
#' is computed.
#'
#' @seealso \code{\link{immunize_duration}}, \code{\link{bond_duration}},
#'   \code{\link{bond_convexity}}
#'
#' @family immunization
#'
#' @examples
#' # Three-asset immunization
#' immunize_duration_convexity(
#'   L = c(5000, 8000, 10000),
#'   t = c(3, 5, 7),
#'   P = c(100, 150, 200),
#'   D = c(2, 5, 8),
#'   C = c(6, 30, 72),
#'   i = 0.05
#' )
#'
#' # Four-asset immunization: minimum-norm solution
#' immunize_duration_convexity(
#'   L = c(5000, 8000, 10000),
#'   t = c(3, 5, 7),
#'   P = c(100, 120, 150, 200),
#'   D = c(2, 4, 6, 8),
#'   C = c(6, 20, 42, 72),
#'   i = 0.05
#' )
#'
#' # Nominal annual interest rate convertible monthly
#' immunize_duration_convexity(
#'   L = c(5000, 8000, 10000),
#'   t = c(3, 5, 7),
#'   P = c(100, 150, 200),
#'   D = c(2, 5, 8),
#'   C = c(6, 30, 72),
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12
#' )
#'
#' @export
immunize_duration_convexity <- function(
    L,
    t,
    P,
    D,
    C,
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

  if (missing(C) || !is.numeric(C)) {
    stop("`C` must be a numeric vector of asset convexities.", call. = FALSE)
  }

  if (length(P) != length(D) || length(P) != length(C)) {
    stop("`P`, `D`, and `C` must have the same length.", call. = FALSE)
  }

  n_assets <- length(P)

  if (n_assets < 3L) {
    stop("At least three assets are required for duration-convexity immunization.",
         call. = FALSE)
  }

  if (any(is.na(P)) || any(!is.finite(P)) || any(P <= 0)) {
    stop("`P` must contain only finite positive values.", call. = FALSE)
  }

  if (any(is.na(D)) || any(!is.finite(D))) {
    stop("`D` must contain only finite numeric values.", call. = FALSE)
  }

  if (any(is.na(C)) || any(!is.finite(C))) {
    stop("`C` must contain only finite numeric values.", call. = FALSE)
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

  # --- Liabilities: PV, duration, convexity ---
  v <- 1 / (1 + i_effective)
  disc <- v^t

  pv_L <- sum(L * disc)

  if (abs(pv_L) < .Machine$double.eps * 100) {
    stop("Present value of liabilities is zero; immunization is not meaningful.",
         call. = FALSE)
  }

  D_L <- sum(t * L * disc) / pv_L
  C_L <- sum(L * t * (t + 1) * v^(t + 2)) / pv_L

  # --- Build linear system A w = b (3 x n_assets) ---
  A <- rbind(
    P,
    P * D,
    P * C
  )
  b <- c(pv_L, pv_L * D_L, pv_L * C_L)

  # --- Solve for weights ---
  if (n_assets == 3L) {
    # 3 x 3 square system
    w <- tryCatch(
      as.numeric(solve(A, b)),
      error = function(e) {
        stop("System matrix is singular; cannot compute unique solution for three assets.",
             call. = FALSE)
      }
    )
  } else {
    # Minimum-norm solution for n_assets >= 4:
    # w = t(A) %*% solve(A %*% t(A), b)
    At <- t(A)
    M  <- A %*% At
    M_solve <- tryCatch(
      solve(M, b),
      error = function(e) {
        stop("System matrix is singular; cannot compute minimum-norm solution.",
             call. = FALSE)
      }
    )
    w <- as.numeric(At %*% M_solve)
  }

  # --- Portfolio measures ---
  pv_A <- sum(w * P)
  D_A  <- sum(w * P * D) / pv_A
  C_A  <- sum(w * P * C) / pv_A

  tibble::tibble(
    w = w,
    PV_L = pv_L,
    D_L = D_L,
    C_L = C_L,
    PV_A = pv_A,
    D_A = D_A,
    C_A = C_A,
    n_assets = n_assets
  )
}
