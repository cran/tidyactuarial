#' Duration and convexity immunization with multiple assets
#'
#' Computes asset weights that immunize a stream of liabilities using three or
#' more assets, enforcing:
#' \enumerate{
#'   \item Present value of assets = PV of liabilities;
#'   \item Macaulay duration of assets = Macaulay duration of liabilities;
#'   \item Convexity of assets = convexity of liabilities.
#' }
#'
#' For exactly three assets, the system is solved directly. For four or more
#' assets, a minimum-norm solution is computed via linear algebra.
#'
#' @param liabilities Numeric vector with liability payments.
#' @param t_liabilities Numeric vector with the times (in periods) at which
#'   each liability payment occurs. Must have the same length as
#'   \code{liabilities}.
#' @param pv_assets Numeric vector with present values (prices) of each asset
#'   evaluated at the same yield rate \code{i}.
#' @param duration_assets Numeric vector with the Macaulay duration of each
#'   asset, in the same time units as \code{t_liabilities}.
#' @param convexity_assets Numeric vector with the discrete convexity of each
#'   asset, in the same time units (periods) as \code{t_liabilities}.
#' @param i Yield rate used to discount the liabilities (effective per period).
#'
#' @return A tibble with:
#' \describe{
#'   \item{weight_asset}{Numeric vector of asset weights (amounts of each asset).}
#'   \item{PV_liabilities}{Present value of the liabilities.}
#'   \item{Duration_liabilities}{Macaulay duration of the liabilities.}
#'   \item{Convexity_liabilities}{Discrete convexity of the liabilities.}
#'   \item{PV_assets}{Present value of the asset portfolio.}
#'   \item{Duration_assets}{Macaulay duration of the asset portfolio.}
#'   \item{Convexity_assets}{Discrete convexity of the asset portfolio.}
#'   \item{n_assets}{Number of assets used.}
#' }
#'
#' @details
#' Let \eqn{PV_L}{PV_L}, \eqn{D_L}{D_L}, and \eqn{C_L}{C_L} be the present
#' value, Macaulay duration, and discrete convexity of the liability stream at
#' yield \eqn{i}. The discrete convexity of the liabilities is computed as:
#' \deqn{C_L = \frac{\sum_t L_t \, t(t+1) \, v^{t+2}}{PV_L}}{C_L = sum(L_t * t*(t+1) * v^(t+2)) / PV_L}
#' where \eqn{v = 1/(1+i)}.
#'
#' The weights \eqn{w_j}{w_j} satisfy the \eqn{3 \times k}{3 x k} system
#' \eqn{A w = b}, where the rows of \eqn{A} are
#' \eqn{(PV_j)}, \eqn{(PV_j D_j)}, \eqn{(PV_j C_j)}, and
#' \eqn{b = (PV_L,\; PV_L D_L,\; PV_L C_L)^T}{b = (PV_L, PV_L*D_L, PV_L*C_L)}.
#'
#' For \eqn{k = 3} assets, the system is square and solved directly.
#' For \eqn{k \geq 4}{k >= 4} assets, the minimum-norm solution
#' \eqn{w = A^T (A A^T)^{-1} b}{w = t(A) solve(A t(A), b)} is computed.
#'
#' @seealso \code{\link{immunize_duration}}, \code{\link{bond_duration}},
#'   \code{\link{bond_convexity}}
#'
#' @family immunization
#'
#' @examples
#' # Three-asset immunization (exact solution)
#' immunize_duration_convexity(
#'   liabilities = c(5000, 8000, 10000),
#'   t_liabilities = c(3, 5, 7),
#'   pv_assets = c(100, 150, 200),
#'   duration_assets = c(2, 5, 8),
#'   convexity_assets = c(6, 30, 72),
#'   i = 0.05
#' )
#'
#' # Four-asset immunization (minimum-norm)
#' immunize_duration_convexity(
#'   liabilities = c(5000, 8000, 10000),
#'   t_liabilities = c(3, 5, 7),
#'   pv_assets = c(100, 120, 150, 200),
#'   duration_assets = c(2, 4, 6, 8),
#'   convexity_assets = c(6, 20, 42, 72),
#'   i = 0.05
#' )
#'
#' @export
immunize_duration_convexity <- function(
    liabilities,
    t_liabilities,
    pv_assets,
    duration_assets,
    convexity_assets,
    i
) {
  # --- Basic checks ---
  if (!is.numeric(liabilities) || !is.numeric(t_liabilities)) {
    stop("'liabilities' and 't_liabilities' must be numeric vectors.", call. = FALSE)
  }
  if (length(liabilities) != length(t_liabilities)) {
    stop("'liabilities' and 't_liabilities' must have the same length.", call. = FALSE)
  }
  if (!is.numeric(pv_assets) || !is.numeric(duration_assets) ||
      !is.numeric(convexity_assets)) {
    stop("'pv_assets', 'duration_assets' and 'convexity_assets' must be numeric vectors.",
         call. = FALSE)
  }
  if (length(pv_assets) != length(duration_assets) ||
      length(pv_assets) != length(convexity_assets)) {
    stop("All asset vectors must have the same length.", call. = FALSE)
  }
  k <- length(pv_assets)
  if (k < 3L) {
    stop("At least three assets are required for duration-convexity immunization.",
         call. = FALSE)
  }
  if (!is.numeric(i) || length(i) != 1L || is.na(i)) {
    stop("'i' must be a numeric scalar (yield rate).", call. = FALSE)
  }

  # --- Liabilities: PV, duration, convexity ---
  v <- 1 / (1 + i)
  disc <- v^t_liabilities

  pv_L <- sum(liabilities * disc)
  if (abs(pv_L) < .Machine$double.eps * 100) {
    stop("Present value of liabilities is zero; immunization is not meaningful.",
         call. = FALSE)
  }

  dur_L  <- sum(t_liabilities * liabilities * disc) / pv_L
  conv_L <- sum(liabilities * t_liabilities * (t_liabilities + 1) * v^(t_liabilities + 2)) / pv_L

  # --- Build linear system A w = b (3 x k) ---
  A <- rbind(
    pv_assets,
    pv_assets * duration_assets,
    pv_assets * convexity_assets
  )
  b <- c(pv_L, pv_L * dur_L, pv_L * conv_L)

  # --- Solve for weights ---
  if (k == 3L) {
    # 3 x 3 square system
    weight_asset <- tryCatch(
      as.numeric(solve(A, b)),
      error = function(e) {
        stop("System matrix is singular; cannot compute unique solution for three assets.",
             call. = FALSE)
      }
    )
  } else {
    # Minimum-norm solution for k >= 4 assets:
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
    weight_asset <- as.numeric(At %*% M_solve)
  }

  # --- Portfolio measures ---
  pv_A   <- sum(weight_asset * pv_assets)
  dur_A  <- sum(weight_asset * pv_assets * duration_assets) / pv_A
  conv_A <- sum(weight_asset * pv_assets * convexity_assets) / pv_A

  tibble::tibble(
    weight_asset          = weight_asset,
    PV_liabilities        = pv_L,
    Duration_liabilities  = dur_L,
    Convexity_liabilities = conv_L,
    PV_assets             = pv_A,
    Duration_assets       = dur_A,
    Convexity_assets      = conv_A,
    n_assets              = k
  )
}
