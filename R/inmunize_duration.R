#' Duration-based immunization with multiple assets
#'
#' Computes asset weights that duration-immunize a stream of liabilities
#' using two or more assets. The method enforces:
#' \enumerate{
#'   \item Present value of assets = PV of liabilities;
#'   \item Macaulay duration of assets = Macaulay duration of liabilities.
#' }
#'
#' For exactly two assets, a closed-form solution is used. For three or more
#' assets, a minimum-norm solution is computed via linear algebra.
#'
#' @param liabilities Numeric vector with liability payments.
#' @param t_liabilities Numeric vector of the same length as \code{liabilities}
#'   with the times (in years or periods) at which each liability payment occurs.
#' @param pv_assets Numeric vector with present values (prices) of each asset
#'   evaluated at the same yield rate \code{i}.
#' @param duration_assets Numeric vector with the Macaulay duration of each
#'   asset, expressed in the same time units as \code{t_liabilities}.
#' @param i Yield rate used to discount the liabilities (effective per period).
#'
#' @return A tibble with:
#' \describe{
#'   \item{weight_asset}{Numeric vector of asset weights (amounts of each asset).}
#'   \item{PV_liabilities}{Present value of the liabilities.}
#'   \item{Duration_liabilities}{Macaulay duration of the liabilities.}
#'   \item{PV_assets}{Present value of the immunized asset portfolio.}
#'   \item{Duration_assets}{Macaulay duration of the asset portfolio.}
#'   \item{n_assets}{Number of assets used.}
#' }
#'
#' @details
#' Let \eqn{PV_L}{PV_L} and \eqn{D_L}{D_L} be the present value and Macaulay
#' duration of the liability stream at yield \eqn{i}. Let \eqn{PV_j}{PV_j} and
#' \eqn{D_j}{D_j} be the present value and duration of asset \eqn{j}. The
#' weights \eqn{w_j}{w_j} are chosen so that:
#'
#' \deqn{\sum_j w_j PV_j = PV_L}{sum(w_j * PV_j) = PV_L}
#' and
#' \deqn{\frac{\sum_j w_j PV_j D_j}{\sum_j w_j PV_j} = D_L}{sum(w_j * PV_j * D_j) / PV_L = D_L}
#'
#' For two assets, the closed-form solution is:
#' \deqn{w_1 = \frac{PV_L (D_L - D_2)}{PV_1 (D_1 - D_2)}, \quad w_2 = \frac{PV_L - w_1 PV_1}{PV_2}}{w1 = PV_L*(D_L - D2) / (PV1*(D1 - D2)), w2 = (PV_L - w1*PV1) / PV2}
#'
#' For \eqn{k \geq 3}{k >= 3} assets, the minimum-norm solution of the linear
#' system \eqn{A w = b} is computed, where \eqn{A} is a \eqn{2 \times k}{2 x k}
#' matrix with rows \eqn{(PV_1, \ldots, PV_k)}{(PV1, ..., PVk)} and
#' \eqn{(PV_1 D_1, \ldots, PV_k D_k)}{(PV1*D1, ..., PVk*Dk)}, and
#' \eqn{b = (PV_L,\; PV_L D_L)^T}{b = (PV_L, PV_L*D_L)}.
#'
#' @seealso \code{\link{immunize_duration_convexity}},
#'   \code{\link{bond_duration}}, \code{\link{bond_convexity}}
#'
#' @family immunization
#'
#' @examples
#' # Two-asset immunization
#' immunize_duration(
#'   liabilities = c(5000, 8000),
#'   t_liabilities = c(3, 7),
#'   pv_assets = c(100, 200),
#'   duration_assets = c(3, 7),
#'   i = 0.05
#' )
#'
#' # Three-asset immunization (minimum-norm)
#' immunize_duration(
#'   liabilities = c(5000, 8000),
#'   t_liabilities = c(3, 7),
#'   pv_assets = c(100, 150, 200),
#'   duration_assets = c(2, 5, 8),
#'   i = 0.05
#' )
#'
#' @export
immunize_duration <- function(
    liabilities,
    t_liabilities,
    pv_assets,
    duration_assets,
    i
) {
  # --- 1. Basic checks ---
  if (!is.numeric(liabilities) || !is.numeric(t_liabilities)) {
    stop("'liabilities' and 't_liabilities' must be numeric vectors.", call. = FALSE)
  }
  if (length(liabilities) != length(t_liabilities)) {
    stop("'liabilities' and 't_liabilities' must have the same length.", call. = FALSE)
  }
  if (!is.numeric(pv_assets) || !is.numeric(duration_assets)) {
    stop("'pv_assets' and 'duration_assets' must be numeric vectors.", call. = FALSE)
  }
  if (length(pv_assets) != length(duration_assets)) {
    stop("'pv_assets' and 'duration_assets' must have the same length.", call. = FALSE)
  }
  k <- length(pv_assets)
  if (k < 2L) {
    stop("At least two assets are required for duration immunization.", call. = FALSE)
  }
  if (!is.numeric(i) || length(i) != 1L || is.na(i)) {
    stop("'i' must be a numeric scalar (yield rate).", call. = FALSE)
  }

  # --- 2. Present value and duration of liabilities ---
  v <- 1 / (1 + i)
  pv_L <- sum(liabilities * v^t_liabilities)

  if (abs(pv_L) < .Machine$double.eps * 100) {
    stop("Present value of liabilities is zero; immunization is not meaningful.", call. = FALSE)
  }

  dur_L <- sum(t_liabilities * liabilities * v^t_liabilities) / pv_L

  # --- 3. Build linear system A w = b ---
  # Row 1: PV constraint
  # Row 2: PV * duration constraint
  A <- rbind(
    pv_assets,
    pv_assets * duration_assets
  )
  b <- c(pv_L, pv_L * dur_L)

  # --- 4. Solve for weights ---
  if (k == 2L) {
    # Closed-form solution for two assets
    pv1 <- pv_assets[1L]
    pv2 <- pv_assets[2L]
    d1  <- duration_assets[1L]
    d2  <- duration_assets[2L]

    den <- pv1 * (d1 - d2)
    if (abs(den) < .Machine$double.eps * 100) {
      stop("Degenerate case: assets have equal duration; cannot solve uniquely.", call. = FALSE)
    }

    w1 <- pv_L * (dur_L - d2) / den
    w2 <- (pv_L - w1 * pv1) / pv2
    weight_asset <- c(w1, w2)
  } else {
    # Minimum-norm solution for k >= 3 assets:
    # w = t(A) %*% solve(A %*% t(A), b)
    At <- t(A)
    M <- A %*% At
    M_solve <- tryCatch(
      solve(M, b),
      error = function(e) {
        stop("System matrix is singular; cannot compute minimum-norm solution.", call. = FALSE)
      }
    )
    weight_asset <- as.numeric(At %*% M_solve)
  }

  # --- 5. Check resulting portfolio ---
  pv_A  <- sum(weight_asset * pv_assets)
  dur_A <- sum(weight_asset * pv_assets * duration_assets) / pv_A

  tibble::tibble(
    weight_asset = weight_asset,
    PV_liabilities = pv_L,
    Duration_liabilities = dur_L,
    PV_assets = pv_A,
    Duration_assets = dur_A,
    n_assets = k
  )
}
