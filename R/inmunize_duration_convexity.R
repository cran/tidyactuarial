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
#' @param C Numeric vector with each asset's relative second derivative
#'   \eqn{P_j''(i_e)/P_j(i_e)} with respect to the annual effective yield
#'   \eqn{i_e}, using the same time units as \code{t}. A coupon-period
#'   convexity divided by squared coupon frequency generally uses a different
#'   yield variable and must be converted before use here.
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
#' It is evaluated by singular value decomposition of the row-scaled
#' constraint matrix, without forming normal equations. The norm is over
#' asset units \code{w}, not market-value fractions; short positions are allowed.
#'
#' All asset sensitivities must use the same annual effective yield as the
#' liabilities, even when \code{i_type} specifies a nominal rate or a force.
#' For coupon frequency \eqn{k}, let \eqn{j=(1+i_e)^{1/k}-1}. If
#' \eqn{C_j=P''(j)/P} and \eqn{D_j^*=-P'(j)/P}, the required convexity is
#' \deqn{C_{i_e}=C_j[j'(i_e)]^2-D_j^*j''(i_e).}
#' The outputs \code{discrete_convexity_periods} from \code{bond_convexity()}
#' and \code{modified_duration_periods_j} from \code{bond_duration()} give
#' \eqn{C_j} and \eqn{D_j^*}. Alternatively, use \code{convexity_cash_flow()}
#' on the asset cash flows at the annual effective yield.
#'
#' Matching convexity makes the surplus's second derivative zero at the
#' base yield. This differs from Redington's strict positive-curvature
#' condition and does not guarantee a nonnegative surplus after rate shifts.
#' Liability values retain the original log accumulation rate; convexity is
#' still with respect to the annual effective yield. Relative moments and row
#' scaling before monetary products prevent avoidable overflow. Scaling uses
#' common monetary and moment magnitudes, preserving the minimum norm over
#' asset units, without rescaling asset columns independently. Time-zero
#' liabilities contribute value but zero duration and convexity. Final values,
#' weights and sensitivities must be finite; nearly cancelled liability values
#' are rejected. Rank and cancellation tests remain subject to double precision;
#' reported portfolio moments should be inspected for ill-conditioned systems.
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
#'   C = c(6, 30, 72) / (1.05)^2,
#'   i = 0.05
#' )
#'
#' # Four-asset immunization: minimum-norm solution
#' immunize_duration_convexity(
#'   L = c(5000, 8000, 10000),
#'   t = c(3, 5, 7),
#'   P = c(100, 120, 150, 200),
#'   D = c(2, 4, 6, 8),
#'   C = c(6, 20, 42, 72) / (1.05)^2,
#'   i = 0.05
#' )
#'
#' # Nominal annual interest rate convertible monthly
#' immunize_duration_convexity(
#'   L = c(5000, 8000, 10000),
#'   t = c(3, 5, 7),
#'   P = c(100, 150, 200),
#'   D = c(2, 5, 8),
#'   C = c(6, 30, 72) / (1 + 0.06 / 12)^(2 * 12),
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

  m <- round(m)
  growth <- .interest_basis(i_type, i, m)$growth
  if (length(growth) != 1L || !is.finite(growth)) {
    stop("The original log accumulation rate must be finite.", call. = FALSE)
  }

  # --- Liabilities: PV, duration, convexity ---
  discounted <- .cash_flow_value(L, -t * growth)
  pv_L <- sum(discounted)
  pv_scale <- sum(abs(discounted))

  if (!is.finite(pv_L) || !is.finite(pv_scale)) {
    stop("Discounting liabilities produced non-finite values.", call. = FALSE)
  }
  if (pv_scale == 0 || abs(pv_L) <= .Machine$double.eps * 100 * pv_scale) {
    stop("Present value of liabilities is zero; immunization is not meaningful.",
         call. = FALSE)
  }

  D_L <- sign(pv_L) * sum(.cash_flow_value(
    discounted, log(t) - log(abs(pv_L))))
  C_L <- sign(pv_L) * sum(.cash_flow_value(
    discounted, log(t) + log1p(t) - 2 * growth - log(abs(pv_L))))
  if (!is.finite(D_L) || !is.finite(C_L)) {
    stop("Liability duration or convexity is non-finite.", call. = FALSE)
  }

  # --- Build linear system A w = b (3 x n_assets) ---
  price_scale <- max(P, abs(pv_L))
  duration_scale <- max(abs(D), abs(D_L))
  convexity_scale <- max(abs(C), abs(C_L))
  if (duration_scale == 0) duration_scale <- 1
  if (convexity_scale == 0) convexity_scale <- 1
  scaled_prices <- P / price_scale
  scaled_durations <- D / duration_scale
  scaled_convexities <- C / convexity_scale
  A <- rbind(scaled_prices, scaled_prices * scaled_durations,
             scaled_prices * scaled_convexities)
  b <- (pv_L / price_scale) * c(1, D_L / duration_scale, C_L / convexity_scale)

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
    w <- .immunization_min_norm(A, b)
  }
  if (any(!is.finite(w))) {
    stop("Immunization weights are non-finite. Check input magnitudes.", call. = FALSE)
  }

  # --- Portfolio measures ---
  market_scaled <- w * scaled_prices
  portfolio_scaled <- sum(market_scaled)
  if (any(!is.finite(market_scaled)) || !is.finite(portfolio_scaled) || portfolio_scaled == 0) {
    stop("Resulting portfolio present value is zero or non-finite.", call. = FALSE)
  }
  pv_A <- .cash_flow_value(portfolio_scaled, log(price_scale))
  D_A <- sign(portfolio_scaled) * .cash_flow_value(
    sum(market_scaled * scaled_durations), log(duration_scale) - log(abs(portfolio_scaled)))
  C_A <- sign(portfolio_scaled) * .cash_flow_value(
    sum(market_scaled * scaled_convexities), log(convexity_scale) - log(abs(portfolio_scaled)))
  if (any(!is.finite(c(pv_A, D_A, C_A)))) {
    stop("Resulting portfolio measures are non-finite.", call. = FALSE)
  }

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
