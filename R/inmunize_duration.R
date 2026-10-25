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
#' The minimum Euclidean norm is taken over asset units \code{w}, not
#' portfolio market-value fractions. It is evaluated by singular value
#' decomposition of the row-scaled constraint matrix, without forming
#' normal equations. Short positions are allowed. Matching present value
#' and duration gives a zero value and first derivative of the surplus at
#' the base yield, but does not ensure nonnegative surplus after rate shifts;
#' the Redington convexity condition must be checked separately.
#' Liability values retain the original log accumulation rate, and relative
#' moments avoid monetary numerator overflow. Payment times include zero:
#' a liability at time zero contributes value but no duration. Constraints are
#' scaled by a common monetary magnitude and by moment magnitudes before
#' products are formed. These are row transformations, preserving the norm over
#' units; asset columns are not rescaled independently. Final values, weights
#' and sensitivities must be finite, and nearly cancelled liability value is
#' rejected. Rank and cancellation tests remain subject to double precision;
#' reported portfolio moments should be inspected for ill-conditioned systems.
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

  m <- round(m)
  growth <- .interest_basis(i_type, i, m)$growth
  if (length(growth) != 1L || !is.finite(growth)) {
    stop("The original log accumulation rate must be finite.", call. = FALSE)
  }

  # --- Present value and duration of liabilities ---
  discounted <- .cash_flow_value(L, -t * growth)
  pv_L <- sum(discounted)
  pv_scale <- sum(abs(discounted))

  if (!is.finite(pv_L) || !is.finite(pv_scale)) {
    stop("Discounting liabilities produced non-finite values.", call. = FALSE)
  }
  if (pv_scale == 0 || abs(pv_L) <= .Machine$double.eps * 100 * pv_scale) {
    stop("Present value of liabilities is zero; immunization is not meaningful.", call. = FALSE)
  }

  dur_L <- sign(pv_L) * sum(.cash_flow_value(
    discounted, log(t) - log(abs(pv_L))))
  if (!is.finite(dur_L)) {
    stop("Liability duration is non-finite.", call. = FALSE)
  }

  # --- Build linear system A w = b ---
  # Row 1: PV constraint
  # Row 2: PV * duration constraint
  price_scale <- max(P, abs(pv_L))
  duration_scale <- max(abs(D), abs(dur_L))
  if (duration_scale == 0) duration_scale <- 1
  scaled_prices <- P / price_scale
  scaled_durations <- D / duration_scale
  A <- rbind(scaled_prices, scaled_prices * scaled_durations)
  b <- (pv_L / price_scale) * c(1, dur_L / duration_scale)

  # --- Solve for weights ---
  if (n_assets == 2L) {
    # Closed-form solution for two assets
    P1 <- P[1L]
    P2 <- P[2L]
    D1 <- scaled_durations[1L]
    D2 <- scaled_durations[2L]

    den <- D1 - D2
    if (abs(den) <= .Machine$double.eps * 100 * max(abs(D1), abs(D2))) {
      stop("Degenerate case: assets have equal duration; cannot solve uniquely.", call. = FALSE)
    }

    fractions <- c((dur_L / duration_scale - D2) / den,
                   (D1 - dur_L / duration_scale) / den)
    w <- .cash_flow_value(sign(fractions) * pv_L,
                          log(abs(fractions)) - log(c(P1, P2)))
  } else {
    w <- .immunization_min_norm(A, b)
  }
  if (any(!is.finite(w))) {
    stop("Immunization weights are non-finite. Check input magnitudes.", call. = FALSE)
  }

  # --- Check resulting portfolio ---
  market_scaled <- w * scaled_prices
  portfolio_scaled <- sum(market_scaled)
  if (any(!is.finite(market_scaled)) || !is.finite(portfolio_scaled) || portfolio_scaled == 0) {
    stop("Resulting portfolio present value is zero or non-finite.", call. = FALSE)
  }
  pv_A <- .cash_flow_value(portfolio_scaled, log(price_scale))
  dur_A <- sign(portfolio_scaled) * .cash_flow_value(
    sum(market_scaled * scaled_durations), log(duration_scale) - log(abs(portfolio_scaled)))
  if (!is.finite(pv_A) || !is.finite(dur_A)) {
    stop("Resulting portfolio measures are non-finite.", call. = FALSE)
  }

  tibble::tibble(
    w = w,
    PV_L = pv_L,
    D_L = dur_L,
    PV_A = pv_A,
    D_A = dur_A,
    n_assets = n_assets
  )
}
