#' Build a multiple decrement table
#'
#' Constructs a multiple decrement table from either independent
#' (absolute) rates or dependent (observed) rates, with automatic
#' conversion between the two under a UDD assumption
#' (Finan, Sections 65-67).
#'
#' @param x Integer vector of ages.
#' @param q_prime Optional numeric matrix of independent (absolute) rates
#'   \eqn{q'^{(j)}_x}, one row per age and one column per cause. Supply
#'   either \code{q_prime} or \code{q_dep} (not both).
#' @param q_dep Optional numeric matrix of dependent rates
#'   \eqn{q^{(j)}_x} in the same format. Supply either \code{q_dep} or
#'   \code{q_prime}.
#' @param radix Starting population \eqn{l^{(\tau)}_a} (default
#'   \code{100000}).
#' @param causes Character vector of cause names. If \code{NULL},
#'   taken from column names of \code{q_prime} or \code{q_dep}, or
#'   defaults to \code{"cause_1"}, \code{"cause_2"}, etc.
#' @param assumption UDD assumption for rate conversion:
#'   \code{"UDD_single"} (default, Finan Sec. 67 second approach) or
#'   \code{"UDD_multi"} (Finan Sec. 67 first approach). Only used
#'   when converting between rate types.
#'
#' @details
#' The function builds the complete multiple decrement table containing:
#' \itemize{
#'   \item \eqn{l^{(\tau)}_x}: survivors under all decrements
#'   \item \eqn{d^{(j)}_x}: expected decrements due to cause \eqn{j}
#'   \item \eqn{q^{(j)}_x}: dependent probability of decrement due to
#'     cause \eqn{j} (Finan, Sec. 66)
#'   \item \eqn{q'^{(j)}_x}: independent (absolute) probability of
#'     decrement (Finan, Sec. 65)
#'   \item \eqn{q^{(\tau)}_x}: total probability of decrement
#'   \item \eqn{p^{(\tau)}_x}: total survival probability
#' }
#'
#' The key relationships (Finan, Sec. 65):
#' \deqn{p^{(\tau)}_x = \prod_{j=1}^{m} p'^{(j)}_x = \prod_{j=1}^{m} (1 - q'^{(j)}_x)}
#' \deqn{q^{(\tau)}_x = \sum_{j=1}^{m} q^{(j)}_x}
#' \deqn{d^{(j)}_x = l^{(\tau)}_x \times q^{(j)}_x}
#'
#' @return A tibble with columns \code{x}, \code{lx_tau},
#'   \code{q_tau}, \code{p_tau}, and for each cause \eqn{j}:
#'   \code{q_dep_j}, \code{q_prime_j}, \code{dx_j}.
#'   The result has class \code{"multi_decrement_table"} in addition
#'   to \code{"tbl_df"}.
#'
#' @export
multi_decrement_table <- function(
    x,
    q_prime = NULL,
    q_dep = NULL,
    radix = 100000,
    causes = NULL,
    assumption = c("UDD_single", "UDD_multi")
) {
  assumption <- match.arg(assumption)

  # --- checks ---
  if (!is.numeric(x) || length(x) < 1) {
    stop("'x' must be a numeric vector of ages.")
  }
  x <- as.integer(round(x))
  n_ages <- length(x)

  if (is.null(q_prime) && is.null(q_dep)) {
    stop("Supply either 'q_prime' or 'q_dep'.")
  }
  if (!is.null(q_prime) && !is.null(q_dep)) {
    stop("Supply only one of 'q_prime' or 'q_dep'.")
  }

  if (!is.numeric(radix) || length(radix) != 1L ||
      radix <= 0) {
    stop("'radix' must be a single positive number.")
  }

  # --- standardize input to matrices ---
  if (!is.null(q_prime)) {
    q_prime <- as.matrix(q_prime)
    if (nrow(q_prime) != n_ages) {
      stop("Number of rows in 'q_prime' must match length of 'x'.")
    }
    n_causes <- ncol(q_prime)
  } else {
    q_dep <- as.matrix(q_dep)
    if (nrow(q_dep) != n_ages) {
      stop("Number of rows in 'q_dep' must match length of 'x'.")
    }
    n_causes <- ncol(q_dep)
  }

  # --- cause names ---
  if (is.null(causes)) {
    src <- if (!is.null(q_prime)) q_prime else q_dep
    if (!is.null(colnames(src))) {
      causes <- colnames(src)
    } else {
      causes <- paste0("cause_", seq_len(n_causes))
    }
  }
  if (length(causes) != n_causes) {
    stop("'causes' length must match number of columns in rates.")
  }

  # --- convert rates ---
  if (!is.null(q_prime)) {
    q_indep <- q_prime
    q_depend <- md_convert_rates(
      q_prime, direction = "indep_to_dep",
      assumption = assumption
    )
  } else {
    q_depend <- q_dep
    q_indep <- md_convert_rates(
      q_dep, direction = "dep_to_indep",
      assumption = assumption
    )
  }

  # --- build table ---
  # Total rates
  q_tau <- rowSums(q_depend)
  p_tau <- 1 - q_tau

  # Survivors
  lx_tau <- numeric(n_ages)
  lx_tau[1] <- radix
  if (n_ages > 1) {
    for (row in 2:n_ages) {
      lx_tau[row] <- lx_tau[row - 1] * p_tau[row - 1]
    }
  }

  # Decrements by cause
  dx_mat <- matrix(NA_real_, nrow = n_ages, ncol = n_causes)
  for (j in seq_len(n_causes)) {
    dx_mat[, j] <- lx_tau * q_depend[, j]
  }

  # --- assemble tibble ---
  out <- tibble::tibble(x = x, lx_tau = lx_tau,
                        q_tau = q_tau, p_tau = p_tau)

  for (j in seq_len(n_causes)) {
    nm_dep   <- paste0("q_dep_", causes[j])
    nm_prime <- paste0("q_prime_", causes[j])
    nm_dx    <- paste0("dx_", causes[j])
    out[[nm_dep]]   <- q_depend[, j]
    out[[nm_prime]] <- q_indep[, j]
    out[[nm_dx]]    <- dx_mat[, j]
  }

  class(out) <- c("multi_decrement_table", class(out))
  attr(out, "causes") <- causes
  attr(out, "radix")  <- radix

  out
}
