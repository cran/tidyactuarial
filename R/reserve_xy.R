#' Benefit reserve schedule for two-life insurance
#'
#' Computes the terminal benefit reserve \eqn{{}_kV} at one or more
#' policy durations for a fully discrete two-life insurance contract
#' (independent lives), using the prospective or recursive method.
#' Generalizes Finan, Sections 47 and 52 to the joint-life (Sec. 58)
#' and last-survivor (Sec. 59) statuses.
#'
#' @param lt Either a single life table data frame (used for both lives) or
#'   a list of two life tables \code{list(lt_x, lt_y)}, one for each life.
#'   Each life table must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for life 1 at issue.
#' @param y Integer actuarial age for life 2 at issue.
#' @param i Annual effective interest rate (must be \code{> -1}).
#' @param type Insurance type: \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param cohort Status cohort: \code{"first"} (joint-life, first death)
#'   or \code{"last"} (last-survivor, second death).
#' @param n Integer insurance term in years. Required for \code{"term"}
#'   and \code{"endowment"}. For \code{"whole"}, determined from the
#'   life table.
#' @param benefit Numeric benefit amount (default \code{1}).
#' @param premium Net premium per payment. If \code{NULL} (default),
#'   computed internally via the equivalence principle using
#'   \code{\link{premium_xy}}.
#' @param h Integer premium-paying term in years. If \code{NULL},
#'   premiums are payable for the full duration of the contract.
#'   Set \code{h < n} for limited-payment policies.
#' @param at Integer vector of policy durations at which to compute
#'   the reserve. Default \code{NULL} computes for all integer
#'   durations \code{0, 1, ..., n}.
#' @param method Computation method: \code{"prospective"} (default)
#'   or \code{"recursive"}.
#' @param tidy Logical. If \code{TRUE} (default), returns a tibble
#'   schedule. If \code{FALSE}, returns a named numeric vector.
#'
#' @details
#' **Prospective method** (generalizing Finan, Sec. 47 to two lives):
#' \deqn{{}_kV = \text{APV(future benefits at ages } x+k, y+k\text{)}
#'   - P \cdot \text{APV(future premiums at ages } x+k, y+k\text{)}}
#'
#' The APV of future benefits uses \code{\link{insurance_xy}} at
#' shifted ages \eqn{(x+k, y+k)}, and the APV of future premiums
#' uses \code{\link{annuity_xy}} on the same status.
#'
#' For endowment insurance at \eqn{k = n}: \eqn{{}_nV = \text{benefit}}.
#'
#' **Recursive method** (generalizing Finan, Sec. 52):
#' \deqn{{}_{k+1}V = \frac{({}_kV + \pi_k)(1+i) - b_{k+1} \cdot
#'   q_{\text{status}}(x+k, y+k)}{p_{\text{status}}(x+k, y+k)}}
#'
#' where \eqn{q_{\text{status}}} and \eqn{p_{\text{status}}} are
#' the one-year death and survival probabilities of the two-life
#' status, computed via \code{\link{t_pxy}}.
#'
#' @return If \code{tidy = TRUE}, a tibble with columns \code{k},
#'   \code{age_x}, \code{age_y}, \code{reserve}, \code{premium_paid},
#'   \code{benefit_due}. If \code{tidy = FALSE}, a named numeric
#'   vector of reserves.
#'
#' @seealso \code{\link{reserve_x}} for single-life reserves,
#'   \code{\link{premium_xy}} for two-life premiums,
#'   \code{\link{insurance_xy}} for two-life insurance APVs,
#'   \code{\link{annuity_xy}} for two-life annuity APVs,
#'   \code{\link{t_pxy}} for two-life survival.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:70,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000,
#'          86000, 81000, 75000, 68000, 60000)
#' )
#'
#' # Joint-life whole life reserve
#' reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'            type = "whole", cohort = "first")
#'
#' # Last-survivor 5-year endowment
#' reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'            type = "endowment", cohort = "last", n = 5)
#'
#' # Joint-life 4-year term
#' reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'            type = "term", cohort = "first", n = 4)
#'
#' # Limited payment: 3-payment, 5-year endowment, joint
#' reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'            type = "endowment", cohort = "first",
#'            n = 5, h = 3)
#'
#' # Verify: prospective = recursive
#' r_pro <- reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'                     type = "endowment", cohort = "first",
#'                     n = 5, method = "prospective")
#' r_rec <- reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'                     type = "endowment", cohort = "first",
#'                     n = 5, method = "recursive")
#' all.equal(r_pro$reserve, r_rec$reserve)
#'
#' # Benefit of $100,000
#' reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'            type = "whole", cohort = "first",
#'            benefit = 100000)
#'
#' # Specific durations
#' reserve_xy(lt, x = 60, y = 62, i = 0.06,
#'            type = "endowment", cohort = "first",
#'            n = 5, at = c(0, 3, 5))
#'
#' @export
reserve_xy <- function(
    lt, x, y, i,
    type = c("whole", "term", "endowment"),
    cohort = c("first", "last"),
    n = NULL,
    benefit = 1,
    premium = NULL,
    h = NULL,
    at = NULL,
    method = c("prospective", "recursive"),
    tidy = TRUE
) {
  type   <- match.arg(type)
  cohort <- match.arg(cohort)
  method <- match.arg(method)

  # --- resolve life table input ---
  if (is.data.frame(lt)) {
    if (!all(c("x", "lx") %in% names(lt))) {
      stop("Life table must contain columns 'x' and 'lx'.")
    }
    lt_use <- lt
    lt_x   <- lt
    lt_y   <- lt
  } else if (is.list(lt) && length(lt) == 2L &&
             all(vapply(lt, is.data.frame, logical(1)))) {
    if (!all(c("x", "lx") %in% names(lt[[1]]))) {
      stop("First life table must contain columns 'x' and 'lx'.")
    }
    if (!all(c("x", "lx") %in% names(lt[[2]]))) {
      stop("Second life table must contain columns 'x' and 'lx'.")
    }
    lt_use <- lt
    lt_x   <- lt[[1]]
    lt_y   <- lt[[2]]
  } else {
    stop("`lt` must be either one life table or a list of two life tables.")
  }

  # --- checks ---
  if (!is.numeric(i) || length(i) != 1L ||
      is.na(i) || i <= -1) {
    stop("'i' must be a single numeric value > -1.")
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      abs(x - round(x)) > 1e-10) {
    stop("'x' must be a single integer age.")
  }
  if (!is.numeric(y) || length(y) != 1L || is.na(y) ||
      abs(y - round(y)) > 1e-10) {
    stop("'y' must be a single integer age.")
  }
  if (!is.numeric(benefit) || length(benefit) != 1L ||
      is.na(benefit) || benefit <= 0) {
    stop("'benefit' must be a single positive number.")
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  omega_x <- max(lt_x$x, na.rm = TRUE)
  omega_y <- max(lt_y$x, na.rm = TRUE)

  # Remaining whole-life horizons from issue
  hx <- max(0L, omega_x - x)
  hy <- max(0L, omega_y - y)

  # Maximum feasible term implied by the tables and the status
  max_n <- if (cohort == "first") min(hx, hy) else max(hx, hy)

  # --- determine n ---
  if (type == "whole") {
    if (is.null(n)) {
      n <- max_n
    } else {
      if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
          n < 0 || abs(n - round(n)) > 1e-10) {
        stop("'n' must be a single nonnegative integer.")
      }
      n <- as.integer(round(n))
      if (n > max_n) {
        stop("'n' exceeds the maximum term implied by the life table(s) for this status.")
      }
    }
    n <- as.integer(round(n))
  } else {
    if (is.null(n)) {
      stop("'n' must be provided for term or endowment.")
    }
    if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
        n < 0 || abs(n - round(n)) > 1e-10) {
      stop("'n' must be a single nonnegative integer.")
    }
    n <- as.integer(round(n))
    if (n > max_n) {
      stop("'n' exceeds the maximum term implied by the life table(s) for this status.")
    }
  }

  # --- determine h (premium payment term) ---
  if (is.null(h)) {
    h <- n
  } else {
    h <- as.integer(round(h))
    if (h < 1 || h > n) stop("'h' must be between 1 and n.")
  }

  # --- compute net premium if not supplied ---
  if (is.null(premium)) {
    premium <- premium_xy(
      lt = lt_use, x = x, y = y, i = i,
      type = type, cohort = cohort,
      benefit = benefit,
      n = if (type == "whole") NULL else n,
      n_prem = h
    )
  }

  # --- determine durations ---
  if (is.null(at)) {
    k_vec <- 0:n
  } else {
    k_vec <- sort(unique(as.integer(round(at))))
    if (any(k_vec < 0) || any(k_vec > n)) {
      stop("'at' durations must be between 0 and n.")
    }
  }

  status <- if (cohort == "first") "joint" else "last"
  d <- i / (1 + i)

  # -------------------------------------------------------
  # PROSPECTIVE METHOD
  # -------------------------------------------------------
  if (method == "prospective") {
    reserves <- vapply(k_vec, function(k) {
      if (k == 0) return(0)
      if (type == "endowment" && k == n) return(benefit)
      if (type == "term" && k >= n) return(0)
      if (k >= n) return(0)

      # APV of future benefits at shifted ages
      apv_ben <- benefit * insurance_xy(
        lt = lt, x = x + k, y = y + k, i = i,
        n = if (type == "whole") NULL else (n - k),
        type = type, cohort = cohort
      )

      # APV of future premiums at shifted ages
      prem_term <- max(0L, h - k)
      if (prem_term == 0L) {
        apv_prem <- 0
      } else {
        apv_prem <- premium * annuity_xy(
          lt = lt, x = x + k, y = y + k, i = i,
          cohort = cohort, n = prem_term,
          timing = "due"
        )
      }

      apv_ben - apv_prem
    }, numeric(1))

    # -------------------------------------------------------
    # RECURSIVE METHOD
    # -------------------------------------------------------
  } else {
    V_full <- numeric(n + 1L)
    V_full[1] <- 0

    for (kk in 0:(n - 1L)) {
      pi_k <- if (kk < h) premium else 0
      b_k1 <- benefit

      # One-year status survival at ages x+kk, y+kk
      p_status <- t_pxy(
        lt = lt, x = x + kk, y = y + kk,
        t = 1, frac = "UDD", status = status
      )
      if (is.na(p_status)) {
        V_full[(kk + 2):(n + 1)] <- 0
        break
      }
      q_status <- 1 - p_status

      if (p_status <= 0) {
        V_full[(kk + 2):(n + 1)] <- 0
        break
      }

      # Forward recursion
      V_next <- ((V_full[kk + 1] + pi_k) * (1 + i) -
                   b_k1 * q_status) / p_status
      V_full[kk + 2] <- V_next
    }

    if (type == "endowment") V_full[n + 1] <- benefit
    if (type == "term") V_full[n + 1] <- 0

    reserves <- V_full[k_vec + 1L]
  }

  # --- output ---
  if (!isTRUE(tidy)) {
    names(reserves) <- paste0("k=", k_vec)
    return(reserves)
  }

  prem_paid <- vapply(k_vec, function(k) {
    if (k >= n) return(0)
    if (k < h) premium else 0
  }, numeric(1))

  ben_due <- vapply(k_vec, function(k) {
    if (k >= n) return(0)
    benefit
  }, numeric(1))

  tibble::tibble(
    k            = k_vec,
    age_x        = x + k_vec,
    age_y        = y + k_vec,
    reserve      = reserves,
    premium_paid = prem_paid,
    benefit_due  = ben_due
  )
}
