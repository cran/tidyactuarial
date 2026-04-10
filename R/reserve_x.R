#' Benefit reserve schedule for single-life insurance
#'
#' Computes the terminal benefit reserve \eqn{{}_kV} at one or more
#' policy durations for a fully discrete single-life insurance contract,
#' using the prospective or recursive method (Finan, Sections 47 and 52).
#'
#' @param lt A life table data frame with columns \code{x} and \code{lx}.
#' @param x Integer actuarial age at issue.
#' @param i Annual effective interest rate (must be \code{> -1}).
#' @param type Insurance type: \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param n Integer insurance term in years. Required for \code{"term"}
#'   and \code{"endowment"}. For \code{"whole"}, determined from the
#'   life table.
#' @param benefit Numeric benefit amount (default \code{1}).
#' @param premium Net premium per payment. If \code{NULL} (default),
#'   computed internally via the equivalence principle using
#'   \code{\link{premium_x}}.
#' @param h Integer premium-paying term in years. If \code{NULL},
#'   premiums are payable for the full duration of the contract
#'   (i.e., \code{h = n} for temporary products, whole life for whole).
#'   Set \code{h < n} for limited-payment policies (Finan, Sec. 47.3).
#' @param at Integer vector of policy durations at which to compute
#'   the reserve. Default \code{NULL} computes for all integer
#'   durations \code{0, 1, ..., n} (or to end of table for whole life).
#' @param method Computation method: \code{"prospective"} (default,
#'   Finan Sec. 47) or \code{"recursive"} (Finan Sec. 52).
#' @param tidy Logical. If \code{TRUE} (default), returns a tibble
#'   schedule. If \code{FALSE}, returns a named numeric vector.
#'
#' @details
#' **Prospective method** (Finan, Sections 47.1--47.3):
#' \deqn{{}_kV = \text{APV(future benefits at } x+k\text{)} -
#'   P \cdot \text{APV(future premiums at } x+k\text{)}}
#'
#' For each product type:
#' \itemize{
#'   \item \strong{Whole life} (Sec. 47.1):
#'     \eqn{{}_kV(A_x) = A_{x+k} - P \, \ddot{a}_{x+k}}
#'   \item \strong{Term} (Sec. 47.2), \eqn{k < n}:
#'     \eqn{{}_kV(A^1_{x:\overline{n}|}) = A^1_{x+k:\overline{n-k}|}
#'     - P \, \ddot{a}_{x+k:\overline{n_p-k}|}}
#'     where \eqn{n_p = \min(n, h)}.
#'   \item \strong{Endowment} (Sec. 47.3), \eqn{k < n}:
#'     \eqn{{}_kV(A_{x:\overline{n}|}) = A_{x+k:\overline{n-k}|}
#'     - P \, \ddot{a}_{x+k:\overline{n_p-k}|}}
#'   \item Endowment at \eqn{k = n}: \eqn{{}_nV = \text{benefit}}.
#' }
#'
#' For limited-payment policies (\code{h < n}), premiums cease after
#' year \eqn{h}. For \eqn{k \ge h}, the premium annuity term is zero.
#'
#' **Recursive method** (Finan, Section 52):
#' \deqn{{}_{k+1}V = \frac{({}_kV + \pi_k)(1+i) - b_{k+1} \cdot
#'   q_{x+k}}{p_{x+k}}}
#' starting from \eqn{{}_0V = 0} (equivalence principle).
#'
#' @return If \code{tidy = TRUE}, a tibble with columns \code{k}
#'   (duration), \code{age} (\eqn{x+k}), \code{reserve},
#'   \code{premium_paid} (premium at start of year \eqn{k+1}, 0
#'   after limited payment), and \code{benefit_due} (death benefit
#'   in year \eqn{k+1}). If \code{tidy = FALSE}, a named numeric
#'   vector of reserves.
#'
#' @seealso \code{\link{premium_x}} for benefit premiums,
#'   \code{\link{insurance_x}} for insurance APVs,
#'   \code{\link{annuity_x}} for annuity APVs,
#'   \code{\link{t_px}} for survival probabilities.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:70,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000,
#'          86000, 81000, 75000, 68000, 60000)
#' )
#'
#' # Whole life reserve schedule (Finan, Sec. 47.1)
#' reserve_x(lt, x = 60, i = 0.06, type = "whole")
#'
#' # 5-year term (Finan, Sec. 47.2)
#' reserve_x(lt, x = 60, i = 0.06, type = "term", n = 5)
#'
#' # 5-year endowment (Finan, Sec. 47.3)
#' reserve_x(lt, x = 60, i = 0.06, type = "endowment", n = 5)
#'
#' # Limited payment: 3-payment, 5-year endowment
#' reserve_x(lt, x = 60, i = 0.06, type = "endowment",
#'           n = 5, h = 3)
#'
#' # Recursive method (Finan, Sec. 52)
#' reserve_x(lt, x = 60, i = 0.06, type = "endowment",
#'           n = 5, method = "recursive")
#'
#' # Verify: prospective = recursive
#' r_pro <- reserve_x(lt, x = 60, i = 0.06, type = "endowment",
#'                    n = 5, method = "prospective")
#' r_rec <- reserve_x(lt, x = 60, i = 0.06, type = "endowment",
#'                    n = 5, method = "recursive")
#' all.equal(r_pro$reserve, r_rec$reserve)
#'
#' # Specific durations only
#' reserve_x(lt, x = 60, i = 0.06, type = "endowment",
#'           n = 5, at = c(0, 3, 5))
#'
#' # Benefit of $100,000
#' reserve_x(lt, x = 60, i = 0.06, type = "whole",
#'           benefit = 100000)
#'
#' # Custom premium (e.g., loaded)
#' reserve_x(lt, x = 60, i = 0.06, type = "endowment",
#'           n = 5, benefit = 100000, premium = 22000)
#'
#' # As named vector
#' reserve_x(lt, x = 60, i = 0.06, type = "endowment",
#'           n = 5, tidy = FALSE)
#'
#' @export
reserve_x <- function(
    lt, x, i,
    type = c("whole", "term", "endowment"),
    n = NULL,
    benefit = 1,
    premium = NULL,
    h = NULL,
    at = NULL,
    method = c("prospective", "recursive"),
    tidy = TRUE
) {
  type   <- match.arg(type)
  method <- match.arg(method)

  # --- checks ---
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) {
    stop("Life table must contain columns 'x' and 'lx'.")
  }
  if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) {
    stop("'i' must be a single numeric value > -1.")
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      abs(x - round(x)) > 1e-10) {
    stop("'x' must be a single integer age.")
  }
  if (!is.numeric(benefit) || length(benefit) != 1L ||
      is.na(benefit) || benefit <= 0) {
    stop("'benefit' must be a single positive number.")
  }

  x <- as.integer(round(x))
  omega <- max(lt$x, na.rm = TRUE)

  # --- determine n ---
  if (type == "whole") {
    max_n <- omega - x
    if (is.null(n)) n <- max_n
    n <- as.integer(round(n))
  } else {
    if (is.null(n)) {
      stop("'n' must be provided for term or endowment.")
    }
    n <- as.integer(round(n))
  }

  # --- determine h (premium payment term) ---
  if (is.null(h)) {
    h <- n
  } else {
    h <- as.integer(round(h))
    if (h < 1 || h > n) {
      stop("'h' must be between 1 and n.")
    }
  }

  # --- compute net premium if not supplied ---
  if (is.null(premium)) {
    premium <- premium_x(
      lt = lt, x = x, i = i,
      product = type, benefit = benefit,
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

  v_fun <- function(tt) (1 + i)^(-tt)
  d <- i / (1 + i)

  # -------------------------------------------------------
  # PROSPECTIVE METHOD (Finan, Sec. 47)
  # -------------------------------------------------------
  if (method == "prospective") {
    reserves <- vapply(k_vec, function(k) {
      # k = 0: reserve is 0 (equivalence principle)
      if (k == 0) return(0)

      # k = n for endowment: maturity value
      if (type == "endowment" && k == n) return(benefit)

      # k = n for term: 0 (coverage expired)
      if (type == "term" && k >= n) return(0)

      # k >= n for whole: not applicable (dead)
      if (k >= n) return(0)

      # APV of future benefits at age x+k
      apv_ben <- if (type == "whole") {
        benefit * insurance_x(lt, x = x + k, i = i, type = "whole")
      } else {
        benefit * insurance_x(
          lt, x = x + k, i = i,
          n = n - k, type = type
        )
      }

      # APV of future premiums at age x+k
      prem_term <- max(0L, h - k)
      if (prem_term == 0L) {
        apv_prem <- 0
      } else {
        apv_prem <- premium * annuity_x(
          lt, x = x + k, i = i,
          n = prem_term, timing = "due"
        )
      }

      apv_ben - apv_prem
    }, numeric(1))

  # -------------------------------------------------------
  # RECURSIVE METHOD (Finan, Sec. 52)
  # -------------------------------------------------------
  } else {
    # Compute full schedule 0..n, then subset
    V_full <- numeric(n + 1L)
    V_full[1] <- 0  # 0V = 0

    for (kk in 0:(n - 1L)) {
      # premium at start of year kk+1
      pi_k <- if (kk < h) premium else 0

      # benefit payable at end of year kk+1
      b_k1 <- benefit

      # mortality at age x + kk
      qxk <- 1 - t_px(lt, x = x + kk, t = 1,
                       frac = "UDD", check = FALSE)
      pxk <- 1 - qxk

      if (pxk <= 0) {
        # everyone dead: remaining reserves are 0
        V_full[(kk + 2):(n + 1)] <- 0
        break
      }

      # Forward recursion (Finan, Sec. 52)
      V_next <- ((V_full[kk + 1] + pi_k) * (1 + i) -
                   b_k1 * qxk) / pxk
      V_full[kk + 2] <- V_next
    }

    # For endowment at k = n: reserve = benefit
    if (type == "endowment") V_full[n + 1] <- benefit

    # For term at k = n: reserve = 0
    if (type == "term") V_full[n + 1] <- 0

    # Subset to requested durations
    reserves <- V_full[k_vec + 1L]
  }

  # --- output ---
  if (!isTRUE(tidy)) {
    names(reserves) <- paste0("k=", k_vec)
    return(reserves)
  }

  # Build schedule
  prem_paid <- vapply(k_vec, function(k) {
    if (k >= n) return(0)
    if (k < h) premium else 0
  }, numeric(1))

  ben_due <- vapply(k_vec, function(k) {
    if (k >= n) return(0)
    benefit
  }, numeric(1))

  tibble::tibble(
    k          = k_vec,
    age        = x + k_vec,
    reserve    = reserves,
    premium_paid = prem_paid,
    benefit_due  = ben_due
  )
}
