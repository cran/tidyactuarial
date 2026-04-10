#' Actuarial present value of a two-life annuity
#'
#' Computes the APV of a discrete annuity contingent on two independent lives
#' aged \code{x} and \code{y}. Supports joint-life, last-survivor, and
#' reversionary-style benefits via state-based weights.
#'
#' @param lt A life table data frame with columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for life 1.
#' @param y Integer actuarial age for life 2.
#' @param i Annual effective interest rate (must be \code{> -1}).
#' @param cohort Status cohort: \code{"first"} (pays while both alive) or
#'   \code{"last"} (pays while at least one alive). Ignored if \code{benefit}
#'   is explicitly supplied.
#' @param benefit Optional list with weights \code{both}, \code{x_only},
#'   \code{y_only}. If supplied, the payment at time \eqn{t} is weighted by
#'   the probability of each state at \eqn{t}.
#' @param n Integer term (years). If \code{NULL}, whole life to end of table.
#' @param m Integer deferment (years). Default \code{0}.
#' @param k Integer payments per year. Default \code{1}.
#' @param timing Payment timing: \code{"immediate"} or \code{"due"}.
#' @param woolhouse Woolhouse approximation for \code{k > 1}:
#'   \code{"none"} (exact UDD k-thly), \code{"first"}, or \code{"second"}.
#' @param frac Fractional-age assumption for \code{k > 1} exact computation:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#'   If not specified and \code{lt} carries a \code{frac} attribute, that
#'   value is used.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#'
#' @details
#' This function assumes independence between lives (Finan, Sections 56--59).
#'
#' **Joint-life annuity-due** (Finan, Section 58):
#' \deqn{\ddot{a}_{xy} = \sum_{k=0}^{\infty} v^k \cdot {}_kp_{xy}.}
#'
#' **Last-survivor annuity-due** (Finan, Section 59):
#' \deqn{\ddot{a}_{\overline{xy}} = \ddot{a}_x + \ddot{a}_y -
#'   \ddot{a}_{xy}.}
#'
#' **State-based benefits**: for reversionary annuities, use
#' \code{benefit = list(both = 1, x_only = alpha, y_only = alpha)} where
#' \eqn{\alpha} is the fraction paid to the survivor.
#'
#' For \code{k > 1} with \code{woolhouse = "none"}, exact k-thly computation
#' is performed under the selected fractional-age assumption using
#' \code{\link{t_px}}.
#'
#' @return A single numeric APV value, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{annuity_x}} for single-life annuity APVs,
#'   \code{\link{insurance_xy}} for two-life insurance,
#'   \code{\link{t_pxy}} for two-life survival,
#'   \code{\link{premium_xy}} for two-life premiums.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Joint-life annuity-due (Finan, Sec. 58)
#' annuity_xy(lt, x = 60, y = 62, i = 0.05, cohort = "first", timing = "due")
#'
#' # Last-survivor annuity-due (Finan, Sec. 59)
#' annuity_xy(lt, x = 60, y = 62, i = 0.05, cohort = "last", timing = "due")
#'
#' # Verify identity: a_{xy-bar} = a_x + a_y - a_{xy}
#' a_joint <- annuity_xy(lt, x = 60, y = 62, i = 0.05,
#'                       cohort = "first", timing = "due")
#' a_last  <- annuity_xy(lt, x = 60, y = 62, i = 0.05,
#'                       cohort = "last", timing = "due")
#' a_x <- annuity_x(lt, x = 60, i = 0.05, timing = "due")
#' a_y <- annuity_x(lt, x = 62, i = 0.05, timing = "due")
#' c(last = a_last, sum_minus_joint = a_x + a_y - a_joint)
#'
#' # Reversionary: full while both alive, 60% to survivor
#' annuity_xy(lt, x = 60, y = 62, i = 0.05,
#'            benefit = list(both = 1, x_only = 0.6, y_only = 0.6),
#'            timing = "due")
#'
#' # 2-year deferred joint-life annuity
#' annuity_xy(lt, x = 60, y = 62, i = 0.05, m = 2,
#'            cohort = "first", timing = "due")
#'
#' # Tidy output
#' annuity_xy(lt, x = 60, y = 62, i = 0.05, cohort = "first",
#'            timing = "due", tidy = TRUE)
#'
#' @export
annuity_xy <- function(
    lt, x, y, i,
    cohort = c("first", "last"),
    benefit = NULL,
    n = NULL,
    m = 0L,
    k = 1L,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    frac,
    tidy = FALSE
) {
  cohort    <- match.arg(cohort)
  timing    <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)

  # --- frac inheritance ---
  if (missing(frac)) {
    lt_frac <- attr(lt, "frac")
    if (!is.null(lt_frac) && lt_frac %in% c("UDD", "CF", "Balducci")) {
      frac <- lt_frac
    } else {
      frac <- "UDD"
    }
  } else {
    frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac == "CML") frac <- "CF"
  }

  # --- checks ---
  if (missing(i) || !is.numeric(i) || length(i) != 1 || i <= -1) {
    stop("'i' must be a single numeric rate > -1.")
  }
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) {
    stop("Life table must contain columns 'x' and 'lx'.")
  }

  x <- as.integer(round(x)); y <- as.integer(round(y))
  m <- as.integer(round(m)); k <- as.integer(round(k))
  if (m < 0) stop("'m' must be nonnegative.")
  if (k < 1) stop("'k' must be >= 1.")

  omega <- max(lt$x, na.rm = TRUE)

  # --- default n: depends on status ---
  if (is.null(n)) {
    if (cohort == "first" && is.null(benefit)) {
      n <- max(0L, omega - max(x, y) - m)
    } else {
      # last survivor or custom benefit: younger life determines horizon
      n <- max(0L, omega - min(x, y) - m)
    }
  } else {
    n <- as.integer(round(n))
    if (n < 0) stop("'n' must be nonnegative or NULL.")
  }

  if (n == 0L) {
    if (isTRUE(tidy)) {
      return(tibble::tibble(x = x, y = y, i = i, n = 0L, m = m,
                            k = k, timing = timing, cohort = cohort, apv = 0))
    }
    return(0)
  }

  v_fun <- function(tt) (1 + i)^(-tt)

  # --- set default benefit weights ---
  if (is.null(benefit)) {
    if (cohort == "first") {
      benefit <- list(both = 1, x_only = 0, y_only = 0)
    } else {
      benefit <- list(both = 1, x_only = 1, y_only = 1)
    }
  } else {
    if (!is.list(benefit) || !all(c("both", "x_only", "y_only") %in% names(benefit))) {
      stop("`benefit` must be a list with names: both, x_only, y_only.")
    }
  }

  # --- expected payment at time tt from issue (using survival from age x, y) ---
  E_pay <- function(tt) {
    px_t <- t_px(lt, x = x, t = tt, frac = frac, check = FALSE)
    py_t <- t_px(lt, x = y, t = tt, frac = frac, check = FALSE)
    if (is.na(px_t) || is.na(py_t)) return(NA_real_)

    p_both   <- px_t * py_t
    p_x_only <- px_t * (1 - py_t)
    p_y_only <- py_t * (1 - px_t)

    benefit$both * p_both + benefit$x_only * p_x_only + benefit$y_only * p_y_only
  }

  # --- exact computation (annual or k-thly) ---
  exact_apv <- function(nn, kk, tim) {
    N <- nn * kk
    if (N == 0L) return(0)

    if (tim == "due") {
      j <- 0:(N - 1L)
    } else {
      j <- 1:N
    }
    u     <- j / kk                     # time offsets from start of annuity
    times <- m + u                       # absolute times from issue
    disc  <- v_fun(times)
    ep    <- vapply(times, E_pay, numeric(1))
    if (anyNA(ep)) return(NA_real_)

    sum((1 / kk) * disc * ep)
  }

  # --- main computation ---
  if (k == 1L || woolhouse == "none") {
    result <- exact_apv(n, k, timing)
  } else {
    # --- Woolhouse approximation on annual annuity-due ---
    adue <- exact_apv(n, 1L, "due")

    # nEx for the status: v^n * E_pay(m+n) / E_pay(m)
    # (ratio of weighted survival at end vs start of annuity)
    ep_start <- E_pay(m)
    ep_end   <- E_pay(m + n)
    if (is.na(ep_start) || ep_start <= 0) ep_start <- 1
    if (is.na(ep_end)) ep_end <- 0
    nEx_status <- v_fun(n) * ep_end / ep_start

    # 2-term Woolhouse (Finan, Problem 38.9a generalized):
    adj1 <- (k - 1) / (2 * k) * (1 - nEx_status)

    if (woolhouse == "first") {
      adue_k <- adue - adj1
    } else {
      # 3-term: need mu of the status at start and end
      delta <- log(1 + i)

      # mu at start: -d/dt ln(E_pay) at t = m
      ep_m1 <- E_pay(m + 1)
      if (is.na(ep_m1) || ep_start <= 0) {
        mu_start <- 0
      } else {
        p_status_1 <- ep_m1 / ep_start
        mu_start <- if (p_status_1 > 0) -log(p_status_1) else 0
      }

      # mu at end: at t = m + n
      ep_n1 <- E_pay(m + n + 1)
      if (is.na(ep_n1) || is.na(ep_end) || ep_end <= 0) {
        mu_end <- 0
      } else {
        p_status_n1 <- ep_n1 / ep_end
        mu_end <- if (p_status_n1 > 0) -log(p_status_n1) else 0
      }

      adj2 <- (k^2 - 1) / (12 * k^2) * (delta + mu_start - nEx_status * (delta + mu_end))
      adue_k <- adue - adj1 - adj2
    }

    # Convert due -> immediate (Finan, Problem 38.6 generalized):
    a_k <- if (timing == "due") adue_k else (adue_k - (1 / k) * (1 - nEx_status))

    result <- a_k
  }

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x, y = y, i = i, n = as.integer(n), m = m, k = k,
      timing = timing, cohort = cohort, woolhouse = woolhouse, apv = result
    ))
  }
  result
}
