#' Actuarial present value of a two-life annuity
#'
#' Computes the APV of a discrete annuity contingent on two independent lives
#' aged \code{x} and \code{y}. Supports joint-life, last-survivor, and
#' reversionary-style benefits via state-based weights.
#'
#' @param lt Either:
#'   \itemize{
#'     \item a single life table data frame, used for both lives; or
#'     \item a list of two life tables \code{list(lt_x, lt_y)}, one for each life.
#'   }
#'   Each life table must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for life 1.
#' @param y Integer actuarial age for life 2.
#' @param i Annual effective interest rate (must be \code{> -1}).
#' @param cohort Status cohort: \code{"first"} (pays while both alive) or
#'   \code{"last"} (pays while at least one alive). Ignored if \code{benefit}
#'   is explicitly supplied.
#' @param benefit Optional list with weights \code{both}, \code{x_only},
#'   \code{y_only}. If supplied, the payment at time \eqn{t} is weighted by
#'   the probability of each state at \eqn{t}.
#' @param n Integer term (years). If \code{NULL}, whole life to the relevant
#'   horizon implied by the supplied table(s) and benefit structure.
#' @param m Integer deferment (years). Default \code{0}.
#' @param k Integer payments per year. Default \code{1}.
#' @param timing Payment timing: \code{"immediate"} or \code{"due"}.
#' @param woolhouse Woolhouse approximation for \code{k > 1}:
#'   \code{"none"} (exact UDD/CF/Balducci k-thly), \code{"first"}, or
#'   \code{"second"}.
#' @param frac Fractional-age assumption for \code{k > 1} exact computation:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#'   If not specified and the supplied life table(s) carry a \code{frac}
#'   attribute, that value is used. If two tables are supplied and their
#'   \code{frac} attributes differ, \code{frac} must be supplied explicitly.
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
#' When \code{lt} is a list of two tables, the first table is used for life
#' \code{x} and the second for life \code{y}. Marginal survival probabilities
#' are always computed via \code{\link{t_px}} on the corresponding life table.
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
#' # Different life tables for the two lives
#' lt_m <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 98500, 96800, 94800, 92400, 89500, 86000)
#' )
#' lt_f <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97800, 96400, 94700, 92700, 90300)
#' )
#' annuity_xy(list(lt_m, lt_f), x = 60, y = 62, i = 0.05,
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

  # --- resolve life table input ---
  if (is.data.frame(lt)) {
    if (!all(c("x", "lx") %in% names(lt))) {
      stop("Life table must contain columns 'x' and 'lx'.")
    }
    lt_x <- lt
    lt_y <- lt
  } else if (is.list(lt) && length(lt) == 2L &&
             all(vapply(lt, is.data.frame, logical(1)))) {
    if (!all(c("x", "lx") %in% names(lt[[1]]))) {
      stop("First life table must contain columns 'x' and 'lx'.")
    }
    if (!all(c("x", "lx") %in% names(lt[[2]]))) {
      stop("Second life table must contain columns 'x' and 'lx'.")
    }
    lt_x <- lt[[1]]
    lt_y <- lt[[2]]
  } else {
    stop("`lt` must be either one life table or a list of two life tables.")
  }

  # --- frac inheritance ---
  if (missing(frac)) {
    frac_x <- attr(lt_x, "frac")
    frac_y <- attr(lt_y, "frac")

    ok_x <- !is.null(frac_x) && frac_x %in% c("UDD", "CF", "Balducci")
    ok_y <- !is.null(frac_y) && frac_y %in% c("UDD", "CF", "Balducci")

    if (ok_x && ok_y) {
      if (!identical(frac_x, frac_y)) {
        stop(
          "The two life tables carry different `frac` attributes. ",
          "Supply `frac` explicitly."
        )
      }
      frac <- frac_x
    } else if (ok_x) {
      frac <- frac_x
    } else if (ok_y) {
      frac <- frac_y
    } else {
      frac <- "UDD"
    }
  } else {
    frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac == "CML") frac <- "CF"
  }

  # --- checks ---
  if (missing(i) || !is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) {
    stop("'i' must be a single numeric rate > -1.")
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      abs(x - round(x)) > 1e-10) {
    stop("'x' must be a single integer age.")
  }
  if (!is.numeric(y) || length(y) != 1L || is.na(y) ||
      abs(y - round(y)) > 1e-10) {
    stop("'y' must be a single integer age.")
  }
  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      abs(m - round(m)) > 1e-10 || m < 0) {
    stop("'m' must be a single nonnegative integer.")
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) ||
      abs(k - round(k)) > 1e-10 || k < 1) {
    stop("'k' must be a single integer >= 1.")
  }
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("'tidy' must be TRUE or FALSE.")
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))
  k <- as.integer(round(k))

  omega_x <- max(lt_x$x, na.rm = TRUE)
  omega_y <- max(lt_y$x, na.rm = TRUE)

  # remaining whole-life horizons from issue
  hx <- max(0L, omega_x - x)
  hy <- max(0L, omega_y - y)

  # --- default n: depends on benefit structure ---
  if (is.null(benefit)) {
    if (cohort == "first") {
      horizon <- min(hx, hy)
    } else {
      horizon <- max(hx, hy)
    }
  } else {
    if (!is.list(benefit) || !all(c("both", "x_only", "y_only") %in% names(benefit))) {
      stop("`benefit` must be a list with names: both, x_only, y_only.")
    }

    h_both   <- if (!is.na(benefit$both)   && benefit$both   != 0) min(hx, hy) else 0L
    h_x_only <- if (!is.na(benefit$x_only) && benefit$x_only != 0) hx else 0L
    h_y_only <- if (!is.na(benefit$y_only) && benefit$y_only != 0) hy else 0L

    horizon <- max(h_both, h_x_only, h_y_only)
  }

  if (is.null(n)) {
    n <- max(0L, horizon - m)
  } else {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
        abs(n - round(n)) > 1e-10 || n < 0) {
      stop("'n' must be a single nonnegative integer or NULL.")
    }
    n <- as.integer(round(n))
  }

  if (n == 0L) {
    if (isTRUE(tidy)) {
      return(tibble::tibble(
        x = x, y = y, i = i, n = 0L, m = m,
        k = k, timing = timing, cohort = cohort,
        woolhouse = woolhouse, apv = 0
      ))
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
  }

  # --- expected payment at time tt from issue ---
  E_pay <- function(tt) {
    px_t <- t_px(lt_x, x = x, t = tt, frac = frac, check = FALSE)
    py_t <- t_px(lt_y, x = y, t = tt, frac = frac, check = FALSE)

    if (is.na(px_t) || is.na(py_t)) return(NA_real_)

    p_both   <- px_t * py_t
    p_x_only <- px_t * (1 - py_t)
    p_y_only <- py_t * (1 - px_t)

    benefit$both * p_both +
      benefit$x_only * p_x_only +
      benefit$y_only * p_y_only
  }

  # --- exact computation (annual or k-thly) ---
  exact_apv <- function(nn, kk, tim) {
    N <- nn * kk
    if (N == 0L) return(0)

    j <- if (tim == "due") 0:(N - 1L) else 1:N

    u     <- j / kk
    times <- m + u
    disc  <- v_fun(times)
    ep    <- vapply(times, E_pay, numeric(1))

    if (anyNA(ep)) return(NA_real_)

    sum((1 / kk) * disc * ep)
  }

  # --- main computation ---
  if (k == 1L || woolhouse == "none") {
    result <- exact_apv(n, k, timing)
  } else {
    # Woolhouse built from annual annuity-due
    adue <- exact_apv(n, 1L, "due")

    # nEx for the status/benefit stream
    ep_start <- E_pay(m)
    ep_end   <- E_pay(m + n)

    if (is.na(ep_start) || ep_start <= 0) ep_start <- 1
    if (is.na(ep_end)) ep_end <- 0

    nEx_status <- v_fun(n) * ep_end / ep_start

    # 2-term Woolhouse adjustment
    adj1 <- (k - 1) / (2 * k) * (1 - nEx_status)

    if (woolhouse == "first") {
      adue_k <- adue - adj1
    } else {
      # 3-term Woolhouse adjustment
      delta <- log(1 + i)

      ep_m1 <- E_pay(m + 1)
      if (is.na(ep_m1) || ep_start <= 0) {
        mu_start <- 0
      } else {
        p_status_1 <- ep_m1 / ep_start
        mu_start <- if (p_status_1 > 0) -log(p_status_1) else 0
      }

      ep_n1 <- E_pay(m + n + 1)
      if (is.na(ep_n1) || is.na(ep_end) || ep_end <= 0) {
        mu_end <- 0
      } else {
        p_status_n1 <- ep_n1 / ep_end
        mu_end <- if (p_status_n1 > 0) -log(p_status_n1) else 0
      }

      adj2 <- (k^2 - 1) / (12 * k^2) *
        (delta + mu_start - nEx_status * (delta + mu_end))

      adue_k <- adue - adj1 - adj2
    }

    # due -> immediate
    result <- if (timing == "due") {
      adue_k
    } else {
      adue_k - (1 / k) * (1 - nEx_status)
    }
  }

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x, y = y, i = i, n = as.integer(n), m = m, k = k,
      timing = timing, cohort = cohort, woolhouse = woolhouse, apv = result
    ))
  }

  result
}
