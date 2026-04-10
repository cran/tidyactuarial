#' Actuarial present value of a life insurance (derived from annuity_x)
#'
#' Computes the actuarial present value (APV) of a discrete life insurance with
#' benefit 1 payable at the end of the year of death. The implementation uses
#' the standard identities that express insurance values as functions of the
#' annuity-due value (Finan, Sections 37 and 41).
#'
#' Supported contracts:
#' \itemize{
#'   \item \strong{Whole life insurance} (Finan, Sec. 37.1):
#'     \eqn{A_x = 1 - d \, \ddot{a}_x}
#'   \item \strong{n-year term insurance} (Finan, Sec. 27):
#'     \eqn{A^1_{x:\overline{n}|} = 1 - d \, \ddot{a}_{x:\overline{n}|}
#'     - v^n \, {}_np_x}
#'   \item \strong{n-year endowment insurance} (Finan, Sec. 26.3.2):
#'     \eqn{A_{x:\overline{n}|} = 1 - d \, \ddot{a}_{x:\overline{n}|}}
#' }
#'
#' The endowment decomposes as (Finan, Example 26.15):
#' \deqn{A_{x:\overline{n}|} = A^1_{x:\overline{n}|} + {}_nE_x}
#'
#' Deferral by \eqn{m} years is handled as:
#' \eqn{v^m \, {}_mp_x \times \text{value at age } x+m}.
#'
#' @param lt A life table data frame. Must contain columns \code{x} and
#'   \code{lx}.
#' @param x Integer actuarial age.
#' @param i Effective annual interest rate (must be \code{> -1}).
#' @param n Integer term in years. Required for \code{type = "term"} or
#'   \code{type = "endowment"}. Ignored for \code{type = "whole"}.
#' @param m Integer deferral in years (default \code{0}).
#' @param type Contract type: \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#'
#' @details
#' This function calls \code{\link{annuity_x}} internally using annual
#' payments (\code{k = 1}) and \code{timing = "due"} (annuity-due).
#' Fractional payments and Woolhouse approximations are not used here
#' because the identities above are stated for annual discrete contracts.
#'
#' The key identity connecting annuities and insurance is
#' (Finan, Sec. 37.1):
#' \deqn{\ddot{a}_x = \frac{1 - A_x}{d}, \quad \text{i.e., } A_x = 1 - d\,\ddot{a}_x}
#' where \eqn{d = i/(1+i)} is the annual effective discount rate.
#'
#' For m-thly and continuous insurance APVs under UDD, use the
#' adjustment factor \eqn{i/i^{(m)}} or \eqn{i/\delta} (Finan, Sec. 30).
#'
#' @return A single numeric APV value, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{annuity_x}} for the annuity-due values used
#'   internally, \code{\link{premium_x}} for benefit premiums
#'   (\eqn{P = A / \ddot{a}}), \code{\link{Var_insurance_x}} for the
#'   variance, \code{\link{t_Ex}} for pure endowments,
#'   \code{\link{insurance_xy}} for two-life insurance.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:65,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000)
#' )
#'
#' # Whole life insurance: A_60 = 1 - d * \ddot{a}_60
#' insurance_x(lt, x = 60, i = 0.06, type = "whole")
#'
#' # 5-year term insurance: A^1_{60:5} (Finan, Sec. 27)
#' insurance_x(lt, x = 60, i = 0.06, n = 5, type = "term")
#'
#' # 5-year endowment insurance: A_{60:5} = A^1 + 5_E_60
#' insurance_x(lt, x = 60, i = 0.06, n = 5, type = "endowment")
#'
#' # Verify endowment decomposition (Finan, Sec. 26.3.2):
#' # A_{x:n} = A^1_{x:n} + nEx
#' A_term  <- insurance_x(lt, x = 60, i = 0.06, n = 5, type = "term")
#' A_endow <- insurance_x(lt, x = 60, i = 0.06, n = 5, type = "endowment")
#' nEx     <- (1.06)^(-5) * lt$lx[lt$x == 65] / lt$lx[lt$x == 60]
#' c(endowment = A_endow, term_plus_nEx = A_term + nEx)  # should match
#'
#' # 2-year deferred whole life
#' insurance_x(lt, x = 60, i = 0.06, m = 2, type = "whole")
#'
#' # Tidy output
#' insurance_x(lt, x = 60, i = 0.06, n = 5, type = "term", tidy = TRUE)
#'
#' @export
insurance_x <- function(
    lt, x, i,
    n = NULL,
    m = 0L,
    type = c("whole", "term", "endowment"),
    tidy = FALSE
) {
  type <- match.arg(type)

  # --- checks ---
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) stop("Life table must contain columns 'x' and 'lx'.")
  if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) stop("'i' must be > -1.")

  if (!is.numeric(x) || length(x) != 1L || is.na(x) || abs(x - round(x)) > 1e-10) {
    stop("'x' must be a single integer age.")
  }
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 0 || abs(m - round(m)) > 1e-10) {
    stop("'m' must be a single nonnegative integer.")
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))

  if (type %in% c("term", "endowment")) {
    if (is.null(n)) stop("'n' must be provided for type 'term' or 'endowment'.")
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 || abs(n - round(n)) > 1e-10) {
      stop("'n' must be a single nonnegative integer.")
    }
    n <- as.integer(round(n))
  }

  lt <- lt[order(lt$x), ]
  if (anyDuplicated(lt$x)) stop("Life table ages 'x' must be unique.")
  if (any(is.na(lt$lx)) || any(lt$lx < 0)) stop("'lx' must be nonnegative and not NA.")

  get_lx <- function(age) {
    idx <- match(age, lt$x)
    if (!is.na(idx)) return(lt$lx[idx])
    if (age == max(lt$x, na.rm = TRUE) + 1L) return(0)
    NA_real_
  }

  t_p_int <- function(age, tt) {
    if (tt == 0) return(1)
    l0 <- get_lx(age)
    l1 <- get_lx(age + tt)
    if (is.na(l0) || is.na(l1) || l0 <= 0) return(NA_real_)
    l1 / l0
  }

  v_fun <- function(tt) (1 + i)^(-tt)
  d <- i / (1 + i)

  # --- deferral factor: v^m * {}_m p_x ---
  defer <- v_fun(m) * t_p_int(x, m)
  if (is.na(defer)) stop("Deferral age x+m is outside the life table (or lx(x)=0).")

  y <- x + m

  # --- whole life: A_x = 1 - d * \ddot{a}_x (Finan, Sec. 37.1) ---
  if (type == "whole") {
    adue_y <- annuity_x(
      lt = lt, x = y, i = i,
      n = NULL, m = 0L, k = 1L,
      timing = "due", woolhouse = "none"
    )
    result <- defer * (1 - d * adue_y)
  } else if (n == 0L) {
    result <- 0
  } else {
    adue_y_n <- annuity_x(
      lt = lt, x = y, i = i,
      n = n, m = 0L, k = 1L,
      timing = "due", woolhouse = "none"
    )

    if (type == "endowment") {
      # A_{x:n} = 1 - d * \ddot{a}_{x:n} (Finan, Sec. 26.3.2 / Example 37.8)
      result <- defer * (1 - d * adue_y_n)
    } else {
      # A^1_{x:n} = 1 - d * \ddot{a}_{x:n} - v^n * n_p_y (Finan, Sec. 27)
      # = A_{x:n} - nEx
      n_p_y <- t_p_int(y, n)
      if (is.na(n_p_y)) stop("Life table does not support {}_n p_{x+m} needed for term insurance.")
      result <- defer * (1 - d * adue_y_n - v_fun(n) * n_p_y)
    }
  }

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x, i = i, n = if (type == "whole") NA_integer_ else n,
      m = m, type = type, apv = result
    ))
  }
  result
}
