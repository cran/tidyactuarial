#' Actuarial present value of a multi-life insurance (N independent lives)
#'
#' Computes the APV of a discrete life insurance contingent on a multi-life
#' status built from multiple independent lives. The benefit is paid at the
#' end of the year of death of the status.
#'
#' @param lt Life table data frame with columns \code{x} and \code{lx}.
#' @param ages Integer vector of issue ages for the lives.
#' @param i Annual effective interest rate (must be \code{> -1}).
#' @param product Insurance type: \code{"whole"}, \code{"term"},
#'   \code{"endowment"}, or \code{"pure_endowment"}.
#' @param cohort Status definition: \code{"first"} (joint-life, first death)
#'   or \code{"last"} (last-survivor, second death).
#' @param benefit Numeric scalar benefit amount (default \code{1}).
#' @param n Integer term in years after deferment. Required for
#'   \code{"term"}, \code{"endowment"}, and \code{"pure_endowment"}.
#' @param m Integer deferment in years (default \code{0}).
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#'
#' @details
#' The implementation uses the standard equivalence between insurance and
#' annuity-due on the same status (generalizing Finan, Sections 27 and 37
#' to N-life statuses from Sections 56--59):
#' \deqn{A_{\text{status}} = 1 - d \, \ddot{a}_{\text{status}},
#'   \quad d = \frac{i}{1+i}.}
#'
#' Supported products:
#' \itemize{
#'   \item \strong{Whole life}: \eqn{A = 1 - d \, \ddot{a}}
#'   \item \strong{n-year term}:
#'     \eqn{A^1 = A_{\text{endow}} - {}_nE_{\text{status}}}
#'   \item \strong{n-year endowment}:
#'     \eqn{A = 1 - d \, \ddot{a}_{\overline{n}|}}
#'   \item \strong{Pure endowment}:
#'     \eqn{{}_nE = v^{m+n} \cdot P(\text{status alive at } m+n)}
#' }
#'
#' Deferral by \eqn{m} years is applied as
#' \eqn{v^m \cdot P(\text{status alive at } m) \times
#' \text{(value at shifted ages)}}.
#'
#' The annuity-due used internally is computed via
#' \code{\link{annuity_multi}}.
#'
#' @return A numeric scalar APV, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{insurance_xy}} for the optimized two-life version,
#'   \code{\link{insurance_x}} for single-life insurance,
#'   \code{\link{annuity_multi}} for multi-life annuity APVs,
#'   \code{\link{Var_insurance_x}} for insurance variance.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:110,
#'   lx = seq(100000, 0, length.out = 51)
#' )
#'
#' # Joint-life (first-death) whole life insurance
#' insurance_multi(lt, ages = c(60, 62), i = 0.05,
#'                 product = "whole", cohort = "first")
#'
#' # Last-survivor 10-year term, deferred 5 years
#' insurance_multi(lt, ages = c(60, 62), i = 0.05,
#'                 product = "term", cohort = "last",
#'                 n = 10, m = 5)
#'
#' # Endowment: verify decomposition A = A^1 + nE
#' A_endow <- insurance_multi(lt, ages = c(60, 62), i = 0.05,
#'                            product = "endowment", cohort = "first",
#'                            n = 10)
#' A_term  <- insurance_multi(lt, ages = c(60, 62), i = 0.05,
#'                            product = "term", cohort = "first",
#'                            n = 10)
#' A_pe    <- insurance_multi(lt, ages = c(60, 62), i = 0.05,
#'                            product = "pure_endowment",
#'                            cohort = "first", n = 10)
#' c(endowment = A_endow, term_plus_pe = A_term + A_pe)
#'
#' # Three lives
#' insurance_multi(lt, ages = c(60, 65, 70), i = 0.05,
#'                 product = "whole", cohort = "first")
#'
#' # Tidy output
#' insurance_multi(lt, ages = c(60, 62), i = 0.05,
#'                 product = "term", cohort = "first",
#'                 n = 10, tidy = TRUE)
#'
#' @export
insurance_multi <- function(
    lt, ages, i,
    product = c("whole", "term", "endowment", "pure_endowment"),
    cohort = c("first", "last"),
    benefit = 1,
    n = NULL,
    m = 0L,
    tidy = FALSE,
    check = TRUE
) {
  product <- match.arg(product)
  cohort  <- match.arg(cohort)

  # --- checks ---
  if (missing(i) || !is.numeric(i) || length(i) != 1L ||
      is.na(i) || i <= -1) {
    stop("'i' must be a single numeric rate > -1.")
  }
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!all(c("x", "lx") %in% names(lt))) {
    stop("Life table must contain columns 'x' and 'lx'.")
  }
  if (!is.numeric(ages) || length(ages) < 1) {
    stop("'ages' must be a non-empty numeric vector.")
  }
  ages <- as.integer(round(ages))
  m <- as.integer(round(m))
  if (m < 0) stop("'m' must be nonnegative.")

  if (!is.numeric(benefit) || length(benefit) != 1L ||
      is.na(benefit) || benefit < 0) {
    stop("'benefit' must be a single nonnegative number.")
  }

  if (product %in% c("term", "endowment", "pure_endowment")) {
    if (is.null(n)) {
      stop("'n' required for term/endowment/pure_endowment.")
    }
    if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
        n < 0 || abs(n - round(n)) > 1e-10) {
      stop("'n' must be a single nonnegative integer.")
    }
    n <- as.integer(round(n))
    if (n == 0L && product == "term") {
      if (isTRUE(tidy)) {
        return(tibble::tibble(
          product = product, cohort = cohort, benefit = benefit,
          m = m, n = 0L, i = i, apv = 0
        ))
      }
      return(0)
    }
  }

  omega <- max(lt$x, na.rm = TRUE)
  v_fun <- function(tt) (1 + i)^(-tt)
  d <- i / (1 + i)

  # --- status survival at integer time via t_px ---
  P_status <- function(tt) {
    p_vec <- vapply(ages, function(a) {
      t_px(lt, x = a, t = tt, frac = "UDD", check = FALSE)
    }, numeric(1))
    if (anyNA(p_vec)) return(NA_real_)
    if (cohort == "first") prod(p_vec) else 1 - prod(1 - p_vec)
  }

  # --- annuity-due at shifted ages (valued at time m) ---
  ann_due_at_m <- function(term_n) {
    shifted <- ages + m
    # max horizon depends on status
    if (cohort == "first") {
      max_n <- max(0L, omega - max(shifted))
    } else {
      max_n <- max(0L, omega - min(shifted))
    }
    if (is.null(term_n)) term_n <- max_n
    if (term_n > max_n) {
      stop("Life table horizon insufficient for requested term.")
    }

    annuity_multi(
      lt = lt, ages = shifted, i = i,
      annuity = "cohort", cohort = cohort,
      n = term_n, m = 0L, k = 1L,
      timing = "due", woolhouse = "none"
    )
  }

  # --- deferral multiplier: v^m * P_status(m) ---
  Pm <- P_status(m)
  if (!is.finite(Pm) || is.na(Pm)) {
    stop("Cannot compute status survival at deferment time m.")
  }
  defer_mult <- v_fun(m) * Pm

  # --- pure endowment: benefit * v^(m+n) * P_status(m+n) ---
  pure_endow_apv <- function(nn) {
    Pmn <- P_status(m + nn)
    if (!is.finite(Pmn) || is.na(Pmn)) {
      stop("Cannot compute status survival at time m+n.")
    }
    benefit * v_fun(m + nn) * Pmn
  }

  # --- compute APV by product ---
  if (product == "pure_endowment") {
    apv <- pure_endow_apv(n)

  } else if (product == "endowment") {
    # A_endow = 1 - d * \ddot{a} (at shifted ages)
    a_due <- ann_due_at_m(n)
    apv <- benefit * defer_mult * (1 - d * a_due)

  } else if (product == "term") {
    # A_term = A_endow - nE_status
    a_due <- ann_due_at_m(n)
    apv_endow <- benefit * defer_mult * (1 - d * a_due)
    apv_pe <- pure_endow_apv(n)
    apv <- apv_endow - apv_pe

  } else {
    # whole life
    a_due <- ann_due_at_m(NULL)
    apv <- benefit * defer_mult * (1 - d * a_due)
  }

  if (!isTRUE(tidy)) return(apv)

  tibble::tibble(
    product = product, cohort = cohort, benefit = benefit,
    m = m,
    n = if (is.null(n)) NA_integer_ else n,
    i = i, apv = apv
  )
}
