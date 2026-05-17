#' Actuarial present value of a two-life insurance
#'
#' Computes the APV of a discrete two-life insurance with benefit 1 payable
#' at the end of the year of the triggering death, assuming independent lives.
#' The implementation uses standard identities that express insurance values
#' in terms of the annuity-due value (generalizing Finan, Sections 27 and 37
#' to two-life statuses from Sections 58--59).
#'
#' Supported contracts:
#' \itemize{
#'   \item \strong{Whole life}: \eqn{A_{xy} = 1 - d \, \ddot{a}_{xy}}
#'   \item \strong{n-year term}:
#'     \eqn{A^1_{xy:\overline{n}|} = 1 - d \, \ddot{a}_{xy:\overline{n}|}
#'     - v^n \, {}_np_{xy}}
#'   \item \strong{n-year endowment}:
#'     \eqn{A_{xy:\overline{n}|} = 1 - d \, \ddot{a}_{xy:\overline{n}|}}
#' }
#'
#' The two-life cohort determines the status:
#' \itemize{
#'   \item \code{cohort = "first"}: joint-life (triggers on first death)
#'   \item \code{cohort = "last"}: last-survivor (triggers on second death)
#' }
#'
#' Deferral by \eqn{m} years is applied as:
#' \eqn{v^m \, {}_mp_{xy} \times \text{(value at ages } x+m, y+m\text{)}}.
#'
#' @param lt Either:
#'   \itemize{
#'     \item a single life table data frame, used for both lives; or
#'     \item a list of two life tables \code{list(lt_x, lt_y)}, one for each life.
#'   }
#'   Each life table must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for life 1 at issue.
#' @param y Integer actuarial age for life 2 at issue.
#' @param i Annual effective interest rate (must be \code{> -1}).
#' @param n Integer term in years. Required for \code{type = "term"} or
#'   \code{type = "endowment"}. Ignored for \code{type = "whole"}.
#' @param m Integer deferral in years (default \code{0}).
#' @param type Insurance type: \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param cohort Status cohort: \code{"first"} (joint-life) or
#'   \code{"last"} (last-survivor).
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#'
#' @details
#' This function calls \code{\link{annuity_xy}} internally using annual
#' payments (\code{k = 1}, \code{timing = "due"}, \code{woolhouse = "none"}).
#'
#' For the term insurance adjustment \eqn{v^n \, {}_np}, the survival
#' probability is computed using \code{\link{t_pxy}} at integer time
#' \code{n} (fractional-age assumption is irrelevant at integer times).
#'
#' The endowment decomposes as:
#' \deqn{A_{xy:\overline{n}|} = A^1_{xy:\overline{n}|} +
#'   {}_nE_{xy}}
#' where \eqn{{}_nE_{xy} = v^n \, {}_np_{xy}} is the two-life pure
#' endowment.
#'
#' @return A single numeric APV value (per unit benefit), or a one-row
#'   tibble if \code{tidy = TRUE}.
#'
#' @seealso \code{\link{annuity_xy}} for two-life annuity APVs,
#'   \code{\link{insurance_x}} for single-life insurance,
#'   \code{insurance_multi} for N-life insurance,
#'   \code{\link{premium_xy}} for two-life benefit premiums,
#'   \code{\link{t_pxy}} for two-life survival.
#'
#' @examples
#' lt <- data.frame(x = 60:110, lx = seq(100000, 0, length.out = 51))
#'
#' # Whole life joint-life insurance (Finan, Sec. 58)
#' insurance_xy(lt, x = 60, y = 62, i = 0.06,
#'              type = "whole", cohort = "first")
#'
#' # 4-year term, joint-life
#' insurance_xy(lt, x = 60, y = 62, i = 0.06,
#'              n = 4, type = "term", cohort = "first")
#'
#' # 4-year endowment, last-survivor
#' insurance_xy(lt, x = 60, y = 62, i = 0.06,
#'              n = 4, type = "endowment", cohort = "last")
#'
#' # Different life tables for the two lives
#' lt_m <- data.frame(x = 60:110, lx = seq(100000, 0, length.out = 51))
#' lt_f <- data.frame(x = 60:110, lx = seq(100000, 1000, length.out = 51))
#' insurance_xy(list(lt_m, lt_f), x = 60, y = 62, i = 0.06,
#'              type = "whole", cohort = "first")
#'
#' @export
insurance_xy <- function(
    lt, x, y, i,
    n = NULL,
    m = 0L,
    type = c("whole", "term", "endowment"),
    cohort = c("first", "last"),
    tidy = FALSE
) {
  type   <- match.arg(type)
  cohort <- match.arg(cohort)

  # --- resolve life table input ---
  if (is.data.frame(lt)) {
    if (!all(c("x", "lx") %in% names(lt))) {
      stop("Life table must contain columns 'x' and 'lx'.")
    }
    lt_use <- lt
  } else if (is.list(lt) && length(lt) == 2L &&
             all(vapply(lt, is.data.frame, logical(1)))) {
    if (!all(c("x", "lx") %in% names(lt[[1]]))) {
      stop("First life table must contain columns 'x' and 'lx'.")
    }
    if (!all(c("x", "lx") %in% names(lt[[2]]))) {
      stop("Second life table must contain columns 'x' and 'lx'.")
    }
    lt_use <- lt
  } else {
    stop("`lt` must be either one life table or a list of two life tables.")
  }

  # --- checks ---
  if (missing(i) || !is.numeric(i) || length(i) != 1L ||
      is.na(i) || i <= -1) {
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
      m < 0 || abs(m - round(m)) > 1e-10) {
    stop("'m' must be a single nonnegative integer.")
  }
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("'tidy' must be TRUE or FALSE.")
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))

  if (type %in% c("term", "endowment")) {
    if (is.null(n)) {
      stop("'n' must be provided for type 'term' or 'endowment'.")
    }
    if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
        n < 0 || abs(n - round(n)) > 1e-10) {
      stop("'n' must be a single nonnegative integer.")
    }
    n <- as.integer(round(n))
  }

  v_fun <- function(tt) (1 + i)^(-tt)
  d <- i / (1 + i)

  # Map cohort -> status for t_pxy()
  status <- if (cohort == "first") "joint" else "last"

  # --- deferral: v^m * m_p_status ---
  p_def <- t_pxy(
    lt = lt_use, x = x, y = y,
    t = m, frac = "UDD", status = status
  )
  if (is.na(p_def)) {
    stop("Cannot compute survival to deferral under the life table.")
  }
  defer <- v_fun(m) * p_def

  # Start ages after deferral
  x1 <- x + m
  y1 <- y + m

  # --- whole life: A = 1 - d * \ddot{a} ---
  if (type == "whole") {
    adue <- annuity_xy(
      lt = lt_use, x = x1, y = y1, i = i,
      cohort = cohort, n = NULL, m = 0L,
      k = 1L, timing = "due", woolhouse = "none"
    )
    result <- defer * (1 - d * adue)

  } else if (n == 0L) {
    result <- 0

  } else {
    adue_n <- annuity_xy(
      lt = lt_use, x = x1, y = y1, i = i,
      cohort = cohort, n = n, m = 0L,
      k = 1L, timing = "due", woolhouse = "none"
    )

    if (type == "endowment") {
      # A_{xy:n} = 1 - d * \ddot{a}_{xy:n}
      result <- defer * (1 - d * adue_n)
    } else {
      # A^1_{xy:n} = 1 - d * \ddot{a}_{xy:n} - v^n * n_p_status
      p_n <- t_pxy(
        lt = lt_use, x = x1, y = y1,
        t = n, frac = "UDD", status = status
      )
      if (is.na(p_n)) {
        stop("Cannot compute {}_n p at deferred ages.")
      }
      result <- defer * (1 - d * adue_n - v_fun(n) * p_n)
    }
  }

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x, y = y, i = i,
      n = if (type == "whole") NA_integer_ else n,
      m = m, type = type, cohort = cohort, apv = result
    ))
  }

  result
}
