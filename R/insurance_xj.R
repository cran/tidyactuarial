#' Cause-specific term/whole-life insurance APV under multiple decrements
#'
#' Computes the actuarial present value (APV) of an annual discrete insurance
#' that pays \code{benefit} at the end of the year of decrement by a specified
#' cause \code{j}, using a multiple decrement table produced by
#' \code{\link{md_table}}.
#'
#' @description
#' This function evaluates cause-specific insurance benefits under a multiple
#' decrement model. It supports whole-life and term insurance, integer
#' deferment, vectorized issue ages and interest-rate assumptions, and compact
#' actuarial notation.
#'
#' @details
#' Let \eqn{q_{x+k}^{(j)}} be the one-year decrement probability for cause
#' \eqn{j} at age \eqn{x+k}, and let \eqn{p_{x+r}^{(\tau)}} be the one-year
#' total survival probability against all decrements at age \eqn{x+r}.
#'
#' For a product with deferment \eqn{h} and term \eqn{n}, with benefit payable
#' at the end of the year of decrement by cause \eqn{j}, the APV is:
#' \deqn{
#' \sum_{k=h}^{h+n-1}
#' v^{k+1}
#' \left(\prod_{r=0}^{k-1} p_{x+r}^{(\tau)}\right)
#' q_{x+k}^{(j)}.
#' }
#'
#' Here \eqn{v = (1+i_e)^{-1}}, where \eqn{i_e} is the annual effective
#' interest rate obtained from \code{i}, \code{i_type}, and \code{m}.
#'
#' If \code{type = "whole"}, the function sets \code{n} to the remaining length
#' of the table after deferment, that is, whole life over the available ages.
#'
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{md} is a multiple decrement table, \code{x} is
#' age at issue, \code{i} is the interest rate, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' \code{n} is the term, and \code{h} is the deferment period.
#'
#' @param md A multiple decrement table produced by \code{\link{md_table}}.
#'   Must contain columns \code{x}, \code{p_total}, and the requested
#'   \code{cause}.
#' @param x Integer age(s) at issue.
#' @param i Annual interest-rate input. Must produce an annual effective rate
#'   greater than \code{-1}.
#' @param cause Character scalar. Name of the cause column in \code{md}, for
#'   example \code{"q_death"}.
#' @param type Character scalar. Insurance type: \code{"whole"} or
#'   \code{"term"}.
#' @param benefit Numeric benefit amount payable at the end of the year of
#'   decrement by the specified cause. Default is \code{1}.
#' @param n Integer term length in years. Required when \code{type = "term"}.
#' @param h Integer deferment period in years. Default is \code{0}.
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param tidy Logical. If \code{TRUE}, returns a tibble with inputs and
#'   \code{insurance_xj}.
#' @param check Logical. If \code{TRUE}, performs input validation.
#' @param tol Numeric tolerance for integer checks.
#'
#' @return Numeric vector of APVs, or a tibble if \code{tidy = TRUE}.
#'
#' @seealso \code{\link{t_qxj}} for cause-specific decrement probabilities,
#'   \code{\link{lt_tau}} to build a single-decrement life table for the total
#'   decrement.
#'
#' @family multiple-decrements
#'
#' @examples
#' qx_df <- tibble::tibble(
#'   x = 30:35,
#'   q_death = c(0.001, 0.0012, 0.0014, 0.0017, 0.0020, 1.0000),
#'   q_disability = c(0.002, 0.0021, 0.0022, 0.0023, 0.0024, 0.0000)
#' )
#' md <- md_table(qx_df, radix = 1e5, close = TRUE)
#'
#' # 5-year term cause-specific insurance for death, i = 5%
#' insurance_xj(
#'   md = md,
#'   x = 30,
#'   i = 0.05,
#'   cause = "q_death",
#'   type = "term",
#'   n = 5
#' )
#'
#' # Whole-life over available ages, 2-year deferred
#' insurance_xj(
#'   md = md,
#'   x = 30,
#'   i = 0.05,
#'   cause = "q_death",
#'   type = "whole",
#'   h = 2,
#'   tidy = TRUE
#' )
#'
#' @export
insurance_xj <- function(
    md,
    x,
    i,
    cause,
    type = c("whole", "term"),
    benefit = 1,
    n = NULL,
    h = 0L,
    i_type = "effective",
    m = 1L,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {
  type <- match.arg(type)

  # --- checks (light, CRAN-friendly) ---
  if (missing(md)) stop("`md` is required.")
  if (!is.data.frame(md)) stop("`md` must be a data frame or tibble.", call. = FALSE)
  if (!all(c("x", "p_total") %in% names(md))) {
    stop("`md` must contain columns `x` and `p_total`.", call. = FALSE)
  }

  if (missing(x) || missing(i) || missing(cause)) {
    stop("`x`, `i`, and `cause` are required.", call. = FALSE)
  }

  if (!is.character(cause) || length(cause) != 1L) {
    stop("`cause` must be a single character string.", call. = FALSE)
  }
  if (!cause %in% names(md)) stop("`cause` not found in `md`.", call. = FALSE)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
  }

  x <- as.numeric(x)
  i <- as.numeric(i)
  h <- as.numeric(h)
  m <- as.numeric(m)
  benefit <- as.numeric(benefit)

  if (check) {
    if (any(!is.finite(x))) stop("`x` must be finite.", call. = FALSE)
    if (any(abs(x - round(x)) > tol)) stop("`x` must be integer ages.", call. = FALSE)

    if (any(!is.finite(i))) stop("`i` must be finite.", call. = FALSE)

    if (!is.character(i_type)) {
      stop("`i_type` must be a character vector.", call. = FALSE)
    }

    valid_i_type <- c(
      "effective",
      "nominal_interest",
      "nominal_discount",
      "force"
    )

    bad_i_type <- !is.na(i_type) & !(i_type %in% valid_i_type)
    if (any(bad_i_type)) {
      stop(
        "`i_type` must contain only: ",
        paste(sprintf("'%s'", valid_i_type), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (any(!is.finite(m))) stop("`m` must be finite.", call. = FALSE)
    if (any(abs(m - round(m)) > tol)) stop("`m` must contain integer values.", call. = FALSE)
    if (any(m < 1)) stop("`m` must contain positive integers.", call. = FALSE)

    if (any(!is.finite(h))) stop("`h` must be finite.", call. = FALSE)
    if (any(abs(h - round(h)) > tol)) stop("`h` must be an integer.", call. = FALSE)
    if (any(h < 0)) stop("`h` must be nonnegative.", call. = FALSE)

    if (!is.null(n)) {
      n <- as.numeric(n)
      if (any(!is.finite(n))) stop("`n` must be finite.", call. = FALSE)
      if (any(abs(n - round(n)) > tol)) stop("`n` must be an integer.", call. = FALSE)
      if (any(n < 0)) stop("`n` must be nonnegative.", call. = FALSE)
    }

    if (any(!is.finite(benefit))) stop("`benefit` must be finite.", call. = FALSE)
    if (any(benefit < 0)) stop("`benefit` must be nonnegative.", call. = FALSE)

    if (type == "term" && is.null(n)) {
      stop("For type = 'term', `n` must be provided.", call. = FALSE)
    }
  }

  x <- as.integer(round(x))
  h <- as.integer(round(h))
  m <- as.integer(round(m))
  if (!is.null(n)) n <- as.integer(round(n))

  md <- dplyr::as_tibble(md) |>
    dplyr::arrange(.data$x)

  ages <- md$x
  p_tot <- md$p_total
  q_j   <- md[[cause]]

  # --- recycle inputs to common length (scalar recycling only) ---
  lengths <- c(length(x), length(i), length(i_type), length(m), length(h), length(benefit))
  if (!is.null(n)) lengths <- c(lengths, length(n))
  L <- max(lengths)

  recycle_ok <- function(z) length(z) %in% c(1L, L)
  if (!all(vapply(list(x, i, i_type, m, h, benefit), recycle_ok, logical(1)))) {
    stop(
      "`x`, `i`, `i_type`, `m`, `h`, and `benefit` must have the same length or length 1.",
      call. = FALSE
    )
  }
  if (!is.null(n) && !recycle_ok(n)) {
    stop("`n` must have the same length as other inputs or length 1.", call. = FALSE)
  }

  if (length(x) == 1L && L > 1L) x <- rep.int(x, L)
  if (length(i) == 1L && L > 1L) i <- rep.int(i, L)
  if (length(i_type) == 1L && L > 1L) i_type <- rep.int(i_type, L)
  if (length(m) == 1L && L > 1L) m <- rep.int(m, L)
  if (length(h) == 1L && L > 1L) h <- rep.int(h, L)
  if (length(benefit) == 1L && L > 1L) benefit <- rep.int(benefit, L)
  if (!is.null(n) && length(n) == 1L && L > 1L) n <- rep.int(n, L)

  i_effective <- standardize_interest(
    type = i_type,
    rate = i,
    m = m
  )

  if (any(is.na(i_effective)) ||
      any(!is.finite(i_effective)) ||
      any(i_effective <= -1)) {
    stop(
      "The standardized annual effective interest rates must be greater than -1.",
      call. = FALSE
    )
  }

  apv <- vapply(seq_len(L), function(idx) {
    x_idx <- x[idx]
    i_idx <- i_effective[idx]
    h_idx <- h[idx]
    benefit_idx <- benefit[idx]

    idx0 <- match(x_idx, ages)
    if (is.na(idx0)) stop("Starting age not found in `md$x`.", call. = FALSE)

    available_years <- length(ages) - idx0 + 1L  # includes last age row
    if (h_idx > available_years) {
      stop("Deferment `h` exceeds the available ages in `md`.", call. = FALSE)
    }

    n_idx <- if (type == "whole") {
      available_years - h_idx
    } else {
      n[idx]
    }

    if (is.na(n_idx)) n_idx <- 0L
    if (n_idx < 0L) stop("`n` must be nonnegative.", call. = FALSE)
    if (n_idx == 0L) return(0)

    K <- h_idx + n_idx  # total years from issue to end of term
    if (K > available_years) {
      stop("Term plus deferment exceeds the available ages in `md`.", call. = FALSE)
    }

    idx_end <- idx0 + K - 1L
    pvec <- p_tot[idx0:idx_end]
    qvec <- q_j[idx0:idx_end]

    v <- (1 + i_idx)^(-1)

    surv_to_start <- cumprod(c(1, pvec))[1:K]  # survival to start of each year
    disc_endyear  <- v^(1:K)                   # discount to end of each year

    idx_sum <- (h_idx + 1L):(h_idx + n_idx)    # years h..h+n-1 (1-based)
    benefit_idx * sum(disc_endyear[idx_sum] * surv_to_start[idx_sum] * qvec[idx_sum])
  }, numeric(1))

  if (!tidy) return(apv)

  tibble::tibble(
    x = x,
    i = i,
    i_type = i_type,
    m = m,
    h = h,
    n = if (type == "whole") NA_integer_ else n,
    type = type,
    cause = cause,
    benefit = benefit,
    insurance_xj = apv
  )
}

#' @rdname insurance_xj
#' @export
A_xj <- insurance_xj
