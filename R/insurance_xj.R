#' Cause-specific term/whole-life insurance APV under multiple decrements: insurance_xj
#'
#' @description
#' Computes the actuarial present value (APV) of an annual (discrete) insurance
#' that pays \code{benefit} at the end of the year of decrement by a specified
#' cause \code{j}, using a multiple decrement table produced by
#' \code{\link{md_table}}.
#'
#' @details
#' Let \eqn{q_{x+k}^{(j)}} be the one-year decrement probability for cause \eqn{j}
#' at age \eqn{x+k}, and \eqn{p_{x+r}^{(\tau)}} be the one-year total survival
#' probability (any cause) at age \eqn{x+r}. For a product with deferment \eqn{m}
#' and term \eqn{n}, with benefit payable at the end of the year of decrement by
#' cause \eqn{j}, the APV is:
#' \deqn{
#' \sum_{k=m}^{m+n-1} v^{k+1}\left(\prod_{r=0}^{k-1} p_{x+r}^{(\tau)}\right) q_{x+k}^{(j)},
#' }
#' where \eqn{v = (1+i)^{-1}} and \eqn{i} is the annual effective interest rate.
#'
#' If \code{product = "whole"}, the function sets \code{n} to the remaining length
#' of the table after deferment (i.e., whole life over the available ages).
#'
#' @param md A multiple decrement table produced by \code{\link{md_table}}.
#'   Must contain columns \code{x}, \code{p_total}, and the requested \code{cause}.
#' @param x Integer age(s) at issue.
#' @param i Annual effective interest rate(s). Must satisfy \code{i > -1}.
#' @param cause Character. Name of the cause column in \code{md}
#'   (e.g., \code{"q_death"}).
#' @param product Character. Insurance type: \code{"whole"} or \code{"term"}.
#' @param benefit Numeric benefit amount payable at the end of the year of
#'   decrement by the specified cause (default \code{1}).
#' @param n Integer term length in years (required when \code{product = "term"}).
#' @param m Integer deferment in years (default \code{0}).
#' @param tidy Logical. If \code{TRUE}, returns a tibble with inputs and \code{insurance_xj}.
#' @param check Logical. If \code{TRUE}, performs input validation.
#' @param tol Numeric tolerance for integer checks.
#'
#' @return Numeric vector of APVs (or a tibble if \code{tidy = TRUE}).
#'
#' @seealso \code{\link{t_qxj}} for cause-specific decrement probabilities,
#'   \code{\link{lt_tau}} to build a single-decrement lifetable for the total decrement.
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
#' insurance_xj(md, x = 30, i = 0.05, cause = "q_death", product = "term", n = 5)
#'
#' # Whole-life (over available ages), 2-year deferred
#' insurance_xj(md, x = 30, i = 0.05, cause = "q_death", product = "whole", m = 2, tidy = TRUE)
#'
#' @export
insurance_xj <- function(
    md,
    x,
    i,
    cause,
    product = c("whole", "term"),
    benefit = 1,
    n = NULL,
    m = 0L,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {
  product <- match.arg(product)

  # --- checks (light, CRAN-friendly) ---
  if (missing(md)) stop("`md` is required.")
  if (!is.data.frame(md)) stop("`md` must be a data.frame/tibble.", call. = FALSE)
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

  x <- as.numeric(x)
  i <- as.numeric(i)
  m <- as.numeric(m)
  benefit <- as.numeric(benefit)

  if (check) {
    if (any(!is.finite(x))) stop("`x` must be finite.", call. = FALSE)
    if (any(abs(x - round(x)) > tol)) stop("`x` must be integer ages.", call. = FALSE)
    if (any(!is.finite(i))) stop("`i` must be finite.", call. = FALSE)
    if (any(i <= -1)) stop("`i` must be greater than -1.", call. = FALSE)
    if (any(!is.finite(m))) stop("`m` must be finite.", call. = FALSE)
    if (any(abs(m - round(m)) > tol)) stop("`m` must be an integer.", call. = FALSE)
    if (any(m < 0)) stop("`m` must be nonnegative.", call. = FALSE)
    if (!is.null(n)) {
      n <- as.numeric(n)
      if (any(!is.finite(n))) stop("`n` must be finite.", call. = FALSE)
      if (any(abs(n - round(n)) > tol)) stop("`n` must be an integer.", call. = FALSE)
      if (any(n < 0)) stop("`n` must be nonnegative.", call. = FALSE)
    }
    if (any(!is.finite(benefit))) stop("`benefit` must be finite.", call. = FALSE)
    if (any(benefit < 0)) stop("`benefit` must be nonnegative.", call. = FALSE)

    if (product == "term" && is.null(n)) {
      stop("For product = 'term', `n` must be provided.", call. = FALSE)
    }
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  if (!is.null(n)) n <- as.integer(round(n))

  md <- dplyr::as_tibble(md) |>
    dplyr::arrange(.data$x)

  ages <- md$x
  p_tot <- md$p_total
  q_j   <- md[[cause]]

  # --- recycle inputs to common length (scalar recycling only) ---
  lengths <- c(length(x), length(i), length(m), length(benefit))
  if (!is.null(n)) lengths <- c(lengths, length(n))
  L <- max(lengths)

  recycle_ok <- function(z) length(z) %in% c(1L, L)
  if (!all(vapply(list(x, i, m, benefit), recycle_ok, logical(1)))) {
    stop("`x`, `i`, `m`, and `benefit` must have the same length or length 1.", call. = FALSE)
  }
  if (!is.null(n) && !recycle_ok(n)) {
    stop("`n` must have the same length as other inputs or length 1.", call. = FALSE)
  }

  if (length(x) == 1L && L > 1L) x <- rep.int(x, L)
  if (length(i) == 1L && L > 1L) i <- rep.int(i, L)
  if (length(m) == 1L && L > 1L) m <- rep.int(m, L)
  if (length(benefit) == 1L && L > 1L) benefit <- rep.int(benefit, L)
  if (!is.null(n) && length(n) == 1L && L > 1L) n <- rep.int(n, L)

  apv <- vapply(seq_len(L), function(k) {
    xk <- x[k]
    ik <- i[k]
    mk <- m[k]
    bk <- benefit[k]

    idx0 <- match(xk, ages)
    if (is.na(idx0)) stop("Starting age not found in `md$x`.", call. = FALSE)

    available_years <- length(ages) - idx0 + 1L  # includes last age row
    if (mk > available_years) stop("Deferment `m` exceeds the available ages in `md`.", call. = FALSE)

    nk <- if (product == "whole") (available_years - mk) else n[k]
    if (is.na(nk)) nk <- 0L
    if (nk < 0L) stop("`n` must be nonnegative.", call. = FALSE)
    if (nk == 0L) return(0)

    K <- mk + nk  # total years from issue to end of term
    if (K > available_years) stop("Term plus deferment exceeds the available ages in `md`.", call. = FALSE)

    idx_end <- idx0 + K - 1L
    pvec <- p_tot[idx0:idx_end]
    qvec <- q_j[idx0:idx_end]

    v <- (1 + ik)^(-1)

    surv_to_start <- cumprod(c(1, pvec))[1:K]  # survival to start of each year
    disc_endyear  <- v^(1:K)                   # discount to end of each year

    idx_sum <- (mk + 1L):(mk + nk)             # years m..m+n-1 (1-based)
    bk * sum(disc_endyear[idx_sum] * surv_to_start[idx_sum] * qvec[idx_sum])
  }, numeric(1))

  if (!tidy) return(apv)

  tibble::tibble(
    x       = x,
    i       = i,
    m       = m,
    n       = if (product == "whole") NA_integer_ else n,
    product = product,
    cause   = cause,
    benefit = benefit,
    insurance_xj = apv
  )
}

#' @rdname insurance_xj
#' @export
A_xj <- insurance_xj
