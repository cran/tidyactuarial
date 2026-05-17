#' t-year probability of decrement by cause j: t_qxj
#'
#' @description
#' Computes \eqn{{}_t q_x^{(j)}}, the probability that a life aged \code{x}
#' decrements by a specific cause \code{j} within \code{t} years, using a multiple
#' decrement table produced by \code{\link{md_table}}.
#'
#' @details
#' Let \eqn{q_x^{(j)}} be the annual decrement probability for cause \eqn{j} and
#' \eqn{q_x^{(\tau)}} be the total decrement probability. For integer \eqn{t},
#' \deqn{
#' {}_t q_x^{(j)} = \sum_{k=0}^{t-1} \left(\prod_{r=0}^{k-1} p_{x+r}^{(\tau)}\right) q_{x+k}^{(j)}.
#' }
#'
#' For non-integer \eqn{t = n + s} with \eqn{n = \lfloor t \rfloor} and
#' \eqn{s \in [0,1)}, this function supports fractional-age assumptions specified
#' by \code{frac} and uses the additional convention that the within-year cause
#' proportions remain constant:
#' \eqn{{}_s q_x^{(j)} = w_j \, {}_s q_x^{(\tau)}} where
#' \eqn{w_j = q_x^{(j)} / q_x^{(\tau)}} (and 0 when \eqn{q_x^{(\tau)}=0}).
#'
#' Supported fractional-age assumptions for the total decrement:
#' \itemize{
#'   \item \code{"UDD"}: \eqn{{}_s q_x^{(\tau)} = s q_x^{(\tau)}}.
#'   \item \code{"CF"}:  \eqn{{}_s p_x^{(\tau)} = (p_x^{(\tau)})^s}.
#'   \item \code{"Balducci"}: \eqn{{}_s q_x^{(\tau)} = \frac{s q_x^{(\tau)}}{1-(1-s)q_x^{(\tau)}}}.
#' }
#'
#' @param md A multiple decrement table produced by \code{\link{md_table}}.
#'   Must contain columns \code{x}, \code{p_total}, and the requested \code{cause}.
#' @param x Numeric vector. Starting age(s) (integer-valued).
#' @param t Numeric vector. Time horizon(s) in years (t >= 0). Can be non-integer
#'   if \code{frac} is supplied.
#' @param cause Character. Name of the cause column in \code{md}
#'   (e.g., \code{"q_death"}).
#' @param frac Character. Fractional-age assumption for non-integer \code{t}:
#'   one of \code{"UDD"}, \code{"CF"}, \code{"Balducci"}, or \code{NULL}.
#'   If \code{NULL} (default), \code{t} must be integer-valued.
#' @param tidy Logical. If \code{TRUE}, returns a tibble with columns
#'   \code{x}, \code{t}, \code{cause}, \code{frac}, and \code{tqxj}.
#' @param check Logical. If \code{TRUE}, performs input validation (default \code{TRUE}).
#' @param tol Numeric tolerance for integer checks (default \code{1e-10}).
#'
#' @return Numeric vector of \eqn{{}_t q_x^{(j)}} (or tibble if \code{tidy=TRUE}).
#'
#' @examples
#' qx_df <- tibble::tibble(
#'   x = 30:35,
#'   q_death = c(0.001, 0.0012, 0.0014, 0.0017, 0.0020, 1.0000),
#'   q_disability = c(0.002, 0.0021, 0.0022, 0.0023, 0.0024, 0.0000)
#' )
#' md <- md_table(qx_df, radix = 1e5, close = TRUE)
#' t_qxj(md, x = 30, t = 5, cause = "q_death")
#' t_qxj(md, x = 30, t = 2.5, cause = "q_death", frac = "CF", tidy = TRUE)
#'
#' @export
t_qxj <- function(
    md,
    x,
    t,
    cause,
    frac = NULL,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {

  if (missing(md)) stop("`md` is required.")
  if (missing(x))  stop("`x` is required.")
  if (missing(t))  stop("`t` is required.")
  if (missing(cause)) stop("`cause` is required.")

  if (check) {
    if (!is.data.frame(md)) stop("`md` must be a data.frame/tibble.", call. = FALSE)
    if (!all(c("x", "p_total", "q_total") %in% names(md))) {
      stop("`md` must contain columns `x`, `p_total`, and `q_total`.", call. = FALSE)
    }
    if (!is.character(cause) || length(cause) != 1L) {
      stop("`cause` must be a single character string.", call. = FALSE)
    }
    if (!cause %in% names(md)) stop("`cause` not found in `md`.", call. = FALSE)
    if (!is.numeric(x) || !is.numeric(t)) stop("`x` and `t` must be numeric.", call. = FALSE)
    if (any(t < -tol, na.rm = TRUE)) stop("`t` must be >= 0.", call. = FALSE)

    if (is.null(frac)) {
      if (any(abs(t - round(t)) > tol, na.rm = TRUE)) {
        stop("When `frac` is NULL, `t` must be integer-valued.", call. = FALSE)
      }
    } else {
      allowed_frac <- c("UDD", "CF", "Balducci")
      if (!is.character(frac) || length(frac) != 1L || !frac %in% allowed_frac) {
        stop("`frac` must be one of: ", paste(allowed_frac, collapse = ", "), ".", call. = FALSE)
      }
    }

    if (any(abs(x - round(x)) > tol, na.rm = TRUE)) {
      stop("`x` must be integer-valued.", call. = FALSE)
    }
  }

  md <- dplyr::as_tibble(md) |>
    dplyr::arrange(.data$x)

  # scalar recycling only
  L <- max(length(x), length(t))
  recycle_ok <- function(z) length(z) %in% c(1L, L)
  if (!recycle_ok(x) || !recycle_ok(t)) {
    stop("`x` and `t` must have the same length or be scalars.", call. = FALSE)
  }
  if (length(x) == 1L) x <- rep.int(x, L)
  if (length(t) == 1L) t <- rep.int(t, L)

  ages <- md$x
  p_tot <- md$p_total
  q_tot <- md$q_total
  q_j <- md[[cause]]

  # helper for fractional one-year cause probability at age index idx and fraction s
  frac_q_j <- function(idx, s) {
    if (s <= 0) return(0)
    qt <- q_tot[idx]
    qj <- q_j[idx]
    if (is.na(qt) || is.na(qj)) return(NA_real_)
    if (qt <= 0) return(0)  # no decrement within the year
    wj <- qj / qt

    if (is.null(frac)) {
      # should not happen if checked
      return(NA_real_)
    }

    if (frac == "UDD") {
      return(s * qj)
    }

    pt <- 1 - qt

    if (frac == "CF") {
      q_tau_s <- 1 - (pt^s)
      return(wj * q_tau_s)
    }

    # Balducci
    denom <- 1 - (1 - s) * qt
    q_tau_s <- (s * qt) / denom
    return(wj * q_tau_s)
  }

  out <- vapply(seq_len(L), function(i) {
    xi <- as.integer(round(x[i]))
    ti <- t[i]

    if (ti <= 0) return(0)

    idx0 <- match(xi, ages)
    if (is.na(idx0)) stop("Starting age not found in `md$x`.", call. = FALSE)

    if (is.null(frac)) {
      n <- as.integer(round(ti))
      idx_end <- idx0 + n - 1L
      if (idx_end > length(ages)) stop("Horizon exceeds the available ages in `md`.", call. = FALSE)

      # integer-horizon formula
      pvec <- p_tot[idx0:idx_end]
      qvec <- q_j[idx0:idx_end]
      surv_to_k <- cumprod(c(1, pvec))[1:n]  # survival to start of year k
      return(sum(surv_to_k * qvec))
    }

    # non-integer allowed when frac provided
    n <- floor(ti)
    s <- ti - n

    # integer part: k = 0..n-1
    ans <- 0
    if (n > 0) {
      idx_end_int <- idx0 + n - 1L
      if (idx_end_int > length(ages)) stop("Horizon exceeds the available ages in `md`.", call. = FALSE)
      pvec <- p_tot[idx0:idx_end_int]
      qvec <- q_j[idx0:idx_end_int]
      surv_to_k <- cumprod(c(1, pvec))[1:n]
      ans <- ans + sum(surv_to_k * qvec)
    }

    # fractional tail at age x+n
    if (s > 0) {
      idx_tail <- idx0 + n
      if (idx_tail > length(ages)) stop("Horizon exceeds the available ages in `md`.", call. = FALSE)

      # survival to time n (start of year n)
      if (n == 0) {
        surv_n <- 1
      } else {
        surv_n <- prod(p_tot[idx0:(idx0 + n - 1L)])
      }

      ans <- ans + surv_n * frac_q_j(idx_tail, s)
    }

    ans
  }, numeric(1))

  if (!tidy) return(out)

  tibble::tibble(
    x     = as.integer(round(x)),
    t     = t,
    cause = cause,
    frac  = if (is.null(frac)) NA_character_ else frac,
    tqxj  = out
  )
}
