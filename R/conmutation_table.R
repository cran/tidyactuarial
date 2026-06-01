#' Build an annual commutation table (discrete ages)
#'
#' Constructs classical annual commutation functions \eqn{D_x}, \eqn{N_x},
#' \eqn{S_x}, \eqn{C_x}, \eqn{M_x}, and \eqn{R_x} from a life table defined
#' at integer ages, using compact actuarial notation.
#'
#' @description
#' The function builds annual commutation columns from \code{x}, \code{lx}, and
#' an interest-rate specification. The interest rate is converted internally to
#' an annual effective rate before constructing the discount factor.
#'
#' @param lt A life table object as produced by \code{\link{lifetable}}. It
#'   must contain columns \code{x} and \code{lx}. If available, \code{qx} or
#'   \code{px} may also be present, but \code{dx} is computed robustly from
#'   successive values of \code{lx}.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#'
#' @return A tibble with columns \code{x}, \code{lx}, \code{dx}, \code{v},
#'   \code{Dx}, \code{Nx}, \code{Sx}, \code{Cx}, \code{Mx}, and \code{Rx}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{i} is the
#' interest-rate input, \code{i_type} is the interest-rate type, and \code{m}
#' is the conversion frequency for nominal rates.
#'
#' The annual effective rate is obtained through \code{\link{standardize_interest}}.
#' If \eqn{i_e} denotes the annual effective rate, then
#' \deqn{v = \frac{1}{1+i_e}.}
#'
#' The annual deaths are computed as
#' \deqn{d_x = l_x - l_{x+1},}
#' closing the table with \eqn{l_{\omega+1}=0}. The main commutation functions
#' are then computed as
#' \deqn{D_x = v^x l_x, \qquad C_x = v^{x+1} d_x,}
#' with reverse cumulative sums used to obtain \eqn{N_x}, \eqn{S_x},
#' \eqn{M_x}, and \eqn{R_x}.
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:65,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000)
#' )
#'
#' commutation_table(
#'   lt = lt,
#'   i = 0.05
#' )
#'
#' commutation_table(
#'   lt = lt,
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12
#' )
#'
#' @export
commutation_table <- function(
    lt,
    i,
    i_type = "effective",
    m = 1L,
    check = TRUE
) {
  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(check)) {
    if (!is.data.frame(lt)) {
      stop("`lt` must be a data.frame or tibble.", call. = FALSE)
    }

    if (!all(c("x", "lx") %in% names(lt))) {
      stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
    }

    if (!is.numeric(lt$x)) {
      stop("Column `x` in `lt` must be numeric.", call. = FALSE)
    }

    if (!is.numeric(lt$lx)) {
      stop("Column `lx` in `lt` must be numeric.", call. = FALSE)
    }

    if (any(is.na(lt$x)) || any(!is.finite(lt$x))) {
      stop("Column `x` in `lt` must contain finite non-missing values.", call. = FALSE)
    }

    if (any(abs(lt$x - round(lt$x)) > 1e-10, na.rm = TRUE)) {
      stop("Column `x` in `lt` must contain integer ages.", call. = FALSE)
    }

    if (any(is.na(lt$lx)) || any(!is.finite(lt$lx)) || any(lt$lx < 0)) {
      stop("Column `lx` in `lt` must contain finite nonnegative values.", call. = FALSE)
    }

    if (!is.numeric(i) || length(i) != 1L || is.na(i) || !is.finite(i)) {
      stop("`i` must be a single finite numeric value.", call. = FALSE)
    }

    if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
      stop("`i_type` must be a single character string.", call. = FALSE)
    }

    valid_i_type <- c(
      "effective",
      "nominal_interest",
      "nominal_discount",
      "force"
    )

    if (!i_type %in% valid_i_type) {
      stop(
        "`i_type` must be one of: ",
        paste(sprintf("'%s'", valid_i_type), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
        !is.finite(m) || m < 1 || abs(m - round(m)) > 1e-10) {
      stop("`m` must be a single positive integer.", call. = FALSE)
    }
  }

  m <- as.integer(round(m))

  i_effective <- standardize_interest(
    type = i_type,
    rate = i,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  lt <- dplyr::arrange(lt, .data$x)

  x  <- lt$x
  lx <- lt$lx

  v <- 1 / (1 + i_effective)

  # Robust annual deaths: dx = lx - l_{x+1}, closing at 0 by default.
  dx <- lx - dplyr::lead(lx, default = 0)

  Dx <- (v^x) * lx
  Cx <- (v^(x + 1)) * dx

  rev_cumsum <- function(z) rev(cumsum(rev(z)))

  Nx <- rev_cumsum(Dx)
  Mx <- rev_cumsum(Cx)
  Sx <- rev_cumsum(Nx)
  Rx <- rev_cumsum(Mx)

  tibble::tibble(
    x = x,
    lx = lx,
    dx = dx,
    v = v,
    Dx = Dx,
    Nx = Nx,
    Sx = Sx,
    Cx = Cx,
    Mx = Mx,
    Rx = Rx
  )
}
