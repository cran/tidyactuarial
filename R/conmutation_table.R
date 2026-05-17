#' Build an annual commutation table (discrete ages)
#'
#' @description
#' Constructs classical annual commutation functions Dx, Nx, Sx, Cx, Mx, Rx
#' from a lifetable defined at integer ages and an annual effective interest rate i.
#'
#' @param lt A lifetable object as produced by \code{\link{lifetable}}.
#'   Must contain columns \code{x} and \code{lx}. If available, \code{qx} or \code{px}
#'   may also be used, but \code{dx} is computed robustly from successive \code{lx}.
#' @param i Numeric. Annual effective interest rate (must satisfy \code{i > -1}).
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#'
#' @return A tibble with columns \code{x}, \code{lx}, \code{dx}, \code{v},
#'   \code{Dx}, \code{Nx}, \code{Sx}, \code{Cx}, \code{Mx}, \code{Rx}.
#'
#' @export
commutation_table <- function(lt, i, check = TRUE) {
  if (check) {
    if (!is.data.frame(lt)) stop("`lt` must be a data.frame/tibble.", call. = FALSE)
    if (!all(c("x", "lx") %in% names(lt))) {
      stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
    }
    if (!is.numeric(i) || length(i) != 1L || is.na(i) || i <= -1) {
      stop("`i` must be a single numeric value with i > -1.", call. = FALSE)
    }
    if (any(abs(lt$x - round(lt$x)) > 0, na.rm = TRUE)) {
      stop("`lt$x` must be integer-valued for annual commutation functions.", call. = FALSE)
    }
  }

  lt <- dplyr::arrange(lt, .data$x)

  x  <- lt$x
  lx <- lt$lx

  v <- 1 / (1 + i)

  # Robust annual deaths: dx = lx - l_{x+1}, closing at 0 by default
  dx <- lx - dplyr::lead(lx, default = 0)

  Dx <- (v^x) * lx
  Cx <- (v^(x + 1)) * dx

  rev_cumsum <- function(z) rev(cumsum(rev(z)))

  Nx <- rev_cumsum(Dx)
  Mx <- rev_cumsum(Cx)
  Sx <- rev_cumsum(Nx)
  Rx <- rev_cumsum(Mx)

  tibble::tibble(
    x  = x,
    lx = lx,
    dx = dx,
    v  = v,
    Dx = Dx,
    Nx = Nx,
    Sx = Sx,
    Cx = Cx,
    Mx = Mx,
    Rx = Rx
  )
}
