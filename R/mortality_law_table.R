#' Generate a tidy life table from a theoretical mortality law
#'
#' @description
#' Creates a tidy life table with one row per integer age from a parametric
#' mortality law. The output follows tidyactuarial conventions and includes
#' columns such as \code{x}, \code{qx}, \code{px}, \code{lx}, \code{dx}, and
#' \code{mx}.
#'
#' @section Supported laws:
#' \itemize{
#'   \item \strong{Exponential}: \eqn{\mu_x = \lambda}
#'   \item \strong{Gompertz}: \eqn{\mu_x = B c^x}
#'   \item \strong{Makeham}: \eqn{\mu_x = A + B c^x}
#'   \item \strong{Weibull}: \eqn{\mu_x = (k/\lambda)(x/\lambda)^{k-1}}
#'   \item \strong{Logistic}: \eqn{\mu_x = (A + B c^x)/(1 + C c^x)}
#'   \item \strong{DeMoivre}: finite lifetime law with \eqn{q_x = 1/(\omega - x)}
#'   for \eqn{x < \omega}
#'   \item \strong{Beta}: scaled lifetime model \eqn{X/\omega \sim Beta(\alpha,\beta)}
#'   \item \strong{HeligmanPollard}: odds model returning \eqn{q_x} directly
#' }
#'
#' @section Converting mu(x) to qx:
#' For laws defined by a force of mortality \eqn{\mu_x}, the one-year death
#' probability \eqn{q_x} is obtained using \code{frac}:
#' \itemize{
#'   \item \code{"CF"}: \eqn{q_x = 1 - \exp(-\mu_x)}
#'   \item \code{"UDD"}: \eqn{q_x = \mu_x}
#'   \item \code{"Balducci"}: \eqn{q_x = \mu_x/(1+\mu_x)}
#' }
#'
#' @section Direct qx laws:
#' \code{"DeMoivre"}, \code{"Beta"}, and \code{"HeligmanPollard"} define
#' \code{qx} directly. For these laws, \code{frac} is not used to derive
#' \code{qx}.
#'
#' @section Closure:
#' If \code{close = TRUE}, the last age is forced to close the table by setting
#' \code{qx[x_max] = 1} and \code{px[x_max] = 0}.
#'
#' @param law Character. Mortality law. One of \code{"Exponential"},
#'   \code{"Gompertz"}, \code{"Makeham"}, \code{"Weibull"}, \code{"Logistic"},
#'   \code{"DeMoivre"}, \code{"Beta"}, or \code{"HeligmanPollard"}.
#' @param x_min Integer. Minimum age, inclusive.
#' @param x_max Integer. Maximum age, inclusive. Must satisfy \code{x_min < x_max}.
#' @param ... Named law parameters. These values override \code{params}. Placing
#'   \code{...} before arguments such as \code{close} and \code{check} avoids
#'   partial-matching conflicts with the Gompertz/Makeham parameter \code{c}.
#' @param params Named list of law parameters, or \code{NULL}. Direct parameters
#'   supplied through \code{...} override values in \code{params}.
#' @param frac Character. Within-year assumption used to convert \eqn{\mu_x} to
#'   \eqn{q_x}. One of \code{"CF"}, \code{"UDD"}, or \code{"Balducci"}.
#' @param radix Numeric. Starting cohort size at age \code{x_min}.
#' @param close Logical. If \code{TRUE}, forces the last age to close.
#' @param ax Numeric. Average fraction of the year lived by those dying.
#' @param check Logical. If \code{TRUE}, performs strict input validation.
#' @param tol Numeric. Tolerance used in checks.
#'
#' @return A tibble with columns \code{x}, \code{law}, \code{frac}, \code{mu_x},
#'   \code{qx}, \code{px}, \code{lx}, \code{dx}, \code{Lx}, \code{Tx},
#'   \code{ex}, and \code{mx}.
#'
#' @examples
#' mortality_law_table("Exponential", 0, 110, lambda = 0.01)
#'
#' mortality_law_table("Gompertz", 0, 110, B = 1e-5, c = 1.08)
#'
#' mortality_law_table("Makeham", 0, 110, A = 5e-4, B = 1e-6, c = 1.10)
#'
#' mortality_law_table("Weibull", 1, 110, k = 2.5, lambda = 90)
#'
#' mortality_law_table("Logistic", 0, 110, A = 1e-4, B = 1e-6, c = 1.10, C = 1e-3)
#'
#' mortality_law_table("DeMoivre", 0, 100, omega = 100)
#'
#' mortality_law_table("Beta", 0, 100, alpha = 2, beta = 5, omega = 101)
#'
#' mortality_law_table(
#'   "HeligmanPollard", 1, 110,
#'   A = 0.0002, B = 0.1, C = 0.03, D = 10, E = 20, F = 0.00005, G = 1.08
#' )
#'
#' @export
mortality_law_table <- function(
    law = c(
      "Exponential", "Gompertz", "Makeham", "Weibull", "Logistic",
      "DeMoivre", "Beta", "HeligmanPollard"
    ),
    x_min,
    x_max,
    ...,
    params = NULL,
    frac = c("CF", "UDD", "Balducci"),
    radix = 1e5,
    close = TRUE,
    ax = 0.5,
    check = TRUE,
    tol = 1e-10
) {
  law <- match.arg(law)
  frac <- match.arg(frac)

  dots <- list(...)

  if (!is.null(params)) {
    if (!is.list(params) || is.null(names(params)) || any(names(params) == "")) {
      stop("`params` must be a named list when provided.", call. = FALSE)
    }
  } else {
    params <- list()
  }

  if (length(dots) > 0 && (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("All parameters in `...` must be named.", call. = FALSE)
  }

  pars <- utils::modifyList(params, dots)

  if (check) {
    if (!is.numeric(x_min) || !is.numeric(x_max) ||
        length(x_min) != 1L || length(x_max) != 1L) {
      stop("`x_min` and `x_max` must be single numeric values.", call. = FALSE)
    }

    if (!is.finite(x_min) || !is.finite(x_max)) {
      stop("`x_min` and `x_max` must be finite.", call. = FALSE)
    }

    if (abs(x_min - round(x_min)) > tol || abs(x_max - round(x_max)) > tol) {
      stop("`x_min` and `x_max` must be integer-valued.", call. = FALSE)
    }

    if (x_min >= x_max) {
      stop("Require `x_min < x_max`.", call. = FALSE)
    }

    if (!is.numeric(radix) || length(radix) != 1L ||
        !is.finite(radix) || radix <= 0) {
      stop("`radix` must be a single positive number.", call. = FALSE)
    }

    if (!is.logical(close) || length(close) != 1L || is.na(close)) {
      stop("`close` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.numeric(ax) || length(ax) != 1L ||
        !is.finite(ax) || ax < 0 || ax > 1) {
      stop("`ax` must be a single number in [0,1].", call. = FALSE)
    }

    if (!is.logical(check) || length(check) != 1L || is.na(check)) {
      stop("`check` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.numeric(tol) || length(tol) != 1L ||
        !is.finite(tol) || tol < 0) {
      stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
    }
  }

  x_min <- as.integer(round(x_min))
  x_max <- as.integer(round(x_max))
  x <- x_min:x_max
  n <- length(x)

  require_names <- function(required) {
    missing <- setdiff(required, names(pars))
    if (length(missing) > 0L) {
      stop(
        "Missing parameter(s) for law = '", law, "': ",
        paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  get_number <- function(name) {
    value <- as.numeric(pars[[name]])
    if (length(value) != 1L) {
      stop("Parameter `", name, "` must be a single numeric value.", call. = FALSE)
    }
    value
  }

  mu_x <- rep(NA_real_, n)
  qx <- rep(NA_real_, n)

  mu_laws <- c("Exponential", "Gompertz", "Makeham", "Weibull", "Logistic")

  if (law == "Exponential") {
    require_names("lambda")

    lambda <- get_number("lambda")

    if (check && (!is.finite(lambda) || lambda <= 0)) {
      stop("Exponential: require lambda > 0.", call. = FALSE)
    }

    mu_x <- rep(lambda, n)
  }

  if (law == "Gompertz") {
    require_names(c("B", "c"))

    B <- get_number("B")
    c0 <- get_number("c")

    if (check) {
      if (!is.finite(B) || B <= 0) stop("Gompertz: require B > 0.", call. = FALSE)
      if (!is.finite(c0) || c0 <= 0) stop("Gompertz: require c > 0.", call. = FALSE)
    }

    mu_x <- B * (c0 ^ x)
  }

  if (law == "Makeham") {
    require_names(c("A", "B", "c"))

    A <- get_number("A")
    B <- get_number("B")
    c0 <- get_number("c")

    if (check) {
      if (!is.finite(A) || A < 0) stop("Makeham: require A >= 0.", call. = FALSE)
      if (!is.finite(B) || B < 0) stop("Makeham: require B >= 0.", call. = FALSE)
      if (!is.finite(c0) || c0 <= 0) stop("Makeham: require c > 0.", call. = FALSE)
      if (A + B == 0) stop("Makeham: at least one of A or B must be positive.", call. = FALSE)
    }

    mu_x <- A + B * (c0 ^ x)
  }

  if (law == "Weibull") {
    require_names(c("k", "lambda"))

    k <- get_number("k")
    lambda <- get_number("lambda")

    if (check) {
      if (!is.finite(k) || k <= 0) stop("Weibull: require k > 0.", call. = FALSE)
      if (!is.finite(lambda) || lambda <= 0) stop("Weibull: require lambda > 0.", call. = FALSE)
      if (any(x < 0)) stop("Weibull: ages must be nonnegative.", call. = FALSE)
      if (any(x == 0L) && k < 1) {
        stop("Weibull: x = 0 with k < 1 gives infinite mu_x; use x_min >= 1.", call. = FALSE)
      }
    }

    mu_x <- (k / lambda) * ((x / lambda) ^ (k - 1))
  }

  if (law == "Logistic") {
    require_names(c("A", "B", "c", "C"))

    A <- get_number("A")
    B <- get_number("B")
    c0 <- get_number("c")
    C <- get_number("C")

    if (check) {
      if (!is.finite(A) || A < 0) stop("Logistic: require A >= 0.", call. = FALSE)
      if (!is.finite(B) || B < 0) stop("Logistic: require B >= 0.", call. = FALSE)
      if (!is.finite(C) || C < 0) stop("Logistic: require C >= 0.", call. = FALSE)
      if (!is.finite(c0) || c0 <= 0) stop("Logistic: require c > 0.", call. = FALSE)
      if (A + B == 0) stop("Logistic: at least one of A or B must be positive.", call. = FALSE)
    }

    cx <- c0 ^ x
    denominator <- 1 + C * cx

    if (check && any(denominator <= 0)) {
      stop("Logistic: denominator must be positive for all ages.", call. = FALSE)
    }

    mu_x <- (A + B * cx) / denominator
  }

  if (law == "DeMoivre") {
    require_names("omega")

    omega <- as.integer(round(get_number("omega")))

    if (check) {
      if (!is.finite(omega)) stop("DeMoivre: omega must be finite.", call. = FALSE)
      if (omega < x_max) stop("DeMoivre: require omega >= x_max.", call. = FALSE)
      if (omega <= x_min) stop("DeMoivre: require omega > x_min.", call. = FALSE)
    }

    qx <- ifelse(x < omega, 1 / (omega - x), 1)
    mu_x <- NA_real_
  }

  if (law == "Beta") {
    require_names(c("alpha", "beta", "omega"))

    alpha <- get_number("alpha")
    beta <- get_number("beta")
    omega <- get_number("omega")

    if (check) {
      if (!is.finite(alpha) || alpha <= 0) stop("Beta: require alpha > 0.", call. = FALSE)
      if (!is.finite(beta) || beta <= 0) stop("Beta: require beta > 0.", call. = FALSE)
      if (!is.finite(omega) || omega <= 0) stop("Beta: require omega > 0.", call. = FALSE)
      if (any(x < 0)) stop("Beta: ages must be nonnegative.", call. = FALSE)
      if (omega < x_max + 1) {
        stop("Beta: require omega >= x_max + 1 to compute qx from S(x+1).", call. = FALSE)
      }
    }

    survival_beta <- function(age) {
      u <- age / omega
      u <- pmin(pmax(u, 0), 1)
      1 - stats::pbeta(u, shape1 = alpha, shape2 = beta)
    }

    sx <- survival_beta(x)
    sx1 <- survival_beta(x + 1)

    qx <- ifelse(sx <= .Machine$double.eps, 1, (sx - sx1) / sx)
    mu_x <- NA_real_
  }

  if (law == "HeligmanPollard") {
    require_names(c("A", "B", "C", "D", "E", "F", "G"))

    A <- get_number("A")
    B <- get_number("B")
    C <- get_number("C")
    D <- get_number("D")
    E <- get_number("E")
    F <- get_number("F")
    G <- get_number("G")

    if (check) {
      if (x_min < 1L) {
        stop("HeligmanPollard: require x_min >= 1 due to log(x).", call. = FALSE)
      }

      if (!is.finite(B)) stop("HeligmanPollard: B must be finite.", call. = FALSE)

      for (nm in c("A", "C", "D", "E", "F")) {
        value <- get(nm)
        if (!is.finite(value) || value < 0) {
          stop("HeligmanPollard: parameters A, C, D, E, and F must be nonnegative.", call. = FALSE)
        }
      }

      if (!is.finite(G) || G <= 0) {
        stop("HeligmanPollard: require G > 0.", call. = FALSE)
      }

      if (A == 0 && C == 0 && F == 0) {
        stop("HeligmanPollard: at least one of A, C, or F must be positive.", call. = FALSE)
      }
    }

    xx <- as.numeric(x)

    odds <- (A ^ (xx + B)) +
      (C * exp(-D * (log(xx) - log(E))^2)) +
      (F * (G ^ xx))

    qx <- odds / (1 + odds)
    mu_x <- NA_real_
  }

  if (law %in% mu_laws) {
    if (check && any(!is.finite(mu_x) | mu_x < 0)) {
      stop("Computed mu_x must be finite and nonnegative.", call. = FALSE)
    }

    qx <- if (frac == "CF") {
      1 - exp(-mu_x)
    } else if (frac == "UDD") {
      mu_x
    } else {
      mu_x / (1 + mu_x)
    }
  }

  if (isTRUE(close)) {
    qx[n] <- 1
  }

  if (check && any(qx < -tol | qx > 1 + tol, na.rm = TRUE)) {
    stop("Derived qx is outside [0,1]. Check parameters and `frac`.", call. = FALSE)
  }

  qx <- pmin(pmax(qx, 0), 1)
  px <- 1 - qx

  lx <- radix * cumprod(c(1, px[1:(n - 1L)]))
  dx <- lx * qx

  lx_next <- c(lx[-1], lx[n] * px[n])

  Lx <- lx_next + ax * dx
  mx <- dx / pmax(Lx, .Machine$double.eps)

  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / pmax(lx, .Machine$double.eps)

  law_out <- law
  frac_out <- if (law_out %in% mu_laws) frac else NA_character_

  tibble::tibble(
    x = x,
    law = law_out,
    frac = frac_out,
    mu_x = mu_x,
    qx = qx,
    px = px,
    lx = lx,
    dx = dx,
    Lx = Lx,
    Tx = Tx,
    ex = ex,
    mx = mx
  )
}
