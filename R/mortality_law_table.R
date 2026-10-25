#' Generate a tidy life table from a theoretical mortality law
#'
#' Creates a tidy life table with one row per integer age from a parametric
#' mortality law, using compact actuarial notation.
#'
#' The output follows tidyactuarial conventions and includes columns such as
#' \code{x}, \code{qx}, \code{px}, \code{lx}, \code{dx}, \code{Lx},
#' \code{Tx}, \code{ex}, and \code{mx}.
#'
#' @section Supported laws:
#' \itemize{
#'   \item \strong{Exponential}: \eqn{\mu_x = \lambda}.
#'   \item \strong{Gompertz}: \eqn{\mu_x = Bc^x}.
#'   \item \strong{Makeham}: \eqn{\mu_x = A + Bc^x}.
#'   \item \strong{Weibull}: \eqn{\mu_x = (shape/scale)(x/scale)^{shape-1}}.
#'   \item \strong{Logistic}: \eqn{\mu_x = (A + Bc^x)/(1 + Cc^x)}.
#'   \item \strong{DeMoivre}: finite lifetime law with
#'   \eqn{q_x = 1/(\omega - x)} for \eqn{x < \omega}.
#'   \item \strong{Beta}: scaled lifetime model
#'   \eqn{X/\omega \sim Beta(\alpha,\beta)}.
#'   \item \strong{HeligmanPollard}: eight-parameter odds model
#'   \eqn{q_x/p_x = A^{(x+B)^C} + D\exp(-E(\log x-\log F)^2) + GH^x},
#'   with \code{F_hp} denoting \eqn{F}.
#' }
#'
#' @section Law parameters:
#' Parameters for the selected law are supplied either through \code{...} or
#' through the named list \code{params}. Direct parameters supplied through
#' \code{...} override values in \code{params}.
#'
#' Required parameters:
#' \itemize{
#'   \item \code{"Exponential"}: \code{lambda}.
#'   \item \code{"Gompertz"}: \code{B}, \code{c}.
#'   \item \code{"Makeham"}: \code{A}, \code{B}, \code{c}.
#'   \item \code{"Weibull"}: \code{shape}, \code{scale}. Transitional aliases
#'   \code{k} and \code{lambda} are accepted.
#'   \item \code{"Logistic"}: \code{A}, \code{B}, \code{c}, \code{C}.
#'   \item \code{"DeMoivre"}: \code{omega}.
#'   \item \code{"Beta"}: \code{alpha}, \code{beta}, \code{omega}.
#'   \item \code{"HeligmanPollard"}: \code{A}, \code{B}, \code{C},
#'   \code{D}, \code{E}, \code{F_hp}, \code{G}, and \code{H}.
#'   \code{A}, \code{B}, \code{D}, and \code{G} must be nonnegative;
#'   \code{C}, \code{E}, \code{F_hp}, and \code{H} must be positive.
#'   At least one of \code{A}, \code{D}, or \code{G} must be positive.
#'   Ages must be at least 1. Transitional alias
#'   \code{F} is accepted and mapped to \code{F_hp}.
#' }
#'
#' @section Converting mu(x) to qx:
#' For laws defined by a force of mortality \eqn{\mu_x}, the one-year death
#' probability is \eqn{q_x = 1 - \exp(-\int_x^{x+1}\mu(s)\,ds)},
#' independently of \code{frac}. The column \code{mu_x} contains the force
#' evaluated at age \eqn{x}, not its annual integral. Exponential, Gompertz,
#' Makeham and Weibull use analytic integrals; Logistic uses numerical
#' integration with relative tolerance \code{1e-10}.
#' The \code{frac} label records a within-year interpolation convention;
#' it does not replace integration of the continuous law. Force products and
#' Logistic ratios use logarithmic evaluation to avoid intermediate power
#' overflow; Weibull annual increments use \code{log1p/expm1} instead of
#' subtracting nearly equal cumulative hazards.
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
#'   \code{"Gompertz"}, \code{"Makeham"}, \code{"Weibull"},
#'   \code{"Logistic"}, \code{"DeMoivre"}, \code{"Beta"}, or
#'   \code{"HeligmanPollard"}.
#' @param x_min Integer. Minimum age, inclusive.
#' @param x_max Integer. Maximum age, inclusive. Must satisfy
#'   \code{x_min < x_max}.
#' @param ... Named law parameters. These values override \code{params}.
#'   Placing \code{...} before arguments such as \code{close} and \code{check}
#'   avoids partial-matching conflicts with the Gompertz and Makeham parameter
#'   \code{c}.
#' @param params Named list of law parameters, or \code{NULL}.
#' @param frac Character. Within-year interpolation label, independent of the
#'   annual probability calculation. One of \code{"CF"}, \code{"UDD"}, \code{"Balducci"}, or
#'   \code{"CML"}. \code{"CML"} is treated as \code{"CF"}.
#' @param l0 Numeric scalar. Initial radix, that is, the starting value
#'   \eqn{l_{x_min}}.
#' @param close Logical. If \code{TRUE}, forces the last age to close.
#' @param a_x Numeric scalar. Average fraction of the year lived by those dying
#'   between age \eqn{x} and \eqn{x+1}. Default is \code{0.5}.
#'   Exposure uses \eqn{L_x = l_{x+1} + a_x d_x}; this is an approximation
#'   controlled by \code{a_x}, not an integration of continuous survival.
#' @param check Logical. If \code{TRUE}, performs strict input validation.
#' @param tol Numeric tolerance used in checks.
#' @param radix Deprecated. Use \code{l0}.
#' @param ax Deprecated. Use \code{a_x}.
#'
#' @return A tibble with columns \code{x}, \code{law}, \code{frac},
#'   \code{mu_x}, \code{qx}, \code{px}, \code{lx}, \code{dx}, \code{Lx},
#'   \code{Tx}, \code{ex}, and \code{mx}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{x} denotes age, \code{l0} denotes the starting
#' radix, \code{a_x} denotes the average fraction of the year lived by those
#' dying during the age interval, \code{qx} denotes the one-year death
#' probability, and \code{px = 1 - qx}.
#' Ratios use their actual denominators: \code{mx = dx / Lx} and
#' \code{ex = Tx / lx}. Positive deaths with zero exposure give infinite
#' \code{mx}; extinct cohorts with zero numerators and denominators give
#' \code{NaN} for these undefined ratios.
#'
#' Heligman-Pollard follows the eight-parameter odds specification in
#' Heligman and Pollard (1980), section 2.1, equation (1),
#' \doi{10.1017/S0020268100040257}. Earlier seven-parameter calls must be
#' respecified: the previous implementation omitted the childhood exponent
#' and assigned different meanings to \code{C} through \code{G}.
#'
#' The parameter \code{F_hp} is used for the Heligman-Pollard law instead of
#' \code{F}, because \code{F} is historically associated with \code{FALSE} in
#' R and should not be promoted in a CRAN-facing API. The old name \code{F} is
#' still accepted as a transitional alias.
#'
#' @seealso \code{\link{lifetable}}, \code{\link{t_px}}, \code{\link{e_xy}}
#'
#' @family life-tables
#'
#' @examples
#' mortality_law_table("Exponential", 0, 110, lambda = 0.01)
#'
#' mortality_law_table("Gompertz", 0, 110, B = 1e-5, c = 1.08)
#'
#' mortality_law_table("Makeham", 0, 110, A = 5e-4, B = 1e-6, c = 1.10)
#'
#' mortality_law_table("Weibull", 1, 110, shape = 2.5, scale = 90)
#'
#' mortality_law_table(
#'   "Logistic",
#'   0, 110,
#'   A = 1e-4,
#'   B = 1e-6,
#'   c = 1.10,
#'   C = 1e-3
#' )
#'
#' mortality_law_table("DeMoivre", 0, 100, omega = 100)
#'
#' mortality_law_table("Beta", 0, 100, alpha = 2, beta = 5, omega = 101)
#'
#' mortality_law_table(
#'   "HeligmanPollard",
#'   1, 110,
#'   A = 0.0002,
#'   B = 0.1,
#'   C = 0.3,
#'   D = 0.03,
#'   E = 10,
#'   F_hp = 20,
#'   G = 0.00005,
#'   H = 1.08
#' )
#'
#' # Transitional compatibility with older names
#' mortality_law_table("Weibull", 1, 110, k = 2.5, lambda = 90)
#'
#' mortality_law_table(
#'   "HeligmanPollard",
#'   1, 110,
#'   A = 0.0002,
#'   B = 0.1,
#'   C = 0.3,
#'   D = 0.03,
#'   E = 10,
#'   F = 20,
#'   G = 0.00005,
#'   H = 1.08
#' )
#'
#' @export
mortality_law_table <- function(
    law = c(
      "Exponential",
      "Gompertz",
      "Makeham",
      "Weibull",
      "Logistic",
      "DeMoivre",
      "Beta",
      "HeligmanPollard"
    ),
    x_min,
    x_max,
    ...,
    params = NULL,
    frac = c("CF", "UDD", "Balducci", "CML"),
    l0 = 1e5,
    close = TRUE,
    a_x = 0.5,
    check = TRUE,
    tol = 1e-10,
    radix = NULL,
    ax = NULL
) {
  law <- match.arg(law)
  frac <- match.arg(frac)

  if (frac == "CML") {
    frac <- "CF"
  }

  dots <- list(...)

  # -------------------------------------------------------------------------
  # Transitional compatibility with older top-level names
  # -------------------------------------------------------------------------

  if (!is.null(radix)) {
    if (!missing(l0)) {
      stop("Provide only one of `l0` or deprecated `radix`.", call. = FALSE)
    }

    l0 <- radix
  }

  if (!is.null(ax)) {
    if (!missing(a_x)) {
      stop("Provide only one of `a_x` or deprecated `ax`.", call. = FALSE)
    }

    a_x <- ax
  }

  # -------------------------------------------------------------------------
  # Parameter handling
  # -------------------------------------------------------------------------

  if (!is.null(params)) {
    if (!is.list(params) || is.null(names(params)) || any(names(params) == "")) {
      stop("`params` must be a named list when provided.", call. = FALSE)
    }
  } else {
    params <- list()
  }

  if (length(dots) > 0L &&
      (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("All parameters in `...` must be named.", call. = FALSE)
  }

  pars <- utils::modifyList(params, dots)

  resolve_alias <- function(old, new) {
    has_old <- old %in% names(pars)
    has_new <- new %in% names(pars)

    if (has_old && has_new) {
      stop(
        "Provide only one of `", new, "` or deprecated `", old, "`.",
        call. = FALSE
      )
    }

    if (has_old) {
      pars[[new]] <<- pars[[old]]
    }

    invisible(TRUE)
  }

  if (law == "Weibull") {
    resolve_alias("k", "shape")
    resolve_alias("lambda", "scale")
  }

  if (law == "HeligmanPollard") {
    resolve_alias("F", "F_hp")
  }

  # -------------------------------------------------------------------------
  # Basic validation
  # -------------------------------------------------------------------------

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(check)) {
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
      stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
    }
    if (!is.numeric(x_min) ||
        !is.numeric(x_max) ||
        length(x_min) != 1L ||
        length(x_max) != 1L) {
      stop("`x_min` and `x_max` must be single numeric values.", call. = FALSE)
    }

    if (!is.finite(x_min) || !is.finite(x_max)) {
      stop("`x_min` and `x_max` must be finite.", call. = FALSE)
    }

    if (abs(x_min - round(x_min)) > tol ||
        abs(x_max - round(x_max)) > tol) {
      stop("`x_min` and `x_max` must be integer-valued.", call. = FALSE)
    }

    if (abs(x_min) > .Machine$integer.max || abs(x_max) > .Machine$integer.max) {
      stop("Ages must be representable in the R integer range.", call. = FALSE)
    }
    if (x_min >= x_max) {
      stop("Require `x_min < x_max`.", call. = FALSE)
    }

    if (!is.numeric(l0) ||
        length(l0) != 1L ||
        !is.finite(l0) ||
        l0 <= 0) {
      stop("`l0` must be a single positive number.", call. = FALSE)
    }

    if (!is.logical(close) || length(close) != 1L || is.na(close)) {
      stop("`close` must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.numeric(a_x) ||
        length(a_x) != 1L ||
        !is.finite(a_x) ||
        a_x < 0 ||
        a_x > 1) {
      stop("`a_x` must be a single number in [0, 1].", call. = FALSE)
    }

  }

  x_min <- as.integer(round(x_min))
  x_max <- as.integer(round(x_max))
  x <- x_min:x_max
  n_age <- length(x)

  require_names <- function(required) {
    missing_params <- setdiff(required, names(pars))

    if (length(missing_params) > 0L) {
      stop(
        "Missing parameter(s) for law = '", law, "': ",
        paste(missing_params, collapse = ", "),
        ".",
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

  mu_x <- rep(NA_real_, n_age)
  qx <- rep(NA_real_, n_age)

  mu_laws <- c("Exponential", "Gompertz", "Makeham", "Weibull", "Logistic")

  # -------------------------------------------------------------------------
  # Force-of-mortality laws
  # -------------------------------------------------------------------------

  if (law == "Exponential") {
    require_names("lambda")

    lambda <- get_number("lambda")

    if (isTRUE(check) && (!is.finite(lambda) || lambda <= 0)) {
      stop("Exponential: require lambda > 0.", call. = FALSE)
    }

    mu_x <- rep(lambda, n_age)
    integrated_mu <- mu_x
  }

  if (law == "Gompertz") {
    require_names(c("B", "c"))

    B <- get_number("B")
    c0 <- get_number("c")

    if (isTRUE(check)) {
      if (!is.finite(B) || B <= 0) {
        stop("Gompertz: require B > 0.", call. = FALSE)
      }

      if (!is.finite(c0) || c0 <= 0) {
        stop("Gompertz: require c > 0.", call. = FALSE)
      }
    }

    log_c <- log(c0)
    log_mu <- log(B) + as.numeric(x) * log_c
    mu_x <- exp(log_mu)
    annual_factor <- if (log_c == 0) 1 else expm1(log_c) / log_c
    integrated_mu <- exp(log_mu + log(annual_factor))
  }

  if (law == "Makeham") {
    require_names(c("A", "B", "c"))

    A <- get_number("A")
    B <- get_number("B")
    c0 <- get_number("c")

    if (isTRUE(check)) {
      if (!is.finite(A) || A < 0) {
        stop("Makeham: require A >= 0.", call. = FALSE)
      }

      if (!is.finite(B) || B < 0) {
        stop("Makeham: require B >= 0.", call. = FALSE)
      }

      if (!is.finite(c0) || c0 <= 0) {
        stop("Makeham: require c > 0.", call. = FALSE)
      }

      if (A + B == 0) {
        stop("Makeham: at least one of A or B must be positive.", call. = FALSE)
      }
    }

    log_c <- log(c0)
    log_gompertz <- log(B) + as.numeric(x) * log_c
    mu_x <- A + exp(log_gompertz)
    annual_factor <- if (log_c == 0) 1 else expm1(log_c) / log_c
    integrated_mu <- A + exp(log_gompertz + log(annual_factor))
  }

  if (law == "Weibull") {
    require_names(c("shape", "scale"))

    shape <- get_number("shape")
    scale <- get_number("scale")

    if (isTRUE(check)) {
      if (!is.finite(shape) || shape <= 0) {
        stop("Weibull: require shape > 0.", call. = FALSE)
      }

      if (!is.finite(scale) || scale <= 0) {
        stop("Weibull: require scale > 0.", call. = FALSE)
      }

      if (any(x < 0)) {
        stop("Weibull: ages must be nonnegative.", call. = FALSE)
      }

      if (any(x == 0L) && shape < 1) {
        stop(
          "Weibull: x = 0 with shape < 1 gives infinite mu_x; use x_min >= 1.",
          call. = FALSE
        )
      }
    }

    xx <- as.numeric(x)
    mu_x <- if (shape == 1) rep(1 / scale, n_age) else
      exp(log(shape) - log(scale) + (shape - 1) * (log(xx) - log(scale)))
    integrated_mu <- numeric(n_age)
    positive <- xx > 0
    increment <- shape * log1p(1 / xx[positive])
    log_increment <- ifelse(increment > log(2),
      increment + log1p(-exp(-increment)), log(expm1(pmin(increment, log(2)))))
    integrated_mu[positive] <- exp(shape * (log(xx[positive]) - log(scale)) + log_increment)
    integrated_mu[!positive] <- exp(-shape * log(scale))
  }

  if (law == "Logistic") {
    require_names(c("A", "B", "c", "C"))

    A <- get_number("A")
    B <- get_number("B")
    c0 <- get_number("c")
    C <- get_number("C")

    if (isTRUE(check)) {
      if (!is.finite(A) || A < 0) {
        stop("Logistic: require A >= 0.", call. = FALSE)
      }

      if (!is.finite(B) || B < 0) {
        stop("Logistic: require B >= 0.", call. = FALSE)
      }

      if (!is.finite(C) || C < 0) {
        stop("Logistic: require C >= 0.", call. = FALSE)
      }

      if (!is.finite(c0) || c0 <= 0) {
        stop("Logistic: require c > 0.", call. = FALSE)
      }

      if (A + B == 0) {
        stop("Logistic: at least one of A or B must be positive.", call. = FALSE)
      }
    }

    # Stable log sums avoid overflowing numerator/denominator powers.
    log_add <- function(a, b) {
      peak <- pmax(a, b)
      result <- peak + log1p(exp(-abs(a - b)))
      result[peak == -Inf] <- -Inf
      result
    }
    force <- function(age) {
      z <- age * log(c0)
      if (C == 0) return(A + exp(log(B) + z))
      # Divide numerator and denominator by exp(z) when it is large.
      # This avoids subtracting nearly equal log sums of size age*log(c).
      positive <- pmax(z, 0)
      negative <- pmin(z, 0)
      exp(log_add(log(A) - positive, log(B) + negative) -
        log_add(-positive, log(C) + negative))
    }
    mu_x <- force(as.numeric(x))
    integrated_mu <- vapply(x, function(age) {
      stats::integrate(force, age, as.numeric(age) + 1, rel.tol = 1e-10)$value
    }, numeric(1))
  }

  # -------------------------------------------------------------------------
  # Direct qx laws
  # -------------------------------------------------------------------------

  if (law == "DeMoivre") {
    require_names("omega")

    omega <- as.integer(round(get_number("omega")))

    if (isTRUE(check)) {
      if (!is.finite(omega)) {
        stop("DeMoivre: omega must be finite.", call. = FALSE)
      }

      if (omega < x_max) {
        stop("DeMoivre: require omega >= x_max.", call. = FALSE)
      }

      if (omega <= x_min) {
        stop("DeMoivre: require omega > x_min.", call. = FALSE)
      }
    }

    qx <- ifelse(x < omega, 1 / (omega - x), 1)
    mu_x <- NA_real_
  }

  if (law == "Beta") {
    require_names(c("alpha", "beta", "omega"))

    alpha <- get_number("alpha")
    beta <- get_number("beta")
    omega <- get_number("omega")

    if (isTRUE(check)) {
      if (!is.finite(alpha) || alpha <= 0) {
        stop("Beta: require alpha > 0.", call. = FALSE)
      }

      if (!is.finite(beta) || beta <= 0) {
        stop("Beta: require beta > 0.", call. = FALSE)
      }

      if (!is.finite(omega) || omega <= 0) {
        stop("Beta: require omega > 0.", call. = FALSE)
      }

      if (any(x < 0)) {
        stop("Beta: ages must be nonnegative.", call. = FALSE)
      }

      if (omega < x_max + 1) {
        stop(
          "Beta: require omega >= x_max + 1 to compute qx from S(x + 1).",
          call. = FALSE
        )
      }
    }

    survival_beta <- function(age) {
      u <- age / omega
      u <- pmin(pmax(u, 0), 1)
      stats::pbeta(u, shape1 = alpha, shape2 = beta,
                   lower.tail = FALSE, log.p = TRUE)
    }

    sx <- survival_beta(x)
    sx1 <- survival_beta(x + 1)

    qx <- -expm1(sx1 - sx)
    mu_x <- NA_real_
  }

  if (law == "HeligmanPollard") {
    require_names(c("A", "B", "C", "D", "E", "F_hp", "G", "H"))

    A <- get_number("A")
    B <- get_number("B")
    C <- get_number("C")
    D <- get_number("D")
    E <- get_number("E")
    F_hp <- get_number("F_hp")
    G <- get_number("G")
    H <- get_number("H")

    if (isTRUE(check)) {
      if (x_min < 1L) {
        stop("HeligmanPollard: require x_min >= 1 due to log(x).", call. = FALSE)
      }

      for (nm in c("A", "B", "D", "G")) {
        value <- get(nm)

        if (!is.finite(value) || value < 0) {
          stop(
            "HeligmanPollard: parameters A, B, D, and G must be nonnegative.",
            call. = FALSE
          )
        }
      }

      for (nm in c("C", "E", "F_hp", "H")) {
        value <- get(nm)
        if (!is.finite(value) || value <= 0) {
          stop("HeligmanPollard: parameters C, E, F_hp, and H must be positive.",
               call. = FALSE)
        }
      }

      if (A == 0 && D == 0 && G == 0) {
        stop(
          "HeligmanPollard: at least one of A, D, or G must be positive.",
          call. = FALSE
        )
      }
    }

    xx <- as.numeric(x)

    # Sum the three odds components in log space to avoid Inf/Inf.
    log_child <- if (A == 1) rep(0, n_age) else log(A) * (xx + B)^C
    log_terms <- cbind(log_child,
                       log(D) - E * (log(xx) - log(F_hp))^2,
                       log(G) + xx * log(H))
    largest <- apply(log_terms, 1, max)
    log_odds <- largest + log(rowSums(exp(log_terms - largest)))
    log_odds[largest == -Inf] <- -Inf
    log_odds[largest == Inf] <- Inf
    qx <- stats::plogis(log_odds)
    mu_x <- NA_real_
  }

  # -------------------------------------------------------------------------
  # Convert force of mortality to qx
  # -------------------------------------------------------------------------

  if (law %in% mu_laws) {
    if (isTRUE(check) && any(!is.finite(mu_x) | mu_x < 0)) {
      stop("Computed mu_x must be finite and nonnegative.", call. = FALSE)
    }

    qx <- -expm1(-integrated_mu)
  }

  if (isTRUE(close)) {
    qx[n_age] <- 1
  }

  if (isTRUE(check) && any(qx < -tol | qx > 1 + tol, na.rm = TRUE)) {
    stop("Derived qx is outside [0, 1]. Check parameters and `frac`.",
         call. = FALSE)
  }

  qx <- pmin(pmax(qx, 0), 1)
  px <- 1 - qx

  lx <- l0 * cumprod(c(1, px[1:(n_age - 1L)]))
  dx <- lx * qx

  lx_next <- c(lx[-1], lx[n_age] * px[n_age])

  Lx <- lx_next + a_x * dx
  # Positive deaths with zero exposure give Inf; 0/0 is undefined (NaN).
  mx <- dx / Lx

  Tx <- rev(cumsum(rev(Lx)))
  ex <- Tx / lx

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
