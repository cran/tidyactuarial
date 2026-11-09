#' Expected future lifetime for two independent lives
#'
#' Computes the expected future lifetime for two independent lives aged
#' \code{x} and \code{y}, for either joint-life (first death) or
#' last-survivor (second death).
#'
#' @param lt A life table data frame shared by both lives, or a list of two
#'   life tables \code{list(lt_x, lt_y)}, one for each life in that order.
#'   Each table must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for life 1.
#' @param y Integer actuarial age for life 2.
#' @param t Optional nonnegative numeric duration(s). If \code{NULL}, uses
#'   the maximum horizon allowed by the tables and selected status.
#' @param type Character: \code{"curtate"} or \code{"complete"}.
#' @param frac Fractional-age assumption for \code{type = "complete"},
#'   as in \code{\link{t_pxy}}: \code{"UDD"}, \code{"CF"},
#'   \code{"CML"}, or \code{"Balducci"}. If omitted, valid \code{frac}
#'   attributes are inherited from the tables, defaulting to \code{"UDD"}.
#'   Inherited \code{"CML"} is normalized to \code{"CF"}; conflicting
#'   attributes require an explicit \code{frac}, as in \code{t_pxy}.
#' @param cohort Legacy survival-state selector: \code{"first"} means
#'   joint-life and \code{"last"} means last survivor. Prefer \code{status}
#'   in new calls. Supply only one of \code{cohort} and \code{status}.
#' @param tidy Logical. If \code{TRUE}, returns a tibble.
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#'   Table consistency and scalar logical controls are always checked.
#' @param tol Numeric tolerance for integer checks.
#' @param status Survival status: \code{"joint"} (first death) or
#'   \code{"last"} (last survivor). Defaults to joint-life when neither
#'   selector is supplied. Appended to preserve existing positional calls.
#'
#' @details
#' \strong{Curtate expectation} (Finan, Section 56.4 / Section 57):
#' \deqn{e_{xy} = \sum_{k=1}^{\infty} {}_kp_{xy}, \quad
#'   e_{\overline{xy}} = \sum_{k=1}^{\infty} {}_kp_{\overline{xy}}.}
#'
#' \strong{Complete expectation} (Finan, Section 56.4):
#' \deqn{\mathring{e}_{xy} = \int_0^{\infty} {}_tp_{xy} \, dt.}
#'
#' The integral is decomposed year-by-year and evaluated analytically under
#' UDD, constant force, or Balducci, retaining small survival in log form.
#' Observed annual \code{lx} ratios determine the fractional survival.
#' Survival is always
#' measured from the initial ages, conditional on both lives being alive
#' at time zero. In particular, last-survivor survival at later times does
#' not imply that both lives remain alive.
#'
#' Each marginal horizon ends at its own table's last recorded age minus
#' its initial age. Joint-life expectation uses the smaller horizon and
#' last-survivor expectation the larger. With a shared table, these are
#' the horizons of the older and younger life, respectively.
#' Survival beyond each marginal horizon is taken as zero, as in
#' \code{t_px}; open tables therefore give truncated expectations.
#' Missing initial ages or
#' zero initial survivors give \code{NA}, including when \code{t = 0}.
#' Tables must have unique consecutive annual integer ages and finite,
#' nonnegative, nonincreasing \code{lx}. Observed \code{px}/\code{qx} must
#' be probabilities consistent with annual \code{lx} ratios within
#' \code{1e-10}; unknown probabilities may be \code{NA}. The age tolerance
#' \code{tol} does not relax this mortality consistency threshold.
#'
#' **Key identity** (Finan, Example 57.4):
#' \deqn{\mathring{e}_{\overline{xy}} = \mathring{e}_x + \mathring{e}_y
#'   - \mathring{e}_{xy}.}
#'
#' This can be used to cross-validate results.
#'
#' @return Numeric vector, or tibble if \code{tidy = TRUE}. The tibble
#'   retains its existing columns, including \code{cohort}: \code{"first"}
#'   for joint-life and \code{"last"} for last survivor, also when the
#'   selection is made using \code{status}.
#'
#' @seealso \code{\link{e_x}} for single-life expectancy,
#'   \code{\link{t_pxy}} for two-life survival probabilities,
#'   \code{\link{annuity_xy}} for two-life annuity APVs.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Curtate joint-life expectancy (Finan, Sec. 56.4)
#' e_xy(lt, x = 60, y = 62, type = "curtate", status = "joint")
#'
#' # Curtate last-survivor expectancy (Finan, Sec. 57)
#' e_xy(lt, x = 60, y = 62, type = "curtate", status = "last")
#'
#' # Verify identity (Finan, Example 57.4):
#' # e_{xy-bar} = e_x + e_y - e_xy
#' e_joint <- e_xy(lt, x = 60, y = 62, type = "curtate", status = "joint")
#' e_last  <- e_xy(lt, x = 60, y = 62, type = "curtate", status = "last")
#' e_x_val <- e_x(lt, x = 60, type = "curtate")
#' e_y_val <- e_x(lt, x = 62, type = "curtate")
#' c(last_surv = e_last, sum_minus_joint = e_x_val + e_y_val - e_joint)
#'
#' # Complete joint-life expectancy under UDD
#' e_xy(lt, x = 60, y = 62, type = "complete", frac = "UDD", status = "joint")
#'
#' # Temporary: 3-year curtate joint-life
#' e_xy(lt, x = 60, y = 62, t = 3, type = "curtate", status = "joint")
#'
#' # Tidy output
#' e_xy(lt, x = 60, y = 62, type = "curtate", status = "joint", tidy = TRUE)
#'
#' # Distinct mortality tables, applied to the corresponding lives
#' lt_x <- data.frame(x = 60:62, lx = c(100, 80, 0))
#' lt_y <- data.frame(x = 60:62, lx = c(100, 50, 0))
#' e_xy(list(lt_x, lt_y), 60, 60, status = "joint") # 0.4
#' e_xy(list(lt_x, lt_y), 60, 60, status = "last")  # 0.9
#'
#' @export
e_xy <- function(
    lt,
    x, y,
    t = NULL,
    type = c("curtate", "complete"),
    frac,
    cohort = c("first", "last"),
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10,
    status = c("joint", "last")
) {
  type   <- match.arg(type)
  if (!missing(status) && !missing(cohort))
    stop("Provide only one of `status` or legacy `cohort`.", call. = FALSE)
  if (missing(status)) {
    cohort <- match.arg(cohort)
    status <- if (cohort == "first") "joint" else "last"
  } else {
    status <- match.arg(status)
    cohort <- if (status == "joint") "first" else "last"
  }

  for (name in c("check", "tidy")) {
    value <- get(name, inherits = FALSE)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("`tol` must be a single finite nonnegative number.", call. = FALSE)

  if (missing(lt)) stop("`lt` is required.")
  if (is.data.frame(lt)) {
    lt_x <- lt_y <- lt
  } else if (is.list(lt) && length(lt) == 2L &&
             all(vapply(lt, is.data.frame, logical(1)))) {
    lt_x <- lt[[1L]]
    lt_y <- lt[[2L]]
  } else {
    stop("`lt` must be either one life table or a list of two life tables.", call. = FALSE)
  }
  .two_life_table(lt_x)
  if (!identical(lt_x, lt_y)) .two_life_table(lt_y)

  # Resolve inherited assumptions exactly as in t_pxy.
  if (missing(frac)) {
    frac_x <- attr(lt_x, "frac")
    frac_y <- attr(lt_y, "frac")
    if (identical(frac_x, "CML")) frac_x <- "CF"
    if (identical(frac_y, "CML")) frac_y <- "CF"
    ok_x <- !is.null(frac_x) && frac_x %in% c("UDD", "CF", "Balducci")
    ok_y <- !is.null(frac_y) && frac_y %in% c("UDD", "CF", "Balducci")
    if (ok_x && ok_y && !identical(frac_x, frac_y))
      stop("The two life tables carry different `frac` attributes. Supply `frac` explicitly.",
           call. = FALSE)
    if (ok_x) {
      frac <- frac_x
    } else if (ok_y) {
      frac <- frac_y
    } else {
      frac <- "UDD"
    }
  } else {
    frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac == "CML") frac <- "CF"
  }

  if (missing(x) || missing(y)) stop("`x` and `y` are required.")

  x <- as.numeric(x); y <- as.numeric(y)

  if (check) {
    if (any(!is.finite(x)) || any(!is.finite(y))) stop("`x` and `y` must be finite.")
    if (any(abs(x) > .Machine$integer.max) || any(abs(y) > .Machine$integer.max))
      stop("`x` and `y` must be representable as integer ages.", call. = FALSE)
    if (any(abs(x - round(x)) > tol) || any(abs(y - round(y)) > tol)) {
      stop("`x` and `y` must be integer ages.")
    }
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))

  # --- duration handling ---
  if (is.null(t)) {
    t <- rep(NA_real_, max(length(x), length(y)))
  } else {
    t <- as.numeric(t)
    if (check) {
      if (any(!is.finite(t))) stop("`t` must be finite when provided.")
      if (any(t < 0)) stop("`t` must be nonnegative.")
    }
  }

  # --- recycle x, y, t ---
  nx <- length(x); ny <- length(y); nt <- length(t)
  nout <- max(nx, ny, nt)
  ok <- (nx == 1L || nx == nout) && (ny == 1L || ny == nout) && (nt == 1L || nt == nout)
  if (!ok) stop("`x`, `y`, and `t` must have the same length or be length 1.")

  if (nx == 1L && nout > 1L) x <- rep(x, nout)
  if (ny == 1L && nout > 1L) y <- rep(y, nout)
  if (nt == 1L && nout > 1L) t <- rep(t, nout)

  omega_x <- max(as.double(lt_x$x))
  omega_y <- max(as.double(lt_y$x))

  ex <- rep(NA_real_, nout)

  for (j in seq_len(nout)) {
    xj <- x[j]; yj <- y[j]; tj <- t[j]

    initial_x <- match(xj, lt_x$x)
    initial_y <- match(yj, lt_y$x)
    if (is.na(initial_x) || is.na(initial_y) ||
        lt_x$lx[initial_x] <= 0 || lt_y$lx[initial_y] <= 0) {
      next
    }

    # Truncate each marginal at its own observed horizon.
    if (status == "joint") {
      max_k <- min(omega_x - xj, omega_y - yj)
    } else {
      max_k <- max(omega_x - xj, omega_y - yj)
    }
    if (max_k <= 0) {
      ex[j] <- 0
      next
    }

    if (is.na(tj)) {
      tj_use <- max_k
    } else {
      tj_use <- min(tj, max_k)
    }

    n_int <- floor(tj_use)
    basis <- .two_life_survival_basis(lt_x, lt_y, xj, yj)

    # --- curtate: e = sum_{k=1}^{n} k_p (Finan, Sec. 56.4) ---
    if (type == "curtate") {
      if (n_int <= 0) {
        ex[j] <- 0
      } else {
        pk <- vapply(1:n_int, function(k) {
          basis$survival(k, frac, status)
        }, numeric(1))
        ex[j] <- if (anyNA(pk)) NA_real_ else sum(pk)
      }
      next
    }

    # Integrate from the original ages, without renewing the two living lives.
    if (tj_use == 0) {
      ex[j] <- 0
      next
    }
    starts <- seq.int(0, ceiling(tj_use) - 1L)
    areas <- vapply(starts, function(k) {
      basis$area(k, min(1, tj_use-k), frac, status)
    }, numeric(1))
    ex[j] <- sum(areas)
  }

  if (isTRUE(tidy)) {
    return(tibble::tibble(
      x = x, y = y, t = t, type = type, frac = frac,
      cohort = cohort, ex = ex
    ))
  }
  ex
}
