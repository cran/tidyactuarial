#' Kaplan--Meier survival curve and a lifetable-style life table
#'
#' Fits the nonparametric Kaplan--Meier estimator \eqn{\hat S(t)} for
#' right-censored time-to-event data, computes Greenwood's variance estimator
#' for \eqn{\hat S(t)}, and constructs a discrete life table by evaluating
#' \eqn{\hat S(t)} at user-provided cut points (\code{breaks}).
#'
#' The resulting life table is intended for *experience-based* (empirical) life
#' tables in actuarial/demographic contexts (e.g., cohort studies, population
#' indicators). It is not a replacement for graduated/regulatory tables when
#' smoothing, extrapolation, or product-specific selection effects are required.
#'
#' @param time Numeric vector. Observed times (event or censoring times).
#' @param status Integer/numeric vector of the same length as \code{time}.
#'   Use \code{1} for event (death), \code{0} for right-censoring.
#' @param entry Optional numeric vector of entry times (left truncation /
#'   delayed entry). If provided, must have the same length as \code{time} and
#'   satisfy \code{entry <= time}. If \code{NULL}, all individuals are assumed
#'   to enter at time 0.
#' @param breaks Optional numeric vector of increasing cut points used to build
#'   the discrete life table (e.g., \code{0:omega}). If \code{NULL}, defaults
#'   to integer ages from \code{0} to \code{ceiling(max(time))}.
#' @param radix Numeric. Life table radix used to scale \eqn{\ell_x}
#'   (default \code{1e5}).
#' @param conf_level Numeric in \code{(0,1)}. Confidence level for pointwise
#'   intervals for \eqn{\hat S(t)} computed via the log(-log) transformation
#'   (default \code{0.95}).
#' @param assumption Character. Fractional-age assumption used to compute
#'   \eqn{L_x} within each interval: \code{"UDD"}, \code{"CF"} (constant
#'   force), or \code{"Balducci"}.
#'
#' @details
#' **Kaplan--Meier estimator.** At each observed event time \eqn{t_j}:
#' \deqn{\hat S(t) = \prod_{t_j \le t} \left(1 - \frac{d_j}{n_j}\right)}
#' where \eqn{n_j} is the risk set size and \eqn{d_j} is the number of events.
#' Greenwood's variance:
#' \deqn{\widehat{\mathrm{Var}}(\hat S(t)) = \hat S(t)^2 \sum_{t_j \le t}
#' \frac{d_j}{n_j(n_j - d_j)}.}
#' Pointwise confidence intervals use the log(-log) transformation.
#'
#' **Life table mapping.** For each interval \eqn{[x, x + \Delta)}:
#' \deqn{\ell_x = \text{radix} \cdot \hat S(x), \quad d_x = \ell_x -
#'   \ell_{x+\Delta}, \quad q_x = d_x / \ell_x.}
#'
#' Exposure \eqn{L_x = \int_x^{x+\Delta} \ell(t)\,dt} is computed using the
#' selected fractional-age assumption (Finan, Section 24):
#' \itemize{
#'   \item UDD (Finan, Sec. 24.1):
#'     \eqn{L_x \approx \tfrac{\ell_x + \ell_{x+\Delta}}{2} \Delta}
#'   \item CF (constant force, Finan, Sec. 24.2):
#'     \eqn{L_x = \Delta \cdot (\ell_x - \ell_{x+\Delta}) /
#'     \ln(\ell_x / \ell_{x+\Delta})}
#'   \item Balducci (Finan, Sec. 24.3):
#'     \eqn{L_x = \Delta \cdot \ell_x \ell_{x+\Delta} /
#'     (\ell_x - \ell_{x+\Delta}) \cdot \ln(\ell_x / \ell_{x+\Delta})}
#' }
#'
#' Additional columns follow Finan, Sections 23.3, 23.8--23.9:
#' \itemize{
#'   \item \eqn{T_x = \sum_{k \ge x} L_k}: total expected years lived after
#'     age \eqn{x} (Finan, Sec. 23.3).
#'   \item \eqn{\mathring{e}_x = T_x / \ell_x}: complete life expectancy
#'     (Finan, Sec. 23.3).
#'   \item \eqn{m_x = d_x / L_x}: central death rate (Finan, Sec. 23.9).
#'   \item \eqn{a_x = (L_x - \ell_{x+\Delta}\Delta) / d_x}: mean time
#'     lived within the interval by those who die, in years (not a fraction
#'     when interval width differs from one). Undefined if there are no deaths.
#' }
#'
#' Nonfinite time/entry rows are omitted with a warning; the remaining sample
#' must be nonempty. Breaks are sorted and deduplicated. Beyond the last
#' observation, the KM step is held constant (no tail extrapolation is fitted).
#' CF/Balducci exposure at an interval ending at zero survival uses its zero
#' limit; positive deaths with zero exposure give infinite central mortality.
#' Confidence intervals remain undefined at survival zero or one.
#' Standard errors and intervals use Greenwood sums directly, preserving
#' small positive survival estimates without a fixed survival floor.
#' For delayed entry the risk set uses the existing inclusive convention
#' \code{entry <= time}; this is a deliberate convention at tied entries/events.
#' The lifetable component may have arbitrary interval widths: annual consumers
#' such as \code{t_px} require conversion to an annual grid with the needed
#' endpoint survivors; its class alone does not certify that conversion.
#'
#' @return A list with two tibbles:
#' \itemize{
#'   \item \code{km}: tibble with columns \code{time}, \code{n_risk},
#'     \code{d}, \code{censored}, \code{S}, \code{varS}, \code{seS},
#'     \code{ci_low}, \code{ci_high}.
#'   \item \code{lifetable}: tibble with columns \code{x}, \code{x_next},
#'     \code{width}, \code{lx}, \code{dx}, \code{qx}, \code{px}, \code{mx},
#'     \code{ax}, \code{Lx}, \code{Tx}, \code{ex}. Carries class
#'     \code{"lifetable"} and standard attributes for compatibility with
#'     downstream functions.
#' }
#'
#' @seealso \code{\link{lifetable}} for building tables from known mortality
#'   inputs, \code{\link{plot_km}} for plotting the KM curve.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' trueT <- rexp(n, rate = 0.08)
#' censT <- rexp(n, rate = 0.04)
#' time <- pmin(trueT, censT)
#' status <- as.integer(trueT <= censT)
#'
#' out <- km_lifetable(time, status, breaks = 0:25, radix = 100000)
#' head(out$km)
#' head(out$lifetable)
#'
#' # Tidy pipeline: filter high-mortality intervals
#' out$lifetable |> dplyr::filter(qx > 0.05)
#'
#' # Compare UDD vs CF assumptions
#' udd <- km_lifetable(time, status, breaks = 0:20, assumption = "UDD")
#' cfm <- km_lifetable(time, status, breaks = 0:20, assumption = "CF")
#' c(ex_udd = udd$lifetable$ex[1], ex_cf = cfm$lifetable$ex[1])
#'
#' # Plot the KM curve with plot_km
#' plot_km(out$km, time_col = "time", surv_col = "S",
#'         lower_col = "ci_low", upper_col = "ci_high")
#'
#' @export
km_lifetable <- function(
    time,
    status,
    entry = NULL,
    breaks = NULL,
    radix = 1e5,
    conf_level = 0.95,
    assumption = c("UDD", "CF", "Balducci")
) {
  assumption <- match.arg(assumption)

  if (!length(time)) stop("`time` must be nonempty.", call. = FALSE)
  if (!is.numeric(radix) || length(radix) != 1L || !is.finite(radix) || radix <= 0)
    stop("`radix` must be positive and finite.", call. = FALSE)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("`conf_level` must be in (0, 1).", call. = FALSE)

  # --- checks ---
  if (length(time) != length(status)) stop("'time' and 'status' must have same length.")
  if (!is.numeric(time)) stop("'time' must be numeric.")
  if (!all(status %in% c(0, 1))) stop("'status' must be 0/1 (0=censor, 1=event).")
  if (any(time < 0, na.rm = TRUE)) stop("'time' must be nonnegative.")

  if (is.null(entry)) {
    entry <- rep(0, length(time))
  } else {
    if (length(entry) != length(time)) stop("'entry' must have same length as 'time'.")
    if (!is.numeric(entry)) stop("'entry' must be numeric.")
    if (any(entry < 0, na.rm = TRUE)) stop("'entry' must be nonnegative.")
  }

  ok <- is.finite(time) & is.finite(status) & is.finite(entry)
  if (any(!ok)) warning("Omitting ", sum(!ok), " rows with nonfinite time/status/entry.", call. = FALSE)
  time <- time[ok]; status <- status[ok]; entry <- entry[ok]
  if (!length(time)) stop("Complete observations must be nonempty.", call. = FALSE)
  if (any(entry > time)) stop("Each entry time must satisfy entry <= time.")

  # --- sort by observed time ---
  o <- order(time)
  time <- time[o]; status <- status[o]; entry <- entry[o]

  # --- tabulate events and censoring at each unique time ---
  obs_df <- tibble::tibble(time = time, status = status, entry = entry)
  tab <- obs_df |>
    dplyr::group_by(time) |>
    dplyr::summarise(
      d        = sum(status == 1L),
      censored = sum(status == 0L),
      .groups  = "drop"
    ) |>
    dplyr::arrange(time)

  ut   <- tab$time
  d    <- tab$d
  cens <- tab$censored
  n_ut <- length(ut)

  # --- risk set with delayed entry ---
  # n_risk at time t: count of individuals with entry <= t and time >= t
  # Vectorized via outer comparison (efficient for moderate n_ut)
  if (n_ut > 0) {
    n_risk <- vapply(ut, function(t) sum(entry <= t & time >= t), numeric(1))
  } else {
    n_risk <- numeric(0)
  }

  # --- KM survival ---
  one_minus_h <- ifelse(n_risk > 0, 1 - d / n_risk, 1)
  S <- cumprod(one_minus_h)

  # --- Greenwood variance ---
  g_term <- ifelse(n_risk > d & n_risk > 0, d / (n_risk * (n_risk - d)), 0)
  g_cum  <- cumsum(g_term)
  varS   <- (S^2) * g_cum
  seS    <- S * sqrt(pmax(g_cum, 0))

  # --- log(-log) CI ---
  alpha <- 1 - conf_level
  z     <- stats::qnorm(1 - alpha / 2)
  interior_S <- S > 0 & S < 1
  ci_low <- ci_high <- rep(NA_real_, length(S))
  loglog <- log(-log(S[interior_S]))
  se_loglog <- sqrt(g_cum[interior_S]) / abs(log(S[interior_S]))
  ci_low[interior_S] <- exp(-exp(loglog + z * se_loglog))
  ci_high[interior_S] <- exp(-exp(loglog - z * se_loglog))

  km_tbl <- tibble::tibble(
    time     = ut,
    n_risk   = n_risk,
    d        = d,
    censored = cens,
    S        = S,
    varS     = varS,
    seS      = seS,
    ci_low   = ci_low,
    ci_high  = ci_high
  )

  # --- helper: evaluate right-continuous step S_hat(t) vectorized ---
  S_at <- function(t) {
    if (n_ut == 0L) return(rep(1, length(t)))
    idx <- findInterval(t, ut)
    c(1, S)[idx + 1L]
  }

  # --- breaks default ---
  if (is.null(breaks)) {
    mx_time <- max(time)
    breaks  <- 0:ceiling(mx_time)
  }
  if (!is.numeric(breaks) || any(!is.finite(breaks)) || any(breaks < 0))
    stop("`breaks` must be finite nonnegative numeric values.", call. = FALSE)
  breaks <- sort(unique(breaks))
  if (length(breaks) < 2) stop("'breaks' must have at least two points.")

  x      <- breaks[-length(breaks)]
  x_next <- breaks[-1]
  width  <- x_next - x

  Sx  <- S_at(x)
  Sx1 <- S_at(x_next)

  lx  <- radix * Sx
  lx1 <- radix * Sx1
  qx <- ifelse(Sx > 0, (Sx - Sx1) / Sx, NA_real_)
  px <- ifelse(Sx > 0, Sx1 / Sx, NA_real_)
  dx <- ifelse(lx == 0, 0, lx * qx)

  # --- Lx under chosen fractional-age assumption (Finan, Section 24) ---
  # Integrate normalized survival, preserving the scale of the cohort.
  factor <- rep(NA_real_, length(px))
  same <- !is.na(px) & px == 1
  zero <- !is.na(px) & px == 0
  interior <- !is.na(px) & px > 0 & px < 1
  factor[same] <- 1
  factor[zero] <- if (assumption == "UDD") .5 else 0
  if (assumption == "UDD") factor[interior] <- 1 - qx[interior] / 2
  log_p <- ifelse(qx[interior] < .5, log1p(-qx[interior]), log(px[interior]))
  if (assumption == "CF") factor[interior] <- qx[interior] / (-log_p)
  if (assumption == "Balducci") factor[interior] <-
    px[interior] * (-log_p) / qx[interior]
  factor[lx == 0] <- 0
  Lx <- .cash_flow_value(lx, log(width) + log(factor))
  if (any(!is.finite(Lx))) stop("Exposure must be representable; rescale the radix.", call. = FALSE)

  # --- Tx, ex (Finan, Sections 23.3 and 23.8) ---
  Tx <- rev(cumsum(rev(Lx)))
  ex <- ifelse(lx > 0, Tx / lx, NA_real_)

  # --- extra empirical columns (Finan, Section 23.9) ---
  mx_col <- ifelse(lx > 0, dx / Lx, NA_real_)
  ax_col <- ifelse(dx > 0, width * (factor - px) / qx, NA_real_)
  # Stable small-decrement limits of the mean time lived by decedents.
  small <- !is.na(qx) & qx > 0 & qx < 1e-4
  if (assumption == "CF") {
    mu <- -log1p(-qx[small])
    ax_col[small] <- width[small] * (.5 - mu / 12 + mu^3 / 720)
  } else if (assumption == "Balducci") {
    q <- qx[small]
    ax_col[small] <- width[small] * (.5 - q/6 - q^2/12 - q^3/20 - q^4/30)
  }
  if (assumption == "UDD") ax_col[!is.na(qx) & qx > 0] <- width[!is.na(qx) & qx > 0] / 2

  life_tbl <- tibble::tibble(
    x      = x,
    x_next = x_next,
    width  = width,
    lx     = lx,
    dx     = dx,
    qx     = qx,
    px     = px,
    mx     = mx_col,
    ax     = ax_col,
    Lx     = Lx,
    Tx     = Tx,
    ex     = ex
  )

  # --- assign lifetable class and standard attributes ---
  class(life_tbl) <- c("lifetable", class(life_tbl))

  attr(life_tbl, "radix")      <- radix
  attr(life_tbl, "omega")      <- max(x)
  attr(life_tbl, "type")       <- "ultimate"
  attr(life_tbl, "frac")       <- assumption
  attr(life_tbl, "closed")     <- FALSE
  attr(life_tbl, "ax")         <- NA_real_
  attr(life_tbl, "source")     <- "Kaplan-Meier"
  attr(life_tbl, "conf_level") <- conf_level

  list(
    km        = km_tbl,
    lifetable = life_tbl
  )
}
