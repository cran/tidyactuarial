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
#'   \item \eqn{a_x = (\ell_x \Delta - L_x) / d_x}: average fraction of the
#'     interval lived by those who die.
#' }
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
  time <- time[ok]; status <- status[ok]; entry <- entry[ok]
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
  seS    <- sqrt(pmax(varS, 0))

  # --- log(-log) CI ---
  alpha <- 1 - conf_level
  z     <- stats::qnorm(1 - alpha / 2)
  eps   <- 1e-16
  loglog    <- log(-log(pmax(S, eps)))
  se_loglog <- ifelse(
    S > 0 & S < 1,
    seS / (abs(log(pmax(S, eps))) * pmax(S, eps)),
    NA_real_
  )
  ci_low  <- exp(-exp(loglog + z * se_loglog))
  ci_high <- exp(-exp(loglog - z * se_loglog))

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
    ifelse(idx == 0L, 1, S[idx])
  }

  # --- breaks default ---
  if (is.null(breaks)) {
    mx_time <- max(time)
    breaks  <- 0:ceiling(mx_time)
  }
  breaks <- sort(unique(breaks))
  if (length(breaks) < 2) stop("'breaks' must have at least two points.")

  x      <- breaks[-length(breaks)]
  x_next <- breaks[-1]
  width  <- x_next - x

  Sx  <- S_at(x)
  Sx1 <- S_at(x_next)

  lx  <- radix * Sx
  lx1 <- radix * Sx1
  dx  <- pmax(lx - lx1, 0)

  qx <- ifelse(lx > 0, dx / lx, NA_real_)
  px <- 1 - qx

  # --- Lx under chosen fractional-age assumption (Finan, Section 24) ---
  same      <- abs(lx - lx1) < 1e-12
  ratio     <- pmax(lx / pmax(lx1, eps), eps)
  log_ratio <- log(ratio)

  Lx <- switch(
    assumption,
    "UDD"     = 0.5 * (lx + lx1) * width,
    "CF"      = ifelse(same, lx * width, width * (lx - lx1) / log_ratio),
    "Balducci" = ifelse(same, lx * width, width * (lx * lx1) / (lx - lx1) * log_ratio)
  )

  # --- Tx, ex (Finan, Sections 23.3 and 23.8) ---
  Tx <- rev(cumsum(rev(Lx)))
  ex <- ifelse(lx > 0, Tx / lx, NA_real_)

  # --- extra empirical columns (Finan, Section 23.9) ---
  mx_col <- ifelse(Lx > 0, dx / Lx, NA_real_)
  ax_col <- ifelse(dx > 0, (lx * width - Lx) / dx, NA_real_)

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
