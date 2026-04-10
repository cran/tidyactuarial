#' Plot a Kaplan--Meier survival curve
#'
#' Creates a step-function plot of the Kaplan--Meier survival estimate
#' \eqn{\hat S(t)} with optional pointwise confidence bands. Designed to
#' work directly with the output of \code{\link{km_lifetable}}.
#'
#' @param km A data frame or tibble with at least columns for time and
#'   survival. Can also be the full list returned by \code{\link{km_lifetable}},
#'   in which case the \code{$km} component is extracted automatically.
#' @param time_col Character. Name of the time column. Default \code{"time"}.
#' @param surv_col Character. Name of the survival column. Default \code{"S"}
#'   (matching \code{\link{km_lifetable}} output).
#' @param lower_col Character. Name of the lower CI column. Default
#'   \code{"ci_low"} (matching \code{\link{km_lifetable}} output).
#' @param upper_col Character. Name of the upper CI column. Default
#'   \code{"ci_high"} (matching \code{\link{km_lifetable}} output).
#' @param conf_int Logical. If \code{TRUE} (default) and CI columns exist in
#'   \code{km}, plot a step-wise confidence ribbon.
#' @param title Optional character string for the plot title.
#'
#' @details
#' Both the survival curve and the confidence band are rendered as step
#' functions (using \code{\link[ggplot2]{geom_step}}), which is the correct
#' representation for the KM estimator - a right-continuous step function
#' that drops at each observed event time.
#'
#' The confidence band uses \code{geom_stepribbon} logic: the data is
#' internally expanded so that a ribbon-fill follows the step pattern rather
#' than interpolating linearly between event times.
#'
#' @return A \code{ggplot} object that can be further customised with
#'   additional ggplot2 layers.
#'
#' @seealso \code{\link{km_lifetable}} for fitting the KM estimator and
#'   building the empirical life table.
#'
#' @examples
#' set.seed(42)
#' n <- 150
#' time   <- rexp(n, rate = 0.05)
#' status <- rbinom(n, 1, prob = 0.7)
#'
#' # Fit KM and plot directly
#' out <- km_lifetable(time, status, breaks = 0:30)
#'
#' # Pass the full list - $km is extracted automatically
#' plot_km(out)
#'
#' # Or pass just the km tibble
#' plot_km(out$km)
#'
#' # Without confidence band
#' plot_km(out, conf_int = FALSE, title = "KM Survival Curve")
#'
#' # Customise with ggplot2 layers
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_km(out) +
#'     ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") +
#'     ggplot2::labs(subtitle = "Dashed line = median survival")
#' }
#'
#' @export
plot_km <- function(
    km,
    time_col  = "time",
    surv_col  = "S",
    lower_col = "ci_low",
    upper_col = "ci_high",
    conf_int  = TRUE,
    title     = NULL
) {
  # --- auto-extract $km from km_lifetable() list output ---
  if (is.list(km) && !is.data.frame(km) && "km" %in% names(km)) {
    km <- km[["km"]]
  }
  if (!is.data.frame(km)) stop("'km' must be a data.frame or tibble.")

  # --- validate required columns ---
  if (!time_col %in% names(km)) {
    stop("Column '", time_col, "' not found in 'km'.")
  }
  if (!surv_col %in% names(km)) {
    stop("Column '", surv_col, "' not found in 'km'.")
  }

  # --- step-wise ribbon helper ---
  # Duplicates rows to create the stair pattern needed for a filled ribbon
  # that follows step-function jumps instead of interpolating linearly.
  make_step_ribbon_df <- function(df, t_col, lo_col, hi_col) {
    tt <- df[[t_col]]
    lo <- df[[lo_col]]
    hi <- df[[hi_col]]
    n  <- length(tt)
    if (n < 2) return(df)

    # For each interior point, insert a duplicate at the next time
    # with the *previous* lo/hi values (the left-continuous step)
    t_out  <- numeric(2 * n - 1)
    lo_out <- numeric(2 * n - 1)
    hi_out <- numeric(2 * n - 1)

    for (i in seq_len(n - 1)) {
      j <- 2 * i - 1
      t_out[j]      <- tt[i]
      lo_out[j]     <- lo[i]
      hi_out[j]     <- hi[i]
      t_out[j + 1]  <- tt[i + 1]
      lo_out[j + 1] <- lo[i]
      hi_out[j + 1] <- hi[i]
    }
    t_out[2 * n - 1]  <- tt[n]
    lo_out[2 * n - 1] <- lo[n]
    hi_out[2 * n - 1] <- hi[n]

    tibble::tibble(
      !!t_col  := t_out,
      !!lo_col := lo_out,
      !!hi_col := hi_out
    )
  }

  # --- base plot: step survival curve ---
  p <- ggplot2::ggplot(
    km,
    ggplot2::aes(x = .data[[time_col]], y = .data[[surv_col]])
  ) +
    ggplot2::geom_step(linewidth = 0.7) +
    ggplot2::labs(
      title = title,
      x     = "Time",
      y     = expression(hat(S)(t))
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme_minimal()

  # --- step-wise confidence ribbon ---
  has_ci <- all(c(lower_col, upper_col) %in% names(km))
  if (conf_int && has_ci) {
    ribbon_df <- make_step_ribbon_df(km, time_col, lower_col, upper_col)

    p <- p +
      ggplot2::geom_ribbon(
        data    = ribbon_df,
        mapping = ggplot2::aes(
          x    = .data[[time_col]],
          ymin = .data[[lower_col]],
          ymax = .data[[upper_col]]
        ),
        inherit.aes = FALSE,
        alpha       = 0.2
      )
  }

  p
}
