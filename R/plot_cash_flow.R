#' Plot a cash flow diagram
#'
#' Creates a visual cash-flow diagram with arrows representing inflows
#' (positive payments, upward) and outflows (negative payments, downward).
#'
#' @param payment Numeric vector of cash flows.
#' @param time Optional numeric vector of times.
#' @param date Optional vector of dates. If supplied, the earliest date is
#'   treated as time 0 for present-value calculations.
#' @param i Optional annual effective interest rate used to compute PV if
#'   \code{pv} is \code{NULL}.
#' @param pv Optional numeric present value to display.
#' @param title Optional character title for the plot.
#' @param financial Logical. If \code{TRUE} (default), uses financial
#'   convention: inflows point up and outflows point down. If \code{FALSE},
#'   all arrows point up with height proportional to absolute value.
#' @param col_inflow Character color for inflow arrows.
#' @param col_outflow Character color for outflow arrows.
#' @param day_count Day-count convention used when \code{date} is supplied.
#'   One of \code{"act/365"} or \code{"act/360"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' The vertical height of each arrow is proportional to the absolute
#' payment amount, normalized so the largest payment reaches unit height.
#' Each arrow is labeled with the formatted payment amount.
#'
#' If an interest rate \code{i} is supplied and \code{pv} is \code{NULL},
#' the present value is computed as
#' \eqn{PV = \sum_k C_k (1+i)^{-t_k}}{PV = sum(C_k * (1+i)^(-t_k))}
#' and displayed in the subtitle.
#'
#' @seealso \code{\link{pv_flow}}, \code{\link{fv_flow}},
#'   \code{\link{irr_flow}}
#'
#' @family time-value
#'
#' @examples
#' # Time-based diagram
#' plot_cash_flow(
#'   payment = c(-1000, 300, 400, 500),
#'   time = c(0, 1, 2, 3),
#'   i = 0.08
#' )
#'
#' # Date-based diagram
#' plot_cash_flow(
#'   payment = c(-1000, 300, 400, 500),
#'   date = as.Date(c("2026-01-01", "2027-01-01", "2028-01-01", "2029-01-01"))
#' )
#'
#' @export
plot_cash_flow <- function(
    payment,
    time = NULL,
    date = NULL,
    i = NULL,
    pv = NULL,
    title = NULL,
    financial = TRUE,
    col_inflow = "#1B9E77",
    col_outflow = "#D95F02",
    day_count = c("act/365", "act/360")
) {
  day_count <- match.arg(day_count)

  if (missing(payment)) {
    stop("`payment` must be provided.", call. = FALSE)
  }
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }
  if (length(payment) == 0L) {
    stop("`payment` must not be empty.", call. = FALSE)
  }
  if (any(is.na(payment)) || any(!is.finite(payment))) {
    stop("`payment` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.null(time) && !is.null(date)) {
    stop("Provide only one of `time` or `date`, not both.", call. = FALSE)
  }
  if (is.null(time) && is.null(date)) {
    stop("You must provide either `time` or `date`.", call. = FALSE)
  }

  if (!is.null(i)) {
    if (!is.numeric(i) || length(i) != 1L || is.na(i) || !is.finite(i) || i <= -1) {
      stop("`i` must be a single finite numeric value greater than -1.", call. = FALSE)
    }
  }

  if (!is.null(pv)) {
    if (!is.numeric(pv) || length(pv) != 1L || is.na(pv) || !is.finite(pv)) {
      stop("`pv` must be a single finite numeric value.", call. = FALSE)
    }
  }

  # ---- Build plotting data ----
  if (!is.null(time)) {
    if (!is.numeric(time)) {
      stop("`time` must be a numeric vector.", call. = FALSE)
    }
    if (length(time) != length(payment)) {
      stop("`time` and `payment` must have the same length.", call. = FALSE)
    }
    if (any(is.na(time)) || any(!is.finite(time))) {
      stop("`time` must contain only finite numeric values.", call. = FALSE)
    }

    x_lab <- "Time"
    t <- as.numeric(time)

    data <- data.frame(
      x = t,
      payment = as.numeric(payment),
      stringsAsFactors = FALSE
    )

  } else {
    date <- as.Date(date)
    if (any(is.na(date))) {
      stop("`date` must contain valid dates.", call. = FALSE)
    }
    if (length(date) != length(payment)) {
      stop("`date` and `payment` must have the same length.", call. = FALSE)
    }

    origin <- min(date)
    denom <- if (day_count == "act/365") 365 else 360
    t <- as.numeric(date - origin) / denom

    x_lab <- "Date"

    data <- data.frame(
      x = date,
      payment = as.numeric(payment),
      time = t,
      stringsAsFactors = FALSE
    )
  }

  # ---- Optional PV ----
  if (is.null(pv) && !is.null(i)) {
    pv <- sum(payment / (1 + i)^t)
  }

  max_abs <- max(abs(data$payment))
  data$height <- if (max_abs == 0) 0 else abs(data$payment) / max_abs

  # Financial convention: inflows up, outflows down
  # Non-financial: all arrows point up (absolute value)
  if (financial) {
    data$y1 <- ifelse(data$payment >= 0, data$height, -data$height)
  } else {
    data$y1 <- data$height
  }

  data$y0 <- 0
  offset <- 0.12
  data$ylab <- ifelse(data$y1 >= 0, data$y1 + offset, data$y1 - offset)

  data$type <- ifelse(data$payment >= 0, "inflow", "outflow")

  data$label <- vapply(
    data$payment,
    function(x) formatC(x, format = "f", digits = 2, big.mark = ","),
    character(1L)
  )

  subtitle_parts <- character(0)

  if (!is.null(i)) {
    subtitle_parts <- c(
      subtitle_parts,
      paste0("i = ", formatC(100 * i, format = "f", digits = 2), "%")
    )
  }

  if (!is.null(pv)) {
    subtitle_parts <- c(
      subtitle_parts,
      paste0("PV = ", formatC(pv, format = "f", digits = 2, big.mark = ","))
    )
  }

  subtitle <- if (length(subtitle_parts)) {
    paste(subtitle_parts, collapse = "   |   ")
  } else {
    NULL
  }

  if (is.null(title)) {
    title <- "Cash flow diagram"
  }

  # Use .data pronoun to avoid NSE NOTEs in R CMD check
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data$x,
        y = .data$y0,
        yend = .data$y1,
        color = .data$type
      ),
      linewidth = 1,
      arrow = ggplot2::arrow(length = grid::unit(0.17, "cm"))
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = .data$ylab, label = .data$label),
      vjust = 0.5
    ) +
    ggplot2::scale_color_manual(
      values = c(
        inflow = col_inflow,
        outflow = col_outflow
      ),
      guide = "none"
    ) +
    ggplot2::scale_y_continuous(limits = c(-1.6, 1.6), breaks = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_lab,
      y = NULL
    )

  if (!is.null(time)) {
    p <- p + ggplot2::scale_x_continuous(breaks = sort(unique(data$x)))
  } else {
    p <- p + ggplot2::scale_x_date(date_labels = "%Y-%m-%d")
  }

  p
}
