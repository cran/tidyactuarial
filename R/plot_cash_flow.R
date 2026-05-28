#' Plot a cash-flow diagram
#'
#' Creates a professional cash-flow diagram with arrows representing inflows
#' and outflows over time.
#'
#' The function supports both classical vector input and tidyverse-style data
#' input through the pipe. It also supports either numeric times or calendar
#' dates.
#'
#' @param .data Optional data.frame or tibble containing cash-flow columns.
#'   If supplied, `payment`, `time`, and `date` can be passed as unquoted column
#'   names or as strings.
#' @param payment Numeric vector of cash flows, or a column name when `.data`
#'   is supplied.
#' @param time Optional numeric vector of times, or a column name when `.data`
#'   is supplied.
#' @param date Optional date vector, or a column name when `.data` is supplied.
#'   If supplied, the earliest date is treated as time 0 for present-value
#'   calculations.
#' @param rate Optional annual effective interest rate used to compute present
#'   value if `pv` is `NULL`.
#' @param i Deprecated alias for `rate`. Use `rate` in new code.
#' @param pv Optional numeric present value to display.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle. If `NULL`, a subtitle is built from
#'   `rate` and `pv` when available.
#' @param x_label Optional x-axis label.
#' @param amount_label Optional y-axis label when `normalize = FALSE`.
#' @param financial Logical. If `TRUE`, positive payments point upward and
#'   negative payments point downward. If `FALSE`, all arrows point upward with
#'   height proportional to absolute value.
#' @param normalize Logical. If `TRUE`, arrow heights are normalized by the
#'   largest absolute payment. If `FALSE`, the vertical axis uses real payment
#'   magnitudes.
#' @param aggregate Logical. If `TRUE`, cash flows occurring at the same time or
#'   date are summed before plotting.
#' @param show_labels Logical. If `TRUE`, labels are shown next to cash-flow
#'   arrows.
#' @param label_size Numeric text size for labels.
#' @param arrow_size Numeric line width for arrows.
#' @param timeline_size Numeric line width for the time axis.
#' @param label_digits Integer number of decimal digits for payment labels.
#' @param currency Optional currency or unit prefix, such as `"$"`.
#' @param col_inflow Character color for inflow arrows.
#' @param col_outflow Character color for outflow arrows.
#' @param date_labels Character date-label format passed to
#'   `ggplot2::scale_x_date()`.
#' @param day_count Day-count convention used when `date` is supplied.
#'   One of `"act/365"` or `"act/360"`.
#'
#' @return A `ggplot2` object.
#'
#' @details
#' By default, the plot uses real cash-flow magnitudes on the vertical axis.
#' This is preferable when the diagram is intended to communicate financial
#' scale. Use `normalize = TRUE` only when an intentionally schematic diagram
#' is desired.
#'
#' If an interest rate is supplied and `pv` is `NULL`, the present value is
#' computed as
#' \deqn{PV = \sum_k C_k (1+i)^{-t_k}}{PV = sum(C_k * (1+i)^(-t_k))}
#' and displayed in the subtitle.
#'
#' @seealso [pv_flow()], [fv_flow()], [irr_flow()]
#'
#' @family time-value
#'
#' @examples
#' plot_cash_flow(
#'   payment = c(-1000, 300, 400, 500),
#'   time = c(0, 1, 2, 3),
#'   rate = 0.08,
#'   currency = "$"
#' )
#'
#' cashflows <- tibble::tibble(
#'   time = c(0, 1, 2, 3),
#'   payment = c(-1000, 300, 400, 500)
#' )
#'
#' cashflows |>
#'   plot_cash_flow(payment = payment, time = time, rate = 0.08)
#'
#' dated_flows <- tibble::tibble(
#'   date = as.Date(c("2026-01-01", "2026-07-01", "2027-01-01")),
#'   payment = c(-1000, 450, 700)
#' )
#'
#' dated_flows |>
#'   plot_cash_flow(payment = payment, date = date, rate = 0.08)
#'
#' @export
plot_cash_flow <- function(
    .data = NULL,
    payment,
    time = NULL,
    date = NULL,
    rate = NULL,
    i = NULL,
    pv = NULL,
    title = NULL,
    subtitle = NULL,
    x_label = NULL,
    amount_label = "Cash flow",
    financial = TRUE,
    normalize = FALSE,
    aggregate = TRUE,
    show_labels = TRUE,
    label_size = 3.5,
    arrow_size = 0.8,
    timeline_size = 0.7,
    label_digits = 2L,
    currency = "",
    col_inflow = "#1B9E77",
    col_outflow = "#D95F02",
    date_labels = "%Y-%m-%d",
    day_count = c("act/365", "act/360")
) {
  day_count <- match.arg(day_count)

  if (!is.null(.data) && !inherits(.data, "data.frame")) {
    if (missing(payment)) {
      payment <- .data
      .data <- NULL
    } else {
      stop(
        "When `.data` is not a data.frame, `payment` must not also be supplied.",
        call. = FALSE
      )
    }
  }

  if (!is.null(i)) {
    if (!is.null(rate)) {
      stop("Use only one of `rate` or `i`, not both.", call. = FALSE)
    }
    rate <- i
  }

  if (!is.null(rate) &&
      (!is.numeric(rate) || length(rate) != 1L ||
       is.na(rate) || !is.finite(rate) || rate <= -1)) {
    stop("`rate` must be a single finite numeric value greater than -1.", call. = FALSE)
  }

  if (!is.null(pv) &&
      (!is.numeric(pv) || length(pv) != 1L ||
       is.na(pv) || !is.finite(pv))) {
    stop("`pv` must be a single finite numeric value.", call. = FALSE)
  }

  for (arg in c("financial", "normalize", "aggregate", "show_labels")) {
    val <- get(arg, inherits = FALSE)
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  for (arg in c("label_size", "arrow_size", "timeline_size")) {
    val <- get(arg, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1L ||
        is.na(val) || !is.finite(val) || val <= 0) {
      stop("`", arg, "` must be a single positive finite number.", call. = FALSE)
    }
  }

  if (!is.numeric(label_digits) || length(label_digits) != 1L ||
      is.na(label_digits) || !is.finite(label_digits) ||
      label_digits < 0 || abs(label_digits - round(label_digits)) > 1e-10) {
    stop("`label_digits` must be a single nonnegative integer.", call. = FALSE)
  }

  label_digits <- as.integer(round(label_digits))

  get_column_name <- function(expr) {
    if (is.character(expr) && length(expr) == 1L) {
      return(expr)
    }

    if (is.name(expr)) {
      return(as.character(expr))
    }

    NULL
  }

  format_amount <- function(x) {
    abs_x <- abs(x)
    formatted <- formatC(abs_x, format = "f", digits = label_digits, big.mark = ",")

    ifelse(
      x < 0,
      paste0("-", currency, formatted),
      paste0(currency, formatted)
    )
  }

  y_axis_labels <- function(x) {
    if (normalize) {
      return(formatC(x, format = "f", digits = 2))
    }

    format_amount(x)
  }

  if (is.null(.data)) {
    if (missing(payment)) {
      stop("`payment` must be provided.", call. = FALSE)
    }

    payment_vec <- payment
    time_vec <- time
    date_vec <- date
  } else {
    if (!inherits(.data, "data.frame")) {
      stop("`.data` must be a data.frame or tibble.", call. = FALSE)
    }

    data_in <- tibble::as_tibble(.data)

    if (missing(payment)) {
      stop("`payment` must be supplied when `.data` is used.", call. = FALSE)
    }

    payment_name <- get_column_name(substitute(payment))

    if (is.null(payment_name) || !payment_name %in% names(data_in)) {
      stop("`payment` must identify a column in `.data`.", call. = FALSE)
    }

    payment_vec <- data_in[[payment_name]]

    time_expr <- substitute(time)
    date_expr <- substitute(date)

    time_is_null <- identical(time_expr, quote(NULL))
    date_is_null <- identical(date_expr, quote(NULL))

    if (!time_is_null && !date_is_null) {
      stop("Provide only one of `time` or `date`, not both.", call. = FALSE)
    }

    if (time_is_null && date_is_null) {
      stop("You must provide either `time` or `date`.", call. = FALSE)
    }

    if (!time_is_null) {
      time_name <- get_column_name(time_expr)

      if (is.null(time_name) || !time_name %in% names(data_in)) {
        stop("`time` must identify a column in `.data`.", call. = FALSE)
      }

      time_vec <- data_in[[time_name]]
      date_vec <- NULL
    } else {
      date_name <- get_column_name(date_expr)

      if (is.null(date_name) || !date_name %in% names(data_in)) {
        stop("`date` must identify a column in `.data`.", call. = FALSE)
      }

      date_vec <- data_in[[date_name]]
      time_vec <- NULL
    }
  }

  if (!is.numeric(payment_vec)) {
    stop("`payment` must be numeric.", call. = FALSE)
  }

  if (length(payment_vec) == 0L) {
    stop("`payment` must not be empty.", call. = FALSE)
  }

  if (any(is.na(payment_vec)) || any(!is.finite(payment_vec))) {
    stop("`payment` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.null(time_vec) && !is.null(date_vec)) {
    stop("Provide only one of `time` or `date`, not both.", call. = FALSE)
  }

  if (is.null(time_vec) && is.null(date_vec)) {
    stop("You must provide either `time` or `date`.", call. = FALSE)
  }

  if (!is.null(time_vec)) {
    if (!is.numeric(time_vec)) {
      stop("`time` must be numeric.", call. = FALSE)
    }

    if (length(time_vec) != length(payment_vec)) {
      stop("`time` and `payment` must have the same length.", call. = FALSE)
    }

    if (any(is.na(time_vec)) || any(!is.finite(time_vec))) {
      stop("`time` must contain only finite numeric values.", call. = FALSE)
    }

    x_values <- as.numeric(time_vec)
    time_for_pv <- as.numeric(time_vec)
    x_is_date <- FALSE

    if (is.null(x_label)) {
      x_label <- "Time"
    }
  } else {
    date_vec <- as.Date(date_vec)

    if (length(date_vec) != length(payment_vec)) {
      stop("`date` and `payment` must have the same length.", call. = FALSE)
    }

    if (any(is.na(date_vec))) {
      stop("`date` must contain valid dates.", call. = FALSE)
    }

    origin <- min(date_vec)
    denom <- if (day_count == "act/365") 365 else 360

    x_values <- date_vec
    time_for_pv <- as.numeric(date_vec - origin) / denom
    x_is_date <- TRUE

    if (is.null(x_label)) {
      x_label <- "Date"
    }
  }

  if (is.null(pv) && !is.null(rate)) {
    pv <- sum(payment_vec / (1 + rate)^time_for_pv)
  }

  plot_data <- data.frame(
    x = x_values,
    time_for_pv = time_for_pv,
    payment = as.numeric(payment_vec)
  )

  if (isTRUE(aggregate)) {
    plot_data <- stats::aggregate(
      payment ~ x + time_for_pv,
      data = plot_data,
      FUN = sum
    )
  }

  plot_data <- plot_data[order(plot_data$x), , drop = FALSE]
  plot_data$type <- ifelse(plot_data$payment >= 0, "inflow", "outflow")

  max_abs_payment <- max(abs(plot_data$payment))

  if (max_abs_payment == 0) {
    plot_data$height <- 0
  } else if (normalize) {
    plot_data$height <- abs(plot_data$payment) / max_abs_payment
  } else {
    plot_data$height <- abs(plot_data$payment)
  }

  plot_data$yend <- if (financial) {
    ifelse(plot_data$payment >= 0, plot_data$height, -plot_data$height)
  } else {
    plot_data$height
  }

  plot_data$y <- 0

  y_range <- max(abs(plot_data$yend))
  if (y_range == 0) y_range <- 1

  label_offset <- 0.06 * y_range

  plot_data$label_y <- ifelse(
    plot_data$yend >= 0,
    plot_data$yend + label_offset,
    plot_data$yend - label_offset
  )

  plot_data$label <- format_amount(plot_data$payment)

  subtitle_parts <- character(0)

  if (!is.null(rate)) {
    subtitle_parts <- c(
      subtitle_parts,
      paste0("rate = ", formatC(100 * rate, format = "f", digits = 2), "%")
    )
  }

  if (!is.null(pv)) {
    subtitle_parts <- c(subtitle_parts, paste0("PV = ", format_amount(pv)))
  }

  if (is.null(subtitle) && length(subtitle_parts) > 0L) {
    subtitle <- paste(subtitle_parts, collapse = "   |   ")
  }

  if (is.null(title)) {
    title <- "Cash-flow diagram"
  }

  y_limit <- y_range * 1.25

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = timeline_size,
      color = "grey35"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data$x,
        y = .data$y,
        yend = .data$yend,
        color = .data$type
      ),
      linewidth = arrow_size,
      arrow = grid::arrow(
        length = grid::unit(0.18, "cm"),
        type = "closed"
      ),
      lineend = "round"
    ) +
    ggplot2::scale_color_manual(
      values = c(inflow = col_inflow, outflow = col_outflow),
      labels = c(inflow = "Inflow", outflow = "Outflow"),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-y_limit, y_limit),
      labels = y_axis_labels
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = if (normalize) "Relative magnitude" else amount_label
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8))
    )

  if (show_labels) {
    p <- p +
      ggplot2::geom_label(
        ggplot2::aes(
          y = .data$label_y,
          label = .data$label,
          fill = .data$type
        ),
        size = label_size,
        label.size = 0,
        alpha = 0.92,
        color = "white",
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = c(inflow = col_inflow, outflow = col_outflow),
        guide = "none"
      )
  }

  if (!x_is_date) {
    unique_x <- sort(unique(plot_data$x))

    if (length(unique_x) <= 15L) {
      p <- p +
        ggplot2::scale_x_continuous(
          breaks = unique_x,
          expand = ggplot2::expansion(mult = 0.05)
        )
    } else {
      p <- p +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = 0.05)
        )
    }
  } else {
    unique_dates <- sort(unique(plot_data$x))

    if (length(unique_dates) <= 8L) {
      p <- p +
        ggplot2::scale_x_date(
          breaks = unique_dates,
          date_labels = date_labels,
          expand = ggplot2::expansion(mult = 0.05)
        )
    } else {
      p <- p +
        ggplot2::scale_x_date(
          date_labels = date_labels,
          expand = ggplot2::expansion(mult = 0.05)
        )
    }
  }

  p
}


