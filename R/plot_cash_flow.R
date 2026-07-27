#' Plot a cash-flow diagram
#'
#' Creates a professional cash-flow diagram with arrows representing inflows
#' and outflows over time, using compact actuarial notation.
#'
#' @param .data Optional data.frame or tibble containing cash-flow columns.
#' @param C Numeric vector of cash flows, or a column name when `.data` is supplied.
#' @param t Optional numeric vector of times in years, or a column name when `.data` is supplied.
#' @param date Optional date vector, or a column name when `.data` is supplied.
#' @param i Optional annual interest-rate input used to compute present value if `PV` is `NULL`.
#' @param i_type Character string indicating the interest-rate type.
#' @param m Positive integer. Conversion frequency for nominal rates.
#' @param PV Optional numeric present value to display.
#' @param payment Deprecated. Use `C`.
#' @param time Deprecated. Use `t`. Kept explicit to avoid partial matching with `timeline_size`.
#' @param rate Deprecated. Use `i`.
#' @param pv Deprecated. Use `PV`.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @param x_label Optional x-axis label.
#' @param amount_label Optional y-axis label when `normalize = FALSE`.
#' @param financial Logical. If `TRUE`, positive cash flows point upward and negative cash flows point downward.
#' @param normalize Logical. If `TRUE`, arrow heights are normalized.
#' @param aggregate Logical. If `TRUE`, cash flows occurring at the same time or date are summed before plotting.
#' @param show_labels Logical. If `TRUE`, labels are shown next to cash-flow arrows.
#' @param label_size Numeric text size for labels.
#' @param arrow_size Numeric line width for arrows.
#' @param timeline_size Numeric line width for the time axis.
#' @param label_digits Integer number of decimal digits for cash-flow labels.
#' @param label_format Character string controlling label notation. Use
#'   \code{"auto"} to abbreviate values from one million onward,
#'   \code{"full"} for unabridged values, or \code{"compact"} for
#'   abbreviated values such as \code{1.25M}.
#' @param currency Optional currency or unit prefix, such as `$`.
#' @param col_inflow Character color for inflow arrows.
#' @param col_outflow Character color for outflow arrows.
#' @param date_labels Character date-label format passed to `ggplot2::scale_x_date()`.
#' @param day_count Day-count convention used when `date` is supplied.
#' @param ... Reserved for future extensions.
#'
#' @return A `ggplot2` object.
#'
#' @details
#' This function uses compact actuarial notation: `C` denotes cash-flow amounts,
#' `t` denotes times in years, `i` denotes the interest-rate input, and `PV`
#' denotes present value.
#'
#' @seealso `pv_flow`, `fv_flow`, `irr_flow`, `standardize_interest`
#'
#' @family time-value
#'
#' @examples
#' plot_cash_flow(
#'   C = c(-1000, 300, 400, 500),
#'   t = c(0, 1, 2, 3),
#'   i = 0.08,
#'   currency = "$"
#' )
#'
#' cashflows <- tibble::tibble(
#'   t = c(0, 1, 2, 3),
#'   C = c(-1000, 300, 400, 500)
#' )
#'
#' cashflows |>
#'   plot_cash_flow(C = C, t = t, i = 0.08)
#'
#' dated_flows <- tibble::tibble(
#'   date = as.Date(c("2026-01-01", "2026-07-01", "2027-01-01")),
#'   C = c(-1000, 450, 700)
#' )
#'
#' dated_flows |>
#'   plot_cash_flow(C = C, date = date, i = 0.08)
#'
#' plot_cash_flow(
#'   payment = c(-1000, 300, 400, 500),
#'   time = c(0, 1, 2, 3),
#'   rate = 0.08,
#'   currency = "$"
#' )
#'
#' @export
plot_cash_flow <- function(
    .data = NULL,
    C,
    t = NULL,
    date = NULL,
    i = NULL,
    i_type = "effective",
    m = 1L,
    PV = NULL,
    payment = NULL,
    time = NULL,
    rate = NULL,
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
    label_format = c("auto", "full", "compact"),
    currency = "",
    col_inflow = "#1B9E77",
    col_outflow = "#D95F02",
    date_labels = "%Y-%m-%d",
    day_count = c("act/365", "act/360"),
    ...
) {
  day_count <- match.arg(day_count)
  label_format <- match.arg(label_format)
  dots <- list(...)

  if (length(dots) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", names(dots)), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(rate)) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }
    i <- rate
  }

  if (!is.null(pv)) {
    if (!is.null(PV)) {
      stop("Provide only one of `PV` or deprecated `pv`.", call. = FALSE)
    }
    PV <- pv
  }

  if (!missing(C) && !is.null(payment)) {
    stop("Provide only one of `C` or deprecated `payment`.", call. = FALSE)
  }

  if (missing(C) && !is.null(payment)) {
    C <- payment
  }

  if (!missing(t) && !is.null(time)) {
    stop("Provide only one of `t` or deprecated `time`.", call. = FALSE)
  }

  if (missing(t) && !is.null(time)) {
    t <- time
  }

  C_expr <- if (missing(C)) quote(expr = ) else substitute(C)
  t_expr <- if (missing(t)) quote(NULL) else substitute(t)
  date_expr <- if (missing(date)) quote(NULL) else substitute(date)

  if (!is.null(.data) && !inherits(.data, "data.frame")) {
    if (missing(C)) {
      C <- .data
      .data <- NULL
      C_expr <- substitute(C)
    } else {
      stop(
        "When `.data` is not a data.frame, `C` must not also be supplied.",
        call. = FALSE
      )
    }
  }

  if (!is.null(i)) {
    if (!is.numeric(i) || length(i) != 1L || is.na(i) || !is.finite(i)) {
      stop("`i` must be a single finite numeric value.", call. = FALSE)
    }

    if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
      stop("`i_type` must be a single character string.", call. = FALSE)
    }

    valid_i_type <- c("effective", "nominal_interest", "nominal_discount", "force")

    if (!i_type %in% valid_i_type) {
      stop(
        "`i_type` must be one of: ",
        paste(sprintf("'%s'", valid_i_type), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
        m < 1 || abs(m - round(m)) > 1e-10) {
      stop("`m` must be a single positive integer.", call. = FALSE)
    }

    m <- as.integer(round(m))
    i_effective <- standardize_interest(i_type = i_type, i = i, m = m)

    if (!is.numeric(i_effective) || length(i_effective) != 1L ||
        is.na(i_effective) || !is.finite(i_effective) || i_effective <= -1) {
      stop(
        "The standardized annual effective interest rate must be greater than -1.",
        call. = FALSE
      )
    }
  } else {
    i_effective <- NULL
  }

  if (!is.null(PV) &&
      (!is.numeric(PV) || length(PV) != 1L || is.na(PV) || !is.finite(PV))) {
    stop("`PV` must be a single finite numeric value.", call. = FALSE)
  }

  for (arg in c("financial", "normalize", "aggregate", "show_labels")) {
    val <- get(arg, inherits = FALSE)
    if (!is.logical(val) || length(val) != 1L || is.na(val)) {
      stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }

  for (arg in c("label_size", "arrow_size", "timeline_size")) {
    val <- get(arg, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        !is.finite(val) || val <= 0) {
      stop("`", arg, "` must be a single positive finite number.", call. = FALSE)
    }
  }

  if (!is.numeric(label_digits) || length(label_digits) != 1L ||
      is.na(label_digits) || !is.finite(label_digits) || label_digits < 0 ||
      abs(label_digits - round(label_digits)) > 1e-10) {
    stop("`label_digits` must be a single nonnegative integer.", call. = FALSE)
  }

  label_digits <- as.integer(round(label_digits))

  get_column_name <- function(expr) {
    if (is.character(expr) && length(expr) == 1L) return(expr)
    if (is.name(expr)) return(as.character(expr))
    NULL
  }

  trim_trailing_zeros <- function(x) {
    x <- sub("(\\.[0-9]*?)0+$", "\\1", x)
    sub("\\.$", "", x)
  }

  format_compact <- function(x) {
    abs_x <- abs(x)
    units <- c("", "K", "M", "B", "T")

    exponent <- ifelse(
      abs_x == 0,
      0,
      pmin(4, pmax(0, floor(log10(abs_x) / 3)))
    )

    scaled <- abs_x / (1000^exponent)
    formatted <- formatC(
      scaled,
      format = "f",
      digits = label_digits,
      big.mark = ""
    )

    paste0(trim_trailing_zeros(formatted), units[exponent + 1L])
  }

  format_amount <- function(x) {
    use_compact <- label_format == "compact" ||
      (label_format == "auto" && max(abs(x), na.rm = TRUE) >= 1e6)

    formatted <- if (use_compact) {
      format_compact(x)
    } else {
      formatC(
        abs(x),
        format = "f",
        digits = label_digits,
        big.mark = ""
      )
    }

    ifelse(
      x < 0,
      paste0("-", currency, formatted),
      paste0(currency, formatted)
    )
  }

  y_axis_labels <- function(x) {
    if (normalize) return(formatC(x, format = "f", digits = 2))
    format_amount(x)
  }

  if (is.null(.data)) {
    if (missing(C)) stop("`C` must be provided.", call. = FALSE)
    C_vec <- C
    t_vec <- t
    date_vec <- date
  } else {
    data_in <- tibble::as_tibble(.data)

    if (missing(C)) {
      stop("`C` must be supplied when `.data` is used.", call. = FALSE)
    }

    C_name <- get_column_name(C_expr)
    if (is.null(C_name) || !C_name %in% names(data_in)) {
      stop("`C` must identify a column in `.data`.", call. = FALSE)
    }
    C_vec <- data_in[[C_name]]

    t_is_null <- identical(t_expr, quote(NULL))
    date_is_null <- identical(date_expr, quote(NULL))

    if (!t_is_null && !date_is_null) {
      stop("Provide only one of `t` or `date`, not both.", call. = FALSE)
    }
    if (t_is_null && date_is_null) {
      stop("You must provide either `t` or `date`.", call. = FALSE)
    }

    if (!t_is_null) {
      t_name <- get_column_name(t_expr)
      if (is.null(t_name) || !t_name %in% names(data_in)) {
        stop("`t` must identify a column in `.data`.", call. = FALSE)
      }
      t_vec <- data_in[[t_name]]
      date_vec <- NULL
    } else {
      date_name <- get_column_name(date_expr)
      if (is.null(date_name) || !date_name %in% names(data_in)) {
        stop("`date` must identify a column in `.data`.", call. = FALSE)
      }
      date_vec <- data_in[[date_name]]
      t_vec <- NULL
    }
  }

  if (!is.numeric(C_vec)) stop("`C` must be numeric.", call. = FALSE)
  if (length(C_vec) == 0L) stop("`C` must not be empty.", call. = FALSE)
  if (any(is.na(C_vec)) || any(!is.finite(C_vec))) {
    stop("`C` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.null(t_vec) && !is.null(date_vec)) {
    stop("Provide only one of `t` or `date`, not both.", call. = FALSE)
  }
  if (is.null(t_vec) && is.null(date_vec)) {
    stop("You must provide either `t` or `date`.", call. = FALSE)
  }

  if (!is.null(t_vec)) {
    if (!is.numeric(t_vec)) stop("`t` must be numeric.", call. = FALSE)
    if (length(t_vec) != length(C_vec)) {
      stop("`t` and `C` must have the same length.", call. = FALSE)
    }
    if (any(is.na(t_vec)) || any(!is.finite(t_vec))) {
      stop("`t` must contain only finite numeric values.", call. = FALSE)
    }
    x_values <- as.numeric(t_vec)
    time_for_pv <- as.numeric(t_vec)
    x_is_date <- FALSE
    if (is.null(x_label)) x_label <- "Time"
  } else {
    date_vec <- as.Date(date_vec)
    if (length(date_vec) != length(C_vec)) {
      stop("`date` and `C` must have the same length.", call. = FALSE)
    }
    if (any(is.na(date_vec))) stop("`date` must contain valid dates.", call. = FALSE)
    origin <- min(date_vec)
    denom <- if (day_count == "act/365") 365 else 360
    x_values <- date_vec
    time_for_pv <- as.numeric(date_vec - origin) / denom
    x_is_date <- TRUE
    if (is.null(x_label)) x_label <- "Date"
  }

  if (is.null(PV) && !is.null(i_effective)) {
    PV <- sum(C_vec / (1 + i_effective)^time_for_pv)
  }

  plot_data <- data.frame(x = x_values, time_for_pv = time_for_pv, C = as.numeric(C_vec))

  if (isTRUE(aggregate)) {
    plot_data <- stats::aggregate(C ~ x + time_for_pv, data = plot_data, FUN = sum)
  }

  plot_data <- plot_data[order(plot_data$x), , drop = FALSE]
  plot_data$type <- ifelse(plot_data$C >= 0, "inflow", "outflow")
  max_abs_C <- max(abs(plot_data$C))

  if (max_abs_C == 0) {
    plot_data$height <- 0
  } else if (normalize) {
    plot_data$height <- abs(plot_data$C) / max_abs_C
  } else {
    plot_data$height <- abs(plot_data$C)
  }

  plot_data$yend <- if (financial) {
    ifelse(plot_data$C >= 0, plot_data$height, -plot_data$height)
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
  plot_data$label <- format_amount(plot_data$C)

  x_min <- min(plot_data$x)
  x_max <- max(plot_data$x)

  plot_data$label_hjust <- dplyr::case_when(
    plot_data$x == x_min & plot_data$x != x_max ~ 0,
    plot_data$x == x_max & plot_data$x != x_min ~ 1,
    TRUE ~ 0.5
  )

  has_positive <- any(plot_data$yend > 0)
  has_negative <- any(plot_data$yend < 0)

  if (!has_positive && !has_negative) {
    y_lower <- -0.15
    y_upper <- 0.15
  } else {
    positive_top <- if (has_positive) {
      max(plot_data$label_y[plot_data$yend > 0])
    } else {
      0
    }

    negative_bottom <- if (has_negative) {
      min(plot_data$label_y[plot_data$yend < 0])
    } else {
      0
    }

    data_span <- positive_top - negative_bottom
    if (!is.finite(data_span) || data_span <= 0) data_span <- 1

    axis_padding <- 0.10 * data_span

    y_lower <- if (has_negative) {
      negative_bottom - axis_padding
    } else {
      0
    }

    y_upper <- if (has_positive) {
      positive_top + axis_padding
    } else {
      0
    }
  }

  show_flow_legend <- dplyr::n_distinct(plot_data$type) > 1L

  subtitle_parts <- character(0)
  if (!is.null(i)) {
    rate_label <- if (identical(i_type, "effective")) {
      paste0("i = ", formatC(100 * i, format = "f", digits = 2), "%")
    } else {
      paste0(
        "i = ", formatC(100 * i, format = "f", digits = 2), "% ",
        "(", i_type, ", m = ", m, "; ",
        "i_eff = ", formatC(100 * i_effective, format = "f", digits = 2), "%)"
      )
    }
    subtitle_parts <- c(subtitle_parts, rate_label)
  }

  if (!is.null(PV)) subtitle_parts <- c(subtitle_parts, paste0("PV = ", format_amount(PV)))
  if (is.null(subtitle) && length(subtitle_parts) > 0L) subtitle <- paste(subtitle_parts, collapse = "   |   ")
  if (is.null(title)) title <- "Cash-flow diagram"

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["x"]])) +
    ggplot2::geom_hline(yintercept = 0, linewidth = timeline_size, color = "grey35") +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data[["x"]],
        y = .data[["y"]],
        yend = .data[["yend"]],
        color = .data[["type"]]
      ),
      linewidth = arrow_size,
      arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
      lineend = "round"
    ) +
    ggplot2::scale_color_manual(
      values = c(inflow = col_inflow, outflow = col_outflow),
      labels = c(inflow = "Inflow", outflow = "Outflow"),
      name = NULL,
      drop = TRUE
    ) +
    ggplot2::scale_y_continuous(
      labels = y_axis_labels,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_cartesian(
      ylim = c(y_lower, y_upper),
      clip = "off"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = if (normalize) "Relative magnitude" else amount_label
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = if (show_flow_legend) "bottom" else "none",
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
          y = .data[["label_y"]],
          label = .data[["label"]],
          fill = .data[["type"]],
          hjust = .data[["label_hjust"]]
        ),
        size = label_size,
        label.size = 0,
        alpha = 0.92,
        color = "white",
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = c(inflow = col_inflow, outflow = col_outflow),
        guide = "none",
        drop = TRUE
      )
  }

  if (!x_is_date) {
    unique_x <- sort(unique(plot_data$x))
    if (length(unique_x) <= 15L) {
      p <- p + ggplot2::scale_x_continuous(breaks = unique_x, expand = ggplot2::expansion(mult = 0.08))
    } else {
      p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08))
    }
  } else {
    unique_dates <- sort(unique(plot_data$x))
    if (length(unique_dates) <= 8L) {
      p <- p + ggplot2::scale_x_date(breaks = unique_dates, date_labels = date_labels, expand = ggplot2::expansion(mult = 0.08))
    } else {
      p <- p + ggplot2::scale_x_date(date_labels = date_labels, expand = ggplot2::expansion(mult = 0.08))
    }
  }

  p
}
