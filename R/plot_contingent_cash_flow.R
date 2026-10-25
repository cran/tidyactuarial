#' Plot discrete life-contingent cash flows
#'
#' Separates contractual amounts, unconditional payment probabilities, expected
#' cash flows and actuarial present-value contributions in editable panels.
#' Consumes the payment-level output of [apv_life_flow()] without recomputing
#' mortality, dates or interest, or an explicitly specified event table.
#'
#' @param .data A nonempty data frame or tibble. Supply an [apv_life_flow()]
#'   result, or columns `t`, `cf`, `probability`, `discount` and `event`.
#'   Optional `date` supplies the horizontal calendar axis. Optional
#'   `expected_cf` and `pv` are reused after consistency checks.
#' @param view One of `"all"`, `"contractual"`, `"probability"`, `"expected"`,
#'   or `"pv"`. The default shows four panels with separate vertical scales.
#' @param title,subtitle Optional plot title and subtitle.
#' @param x_label Optional horizontal axis label.
#' @param show_labels Logical. Show amount/probability labels, omitting text
#'   that overlaps at rendering time. All payment rows remain in the plot.
#' @param label_size Positive text size.
#' @param label_digits Integer between 0 and 50 for monetary labels.
#' @param label_format Monetary notation: `"auto"`, `"full"`, or `"compact"`,
#'   as in [plot_cash_flow()]. Probabilities/weights use separate unitless labels.
#' @param currency Monetary prefix; never applied to probabilities or weights.
#' @param col_inflow,col_outflow Colors for positive and negative monetary flows.
#' @param date_labels Calendar label format for [ggplot2::scale_x_date()].
#'
#' @details
#' For a binary payment event \eqn{E_j}, contractual signed amount \eqn{C_j},
#' payment date \eqn{t_j} and deterministic discount \eqn{d_j}, the panels show
#' \deqn{C_j,\quad p_j=P(E_j),\quad C_jp_j,\quad C_jp_jd_j.}
#' Probabilities are unconditional relative to the common initial actuarial
#' state, not conditional annual mortality rates. For survival payments,
#' \eqn{p_j={}_{t_j}p_x}. For death in \eqn{(a,b]},
#' \eqn{p_j={}_ap_x-{}_bp_x={}_ap_x\,{}_{b-a}q_{x+a}}. A benefit payable at
#' the end of that interval is discounted to its payment date, not the exact
#' time of death. Continuous death-time benefits require a different model.
#'
#' In explicit event tables, `event` contains `"survival"`, `"death_interval"`
#' or `"certain"` (the last requires probability one). Interval probabilities
#' must be supplied; the plot never differences arbitrary rows of a stream.
#' Rows may represent different events at the same date; probabilities are
#' never added or netted. Coincident stems can overlap; points identify their
#' endpoints and no payment date is moved. Negative amounts denote outflows.
#'
#' For new [apv_life_flow()] results, contingency metadata distinguishes binary
#' survival from a reversionary expected benefit fraction. For the latter,
#' `surv_prob` is \eqn{P(\text{all alive})+\alpha P(\text{some but not all alive})},
#' not generally a probability: the second panel is labeled **Expected benefit
#' fraction**. Results lacking metadata use the conservative label **Payment
#' weight**; no event is inferred from column names alone. Explicit tables can
#' use `event = "weighted"` for such finite multipliers as well.
#'
#' `life_contract()` stores a specification, and `insurance_x()`/`annuity_x()`
#' return aggregate values rather than payment schedules. These objects cannot
#' be plotted directly: supply an explicit payment table instead. The plot
#' does not expand contracts or reevaluate actuarial assumptions.
#'
#' All required numeric values must be finite, with nonnegative `t` and
#' `discount`, and probabilities in \eqn{[0,1]} for binary events. Missing or
#' inconsistent rows are rejected, never silently removed. Supplied expected
#' amounts and PVs must agree with the identities to relative tolerance 1e-10;
#' they are retained, not replaced. The plot ignores cached APV/cumulative
#' attributes and sums the displayed contributions. Calendar dates and stored
#' `t` are retained; in particular [apv_life_flow()] uses actual days/365.25.
#' All rows must share a valuation origin and monetary unit.
#'
#' @return An editable `ggplot2` object. Its tibble `data` has one row per
#'   payment per selected measure: `row`, `t`, `x`, `event`, `measure`, `value`,
#'   `plot_value` and `label`. `value` retains original units; extreme monetary
#'   panels use scaled `plot_value` with the power of ten in the panel title.
#' @seealso [plot_cash_flow()], [apv_life_flow()], [life_contract()],
#'   [insurance_x()], [annuity_x()]
#' @examples
#' lt <- tibble::tibble(x = 60:64, lx = c(1000, 800, 600, 300, 0))
#' apv_life_flow(lt, ages = 60, t = 1:3, cf = rep(1200, 3), i = 0.04) |>
#'   plot_contingent_cash_flow(currency = "$")
#'
#' # Three-year term insurance, paid at each year-end of death.
#' # P(j-1 < T <= j) = (l[x+j-1] - l[x+j]) / l[x], not q[x+j-1].
#' death_flows <- tibble::tibble(
#'   t = 1:3, cf = rep(10000, 3),
#'   probability = c(200, 200, 300) / 1000,
#'   discount = 1.04^-(1:3), event = "death_interval"
#' )
#' death_flows |> plot_contingent_cash_flow(currency = "$")
#' death_flows |> plot_contingent_cash_flow(view = "pv") +
#'   ggplot2::labs(title = "Term insurance: contributions to APV")
#' @md
#' @export
plot_contingent_cash_flow <- function(
    .data,
    view = c("all", "contractual", "probability", "expected", "pv"),
    title = NULL, subtitle = NULL, x_label = NULL,
    show_labels = TRUE, label_size = 3.2, label_digits = 2L,
    label_format = c("auto", "full", "compact"), currency = "",
    col_inflow = "#1B9E77", col_outflow = "#D95F02",
    date_labels = "%Y-%m-%d"
) {
  view <- match.arg(view)
  label_format <- match.arg(label_format)
  if (!is.data.frame(.data) || nrow(.data) == 0L) {
    stop("`.data` must be a nonempty payment-level data frame, not a contract or APV scalar.",
         call. = FALSE)
  }
  required <- c("t", "cf", "discount")
  if (!all(required %in% names(.data))) {
    stop("`.data` must contain payment-level `t`, `cf` and `discount` columns.", call. = FALSE)
  }
  if (!is.logical(show_labels) || length(show_labels) != 1L || is.na(show_labels)) {
    stop("`show_labels` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(label_size) || length(label_size) != 1L ||
      !is.finite(label_size) || label_size <= 0) {
    stop("`label_size` must be positive and finite.", call. = FALSE)
  }
  if (!is.numeric(label_digits) || length(label_digits) != 1L ||
      !is.finite(label_digits) || label_digits < 0 || label_digits > 50 ||
      label_digits != round(label_digits)) {
    stop("`label_digits` must be an integer between 0 and 50.", call. = FALSE)
  }
  validate_numeric <- function(value, name) {
    if (!is.numeric(value) || any(!is.finite(value))) {
      stop("`", name, "` must contain only finite numeric values.", call. = FALSE)
    }
  }
  for (name in required) validate_numeric(.data[[name]], name)
  if (any(.data$t < 0) || any(.data$discount < 0)) {
    stop("`t` and `discount` must be nonnegative.", call. = FALSE)
  }
  if ("probability" %in% names(.data)) {
    if (!"event" %in% names(.data)) {
      stop("Explicit `probability` requires an `event` column.", call. = FALSE)
    }
    weight <- .data$probability
    event <- as.character(.data$event)
  } else if ("surv_prob" %in% names(.data)) {
    weight <- .data$surv_prob
    kind <- attr(.data, "contingency", exact = TRUE)
    if (is.null(kind)) kind <- "weighted"
    event <- rep(kind, nrow(.data))
  } else {
    stop("Supply `probability` and `event`, or an `apv_life_flow()` result.", call. = FALSE)
  }
  validate_numeric(weight, "probability/payment weight")
  if (length(event) != nrow(.data) || anyNA(event) ||
      any(!event %in% c("survival", "death_interval", "certain", "weighted"))) {
    stop("`event` must be survival, death_interval, certain or weighted.", call. = FALSE)
  }
  binary <- event != "weighted"
  if (any(weight[binary] < 0 | weight[binary] > 1) ||
      any(weight[event == "certain"] != 1)) {
    stop("Binary event probabilities must be in [0,1]; certain payments require 1.", call. = FALSE)
  }
  # Retain supplied values, checking identities without a monetary absolute floor.
  expected <- .cash_flow_value(.data$cf, log(abs(weight))) * sign(weight)
  computed_pv <- .cash_flow_value(.data$cf,
    log(abs(weight)) + log(.data$discount)) * sign(weight)
  computed_pv[weight == 0 | .data$discount == 0 | .data$cf == 0] <- 0
  checked <- function(name, computed) {
    validate_numeric(computed, paste0("computed ", name))
    if (!name %in% names(.data)) return(computed)
    value <- .data[[name]]
    validate_numeric(value, name)
    scale <- pmax(abs(value), abs(computed))
    keep <- scale > 0
    if (any(abs(value[keep] / scale[keep] - computed[keep] / scale[keep]) > 1e-10)) {
      stop("`", name, "` is inconsistent with contractual amounts, weights and discounts.", call. = FALSE)
    }
    value
  }
  expected <- checked("expected_cf", expected)
  pv <- checked("pv", computed_pv)
  total <- sum(pv)
  if (!is.finite(total)) stop("Total APV must be finite.", call. = FALSE)
  dated <- "date" %in% names(.data)
  x <- .data$t
  if (dated) {
    if (!inherits(.data$date, "Date") || any(!is.finite(as.numeric(.data$date)))) {
      stop("`date` must contain finite Date values.", call. = FALSE)
    }
    x <- .data$date
  }
  weight_label <- if (any(event == "weighted")) {
    if (identical(attr(.data, "contingency", exact = TRUE), "weighted"))
      "Expected benefit fraction" else "Payment weight"
  } else "Unconditional payment probability"
  panel_labels <- c(contractual = "Contractual amount", probability = weight_label,
             expected = "Expected cash flow", pv = "APV contribution")
  values <- list(contractual = .data$cf, probability = weight, expected = expected, pv = pv)
  selected <- if (view == "all") names(values) else view
  rows <- seq_len(nrow(.data))
  payment_times <- .data$t
  panels <- lapply(selected, function(key) {
    value <- values[[key]]
    monetary <- key != "probability"
    display <- if (monetary) .cash_flow_display_scale(value) else list(divisor = 1, exponent = 0)
    panel <- panel_labels[[key]]
    if (display$exponent != 0) panel <- paste0(panel, " (x 10^", display$exponent, ")")
    tibble::tibble(row = rows, t = payment_times, x = x, event = event,
      measure = panel, value = value, plot_value = value / display$divisor,
      label = if (monetary) .cash_flow_label(value, label_digits, label_format, currency) else
        formatC(value, format = "g", digits = 3))
  })
  data <- dplyr::bind_rows(panels)
  data$measure <- factor(data$measure, levels = unique(data$measure))
  # Unitless weights have their own color; monetary signs retain family colors.
  colors <- c(inflow = col_inflow, outflow = col_outflow, probability = "#355C7D")
  weight_panel <- levels(data$measure)[grepl("probability|weight|fraction", levels(data$measure))]
  if (is.null(subtitle)) subtitle <- paste0("Total APV = ",
    .cash_flow_label(total, label_digits, label_format, currency))
  if (is.null(title)) title <- "Life-contingent cash flows"
  if (is.null(x_label)) x_label <- if (dated) "Payment date" else "Payment time (years)"
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[["x"]], y = .data[["plot_value"]])) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data[["x"]], y = 0,
      yend = .data[["plot_value"]], color = ifelse(.data[["measure"]] %in% weight_panel,
        "probability", ifelse(.data[["value"]] < 0, "outflow", "inflow"))), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(shape = .data[["event"]]), color = "#263746", size = 2) +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::scale_shape_manual(values = c(survival = 16, death_interval = 17,
      certain = 15, weighted = 18), labels = c(survival = "Survival",
      death_interval = "Death during interval", certain = "Certain",
      weighted = "Weighted benefit"), name = "Payment event") +
    ggplot2::facet_wrap(ggplot2::vars(.data[["measure"]]), scales = "free_y",
      ncol = if (view == "all") 2L else 1L) +
    ggplot2::scale_y_continuous(labels = function(x) .cash_flow_label(x, 2L, "auto"),
      expand = ggplot2::expansion(mult = c(.15, .23))) +
    ggplot2::labs(title = title, subtitle = subtitle, x = x_label, y = NULL,
      caption = "Discrete payments at stated dates. Panels use separate vertical scales.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(), legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      panel.spacing = grid::unit(1.3, "lines"))
  if (show_labels) p <- p + ggplot2::geom_text(ggplot2::aes(label = .data[["label"]],
    vjust = ifelse(.data[["plot_value"]] < 0, 1.6, -.7)),
    size = label_size, check_overlap = TRUE, color = "#263746")
  if (dated) {
    p <- p + ggplot2::scale_x_date(date_labels = date_labels,
      expand = ggplot2::expansion(mult = .12))
  } else {
    times <- sort(unique(data$x))
    p <- p + ggplot2::scale_x_continuous(breaks = if (length(times) <= 12L) times else
      scales::breaks_pretty(), expand = ggplot2::expansion(mult = .12))
  }
  p
}
