#' Actuarial present value of a payment stream under mortality
#'
#' Computes the actuarial present value (APV) of a cash-flow stream contingent
#' on survival. The life table is supplied as the first argument
#' (pipe-friendly). Payments may be specified by numeric times or by calendar
#' dates.
#'
#' Multiple lives are supported under an independence assumption, through
#' common statuses: single-life, first-death (all alive), last-survivor
#' (any alive), and reversionary (joint-and-survivor) with fraction
#' \code{alpha}.
#'
#' @param lt A life table data frame with column \code{x} and at least one of
#'   \code{lx}, \code{px}, or \code{qx}.
#' @param ages Integer vector of actuarial ages. Use length 1 for a single life
#'   and length 2 or more for multiple lives.
#' @param t Numeric vector of payment times in years, measured from time 0.
#'   Provide either \code{t} or \code{date}.
#' @param date Optional vector of \code{Date} payment dates. Provide either
#'   \code{t} or \code{date}.
#' @param date0 Optional \code{Date} used as time 0 when \code{date} is
#'   provided. If missing, the minimum of \code{date} is used.
#' @param cf Numeric vector of cash flows. Must have the same length as
#'   \code{t} or \code{date}.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param status Survival status: \code{"single"}, \code{"first"},
#'   \code{"last"}, or \code{"reversionary"}.
#' @param alpha Reversionary fraction for \code{status = "reversionary"}.
#'   While all lives are alive, full benefit is paid; while at least one but
#'   not all are alive, \code{alpha} times the benefit is paid.
#' @param plot Logical. If \code{TRUE}, attaches a ggplot object in
#'   \code{attr(result, "plot")} showing cumulative APV over time.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{t} denotes payment time, \code{cf} denotes cash
#' flows, \code{i} denotes the interest rate, \code{i_type} denotes the
#' interest-rate type, and \code{m} denotes the conversion frequency for nominal
#' rates.
#'
#' For each payment at time \eqn{t}, the APV contribution is
#' \deqn{PV(t) = C(t) v^t P(\text{status alive at } t).}
#'
#' The survival probability depends on the \code{status}:
#' \itemize{
#'   \item \code{"single"}: \eqn{{}_t p_x} for a single life.
#'   \item \code{"first"}: \eqn{{}_t p_{x_1} \cdot {}_t p_{x_2} \cdots},
#'     so all lives must be alive.
#'   \item \code{"last"}: \eqn{1 - \prod_j (1 - {}_t p_{x_j})},
#'     so at least one life must be alive.
#'   \item \code{"reversionary"}: full benefit while all lives are alive, and
#'     fraction \eqn{\alpha} while at least one but not all lives are alive.
#' }
#'
#' Fractional-year survival is computed under UDD within each year.
#' Discounting uses the log accumulation rate from the original interest
#' specification. A displayed discount factor may be zero or infinite even
#' when its monetary APV contribution is representable.
#'
#' Unknown mortality intervals produce missing contributions, unless survival
#' has already reached zero. Missing numeric payment times and cash flows are
#' retained as missing rows. Cumulative APV follows input row order and remains
#' missing after a missing contribution. For compatibility, the `apv` attribute
#' is the sum of nonmissing contributions (zero if none are known), so it is a
#' subtotal when any contribution is missing, not a complete contract valuation.
#' Calendar times use actual elapsed days divided by 365.25; all dates must be
#' finite and nonmissing. Payments at the same time remain separate rows.
#'
#' @return A tibble with one row per payment and columns:
#'   \code{t}, \code{cf}, \code{surv_prob}, \code{discount},
#'   \code{expected_cf}, \code{pv}, and \code{pv_cum}. If \code{date} was
#'   provided, a \code{date} column is included. The total APV is stored as
#'   \code{attr(result, "apv")}.
#'   The `contingency` attribute identifies `"survival"` or `"weighted"`
#'   (reversionary benefit fraction) for [plot_contingent_cash_flow()].
#'
#' @examples
#' lt <- data.frame(
#'   x = 40:100,
#'   lx = seq(100000, 0, length.out = 61)
#' )
#'
#' apv_life_flow(
#'   lt = lt,
#'   ages = 40,
#'   t = c(1, 2, 3),
#'   cf = c(100, 100, 100),
#'   i = 0.05
#' )
#'
#' apv_life_flow(
#'   lt = lt,
#'   ages = c(60, 58),
#'   t = c(1, 2, 3),
#'   cf = c(100, 100, 100),
#'   i = 0.05,
#'   status = "first"
#' )
#'
#' @md
#' @export
apv_life_flow <- function(
    lt,
    ages,
    t = NULL,
    date = NULL,
    date0 = NULL,
    cf,
    i,
    i_type = "effective",
    m = 1L,
    status = c("single", "first", "last", "reversionary"),
    alpha = NULL,
    plot = FALSE
) {
  status <- match.arg(status)

  if (missing(i) || !is.numeric(i) || length(i) != 1L || is.na(i) ||
      !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
    stop("`i_type` must be a single character string.", call. = FALSE)
  }

  valid_i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

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
  m <- round(m)

  if (!is.data.frame(lt)) {
    stop("`lt` must be a data frame.", call. = FALSE)
  }

  if (!("x" %in% names(lt))) {
    stop("`lt` must contain column `x`.", call. = FALSE)
  }

  if (!("lx" %in% names(lt)) &&
      !("px" %in% names(lt)) &&
      !("qx" %in% names(lt))) {
    stop("`lt` must contain `lx`, `px`, or `qx`.", call. = FALSE)
  }

  if (!is.numeric(lt$x) || !nrow(lt) || any(!is.finite(lt$x)) ||
      any(lt$x != round(lt$x)) || anyDuplicated(lt$x) ||
      any(abs(lt$x) > .Machine$integer.max)) {
    stop("`lt$x` must contain unique finite integer ages.", call. = FALSE)
  }
  for (column in intersect(c("lx", "px", "qx"), names(lt))) {
    values <- lt[[column]]
    if (!is.numeric(values) || any(is.nan(values)) ||
        any(!is.finite(values) & !is.na(values)) ||
        any(values < 0, na.rm = TRUE) ||
        (column != "lx" && any(values > 1, na.rm = TRUE))) {
      stop("`lt$", column, "` contains invalid mortality values.", call. = FALSE)
    }
  }

  if (is.null(t) && is.null(date)) {
    stop("Provide either `t` or `date`.", call. = FALSE)
  }

  if (!is.null(t) && !is.null(date)) {
    stop("Provide only one of `t` or `date`.", call. = FALSE)
  }

  if (missing(cf) || !is.numeric(cf)) {
    stop("`cf` must be a numeric vector.", call. = FALSE)
  }

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be a logical scalar.", call. = FALSE)
  }

  # Build t from dates if needed
  date_out <- NULL
  if (!is.null(date)) {
    date_out <- as.Date(date)

    if (any(!is.finite(as.numeric(date_out)))) {
      stop("`date` must be coercible to finite non-missing Date values.", call. = FALSE)
    }

    if (is.null(date0)) {
      date0 <- min(date_out, na.rm = TRUE)
    }

    date0 <- as.Date(date0)

    if (length(date0) != 1L || !is.finite(as.numeric(date0))) {
      stop("`date0` must be a single finite non-missing Date.", call. = FALSE)
    }

    t <- as.numeric(date_out - date0) / 365.25
  }

  if (!is.numeric(t)) {
    stop("`t` must be numeric, in years.", call. = FALSE)
  }

  if (length(t) != length(cf)) {
    stop("`t`/`date` and `cf` must have the same length.", call. = FALSE)
  }

  if (any(t < 0, na.rm = TRUE)) {
    stop("All payment times `t` must be >= 0.", call. = FALSE)
  }

  if (any(!is.finite(t) & !is.na(t))) {
    stop("`t` must contain finite values or NA.", call. = FALSE)
  }

  if (any(!is.finite(cf) & !is.na(cf))) {
    stop("`cf` must contain finite values or NA.", call. = FALSE)
  }

  if (!is.numeric(ages) || length(ages) < 1L) {
    stop("`ages` must be a non-empty numeric vector.", call. = FALSE)
  }

  if (any(is.na(ages)) || any(!is.finite(ages))) {
    stop("`ages` must contain finite non-missing values.", call. = FALSE)
  }

  if (any(abs(ages - round(ages)) > 1e-10) ||
      any(abs(ages) > .Machine$integer.max)) {
    stop("`ages` must contain integer ages.", call. = FALSE)
  }

  ages <- round(ages)

  if (status == "reversionary") {
    if (is.null(alpha) || !is.numeric(alpha) || length(alpha) != 1L ||
        is.na(alpha) || !is.finite(alpha)) {
      stop("`alpha` must be a single finite numeric value.", call. = FALSE)
    }
  }

  basis <- .interest_basis(
    i_type = i_type,
    i = i,
    m = m
  )

  lt <- lt[order(lt$x), , drop = FALSE]
  if ("lx" %in% names(lt) &&
      any(diff(lt$lx[!is.na(lt$lx)]) > 0)) {
    stop("`lt$lx` must be nonincreasing.", call. = FALSE)
  }
  if (all(c("px", "qx") %in% names(lt)) &&
      any(abs(lt$px + lt$qx - 1) > 1e-10, na.rm = TRUE)) {
    stop("`lt$px` and `lt$qx` must be complementary.", call. = FALSE)
  }

  # --- one-year px at integer age ---
  get_px <- function(age) {
    if ("px" %in% names(lt)) {
      return(lt$px[match(age, lt$x)])
    }

    if ("qx" %in% names(lt)) {
      q <- lt$qx[match(age, lt$x)]
      return(1 - q)
    }

    l0 <- lt$lx[match(age, lt$x)]
    l1 <- lt$lx[match(age + 1L, lt$x)]

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    l1 / l0
  }

  get_qx <- function(age) {
    if ("qx" %in% names(lt)) {
      return(lt$qx[match(age, lt$x)])
    }

    px <- get_px(age)

    if (is.na(px)) {
      return(NA_real_)
    }

    1 - px
  }

  # --- {}_t p_x using UDD for fractional years ---
  t_px_udd <- function(x_age, tt) {
    if (is.na(tt)) return(NA_real_)
    if (tt <= 0) {
      return(1)
    }

    t_int <- floor(tt)
    u <- tt - t_int

    p_k <- 1
    # Never allocate a sequence proportional to an unsupported payment horizon.
    # Extinction is absorbing; an unknown earlier interval remains unknown.
    if (t_int > 0) {
      for (j in seq_len(min(t_int, nrow(lt)))) {
        px <- get_px(x_age + j - 1)
        if (is.na(px)) return(NA_real_)
        if (px == 0) return(0)
        p_k <- p_k * px
      }
      if (t_int > nrow(lt)) return(NA_real_)
    }

    if (u == 0) {
      return(p_k)
    }

    q <- get_qx(x_age + t_int)

    if (is.na(q)) {
      return(NA_real_)
    }

    p_k * (1 - u * q)
  }

  # --- survival probability by status ---
  p_status <- function(tt) {
    p_vec <- vapply(ages, function(a) t_px_udd(a, tt), numeric(1L))

    if (anyNA(p_vec)) {
      return(NA_real_)
    }

    if (length(ages) == 1L || status == "single") {
      return(p_vec[[1L]])
    }

    p_all <- prod(p_vec)
    p_any <- -expm1(sum(log1p(-p_vec)))

    if (status == "first") {
      return(p_all)
    }

    if (status == "last") {
      return(p_any)
    }

    # reversionary: full while all alive; alpha while partially alive
    p_all + alpha * (p_any - p_all)
  }

  surv   <- vapply(t, p_status, numeric(1L))
  log_disc <- -basis$growth * t
  disc   <- exp(log_disc)
  exp_cf <- cf * surv
  # Combine all three factors before rounding intermediate monetary products.
  pv <- .cash_flow_value(cf * sign(surv), log(abs(surv)) + log_disc)
  known_zero <- !is.na(cf) & !is.na(surv) & !is.na(t) &
    (cf == 0 | surv == 0)
  pv[known_zero] <- 0

  out <- tibble::tibble(
    t           = as.numeric(t),
    cf          = as.numeric(cf),
    surv_prob   = as.numeric(surv),
    discount    = as.numeric(disc),
    expected_cf = as.numeric(exp_cf),
    pv          = as.numeric(pv),
    pv_cum      = cumsum(as.numeric(pv))
  )

  # Insert date column if dates were provided
  if (!is.null(date_out)) {
    out <- tibble::tibble(date = date_out, out)
  }

  attr(out, "apv") <- sum(pv, na.rm = TRUE)
  attr(out, "contingency") <- if (status == "reversionary" && length(ages) > 1L)
    "weighted" else "survival"

  if (isTRUE(plot)) {
    plot_df <- out

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data[["t"]], y = .data[["pv_cum"]]
    )) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::labs(
        x     = "Time (years)",
        y     = "Cumulative APV",
        title = "Actuarial Present Value (cumulative)"
      ) +
      ggplot2::theme_minimal()

    attr(out, "plot") <- p
  }

  out
}
