#' Actuarial present value of a payment stream under mortality
#'
#' Computes the actuarial present value (APV) of a cash-flow stream contingent
#' on survival. The life table is supplied as the first argument
#' (pipe-friendly). Payments may be specified by numeric times (years from 0)
#' or by calendar dates.
#'
#' Multiple lives are supported under an independence assumption, through
#' common statuses: single-life, first-death (all alive), last-survivor
#' (any alive), and reversionary (joint-and-survivor) with fraction
#' \code{alpha}.
#'
#' @param lt A life table data frame with column \code{x} and at least one of
#'   \code{lx}, \code{px}, or \code{qx}.
#' @param ages Integer vector of actuarial ages (length 1 for single life,
#'   length 2+ for multiple lives).
#' @param time Numeric vector of payment times in years (\code{>= 0}). Provide
#'   either \code{time} or \code{date}.
#' @param date Optional vector of \code{Date} payment dates. Provide either
#'   \code{time} or \code{date}.
#' @param start_date Optional \code{Date} used as time 0 when \code{date} is
#'   provided. If missing, the minimum of \code{date} is used.
#' @param cf Numeric vector of cash flows (same length as \code{time} or
#'   \code{date}).
#' @param i Annual effective interest rate (single numeric value).
#' @param status Survival status: \code{"single"}, \code{"first"},
#'   \code{"last"}, or \code{"reversionary"}.
#' @param alpha Reversionary fraction for \code{status = "reversionary"}
#'   (single numeric value). While all lives are alive, full benefit is paid;
#'   while at least one but not all are alive, \code{alpha} times the benefit
#'   is paid.
#' @param plot Logical; if \code{TRUE}, attaches a ggplot object in
#'   \code{attr(result, "plot")} showing cumulative APV over time.
#'
#' @details
#' For each payment at time \eqn{t}, the contribution to the APV is
#' (Finan, Sections 33 and 37):
#' \deqn{PV(t) = C(t) \times v^t \times P(StatusAliveMatT)}
#'
#' The survival probability depends on the \code{status}:
#' \itemize{
#'   \item \code{"single"}: \eqn{{}_t p_x} (single life).
#'   \item \code{"first"}: \eqn{{}_t p_{x_1} \cdot {}_t p_{x_2} \dots}
#'     (all must be alive - joint life).
#'   \item \code{"last"}: \eqn{1 - \prod (1 - {}_t p_{x_j})}
#'     (at least one alive - last survivor).
#'   \item \code{"reversionary"}: full benefit while all alive, fraction
#'     \eqn{\alpha} while partially alive.
#' }
#'
#' Fractional-year survival is computed under UDD within each year
#' (Finan, Section 24.1).
#'
#' @return A tibble with one row per payment and columns:
#'   \code{time}, \code{cf}, \code{surv_prob}, \code{discount},
#'   \code{expected_cf}, \code{pv}, \code{pv_cum}.
#'   If \code{date} was provided, a \code{date} column is included.
#'   The total APV is stored as \code{attr(result, "apv")}.
#' @export
apv_life_flow <- function(
    lt,
    ages,
    time = NULL,
    date = NULL,
    start_date = NULL,
    cf,
    i,
    status = c("single", "first", "last", "reversionary"),
    alpha = NULL,
    plot = FALSE
) {
  status <- match.arg(status)

  if (missing(i) || !is.numeric(i) || length(i) != 1) stop("'i' must be a single numeric rate.")
  if (!is.data.frame(lt)) stop("'lt' must be a data.frame.")
  if (!("x" %in% names(lt))) stop("Life table must contain column 'x'.")
  if (!("lx" %in% names(lt)) && !("px" %in% names(lt)) && !("qx" %in% names(lt))) {
    stop("Life table must contain 'lx', 'px', or 'qx'.")
  }

  if (is.null(time) && is.null(date)) stop("Provide either 'time' or 'date'.")
  if (!is.null(time) && !is.null(date)) stop("Provide only one of 'time' or 'date'.")
  if (missing(cf) || !is.numeric(cf)) stop("'cf' must be a numeric vector.")

  # Build time from dates if needed
  date_out <- NULL
  if (!is.null(date)) {
    date_out <- as.Date(date)
    if (is.null(start_date)) start_date <- min(date_out, na.rm = TRUE)
    start_date <- as.Date(start_date)
    time <- as.numeric(date_out - start_date) / 365.25
  }

  if (!is.numeric(time)) stop("'time' must be numeric (in years).")
  if (length(time) != length(cf)) stop("'time'/'date' and 'cf' must have the same length.")
  if (any(time < 0, na.rm = TRUE)) stop("All times must be >= 0.")

  if (!is.numeric(ages) || length(ages) < 1) stop("'ages' must be a non-empty numeric vector.")
  ages <- as.integer(ages)

  if (status == "reversionary") {
    if (is.null(alpha) || !is.numeric(alpha) || length(alpha) != 1) {
      stop("'alpha' must be a single numeric value.")
    }
  }

  v_fun <- function(tt) (1 + i)^(-tt)

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
    l1 <- lt$lx[match(age + 1, lt$x)]
    if (is.na(l0) || is.na(l1) || l0 <= 0) return(NA_real_)
    l1 / l0
  }

  get_qx <- function(age) {
    if ("qx" %in% names(lt)) return(lt$qx[match(age, lt$x)])
    px <- get_px(age)
    if (is.na(px)) return(NA_real_)
    1 - px
  }

  # --- {}_t p_x using UDD for fractional years (Finan, Sec. 24.1) ---
  t_px_udd <- function(x_age, tt) {
    if (tt <= 0) return(1)
    t_int <- floor(tt)
    u <- tt - t_int

    if (t_int > 0) {
      px_vec <- vapply(x_age:(x_age + t_int - 1), get_px, numeric(1))
      if (anyNA(px_vec)) return(NA_real_)
      p_k <- prod(px_vec)
    } else {
      p_k <- 1
    }

    if (u == 0) return(p_k)

    q <- get_qx(x_age + t_int)
    if (is.na(q)) return(NA_real_)
    p_k * (1 - u * q)
  }

  # --- survival probability by status ---
  P_status <- function(tt) {
    p_vec <- vapply(ages, function(a) t_px_udd(a, tt), numeric(1))
    if (anyNA(p_vec)) return(NA_real_)

    if (length(ages) == 1L || status == "single") return(p_vec[1])

    p_all <- prod(p_vec)
    p_any <- 1 - prod(1 - p_vec)

    if (status == "first") return(p_all)
    if (status == "last")  return(p_any)

    # reversionary: full while all alive; alpha while partially alive
    p_all + alpha * (p_any - p_all)
  }

  surv   <- vapply(time, P_status, numeric(1))
  disc   <- v_fun(time)
  exp_cf <- cf * surv
  pv     <- exp_cf * disc

  out <- tibble::tibble(
    time        = as.numeric(time),
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

  if (isTRUE(plot)) {
    plot_df <- out
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(
      x = .data[["time"]], y = .data[["pv_cum"]]
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
