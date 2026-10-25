#' Accrued coupon interest and clean/dirty bond prices
#'
#' @md
#' @description
#' Computes the coupon interest accrued between the previous coupon date and
#' settlement, and optionally reconciles a quoted clean or dirty bond price.
#' This is a day-count calculation, not a yield-based bond valuation.
#'
#' @param .data Optional data frame containing one bond per row. Required
#'   columns are `face`, `c`, `settlement`, `previous_coupon`, and
#'   `next_coupon`. Optional columns are `k`, `day_count`, `clean_price`, and
#'   `dirty_price`. Other columns are preserved. When `.data` is supplied,
#'   `face`, `c`, and the date arguments must come from its columns.
#' @param face Positive face value of the bond. Numeric scalar or vector in
#'   direct-input mode; otherwise the corresponding column of `.data`.
#' @param c Nonnegative annual nominal coupon rate, as a proportion (0.06 for
#'   6 percent). Numeric scalar or vector.
#' @param settlement `Date` of settlement. The accrued-interest convention
#'   excludes this date from the elapsed interval.
#' @param previous_coupon `Date` of the immediately preceding coupon date.
#' @param next_coupon `Date` of the immediately following coupon date.
#'   Settlement must satisfy `previous_coupon <= settlement < next_coupon`.
#'   At a coupon date, supply that date as `previous_coupon` to obtain zero
#'   accrued interest.
#' @param k Number of regular coupon payments per year, a positive integer;
#'   defaults to 1. For data-frame input, the `k` column, if present, takes
#'   precedence over this argument.
#' @param day_count Day-count convention, one of `"act/act"`, `"act/360"`,
#'   `"act/365"`, or `"30e/360"`. The default is `"act/act"`, defined as
#'   actual days divided by actual days in the coupon period, with annual
#'   coupon fraction `1/k`. For data-frame input, a `day_count` column, if
#'   present, takes precedence over this argument. `"30e/360"` denotes the
#'   European 30E/360 convention with both day-of-month values capped at 30.
#' @param clean_price Optional quoted clean price, in the same currency units
#'   and for the same face value as `face`.
#' @param dirty_price Optional dirty price, on the same basis. If either price
#'   is provided, the other is determined by the clean/dirty price identity.
#'   If both are supplied, their consistency is checked.
#' @param tol Positive scalar tolerance for reconciling prices when both are
#'   supplied; finite and used as a scaled absolute tolerance. The bound is
#'   `tol * max(1, abs(clean_price), abs(dirty_price), abs(accrued_interest))`
#'   for each row. The unit floor is part of this price-reconciliation rule;
#'   it is not a yield-solver tolerance.
#'
#' @details
#' Let `F` be face value, `c` the nominal annual coupon rate, and `k`
#' the coupon frequency. The regular coupon is `C = F * c / k`.
#'
#' For `"act/act"`, with `d` actual elapsed calendar days and `D` actual
#' days in the full coupon period, the year fraction and accrued interest are
#' `alpha = d / (k * D)` and `AI = F * c * alpha = C * d / D`.
#'
#' For `"act/360"` and `"act/365"`, `alpha` equals actual elapsed calendar
#' days divided by 360 or 365 respectively; for `"30e/360"`, it equals the
#' number of 30E/360 days divided by 360. These conventions accrue at the
#' annual coupon rate and therefore need not reach exactly one regular coupon
#' at the end of every actual calendar coupon period.
#'
#' The price identity is `dirty_price = clean_price + accrued_interest`.
#' A price supplied to this function is *not* computed from a yield curve.
#' Calculations assume a fixed-rate bond with regular coupon periods, no
#' stubs, no ex-coupon adjustment, no withholding tax, no amortizing principal,
#' and no transaction costs. The caller supplies the correct bracketing coupon
#' dates: a full coupon schedule is not inferred.
#'
#' All direct-input arguments of length greater than one must have the same
#' length; length-one values are recycled. Calendar dates must have class
#' `Date`, not bare character strings or numeric date serials.
#' Finite coupon and accrued amounts are recovered from overflowing intermediate
#' products when possible; nonrepresentable final amounts cause an error.
#'
#' @return A `tibble`. For direct input, it contains the bond parameters,
#'   settlement/coupon dates, the chosen convention, `coupon_per_period`,
#'   `days_elapsed`, `days_in_period`, `year_fraction`, `coupon_fraction`,
#'   `accrued_interest`, `clean_price`, and `dirty_price`. For data-frame
#'   input, original columns are retained and the audit columns added.
#'   When neither price is supplied, both price columns contain `NA_real_`.
#'
#' @seealso [bond_cash_flows()], [bond_price()], [bond_ytm()]
#' @family bonds
#'
#' @examples
#' accrued_interest(
#'   face = 1000, c = 0.06, k = 2,
#'   previous_coupon = as.Date("2026-01-01"),
#'   settlement = as.Date("2026-04-01"),
#'   next_coupon = as.Date("2026-07-01"),
#'   clean_price = 1005
#' )
#'
#' # The same position under two market day-count conventions.
#' accrued_interest(
#'   face = 1000, c = 0.06, k = 2,
#'   previous_coupon = as.Date("2026-01-01"),
#'   settlement = as.Date("2026-04-01"),
#'   next_coupon = as.Date("2026-07-01"),
#'   day_count = "30e/360"
#' )
#'
#' # Pipe-friendly portfolio input: each row is one bond position.
#' bonds <- tibble::tibble(
#'   bond = c("A", "B"), face = c(1000, 500), c = c(0.06, 0.04),
#'   k = c(2L, 2L),
#'   previous_coupon = as.Date(c("2026-01-01", "2026-01-01")),
#'   settlement = as.Date(c("2026-04-01", "2026-04-01")),
#'   next_coupon = as.Date(c("2026-07-01", "2026-07-01"))
#' )
#' bonds |> accrued_interest()
#'
#' @export
accrued_interest <- function(
    .data = NULL,
    face = NULL,
    c = NULL,
    settlement = NULL,
    previous_coupon = NULL,
    next_coupon = NULL,
    k = 1L,
    day_count = "act/act",
    clean_price = NULL,
    dirty_price = NULL,
    tol = 1e-10
) {
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a finite positive numeric scalar.", call. = FALSE)
  }

  if (is.null(.data)) {
    if (is.null(face) || is.null(c) || is.null(settlement) ||
        is.null(previous_coupon) || is.null(next_coupon)) {
      stop("`face`, `c`, `settlement`, `previous_coupon`, and `next_coupon` are required.",
           call. = FALSE)
    }
    original <- NULL
  } else {
    if (!is.data.frame(.data)) {
      stop("`.data` must be a data frame or NULL.", call. = FALSE)
    }
    if (nrow(.data) == 0L) {
      stop("`.data` must have at least one row.", call. = FALSE)
    }
    required <- c("face", "c", "settlement", "previous_coupon", "next_coupon")
    absent <- setdiff(required, names(.data))
    if (length(absent)) {
      stop("Missing required column(s) in `.data`: ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
    if (!is.null(face) || !is.null(c) || !is.null(settlement) ||
        !is.null(previous_coupon) || !is.null(next_coupon)) {
      stop("With `.data`, read `face`, `c`, and coupon dates from the data columns.",
           call. = FALSE)
    }
    original <- tibble::as_tibble(.data)
    face <- .data[["face"]]
    c <- .data[["c"]]
    settlement <- .data[["settlement"]]
    previous_coupon <- .data[["previous_coupon"]]
    next_coupon <- .data[["next_coupon"]]
    if ("k" %in% names(.data)) k <- .data[["k"]]
    if ("day_count" %in% names(.data)) day_count <- .data[["day_count"]]
    if ("clean_price" %in% names(.data)) {
      if (!is.null(clean_price)) {
        stop("Supply `clean_price` in `.data` or as an argument, not both.", call. = FALSE)
      }
      clean_price <- .data[["clean_price"]]
    }
    if ("dirty_price" %in% names(.data)) {
      if (!is.null(dirty_price)) {
        stop("Supply `dirty_price` in `.data` or as an argument, not both.", call. = FALSE)
      }
      dirty_price <- .data[["dirty_price"]]
    }
    audit_names <- c("coupon_per_period", "days_elapsed", "days_in_period",
                     "year_fraction", "coupon_fraction", "accrued_interest")
    if (any(audit_names %in% names(.data))) {
      stop("`.data` contains existing audit columns; refusing to overwrite them.",
           call. = FALSE)
    }
  }

  required_inputs <- list(face = face, c = c, k = k, settlement = settlement,
                          previous_coupon = previous_coupon, next_coupon = next_coupon,
                          day_count = day_count)
  if (!is.null(clean_price)) required_inputs$clean_price <- clean_price
  if (!is.null(dirty_price)) required_inputs$dirty_price <- dirty_price
  sizes <- vapply(required_inputs, length, integer(1))
  if (any(sizes == 0L)) {
    stop("Required arguments must have positive length.", call. = FALSE)
  }
  n <- if (is.null(original)) max(sizes) else nrow(original)
  if (any(!(sizes %in% c(1L, n)))) {
    stop("Inputs must have length one or exactly one value per bond.", call. = FALSE)
  }
  expand <- function(x) rep(x, length.out = n)

  face <- expand(face)
  c <- expand(c)
  k <- expand(k)
  settlement <- expand(settlement)
  previous_coupon <- expand(previous_coupon)
  next_coupon <- expand(next_coupon)
  day_count <- expand(day_count)
  clean_price <- if (is.null(clean_price)) rep(NA_real_, n) else expand(clean_price)
  dirty_price <- if (is.null(dirty_price)) rep(NA_real_, n) else expand(dirty_price)

  valid_numeric <- function(x, lower = -Inf, strict = FALSE) {
    is.numeric(x) && !is.object(x) &&
      !anyNA(x) && all(is.finite(x)) &&
      (if (strict) all(x > lower) else all(x >= lower))
  }
  if (!valid_numeric(face, 0, strict = TRUE)) {
    stop("`face` must contain finite strictly positive numbers.", call. = FALSE)
  }
  if (!valid_numeric(c, 0)) {
    stop("`c` must contain finite nonnegative coupon rates.", call. = FALSE)
  }
  if (!valid_numeric(k, 0, strict = TRUE) || any(k != floor(k)) ||
      any(k > .Machine$integer.max)) {
    stop("`k` must contain positive integers.", call. = FALSE)
  }
  k <- as.integer(k)
  if (!is.character(day_count) || anyNA(day_count) ||
      any(!(day_count %in% c("act/act", "act/360", "act/365", "30e/360")))) {
    stop("`day_count` must be one of 'act/act', 'act/360', 'act/365', or '30e/360'.",
         call. = FALSE)
  }
  for (nm in c("settlement", "previous_coupon", "next_coupon")) {
    if (!inherits(get(nm), "Date") || anyNA(get(nm))) {
      stop("`", nm, "` must contain valid Date values.", call. = FALSE)
    }
  }
  if (any(previous_coupon > settlement) || any(settlement >= next_coupon)) {
    stop("Dates must satisfy previous_coupon <= settlement < next_coupon.",
         call. = FALSE)
  }
  for (nm in c("clean_price", "dirty_price")) {
    x <- get(nm)
    if (!is.numeric(x) || any(!is.na(x) & !is.finite(x))) {
      stop("`", nm, "` must contain finite numeric prices or NA.",
           call. = FALSE)
    }
  }

  actual_elapsed <- as.numeric(settlement - previous_coupon)
  actual_period <- as.numeric(next_coupon - previous_coupon)
  day_30e_360 <- function(from, to) {
    y1 <- as.integer(format(from, "%Y"))
    m1 <- as.integer(format(from, "%m"))
    d1 <- pmin(as.integer(format(from, "%d")), 30L)
    y2 <- as.integer(format(to, "%Y"))
    m2 <- as.integer(format(to, "%m"))
    d2 <- pmin(as.integer(format(to, "%d")), 30L)
    (y2 - y1) * 360 + (m2 - m1) * 30 + d2 - d1
  }
  elapsed_30e <- day_30e_360(previous_coupon, settlement)
  year_fraction <- numeric(n)
  use_act_act <- day_count == "act/act"
  use_act_360 <- day_count == "act/360"
  use_act_365 <- day_count == "act/365"
  use_30e <- day_count == "30e/360"
  year_fraction[use_act_act] <- actual_elapsed[use_act_act] /
    (k[use_act_act] * actual_period[use_act_act])
  year_fraction[use_act_360] <- actual_elapsed[use_act_360] / 360
  year_fraction[use_act_365] <- actual_elapsed[use_act_365] / 365
  year_fraction[use_30e] <- elapsed_30e[use_30e] / 360

  if (any(!is.finite(year_fraction)) || any(year_fraction < 0)) {
    stop("Accrued day-count fractions must be finite and nonnegative.",
         call. = FALSE)
  }
  coupon_per_period <- face * c / k
  accrued <- face * c * year_fraction
  recover_coupon <- !is.finite(coupon_per_period) |
    (coupon_per_period == 0 & c > 0)
  coupon_per_period[recover_coupon] <- .cash_flow_value(
    face[recover_coupon], log(c[recover_coupon]) - log(k[recover_coupon]))
  recover_accrued <- !is.finite(accrued) |
    (accrued == 0 & c > 0 & year_fraction > 0)
  accrued[recover_accrued] <- .cash_flow_value(
    face[recover_accrued], log(c[recover_accrued]) + log(year_fraction[recover_accrued]))
  if (any(!is.finite(coupon_per_period)) || any(!is.finite(accrued))) {
    stop("Coupon/accrued interest overflowed; check monetary inputs.",
         call. = FALSE)
  }

  both <- !is.na(clean_price) & !is.na(dirty_price)
  if (any(both)) {
    residual <- dirty_price[both] - clean_price[both] - accrued[both]
    scale <- pmax(1, abs(clean_price[both]), abs(dirty_price[both]),
                  abs(accrued[both]))
    if (any(abs(residual) > tol * scale)) {
      stop("Supplied `clean_price` and `dirty_price` are inconsistent with accrued interest.",
           call. = FALSE)
    }
  }
  missing_dirty <- !is.na(clean_price) & is.na(dirty_price)
  dirty_price[missing_dirty] <- clean_price[missing_dirty] + accrued[missing_dirty]
  missing_clean <- is.na(clean_price) & !is.na(dirty_price)
  clean_price[missing_clean] <- dirty_price[missing_clean] - accrued[missing_clean]
  if (any(!is.finite(clean_price[!is.na(clean_price)])) ||
      any(!is.finite(dirty_price[!is.na(dirty_price)]))) {
    stop("The implied price is not finite.", call. = FALSE)
  }

  coupon_fraction <- k * year_fraction
  audit <- tibble::tibble(
    coupon_per_period = coupon_per_period,
    days_elapsed = actual_elapsed,
    days_in_period = actual_period,
    year_fraction = year_fraction,
    coupon_fraction = coupon_fraction,
    accrued_interest = accrued,
    clean_price = clean_price,
    dirty_price = dirty_price
  )
  if (is.null(original)) {
    return(tibble::tibble(
      face = face,
      c = c,
      k = k,
      settlement = settlement,
      previous_coupon = previous_coupon,
      next_coupon = next_coupon,
      day_count = day_count,
      coupon_per_period = audit$coupon_per_period,
      days_elapsed = audit$days_elapsed,
      days_in_period = audit$days_in_period,
      year_fraction = audit$year_fraction,
      coupon_fraction = audit$coupon_fraction,
      accrued_interest = audit$accrued_interest,
      clean_price = audit$clean_price,
      dirty_price = audit$dirty_price
    ))
  }
  # Preserve the original bond identifiers and columns, adding audit fields.
  original[["coupon_per_period"]] <- audit$coupon_per_period
  original[["days_elapsed"]] <- audit$days_elapsed
  original[["days_in_period"]] <- audit$days_in_period
  original[["year_fraction"]] <- audit$year_fraction
  original[["coupon_fraction"]] <- audit$coupon_fraction
  original[["accrued_interest"]] <- audit$accrued_interest
  original[["clean_price"]] <- audit$clean_price
  original[["dirty_price"]] <- audit$dirty_price
  original
}
