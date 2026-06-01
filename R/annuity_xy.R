#' Actuarial present value of a two-life annuity
#'
#' Computes the actuarial present value of a discrete annuity contingent on two
#' independent lives, using compact actuarial notation.
#'
#' The function supports joint-life, last-survivor, and state-based
#' reversionary-style payments. The life table input may be either one common
#' table for both lives or a list of two life tables, one for each life.
#'
#' @param lt A life table, a list of two life tables \code{list(lt_x, lt_y)}, or
#'   a \code{tidyact_life_contract} object created by \code{\link{life_contract}}.
#'   Each life table must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age for the first life.
#' @param y Integer actuarial age for the second life.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param status Character string. Use \code{"joint"} for joint-life payments
#'   while both lives are alive, or \code{"last"} for last-survivor payments
#'   while at least one life is alive. Ignored when \code{benefit} is supplied.
#' @param benefit Optional list with numeric scalar weights \code{both},
#'   \code{x_only}, and \code{y_only}. These weights define the payment rate
#'   according to the survival state of the two lives.
#' @param n Term in years. Use \code{Inf} for the maximum horizon allowed by
#'   the life tables.
#' @param h Nonnegative integer deferment period in years.
#' @param k Positive integer. Number of payments per year.
#' @param timing Payment timing. Use \code{"immediate"} or \code{"due"}.
#' @param woolhouse Woolhouse approximation for \code{k > 1}:
#'   \code{"none"}, \code{"first"}, or \code{"second"}.
#' @param frac Fractional-age assumption for exact k-thly computation:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric APV. If
#'   \code{TRUE}, returns a one-row tibble.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age_x}, \code{age_y}, \code{rate},
#'   \code{rate_type}, \code{cohort}, \code{term_years},
#'   \code{deferment_years}, \code{payments_per_year}, and \code{output}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table input, \code{x} and
#' \code{y} are the two actuarial ages, \code{i} is the interest-rate input,
#' \code{i_type} is the interest-rate type, \code{m} is the conversion
#' frequency for nominal rates, \code{n} is the term, \code{h} is the
#' deferment period, and \code{k} is the payment frequency.
#'
#' The function assumes independent future lifetimes. For state-based benefits,
#' the expected payment at time \eqn{t} is
#' \deqn{
#' b_{\text{both}}\,{}_tp_x{}_tp_y
#' + b_{x\text{ only}}\,{}_tp_x(1-{}_tp_y)
#' + b_{y\text{ only}}\,{}_tp_y(1-{}_tp_x).
#' }
#'
#' When \code{benefit = NULL}, \code{status = "joint"} uses
#' \code{benefit = list(both = 1, x_only = 0, y_only = 0)}, while
#' \code{status = "last"} uses
#' \code{benefit = list(both = 1, x_only = 1, y_only = 1)}.
#'
#' For \code{k}-thly payments, the function values a payment rate of 1 per year,
#' so each payment has size \eqn{1/k}.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric scalar.
#'
#' If \code{tidy = TRUE}, a one-row tibble with input values, standardized
#' interest rate, term used, and APV.
#'
#' @seealso \code{\link{annuity_x}}, \code{\link{insurance_xy}},
#'   \code{\link{premium_xy}}, \code{\link{t_pxy}}, \code{\link{t_px}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Joint-life annuity-due
#' annuity_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   status = "joint",
#'   timing = "due"
#' )
#'
#' # Last-survivor annuity-due
#' annuity_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   status = "last",
#'   timing = "due"
#' )
#'
#' # Different life tables for the two lives
#' lt_m <- lt
#' lt_f <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99200, 98100, 96500, 94500, 92000, 89000)
#' )
#'
#' annuity_xy(
#'   lt = list(lt_m, lt_f),
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   status = "joint",
#'   timing = "due"
#' )
#'
#' # State-based reversionary-style payments
#' annuity_xy(
#'   lt = list(lt_m, lt_f),
#'   x = 60,
#'   y = 62,
#'   i = 0.05,
#'   benefit = list(both = 0, x_only = 1, y_only = 0),
#'   timing = "due"
#' )
#'
#' @export
annuity_xy <- function(
    lt,
    x = NULL,
    y = NULL,
    i = NULL,
    i_type = "effective",
    m = 1L,
    status = c("joint", "last"),
    benefit = NULL,
    n = Inf,
    h = 0L,
    k = 1L,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    frac,
    tidy = FALSE,
    ...
) {
  dots <- list(...)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "mortality_table",
    "age_x",
    "age_y",
    "rate",
    "rate_type",
    "cohort",
    "term_years",
    "deferment_years",
    "payments_per_year",
    "output"
  )

  bad_dots <- setdiff(names(dots), allowed_old)

  if (length(bad_dots) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad_dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(dots$mortality_table)) {
    if (!missing(lt)) {
      stop("Provide only one of `lt` or deprecated `mortality_table`.", call. = FALSE)
    }
    lt <- dots$mortality_table
  }

  if (!is.null(dots$age_x)) {
    if (!is.null(x)) {
      stop("Provide only one of `x` or deprecated `age_x`.", call. = FALSE)
    }
    x <- dots$age_x
  }

  if (!is.null(dots$age_y)) {
    if (!is.null(y)) {
      stop("Provide only one of `y` or deprecated `age_y`.", call. = FALSE)
    }
    y <- dots$age_y
  }

  if (!is.null(dots$rate)) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }
    i <- dots$rate
  }

  if (!is.null(dots$rate_type)) {
    if (!identical(i_type, "effective")) {
      stop("Provide only one of `i_type` or deprecated `rate_type`.", call. = FALSE)
    }
    i_type <- dots$rate_type
  }

  if (!is.null(dots$cohort)) {
    if (!missing(status)) {
      stop("Provide only one of `status` or deprecated `cohort`.", call. = FALSE)
    }

    old_cohort <- match.arg(dots$cohort, c("first", "last"))
    status <- if (old_cohort == "first") "joint" else "last"
  }

  if (!is.null(dots$term_years)) {
    if (!is.infinite(n)) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }
    n <- dots$term_years
  }

  if (!is.null(dots$deferment_years)) {
    if (!identical(h, 0L) && !identical(h, 0)) {
      stop("Provide only one of `h` or deprecated `deferment_years`.", call. = FALSE)
    }
    h <- dots$deferment_years
  }

  if (!is.null(dots$payments_per_year)) {
    if (!identical(k, 1L) && !identical(k, 1)) {
      stop("Provide only one of `k` or deprecated `payments_per_year`.", call. = FALSE)
    }
    k <- dots$payments_per_year
  }

  if (!is.null(dots$output)) {
    if (!identical(tidy, FALSE)) {
      stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
    }

    output <- match.arg(dots$output, c("value", "table"))
    tidy <- identical(output, "table")
  }

  status <- match.arg(status)
  timing <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (!missing(lt) && inherits(lt, "tidyact_life_contract")) {
    contract <- lt

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop("`annuity_xy()` requires a two-life `life_contract()` object.", call. = FALSE)
    }

    lt <- contract$mortality_table

    if (is.null(x)) {
      x <- contract$x %||% contract$age_x
    }

    if (is.null(y)) {
      y <- contract$y %||% contract$age_y
    }

    if (is.null(i)) {
      i <- contract$i %||% contract$rate
    }

    if (identical(i_type, "effective")) {
      i_type <- contract$i_type %||% contract$rate_type %||% i_type
    }

    if (identical(m, 1L) || identical(m, 1)) {
      m <- contract$m %||% m
    }

    if (contract$lives == "last_survivor") {
      status <- "last"
    }
  }

  if (is.null(i_type)) {
    i_type <- "effective"
  }

  if (is.null(m)) {
    m <- 1L
  }

  # -------------------------------------------------------------------------
  # Life table handling
  # -------------------------------------------------------------------------

  validate_lifetable <- function(tbl, arg) {
    if (!is.data.frame(tbl)) {
      stop("`", arg, "` must be a data frame or tibble.", call. = FALSE)
    }

    if (!all(c("x", "lx") %in% names(tbl))) {
      stop("`", arg, "` must contain columns `x` and `lx`.", call. = FALSE)
    }

    if (!is.numeric(tbl$x) || !is.numeric(tbl$lx)) {
      stop("Columns `x` and `lx` in `", arg, "` must be numeric.", call. = FALSE)
    }

    if (any(is.na(tbl$x)) || any(!is.finite(tbl$x)) ||
        any(abs(tbl$x - round(tbl$x)) > 1e-10)) {
      stop("Column `x` in `", arg, "` must contain finite integer ages.", call. = FALSE)
    }

    if (anyDuplicated(tbl$x)) {
      stop("Column `x` in `", arg, "` must not contain duplicated ages.", call. = FALSE)
    }

    if (any(is.na(tbl$lx)) || any(!is.finite(tbl$lx)) || any(tbl$lx < 0)) {
      stop("Column `lx` in `", arg, "` must contain finite nonnegative values.", call. = FALSE)
    }

    tbl[order(tbl$x), , drop = FALSE]
  }

  if (missing(lt)) {
    stop("`lt` must be provided.", call. = FALSE)
  }

  if (is.data.frame(lt)) {
    lt_x <- validate_lifetable(lt, "lt")
    lt_y <- lt_x
  } else if (
    is.list(lt) &&
    length(lt) == 2L &&
    all(vapply(lt, is.data.frame, logical(1L)))
  ) {
    lt_x <- validate_lifetable(lt[[1L]], "lt[[1]]")
    lt_y <- validate_lifetable(lt[[2L]], "lt[[2]]")
  } else {
    stop(
      "`lt` must be one life table, a list of two life tables, ",
      "or a `tidyact_life_contract` object.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Fractional-age assumption
  # -------------------------------------------------------------------------

  if (missing(frac)) {
    frac_x <- attr(lt_x, "frac")
    frac_y <- attr(lt_y, "frac")
    ok_x <- !is.null(frac_x) && frac_x %in% c("UDD", "CF", "Balducci")
    ok_y <- !is.null(frac_y) && frac_y %in% c("UDD", "CF", "Balducci")

    if (ok_x && ok_y) {
      if (!identical(frac_x, frac_y)) {
        stop(
          "The two life tables carry different `frac` attributes. ",
          "Supply `frac` explicitly.",
          call. = FALSE
        )
      }

      frac <- frac_x
    } else if (ok_x) {
      frac <- frac_x
    } else if (ok_y) {
      frac <- frac_y
    } else {
      frac <- "UDD"
    }
  } else {
    frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))

    if (frac == "CML") {
      frac <- "CF"
    }
  }

  # -------------------------------------------------------------------------
  # Basic validation
  # -------------------------------------------------------------------------

  if (is.null(i) ||
      !is.numeric(i) ||
      length(i) != 1L ||
      is.na(i) ||
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

  if (is.null(x) ||
      !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      abs(x - round(x)) > 1e-10) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (is.null(y) ||
      !is.numeric(y) ||
      length(y) != 1L ||
      is.na(y) ||
      !is.finite(y) ||
      abs(y - round(y)) > 1e-10) {
    stop("`y` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      !is.finite(m) || m < 1 || abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(h) || length(h) != 1L ||
      is.na(h) || !is.finite(h) ||
      h < 0 || abs(h - round(h)) > 1e-10) {
    stop("`h` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1L ||
      is.na(k) || !is.finite(k) ||
      k < 1 || abs(k - round(k)) > 1e-10) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L ||
      is.na(n) || n < 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) || abs(n - round(n)) > 1e-10))) {
    stop("`n` must be `Inf` or a single nonnegative integer.", call. = FALSE)
  }

  if (!is.null(benefit)) {
    if (!is.list(benefit) || !all(c("both", "x_only", "y_only") %in% names(benefit))) {
      stop("`benefit` must be a list with names `both`, `x_only`, and `y_only`.",
           call. = FALSE)
    }

    benefit_vals <- benefit[c("both", "x_only", "y_only")]

    ok_benefit <- vapply(
      benefit_vals,
      function(z) {
        is.numeric(z) && length(z) == 1L && !is.na(z) && is.finite(z)
      },
      logical(1L)
    )

    if (!all(ok_benefit)) {
      stop(
        "`benefit$both`, `benefit$x_only`, and `benefit$y_only` ",
        "must be finite numeric scalars.",
        call. = FALSE
      )
    }
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (!x %in% as.integer(round(lt_x$x))) {
    stop("`x` must be present in the first life table.", call. = FALSE)
  }

  if (!y %in% as.integer(round(lt_y$x))) {
    stop("`y` must be present in the second life table.", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Interest conversion
  # -------------------------------------------------------------------------

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  v_fun <- function(tt) {
    (1 + i_effective)^(-tt)
  }

  # -------------------------------------------------------------------------
  # Horizon and benefit state weights
  # -------------------------------------------------------------------------

  omega_x <- max(lt_x$x, na.rm = TRUE)
  omega_y <- max(lt_y$x, na.rm = TRUE)

  horizon_x <- max(0L, as.integer(round(omega_x - x)))
  horizon_y <- max(0L, as.integer(round(omega_y - y)))

  if (is.null(benefit)) {
    horizon <- if (status == "joint") {
      min(horizon_x, horizon_y)
    } else {
      max(horizon_x, horizon_y)
    }
  } else {
    h_both <- if (benefit$both != 0) min(horizon_x, horizon_y) else 0L
    h_x_only <- if (benefit$x_only != 0) horizon_x else 0L
    h_y_only <- if (benefit$y_only != 0) horizon_y else 0L
    horizon <- max(h_both, h_x_only, h_y_only)
  }

  term_used <- if (is.infinite(n)) {
    max(0L, horizon - h)
  } else {
    n
  }

  if (term_used == 0L) {
    result <- 0

    if (!tidy) {
      return(result)
    }

    return(tibble::tibble(
      x = x,
      y = y,
      age_x = x,
      age_y = y,
      i = i,
      rate = i,
      i_type = i_type,
      rate_type = i_type,
      m = m,
      i_effective = i_effective,
      n = n,
      term_years = n,
      term_used = term_used,
      h = h,
      deferment_years = h,
      k = k,
      payments_per_year = k,
      timing = timing,
      status = status,
      cohort = if (status == "joint") "first" else "last",
      woolhouse = woolhouse,
      frac = frac,
      apv = result
    ))
  }

  if (is.null(benefit)) {
    benefit <- if (status == "joint") {
      list(both = 1, x_only = 0, y_only = 0)
    } else {
      list(both = 1, x_only = 1, y_only = 1)
    }
  }

  expected_payment <- function(tt) {
    px_t <- t_px(
      lt = lt_x,
      x = x,
      t = tt,
      frac = frac,
      tidy = FALSE,
      check = FALSE
    )

    py_t <- t_px(
      lt = lt_y,
      x = y,
      t = tt,
      frac = frac,
      tidy = FALSE,
      check = FALSE
    )

    if (is.na(px_t) || is.na(py_t)) {
      return(NA_real_)
    }

    p_both <- px_t * py_t
    p_x_only <- px_t * (1 - py_t)
    p_y_only <- py_t * (1 - px_t)

    benefit$both * p_both +
      benefit$x_only * p_x_only +
      benefit$y_only * p_y_only
  }

  exact_apv <- function(nn, kk, tim) {
    N <- nn * kk

    if (N == 0L) {
      return(0)
    }

    j <- if (tim == "due") {
      0:(N - 1L)
    } else {
      1:N
    }

    u <- j / kk
    times <- h + u
    disc <- v_fun(times)
    ep <- vapply(times, expected_payment, numeric(1L))

    if (anyNA(ep)) {
      return(NA_real_)
    }

    sum((1 / kk) * disc * ep)
  }

  # -------------------------------------------------------------------------
  # APV computation
  # -------------------------------------------------------------------------

  if (k == 1L || woolhouse == "none") {
    result <- exact_apv(term_used, k, timing)
  } else {
    annual_due <- exact_apv(term_used, 1L, "due")

    ep_start <- expected_payment(h)
    ep_end <- expected_payment(h + term_used)

    if (is.na(ep_start) || ep_start <= 0) {
      ep_start <- 1
    }

    if (is.na(ep_end)) {
      ep_end <- 0
    }

    nEx_status <- v_fun(term_used) * ep_end / ep_start

    adj1 <- (k - 1) / (2 * k) * (1 - nEx_status)

    if (woolhouse == "first") {
      due_k <- annual_due - adj1
    } else {
      delta <- log1p(i_effective)

      ep_m1 <- expected_payment(h + 1)

      if (is.na(ep_m1) || ep_start <= 0) {
        mu_start <- 0
      } else {
        p_status_1 <- ep_m1 / ep_start
        mu_start <- if (p_status_1 > 0) -log(p_status_1) else 0
      }

      ep_n1 <- expected_payment(h + term_used + 1)

      if (is.na(ep_n1) || is.na(ep_end) || ep_end <= 0) {
        mu_end <- 0
      } else {
        p_status_n1 <- ep_n1 / ep_end
        mu_end <- if (p_status_n1 > 0) -log(p_status_n1) else 0
      }

      adj2 <- (k^2 - 1) /
        (12 * k^2) *
        (delta + mu_start - nEx_status * (delta + mu_end))

      due_k <- annual_due - adj1 - adj2
    }

    result <- if (timing == "due") {
      due_k
    } else {
      due_k - (1 / k) * (1 - nEx_status)
    }
  }

  if (!is.finite(result) && !is.na(result)) {
    stop("The annuity APV calculation produced a non-finite value.", call. = FALSE)
  }

  if (!tidy) {
    return(result)
  }

  tibble::tibble(
    x = x,
    y = y,
    age_x = x,
    age_y = y,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective,
    n = n,
    term_years = n,
    term_used = term_used,
    h = h,
    deferment_years = h,
    k = k,
    payments_per_year = k,
    timing = timing,
    status = status,
    cohort = if (status == "joint") "first" else "last",
    woolhouse = woolhouse,
    frac = frac,
    apv = result
  )
}
