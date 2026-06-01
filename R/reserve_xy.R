#' Benefit reserve schedule for two-life insurance
#'
#' Computes terminal benefit reserves at selected policy durations for a fully
#' discrete two-life insurance contract, assuming independent future lifetimes
#' and using compact actuarial notation.
#'
#' The function supports joint-life and last-survivor statuses, one common life
#' table for both lives, or two different life tables supplied as
#' \code{list(lt_x, lt_y)}.
#'
#' @param lt Either a single life table used for both lives, a list of two life
#'   tables \code{list(lt_x, lt_y)}, or a \code{tidyact_life_contract} object
#'   created by \code{\link{life_contract}}. Each table must contain columns
#'   \code{x} and \code{lx}.
#' @param x Integer actuarial age for the first life at issue.
#' @param y Integer actuarial age for the second life at issue.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param type Insurance type. One of \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param status Two-life status definition. Use \code{"joint"} for the
#'   joint-life status or \code{"last"} for the last-survivor status.
#' @param n Insurance term in years after deferment. Required as finite for
#'   term and endowment insurance. Use \code{Inf} for whole-life insurance.
#' @param h Nonnegative integer deferment period in years.
#' @param benefit Numeric benefit amount.
#' @param P Optional net premium per payment. If \code{NULL}, it is computed
#'   internally using \code{\link{premium_xy}}.
#' @param n_prem Optional premium-paying term in years, counted from
#'   \code{premium_start}. If \code{NULL}, premiums are payable for the
#'   corresponding default term.
#' @param k Positive integer. Number of premium payments per year.
#' @param timing Timing of premium payments. Use \code{"due"} or
#'   \code{"immediate"}. The recursive method currently requires annual due
#'   premiums.
#' @param premium_start Start of premium payments. Use \code{"issue"} for time
#'   0 or \code{"deferred"} for time \code{h}.
#' @param frac Fractional-age assumption used for status survival probabilities:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}.
#' @param t Integer vector of policy durations at which to compute reserves.
#'   If \code{NULL}, reserves are computed for all integer durations from
#'   issue to the contract horizon.
#' @param method Computation method. Use \code{"prospective"} or
#'   \code{"recursive"}.
#' @param tidy Logical scalar. If \code{TRUE}, returns a reserve schedule as a
#'   tibble. If \code{FALSE}, returns a named numeric vector.
#' @param check Logical scalar. If \code{TRUE}, performs basic input checks.
#' @param tol Numeric tolerance used for integer-grid checks.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age_x}, \code{age_y}, \code{rate},
#'   \code{rate_type}, \code{insurance_type}, \code{cohort},
#'   \code{term_years}, \code{deferment_years}, \code{premium},
#'   \code{premium_term_years}, \code{payments_per_year},
#'   \code{premium_timing}, \code{at}, and \code{output}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table input, \code{x} and
#' \code{y} are the two actuarial ages, \code{i} is the interest-rate input,
#' \code{i_type} is the interest-rate type, \code{m} is the conversion
#' frequency for nominal rates, \code{n} is the insurance term, \code{h} is the
#' deferment period, \code{k} is the premium payment frequency, \code{P} is the
#' premium per payment, and \code{t} is the policy duration.
#'
#' The prospective reserve at duration \eqn{t} is computed as
#' \deqn{
#' {}_tV =
#' APV_t(\text{future benefits}) -
#' P\,APV_t(\text{future premiums}).
#' }
#'
#' The recursive method is implemented only for nondeferred contracts with
#' annual due premiums. For deferred or subannual premium structures, use the
#' prospective method.
#'
#' @return
#' If \code{tidy = TRUE}, a tibble with reserve schedule details.
#' If \code{tidy = FALSE}, a named numeric vector.
#'
#' @seealso \code{\link{reserve_x}}, \code{\link{premium_xy}},
#'   \code{\link{insurance_xy}}, \code{\link{annuity_xy}}, \code{\link{t_pxy}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:70,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000,
#'          86000, 81000, 75000, 68000, 60000)
#' )
#'
#' reserve_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.06,
#'   type = "term",
#'   status = "joint",
#'   n = 4
#' )
#'
#' reserve_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.06,
#'   type = "endowment",
#'   status = "last",
#'   n = 5,
#'   benefit = 100000
#' )
#'
#' # Different life tables
#' lt2 <- data.frame(
#'   x  = 60:70,
#'   lx = c(100000, 99200, 98100, 96500, 94500, 92000,
#'          89000, 85000, 80000, 74000, 67000)
#' )
#'
#' reserve_xy(
#'   lt = list(lt, lt2),
#'   x = 60,
#'   y = 62,
#'   i = 0.06,
#'   type = "term",
#'   status = "joint",
#'   n = 4,
#'   tidy = TRUE
#' )
#'
#' @export
reserve_xy <- function(
    lt,
    x = NULL,
    y = NULL,
    i = NULL,
    i_type = "effective",
    m = 1L,
    type = c("whole", "term", "endowment"),
    status = c("joint", "last"),
    n = Inf,
    h = 0L,
    benefit = 1,
    P = NULL,
    n_prem = NULL,
    k = 1L,
    timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    frac = c("UDD", "CF", "CML", "Balducci"),
    t = NULL,
    method = c("prospective", "recursive"),
    tidy = TRUE,
    check = TRUE,
    tol = 1e-10,
    ...
) {
  dots <- list(...)
  status_missing <- missing(status)

  # Use exact matching for deprecated arguments.
  # This avoids partial matching problems such as `premium` matching
  # `premium_term_years`.
  dot_has <- function(nm) {
    nm %in% names(dots)
  }

  dot_get <- function(nm) {
    dots[[nm]]
  }

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "mortality_table",
    "age_x",
    "age_y",
    "rate",
    "rate_type",
    "insurance_type",
    "cohort",
    "term_years",
    "deferment_years",
    "premium",
    "premium_term_years",
    "payments_per_year",
    "premium_timing",
    "at",
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

  if (dot_has("mortality_table")) {
    if (!missing(lt)) {
      stop("Provide only one of `lt` or deprecated `mortality_table`.", call. = FALSE)
    }

    lt <- dot_get("mortality_table")
  }

  if (dot_has("age_x")) {
    if (!is.null(x)) {
      stop("Provide only one of `x` or deprecated `age_x`.", call. = FALSE)
    }

    x <- dot_get("age_x")
  }

  if (dot_has("age_y")) {
    if (!is.null(y)) {
      stop("Provide only one of `y` or deprecated `age_y`.", call. = FALSE)
    }

    y <- dot_get("age_y")
  }

  if (dot_has("rate")) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dot_get("rate")
  }

  if (dot_has("rate_type")) {
    if (!identical(i_type, "effective")) {
      stop("Provide only one of `i_type` or deprecated `rate_type`.", call. = FALSE)
    }

    i_type <- dot_get("rate_type")
  }

  if (dot_has("insurance_type")) {
    type <- dot_get("insurance_type")
  }

  if (dot_has("cohort")) {
    if (!status_missing) {
      stop("Provide only one of `status` or deprecated `cohort`.", call. = FALSE)
    }

    old_cohort <- match.arg(dot_get("cohort"), c("first", "last"))
    status <- if (old_cohort == "first") "joint" else "last"
  }

  if (dot_has("term_years")) {
    if (!is.infinite(n)) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }

    n <- dot_get("term_years")
  }

  if (dot_has("deferment_years")) {
    if (!identical(h, 0L) && !identical(h, 0)) {
      stop("Provide only one of `h` or deprecated `deferment_years`.", call. = FALSE)
    }

    h <- dot_get("deferment_years")
  }

  if (dot_has("premium")) {
    if (!is.null(P)) {
      stop("Provide only one of `P` or deprecated `premium`.", call. = FALSE)
    }

    P <- dot_get("premium")
  }

  if (dot_has("premium_term_years")) {
    if (!is.null(n_prem)) {
      stop("Provide only one of `n_prem` or deprecated `premium_term_years`.", call. = FALSE)
    }

    n_prem <- dot_get("premium_term_years")
  }

  if (dot_has("payments_per_year")) {
    if (!identical(k, 1L) && !identical(k, 1)) {
      stop("Provide only one of `k` or deprecated `payments_per_year`.", call. = FALSE)
    }

    k <- dot_get("payments_per_year")
  }

  if (dot_has("premium_timing")) {
    timing <- dot_get("premium_timing")
  }

  if (dot_has("at")) {
    if (!is.null(t)) {
      stop("Provide only one of `t` or deprecated `at`.", call. = FALSE)
    }

    t <- dot_get("at")
  }

  if (dot_has("output")) {
    if (!identical(tidy, TRUE)) {
      stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
    }

    output <- match.arg(dot_get("output"), c("table", "value"))
    tidy <- identical(output, "table")
  }

  type <- match.arg(type)
  status <- match.arg(status)
  timing <- match.arg(timing)
  premium_start <- match.arg(premium_start)
  frac <- match.arg(frac)
  method <- match.arg(method)

  if (frac == "CML") {
    frac <- "CF"
  }

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  # -------------------------------------------------------------------------
  # Resolve life_contract input
  # -------------------------------------------------------------------------

  if (!missing(lt) && inherits(lt, "tidyact_life_contract")) {
    contract <- lt

    if (!contract$lives %in% c("joint", "last_survivor")) {
      stop("`reserve_xy()` requires a two-life `life_contract()` object.", call. = FALSE)
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

    if (contract$lives == "last_survivor" && status_missing) {
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
        any(abs(tbl$x - round(tbl$x)) > tol)) {
      stop("Column `x` in `", arg, "` must contain finite integer ages.", call. = FALSE)
    }

    if (anyDuplicated(tbl$x)) {
      stop("Column `x` in `", arg, "` must not contain duplicated ages.", call. = FALSE)
    }

    if (any(is.na(tbl$lx)) || any(!is.finite(tbl$lx)) || any(tbl$lx < 0)) {
      stop("Column `lx` in `", arg, "` must contain finite nonnegative values.", call. = FALSE)
    }

    tbl <- tbl[order(tbl$x), , drop = FALSE]
    tbl$x <- as.integer(round(tbl$x))
    tbl
  }

  if (missing(lt)) {
    stop("`lt` must be provided.", call. = FALSE)
  }

  if (is.data.frame(lt)) {
    table_x <- validate_lifetable(lt, "lt")
    table_y <- table_x
    table_use <- table_x
  } else if (
    is.list(lt) &&
    length(lt) == 2L &&
    all(vapply(lt, is.data.frame, logical(1L)))
  ) {
    table_x <- validate_lifetable(lt[[1L]], "lt[[1]]")
    table_y <- validate_lifetable(lt[[2L]], "lt[[2]]")
    table_use <- list(table_x, table_y)
  } else {
    stop(
      "`lt` must be one life table, a list of two life tables, ",
      "or a `tidyact_life_contract` object.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (is.null(i) ||
      !is.numeric(i) ||
      length(i) != 1L ||
      is.na(i) ||
      !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) ||
      length(i_type) != 1L ||
      is.na(i_type)) {
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
      abs(x - round(x)) > tol) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (is.null(y) ||
      !is.numeric(y) ||
      length(y) != 1L ||
      is.na(y) ||
      !is.finite(y) ||
      abs(y - round(y)) > tol) {
    stop("`y` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      is.na(m) ||
      !is.finite(m) ||
      m < 1 ||
      abs(m - round(m)) > tol) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(h) ||
      length(h) != 1L ||
      is.na(h) ||
      !is.finite(h) ||
      h < 0 ||
      abs(h - round(h)) > tol) {
    stop("`h` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k < 1 ||
      abs(k - round(k)) > tol) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit) ||
      benefit <= 0) {
    stop("`benefit` must be a single positive finite number.", call. = FALSE)
  }

  if (!is.null(P) &&
      (!is.numeric(P) ||
       length(P) != 1L ||
       is.na(P) ||
       !is.finite(P) ||
       P < 0)) {
    stop("`P` must be NULL or a single nonnegative finite number.", call. = FALSE)
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n < 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) || abs(n - round(n)) > tol))) {
    stop("`n` must be `Inf` or a single nonnegative integer.", call. = FALSE)
  }

  if (type %in% c("term", "endowment") && is.infinite(n)) {
    stop("`n` must be finite for term and endowment insurance.", call. = FALSE)
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  m <- as.integer(round(m))
  h <- as.integer(round(h))
  k <- as.integer(round(k))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  if (!x %in% table_x$x) {
    stop("`x` must be present in the first life table.", call. = FALSE)
  }

  if (!y %in% table_y$x) {
    stop("`y` must be present in the second life table.", call. = FALSE)
  }

  if (method == "recursive" &&
      (h != 0L || k != 1L || timing != "due" || premium_start != "issue")) {
    stop(
      "The recursive method currently requires `h = 0`, `k = 1`, ",
      "`timing = 'due'`, and `premium_start = 'issue'`.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Interest
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

  v <- 1 / (1 + i_effective)

  # -------------------------------------------------------------------------
  # Contract horizon and premium horizon
  # -------------------------------------------------------------------------

  omega_x <- max(table_x$x, na.rm = TRUE)
  omega_y <- max(table_y$x, na.rm = TRUE)

  horizon_x <- max(0L, as.integer(round(omega_x - x)))
  horizon_y <- max(0L, as.integer(round(omega_y - y)))

  status_horizon <- if (status == "joint") {
    min(horizon_x, horizon_y)
  } else {
    max(horizon_x, horizon_y)
  }

  if (h > status_horizon) {
    stop("`h` exceeds the available two-life status horizon.", call. = FALSE)
  }

  contract_horizon <- if (type == "whole") {
    if (is.infinite(n)) {
      status_horizon
    } else {
      min(status_horizon, h + n)
    }
  } else {
    h + n
  }

  if (contract_horizon > status_horizon) {
    stop(
      "`h + n` exceeds the maximum term implied by the life tables for this status.",
      call. = FALSE
    )
  }

  contract_horizon <- as.integer(round(contract_horizon))

  if (is.null(n_prem)) {
    n_prem <- if (type == "whole") {
      if (premium_start == "issue") {
        contract_horizon
      } else {
        max(0L, contract_horizon - h)
      }
    } else {
      n
    }
  }

  if (!is.numeric(n_prem) ||
      length(n_prem) != 1L ||
      is.na(n_prem) ||
      !is.finite(n_prem) ||
      n_prem <= 0 ||
      abs(n_prem * k - round(n_prem * k)) > tol) {
    stop(
      "`n_prem` must be a positive value satisfying `n_prem * k` integer.",
      call. = FALSE
    )
  }

  n_prem <- round(n_prem * k) / k

  premium_start_time <- if (premium_start == "issue") {
    0
  } else {
    h
  }

  premium_end_time <- premium_start_time + n_prem

  if (premium_end_time > contract_horizon + tol) {
    stop(
      "The premium-paying period must not extend beyond the contract horizon.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Net premium
  # -------------------------------------------------------------------------

  premium_was_computed <- is.null(P)

  if (is.null(P)) {
    P <- premium_xy(
      lt = table_use,
      x = x,
      y = y,
      i = i,
      i_type = i_type,
      m = m,
      type = type,
      status = status,
      n = if (type == "whole") Inf else n,
      h = h,
      benefit = benefit,
      n_prem = n_prem,
      k = k,
      timing = timing,
      premium_start = premium_start,
      frac = frac,
      tidy = FALSE,
      check = check,
      tol = tol
    )
  }

  # -------------------------------------------------------------------------
  # Survival and APV helpers
  # -------------------------------------------------------------------------

  status_survival_from <- function(current_x, current_y, tt) {
    val <- t_pxy(
      lt = table_use,
      x = current_x,
      y = current_y,
      t = tt,
      frac = frac,
      status = status
    )

    if (is.na(val)) {
      stop(
        "Cannot compute two-life status survival at duration ",
        tt,
        ".",
        call. = FALSE
      )
    }

    val
  }

  future_benefit_apv <- function(duration) {
    if (duration >= contract_horizon) {
      if (type == "endowment" && duration == contract_horizon) {
        return(benefit)
      }

      return(0)
    }

    current_x <- x + duration
    current_y <- y + duration

    if (type == "whole") {
      h_remaining <- max(0L, h - duration)

      return(insurance_xy(
        lt = table_use,
        x = current_x,
        y = current_y,
        i = i,
        i_type = i_type,
        m = m,
        type = "whole",
        status = status,
        n = Inf,
        h = h_remaining,
        benefit = benefit,
        frac = frac,
        tidy = FALSE
      ))
    }

    coverage_end <- h + n

    if (duration < h) {
      h_remaining <- h - duration
      n_remaining <- n
    } else {
      h_remaining <- 0L
      n_remaining <- coverage_end - duration
    }

    if (n_remaining <= 0L) {
      if (type == "endowment" && duration == coverage_end) {
        return(benefit)
      }

      return(0)
    }

    insurance_xy(
      lt = table_use,
      x = current_x,
      y = current_y,
      i = i,
      i_type = i_type,
      m = m,
      type = type,
      status = status,
      n = n_remaining,
      h = h_remaining,
      benefit = benefit,
      frac = frac,
      tidy = FALSE
    )
  }

  future_premium_factor <- function(duration) {
    n_payments <- as.integer(round(n_prem * k))

    payment_times <- if (timing == "due") {
      premium_start_time + (0:(n_payments - 1L)) / k
    } else {
      premium_start_time + (1:n_payments) / k
    }

    future_times <- payment_times[payment_times >= duration - tol]

    if (length(future_times) == 0L) {
      return(0)
    }

    current_x <- x + duration
    current_y <- y + duration

    sum(vapply(future_times, function(abs_time) {
      elapsed <- abs_time - duration
      (v^elapsed) * status_survival_from(
        current_x = current_x,
        current_y = current_y,
        tt = elapsed
      )
    }, numeric(1L)))
  }

  # -------------------------------------------------------------------------
  # Durations
  # -------------------------------------------------------------------------

  if (is.null(t)) {
    t_vec <- 0:contract_horizon
  } else {
    if (!is.numeric(t) ||
        any(is.na(t)) ||
        any(!is.finite(t)) ||
        any(abs(t - round(t)) > tol)) {
      stop("`t` must be NULL or an integer-valued numeric vector.", call. = FALSE)
    }

    t_vec <- sort(unique(as.integer(round(t))))

    if (any(t_vec < 0) || any(t_vec > contract_horizon)) {
      stop("`t` durations must be between 0 and the contract horizon.", call. = FALSE)
    }
  }

  # -------------------------------------------------------------------------
  # Prospective or recursive reserves
  # -------------------------------------------------------------------------

  if (method == "prospective") {
    reserves <- vapply(t_vec, function(duration) {
      apv_benefits <- future_benefit_apv(duration)
      apv_premiums <- P * future_premium_factor(duration)

      apv_benefits - apv_premiums
    }, numeric(1L))

    if (premium_was_computed && 0L %in% t_vec) {
      reserves[t_vec == 0L] <- 0
    }
  } else {
    reserves_full <- numeric(contract_horizon + 1L)
    reserves_full[[1L]] <- 0

    for (duration in 0:(contract_horizon - 1L)) {
      premium_due <- if (duration < n_prem) P else 0
      benefit_due <- benefit

      p_status <- status_survival_from(
        current_x = x + duration,
        current_y = y + duration,
        tt = 1
      )

      if (p_status <= 0) {
        if (duration + 2L <= length(reserves_full)) {
          reserves_full[(duration + 2L):length(reserves_full)] <- 0
        }

        break
      }

      q_status <- 1 - p_status

      reserves_full[[duration + 2L]] <- (
        (reserves_full[[duration + 1L]] + premium_due) *
          (1 + i_effective) -
          benefit_due * q_status
      ) / p_status
    }

    if (type == "endowment") {
      reserves_full[[contract_horizon + 1L]] <- benefit
    }

    if (type == "term") {
      reserves_full[[contract_horizon + 1L]] <- 0
    }

    reserves <- reserves_full[t_vec + 1L]
  }

  # -------------------------------------------------------------------------
  # Output
  # -------------------------------------------------------------------------

  if (!tidy) {
    names(reserves) <- paste0("t=", t_vec)
    return(reserves)
  }

  premium_paid <- vapply(t_vec, function(duration) {
    if (duration >= premium_end_time - tol) {
      return(0)
    }

    if (duration >= premium_start_time - tol) {
      P
    } else {
      0
    }
  }, numeric(1L))

  benefit_due <- vapply(t_vec, function(duration) {
    if (duration >= contract_horizon) {
      return(0)
    }

    benefit
  }, numeric(1L))

  cohort_out <- if (identical(status, "joint")) {
    "first"
  } else {
    "last"
  }

  n_out <- if (identical(type, "whole")) {
    Inf
  } else {
    n
  }

  tibble::tibble(
    t = t_vec,
    duration = t_vec,
    x_t = x + t_vec,
    y_t = y + t_vec,
    age_x = x + t_vec,
    age_y = y + t_vec,
    V = reserves,
    reserve = reserves,
    P_paid = premium_paid,
    premium_paid = premium_paid,
    benefit_due = benefit_due,
    type = type,
    insurance_type = type,
    status = status,
    cohort = cohort_out,
    benefit = benefit,
    P = P,
    premium = P,
    n_prem = n_prem,
    premium_term_years = n_prem,
    k = k,
    payments_per_year = k,
    timing = timing,
    premium_timing = timing,
    premium_start = premium_start,
    n = n_out,
    term_years = n_out,
    h = h,
    deferment_years = h,
    contract_horizon = contract_horizon,
    method = method,
    frac = frac,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective
  )
}
