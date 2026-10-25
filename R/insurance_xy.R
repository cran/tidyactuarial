#' Actuarial present value of a two-life insurance
#'
#' Computes the actuarial present value of a discrete two-life insurance with
#' benefit payable at the end of the year of the triggering death, assuming
#' independent future lifetimes and using compact actuarial notation.
#'
#' Supported contracts:
#' \itemize{
#'   \item \code{"whole"}: whole-life two-life insurance.
#'   \item \code{"term"}: n-year two-life term insurance.
#'   \item \code{"endowment"}: n-year two-life endowment insurance.
#' }
#'
#' The \code{status} argument determines the two-life status:
#' \itemize{
#'   \item \code{status = "joint"}: first-death insurance, based on the
#'   joint-life status.
#'   \item \code{status = "last"}: second-death insurance, based on the
#'   last-survivor status.
#' }
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
#' @param type Character string. One of \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param status Character string. Use \code{"joint"} for first-death insurance
#'   or \code{"last"} for second-death insurance.
#' @param n Term in years. Required as finite for term and endowment insurance.
#'   Use \code{Inf} for whole-life insurance.
#' @param h Nonnegative integer deferment period in years.
#' @param benefit Numeric scalar. Insurance benefit.
#' @param frac Fractional-age assumption used in two-life survival
#'   probabilities. One of \code{"UDD"}, \code{"CF"}, \code{"CML"}, or
#'   \code{"Balducci"}.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric APV. If
#'   \code{TRUE}, returns a one-row tibble.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age_x}, \code{age_y}, \code{rate},
#'   \code{rate_type}, \code{insurance_type}, \code{cohort},
#'   \code{term_years}, \code{deferment_years}, and \code{output}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table input, \code{x} and
#' \code{y} are the two actuarial ages, \code{i} is the interest-rate input,
#' \code{i_type} is the interest-rate type, \code{m} is the conversion
#' frequency for nominal rates, \code{n} is the insurance term, and \code{h} is
#' the deferment period.
#'
#' Both lives are alive at issue (state 11). Internal states 10 and 01 mean
#' only x or only y is alive; 00 is absorbing. The API does not take an
#' observed single-survivor initial state. In a joint-life contract, leaving
#' 11 triggers the first-death benefit once. A last-survivor contract continues
#' through 10 or 01 and pays on entering 00. Both deaths within one annual
#' interval trigger one benefit, payable at the end of that year.
#'
#' Valuation sums discounted annual transition masses from the issue ages;
#' it does not subtract rounded annuity values. With marginal survival Sx, Sy,
#' cumulative death Qx, Qy, and interval death masses Dx, Dy, the first-death
#' mass in (t,t+1] is Dx(t) Sy(t) + Sx(t+1) Dy(t); the second-death mass is
#' Dx(t) Qy(t) + Qx(t) Dy(t) + Dx(t) Dy(t). Counts supply the decrements
#' directly. Original log interest and monetary amounts are combined before
#' exponentiating when necessary. The annual payment model is independent of
#' fractional interpolation; frac is retained as an audit assumption.
#'
#' The following fully discrete identities remain checks on the result:
#'
#' For a whole-life contract:
#' \deqn{A = 1 - d \ddot{a}.}
#'
#' For an n-year term insurance:
#' \deqn{A^1_{:\overline{n}|} =
#' 1 - d\ddot{a}_{:\overline{n}|} - v^n\,{}_np.}
#'
#' For an n-year endowment insurance:
#' \deqn{A_{:\overline{n}|} =
#' 1 - d\ddot{a}_{:\overline{n}|}.}
#'
#' Whole-life valuation requires observed mortality through extinction of the
#' selected status: at least one marginal must reach zero within the joint
#' horizon, and both must reach zero for last-survivor insurance. An open
#' table cannot imply a terminal death or maturity benefit. Use an explicit
#' term or endowment for the observed horizon, or extend the table.
#' No missing annual interval of an open marginal table is extrapolated.
#' A zero-term endowment pays benefit at h if the status survives; a zero-term
#' term insurance has no death coverage and pays zero. Annual tables must
#' satisfy the consistency requirements of \code{\link{t_pxy}}, and issue
#' ages must have positive lx. One call values one contract, not scenarios.
#'
#' A deferred insurance retains survival probabilities from the issue ages.
#' Its continuation value conditional on the selected status surviving to
#' \code{h} is multiplied by the deferment factor
#' \deqn{v^h\,{}_hp_{xy}}
#' for the selected two-life status.
#' For last-survivor insurance, this conditioning includes paths with only
#' one survivor at \code{h}; it does not restart both lives at deferred ages.
#' In tidy output, \code{annuity_due_value} is the annuity at time \code{h}
#' conditional on that status surviving to \code{h}, and
#' \code{survival_to_maturity}, when reported, uses the same conditioning.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric scalar.
#'
#' If \code{tidy = TRUE}, a one-row tibble with input values, standardized
#' interest rate, deferment factor, unit APV, and APV.
#'
#' @seealso \code{\link{annuity_xy}}, \code{\link{premium_xy}},
#'   \code{\link{insurance_x}}, \code{\link{t_pxy}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:110,
#'   lx = seq(100000, 0, length.out = 51)
#' )
#'
#' insurance_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.06,
#'   type = "whole",
#'   status = "joint"
#' )
#'
#' insurance_xy(
#'   lt = lt,
#'   x = 60,
#'   y = 62,
#'   i = 0.06,
#'   type = "term",
#'   status = "last",
#'   n = 10,
#'   tidy = TRUE
#' )
#'
#' lt |>
#'   life_contract(lives = "joint", x = 60, y = 62, i = 0.06) |>
#'   insurance_xy(
#'     type = "term",
#'     n = 4,
#'     status = "joint"
#'   )
#'
#' @export
insurance_xy <- function(
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
    frac,
    tidy = FALSE,
    ...
) {
  dots <- list(...)
  status_missing <- missing(status)

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

  if (!is.null(dots$insurance_type)) {
    type <- dots$insurance_type
  }

  if (!is.null(dots$cohort)) {
    if (!status_missing) {
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

  if (!is.null(dots$output)) {
    if (!identical(tidy, FALSE)) {
      stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
    }

    output <- match.arg(dots$output, c("value", "table"))
    tidy <- identical(output, "table")
  }

  type <- match.arg(type)
  status <- match.arg(status)

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
      stop("`insurance_xy()` requires a two-life `life_contract()` object.", call. = FALSE)
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
    table_use <- lt_x
  } else if (
    is.list(lt) &&
    length(lt) == 2L &&
    all(vapply(lt, is.data.frame, logical(1L)))
  ) {
    lt_x <- validate_lifetable(lt[[1L]], "lt[[1]]")
    lt_y <- validate_lifetable(lt[[2L]], "lt[[2]]")
    table_use <- list(lt_x, lt_y)
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
    if (identical(frac_x, "CML")) frac_x <- "CF"
    if (identical(frac_y, "CML")) frac_y <- "CF"
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
      abs(x) > .Machine$integer.max ||
      abs(x - round(x)) > 1e-10) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (is.null(y) ||
      !is.numeric(y) ||
      length(y) != 1L ||
      is.na(y) ||
      !is.finite(y) ||
      abs(y) > .Machine$integer.max ||
      abs(y - round(y)) > 1e-10) {
    stop("`y` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      is.na(m) ||
      !is.finite(m) ||
      m < 1 ||
      abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(h) ||
      length(h) != 1L ||
      is.na(h) ||
      !is.finite(h) ||
      h < 0 ||
      abs(h - round(h)) > 1e-10) {
    stop("`h` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n < 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) || abs(n - round(n)) > 1e-10))) {
    stop("`n` must be `Inf` or a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit) ||
      benefit < 0) {
    stop("`benefit` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (type %in% c("term", "endowment") && is.infinite(n)) {
    stop(
      "`n` must be finite for term and endowment insurance.",
      call. = FALSE
    )
  }

  x <- as.integer(round(x))
  y <- as.integer(round(y))
  store_integer <- function(z) {
    z <- round(z)
    if (abs(z) <= .Machine$integer.max) as.integer(z) else z
  }
  m <- store_integer(m)
  h <- store_integer(h)

  if (!is.infinite(n)) {
    n <- store_integer(n)
  }

  if (!x %in% as.integer(round(lt_x$x))) {
    stop("`x` must be present in the first life table.", call. = FALSE)
  }

  if (!y %in% as.integer(round(lt_y$x))) {
    stop("`y` must be present in the second life table.", call. = FALSE)
  }

  # Annual transition masses and original log interest are the valuation basis.
  basis <- .interest_basis(i_type = i_type, i = i, m = m)
  i_effective <- basis$effective
  growth <- basis$growth
  d_effective <- -expm1(-growth)
  state <- .two_life_insurance_basis(lt_x, lt_y, x, y, status)
  if (type == "whole" && !state$complete)
    stop("Whole-life insurance requires observed mortality through extinction of the selected status; extend the open table or use an explicit term/endowment.", call. = FALSE)
  log_h <- state$survival(h)
  deferment_factor <- if (log_h == -Inf) 0 else exp(log_h - growth * h)
  annuity_due <- NA_real_
  survival_to_maturity <- NA_real_
  unit_apv <- result <- 0

  if (log_h != -Inf) {
    years <- if (type == "whole") max(0, state$horizon - h) else n
    if (h + years > state$horizon)
      stop("`n` exceeds the horizon supported by the two life tables after deferment.", call. = FALSE)
    if (years == 0 && type != "whole") {
      # C14 authorized: an empty endowment matures at the deferred start.
      annuity_due <- 0
      survival_to_maturity <- 1
      if (type == "endowment") {
        unit_apv <- deferment_factor
        result <- .cash_flow_value(benefit, log_h - growth * h)
      }
    } else {
      times <- if (years == 0) numeric() else h + seq_len(years) - 1
      log_s <- vapply(times, state$survival, 0)
      annuity_due <- sum(vapply(seq_along(times), function(j) {
        if (log_s[j] == -Inf) 0 else exp(log_s[j] - log_h - growth * (times[j] - h))
      }, 0))
      log_mass <- vapply(times, state$death, 0)
      log_factors <- log_mass - growth * (times + 1)
      # A zero mass is exactly zero, even if the discount factor overflows.
      log_factors[log_mass == -Inf] <- -Inf
      maturity <- state$survival(h + years)
      if (type == "term") survival_to_maturity <- exp(maturity - log_h)
      if (type == "endowment") {
        terminal <- if (maturity == -Inf) -Inf else maturity - growth * (h + years)
        log_factors <- c(log_factors, terminal)
      }
      unit_apv <- sum(.cash_flow_value(rep(1, length(log_factors)), log_factors))
      result <- sum(.cash_flow_value(rep(benefit, length(log_factors)), log_factors))
    }
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
    d_effective = d_effective,
    type = type,
    insurance_type = type,
    status = status,
    cohort = if (status == "joint") "first" else "last",
    n = if (type == "whole") Inf else n,
    term_years = if (type == "whole") Inf else n,
    h = h,
    deferment_years = h,
    benefit = benefit,
    frac = frac,
    deferment_factor = deferment_factor,
    annuity_due_value = annuity_due,
    survival_to_maturity = survival_to_maturity,
    unit_apv = unit_apv,
    apv = result
  )
}
