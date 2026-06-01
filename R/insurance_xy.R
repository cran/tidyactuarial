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
#' The function uses the standard fully discrete identities that express
#' insurance values through two-life annuity-due values.
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
#' A deferred insurance is valued by multiplying the value at deferred ages
#' \code{x + h} and \code{y + h} by the deferment factor
#' \deqn{v^h\,{}_hp_{xy}}
#' for the selected two-life status.
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
  m <- as.integer(round(m))
  h <- as.integer(round(h))

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

  d_effective <- i_effective / (1 + i_effective)

  # -------------------------------------------------------------------------
  # Deferment
  # -------------------------------------------------------------------------

  p_def <- t_pxy(
    lt = table_use,
    x = x,
    y = y,
    t = h,
    frac = frac,
    status = status
  )

  if (is.na(p_def)) {
    stop("Cannot compute survival to deferment under the life table.", call. = FALSE)
  }

  deferment_factor <- v_fun(h) * p_def

  start_x <- x + h
  start_y <- y + h

  # If the selected status does not survive to the deferred date, APV is zero.
  if (deferment_factor == 0) {
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
      annuity_due_value = NA_real_,
      survival_to_maturity = NA_real_,
      unit_apv = 0,
      apv = result
    ))
  }

  # -------------------------------------------------------------------------
  # APV computation
  # -------------------------------------------------------------------------

  if (type == "whole") {
    annuity_due <- annuity_xy(
      lt = table_use,
      x = start_x,
      y = start_y,
      i = i,
      i_type = i_type,
      m = m,
      status = status,
      n = Inf,
      h = 0L,
      k = 1L,
      timing = "due",
      woolhouse = "none",
      frac = frac,
      tidy = FALSE
    )

    survival_to_maturity <- NA_real_
    unit_apv <- deferment_factor * (1 - d_effective * annuity_due)
  } else if (n == 0L) {
    annuity_due <- 0
    survival_to_maturity <- 1
    unit_apv <- 0
  } else {
    annuity_due <- annuity_xy(
      lt = table_use,
      x = start_x,
      y = start_y,
      i = i,
      i_type = i_type,
      m = m,
      status = status,
      n = n,
      h = 0L,
      k = 1L,
      timing = "due",
      woolhouse = "none",
      frac = frac,
      tidy = FALSE
    )

    if (type == "endowment") {
      survival_to_maturity <- NA_real_
      unit_apv <- deferment_factor * (1 - d_effective * annuity_due)
    } else {
      survival_to_maturity <- t_pxy(
        lt = table_use,
        x = start_x,
        y = start_y,
        t = n,
        frac = frac,
        status = status
      )

      if (is.na(survival_to_maturity)) {
        stop("Cannot compute n-year survival at deferred ages.", call. = FALSE)
      }

      unit_apv <- deferment_factor *
        (1 - d_effective * annuity_due - v_fun(n) * survival_to_maturity)
    }
  }

  result <- benefit * unit_apv

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
