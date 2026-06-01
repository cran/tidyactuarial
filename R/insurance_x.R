#' Actuarial present value of a life insurance
#'
#' Computes the actuarial present value of a discrete single-life insurance
#' using a life table and compact actuarial notation.
#'
#' The benefit is paid at the end of the year of death for whole-life and term
#' insurance. For endowment insurance, the same benefit is paid either at death
#' within the term or at the end of the term if the life survives.
#'
#' Supported contracts:
#' \itemize{
#'   \item \code{"whole"}: whole-life insurance.
#'   \item \code{"term"}: n-year term insurance.
#'   \item \code{"endowment"}: n-year endowment insurance.
#' }
#'
#' @param lt A life table as produced by \code{\link{lifetable}}. It must
#'   contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age at issue.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param n Integer insurance term in years. Required for
#'   \code{type = "term"} and \code{type = "endowment"}. Use \code{Inf} only
#'   for whole-life insurance.
#' @param h Integer deferment period in years.
#' @param type Character string. One of \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param benefit Numeric scalar. Benefit amount.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric APV. If
#'   \code{TRUE}, returns a one-row tibble with intermediate quantities.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{term_years}, \code{deferral_years}, \code{insurance_type}, and
#'   \code{output}.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric scalar containing the actuarial present
#' value.
#'
#' If \code{tidy = TRUE}, a one-row tibble with the main input values,
#' equivalent interest rate, deferral factor, pure endowment factor,
#' annuity-due value used in the standard identities, and APV.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age at
#' issue, \code{i} is the interest-rate input, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' \code{n} is the insurance term, and \code{h} is the deferment period.
#'
#' The function computes APVs directly from \code{lx}. For a deferred
#' insurance, the value at age \code{x + h} is multiplied by
#' \deqn{v^h\,{}_hp_x.}
#'
#' For whole-life insurance, the death benefit APV at the deferred starting age
#' is computed over the available life-table horizon. For term and endowment
#' insurance, the death benefit is computed over the first \code{n} years.
#' Endowment insurance additionally includes the pure endowment component
#' \deqn{v^n\,{}_np_x.}
#'
#' @seealso \code{\link{annuity_x}}, \code{\link{premium_x}},
#'   \code{\link{reserve_x}}, \code{\link{t_Ex}}, \code{\link{insurance_xy}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:65,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000)
#' )
#'
#' # Whole-life insurance
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   type = "whole"
#' )
#'
#' # 5-year term insurance
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 5,
#'   type = "term"
#' )
#'
#' # 5-year endowment insurance
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 5,
#'   type = "endowment"
#' )
#'
#' # Deferred whole-life insurance
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   h = 2,
#'   type = "whole"
#' )
#'
#' # Tidy output
#' insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 5,
#'   type = "term",
#'   tidy = TRUE
#' )
#'
#' @export
insurance_x <- function(
    lt,
    x,
    i,
    i_type = "effective",
    m = 1L,
    n = Inf,
    h = 0L,
    type = c("whole", "term", "endowment"),
    benefit = 1,
    tidy = FALSE,
    ...
) {
  dots <- list(...)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "mortality_table",
    "age",
    "rate",
    "rate_type",
    "term_years",
    "deferral_years",
    "insurance_type",
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

  if (!is.null(dots$age)) {
    if (!missing(x)) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }
    x <- dots$age
  }

  if (!is.null(dots$rate)) {
    if (!missing(i)) {
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

  if (!is.null(dots$term_years)) {
    if (!is.infinite(n)) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }
    n <- dots$term_years
  }

  if (!is.null(dots$deferral_years)) {
    if (!identical(h, 0L) && !identical(h, 0)) {
      stop("Provide only one of `h` or deprecated `deferral_years`.", call. = FALSE)
    }
    h <- dots$deferral_years
  }

  if (!is.null(dots$insurance_type)) {
    type <- dots$insurance_type
  }

  if (!is.null(dots$output)) {
    if (!identical(tidy, FALSE)) {
      stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
    }

    output <- match.arg(dots$output, c("value", "table"))
    tidy <- identical(output, "table")
  }

  type <- match.arg(type)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (!missing(lt) &&
      exists(".as_life_contract", mode = "function") &&
      .as_life_contract(lt)) {
    contract <- lt

    if (!identical(contract$lives, "single")) {
      stop(
        "`insurance_x()` currently supports only single-life `life_contract()` objects.",
        call. = FALSE
      )
    }

    lt <- contract$mortality_table

    if (missing(x) || is.null(x)) {
      x <- contract$age %||% contract$x
    }

    if (missing(i) || is.null(i)) {
      i <- contract$rate %||% contract$i
    }

    if (identical(i_type, "effective")) {
      i_type <- contract$rate_type %||% contract$i_type %||% i_type
    }

    if (identical(m, 1L) || identical(m, 1)) {
      m <- contract$m %||% m
    }
  }

  # -------------------------------------------------------------------------
  # Basic validation
  # -------------------------------------------------------------------------

  if (missing(lt)) {
    stop("`lt` must be provided.", call. = FALSE)
  }

  if (missing(x)) {
    stop("`x` must be provided.", call. = FALSE)
  }

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  if (!is.data.frame(lt)) {
    stop("`lt` must be a data.frame or tibble.", call. = FALSE)
  }

  if (!all(c("x", "lx") %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
  }

  if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      abs(x - round(x)) > 1e-10) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (!is.numeric(i) ||
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

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit)) {
    stop("`benefit` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n < 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) || abs(n - round(n)) > 1e-10))) {
    stop(
      "`n` must be `Inf` or a single nonnegative integer.",
      call. = FALSE
    )
  }

  if (type %in% c("term", "endowment") && is.infinite(n)) {
    stop("`n` must be finite for term and endowment insurance.", call. = FALSE)
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  h <- as.integer(round(h))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
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
  # Life table preparation
  # -------------------------------------------------------------------------

  lt <- lt[order(lt$x), , drop = FALSE]

  if (!is.numeric(lt$x)) {
    stop("Column `x` in `lt` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(lt$lx)) {
    stop("Column `lx` in `lt` must be numeric.", call. = FALSE)
  }

  if (any(is.na(lt$x)) || any(!is.finite(lt$x))) {
    stop("Column `x` must contain finite non-missing values.", call. = FALSE)
  }

  if (any(abs(lt$x - round(lt$x)) > 1e-10)) {
    stop("Column `x` must contain integer ages.", call. = FALSE)
  }

  if (anyDuplicated(lt$x)) {
    stop("Life table ages in column `x` must be unique.", call. = FALSE)
  }

  if (any(is.na(lt$lx)) || any(!is.finite(lt$lx)) || any(lt$lx < 0)) {
    stop("Column `lx` must contain finite nonnegative values.", call. = FALSE)
  }

  ages <- as.integer(round(lt$x))
  lx <- as.numeric(lt$lx)
  omega <- max(ages)

  get_lx <- function(current_age) {
    idx <- match(current_age, ages)

    if (!is.na(idx)) {
      return(lx[[idx]])
    }

    if (current_age == omega + 1L) {
      return(0)
    }

    NA_real_
  }

  t_p_int <- function(current_age, tt) {
    if (tt == 0L) {
      return(1)
    }

    l0 <- get_lx(current_age)
    l1 <- get_lx(current_age + tt)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    l1 / l0
  }

  one_year_qx <- function(current_age) {
    l0 <- get_lx(current_age)
    l1 <- get_lx(current_age + 1L)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    1 - l1 / l0
  }

  annuity_due_from <- function(current_age, years) {
    years <- as.integer(round(years))

    if (years <= 0L) {
      return(0)
    }

    sum(vapply(0:(years - 1L), function(r) {
      px_r <- t_p_int(current_age, r)

      if (is.na(px_r)) {
        stop("The life table does not support the requested annuity APV.", call. = FALSE)
      }

      v_fun(r) * px_r
    }, numeric(1L)))
  }

  death_apv_from <- function(current_age, years) {
    years <- as.integer(round(years))

    if (years <= 0L) {
      return(0)
    }

    sum(vapply(0:(years - 1L), function(r) {
      px_r <- t_p_int(current_age, r)
      qx_r <- one_year_qx(current_age + r)

      if (is.na(px_r) || is.na(qx_r)) {
        stop("The life table does not support the requested insurance APV.", call. = FALSE)
      }

      v_fun(r + 1L) * px_r * qx_r
    }, numeric(1L)))
  }

  # -------------------------------------------------------------------------
  # Deferral
  # -------------------------------------------------------------------------

  deferment_factor <- v_fun(h) * t_p_int(x, h)

  if (is.na(deferment_factor)) {
    stop(
      "The deferral age `x + h` is outside the life table or `lx(x)` is zero.",
      call. = FALSE
    )
  }

  start_x <- x + h

  # If deferment reaches a point where survival is zero, the APV is zero.
  if (deferment_factor == 0) {
    result <- 0

    if (!tidy) {
      return(result)
    }

    return(tibble::tibble(
      x = x,
      age = x,
      i = i,
      rate = i,
      i_type = i_type,
      rate_type = i_type,
      m = m,
      i_effective = i_effective,
      d_effective = d_effective,
      n = n,
      term_years = n,
      h = h,
      deferral_years = h,
      start_x = start_x,
      start_age = start_x,
      type = type,
      insurance_type = type,
      benefit = benefit,
      deferment_factor = deferment_factor,
      annuity_due_value = NA_real_,
      pure_endowment_factor = NA_real_,
      value_at_start = 0,
      apv = result
    ))
  }

  if (start_x > omega) {
    stop(
      "The deferred starting age `x + h` is outside the available life table.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------------
  # Insurance value at the deferred starting age
  # -------------------------------------------------------------------------

  pure_endowment_factor <- NA_real_

  contract_years <- if (type == "whole") {
    omega - start_x + 1L
  } else {
    n
  }

  if (contract_years < 0L) {
    stop("The requested contract term is invalid.", call. = FALSE)
  }

  if (type %in% c("term", "endowment") &&
      start_x + contract_years > omega + 1L) {
    stop(
      "The requested term exceeds the available life table horizon.",
      call. = FALSE
    )
  }

  if (type == "whole" && contract_years == 0L) {
    value_at_start <- 0
    annuity_due_value <- 0
  } else if (type %in% c("term", "endowment") && contract_years == 0L) {
    # Preserve previous behavior for zero-year term/endowment insurance.
    value_at_start <- 0
    annuity_due_value <- 0
    pure_endowment_factor <- 1
  } else {
    death_apv <- death_apv_from(
      current_age = start_x,
      years = contract_years
    )

    annuity_due_value <- annuity_due_from(
      current_age = start_x,
      years = contract_years
    )

    if (type == "endowment") {
      n_p_start <- t_p_int(start_x, contract_years)

      if (is.na(n_p_start)) {
        stop(
          "The life table does not support the survival probability needed ",
          "for this endowment insurance calculation.",
          call. = FALSE
        )
      }

      pure_endowment_factor <- v_fun(contract_years) * n_p_start
      value_at_start <- death_apv + pure_endowment_factor
    } else {
      value_at_start <- death_apv
    }
  }

  result <- benefit * deferment_factor * value_at_start

  if (!tidy) {
    return(result)
  }

  tibble::tibble(
    x = x,
    age = x,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective,
    d_effective = d_effective,
    n = n,
    term_years = n,
    h = h,
    deferral_years = h,
    start_x = start_x,
    start_age = start_x,
    type = type,
    insurance_type = type,
    benefit = benefit,
    deferment_factor = deferment_factor,
    annuity_due_value = annuity_due_value,
    pure_endowment_factor = pure_endowment_factor,
    value_at_start = value_at_start,
    apv = result
  )
}
