#' Actuarial present value of a life annuity
#'
#' Computes the actuarial present value of a discrete life annuity using compact
#' actuarial notation.
#'
#' The function supports:
#' \itemize{
#'   \item whole-life annuities,
#'   \item temporary annuities,
#'   \item integer deferment,
#'   \item annual or k-thly payments,
#'   \item exact fractional survival for k-thly payments,
#'   \item first- and second-order Woolhouse approximations.
#' }
#'
#' @param lt A life table as produced by \code{\link{lifetable}}, or a
#'   \code{tidyact_life_contract} object created by \code{\link{life_contract}}.
#'   A life table must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age, or a vector of scenario ages. Optional when \code{lt} is a single-life
#'   \code{tidyact_life_contract}.
#' @param i Numeric annual interest-rate input, scalar or by scenario. Optional when \code{lt}
#'   is a single-life \code{tidyact_life_contract}.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal interest rates.
#'   Ignored for \code{i_type = "effective"} and \code{i_type = "force"}.
#'   In \code{tidyactuarial}, \code{m} is reserved for interest conversion
#'   frequency, not deferment.
#' @param n Numeric term in years. Use \code{NULL} or \code{Inf} for a
#'   whole-life annuity. For exact k-thly valuation
#'   (\code{woolhouse = "none"}), fractional terms are allowed when
#'   \code{n * k} is an integer. Woolhouse approximations require an
#'   integer number of years.
#' @param h Integer deferment period in years.
#' @param k Positive integer. Number of annuity payments per year. For example,
#'   use \code{k = 12} for monthly payments. The annuity is normalized to
#'   an annual payment rate of 1, so each individual payment has amount
#'   \code{1 / k}.
#' @param timing Character string. Either \code{"immediate"} for payments at
#'   the end of each payment period or \code{"due"} for payments at the
#'   beginning of each payment period.
#' @param woolhouse Character string. For \code{k > 1}, use \code{"none"} for
#'   exact fractional-age computation, \code{"first"} for the first-order
#'   Woolhouse approximation, or \code{"second"} for the second-order
#'   Woolhouse approximation.
#' @param frac Character string. Fractional-age assumption used when
#'   \code{k > 1} and \code{woolhouse = "none"}. Allowed values are
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, and \code{"Balducci"}. If
#'   \code{NULL}, the \code{frac} attribute of \code{lt} is used when available;
#'   otherwise \code{"UDD"} is used.
#' @param tidy Logical. If \code{FALSE}, returns a numeric APV. If \code{TRUE},
#'   returns one tibble row per scenario with intermediate quantities.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of actuarial present values,
#' of length one for a single scenario.
#'
#' If \code{tidy = TRUE}, a tibble with one row per scenario, the main input values,
#' equivalent interest rate, deferment factor, pure endowment factor, and APV.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}:
#'
#' \itemize{
#'   \item \code{lt}: life table;
#'   \item \code{x}: actuarial age;
#'   \item \code{i}: interest rate;
#'   \item \code{i_type}: interest-rate type;
#'   \item \code{m}: interest conversion frequency;
#'   \item \code{n}: annuity term;
#'   \item \code{h}: deferment period;
#'   \item \code{k}: payment frequency.
#' }
#'
#' For annual annuities-due,
#' \deqn{\ddot{a}_{x:\overline{n}|} =
#' \sum_{j=0}^{n-1} v^j\,{}_jp_x.}
#'
#' For annual annuities-immediate,
#' \deqn{a_{x:\overline{n}|} =
#' \sum_{j=1}^{n} v^j\,{}_jp_x.}
#'
#' Deferment is handled through
#' \deqn{v^h\,{}_hp_x,}
#' where \eqn{h} is the deferment period.
#'
#' For k-thly payments with \code{woolhouse = "none"}, fractional survival is
#' computed under the selected fractional-age assumption. The payment stream
#' is normalized to an annual payment rate of 1:
#' \deqn{\ddot{a}_{x:\overline{n}|}^{(k)} =
#' \frac{1}{k}\sum_{j=0}^{kn-1}
#' v^{j/k}\,{}_{j/k}p_x.}
#' Consequently, each k-thly installment has amount \eqn{1/k}. This
#' normalization is essential when the function is used as the premium
#' annuity in the equivalence principle: the resulting premium is annualized,
#' and the amount paid at each installment is the annualized premium divided
#' by \eqn{k}.
#'
#' \strong{Scenarios and tidy workflows.} All arguments except \code{lt} and
#' \code{tidy} may have a common positive length; only length-one arguments
#' recycle. \code{n = NULL} and \code{frac = NULL} keep their scalar default
#' meanings. Scenarios are paired by position, not expanded to a Cartesian
#' product, and each describes one life under the shared table. Invalid or
#' missing inputs are rejected with the scenario number; no partial result is
#' returned. Use \code{dplyr::mutate()} with \code{tidy = FALSE} for portfolios,
#' or \code{tidy = TRUE} for an auditable scenario table. Different mortality
#' tables require separate calls, for example in a list-column workflow.
#'
#' \strong{Model and numerical limits.} Discounting retains the original log
#' accumulation rate; exact payments are valued from the issue age, avoiding
#' products of separately rounded deferment and conditional values. Intermediate
#' diagnostic values can be zero or infinite while the APV remains finite.
#' Tables must have an annual age grid and nonincreasing finite \code{lx},
#' with a positive count at issue. For compatibility, a count of zero is assumed
#' at one year beyond the last recorded age, even for an open input table.
#' Known extinct tails contribute zero. Woolhouse is an approximation using
#' annual survival estimates of boundary forces, not an exact integral; its
#' second order requires positive survival where a boundary force is needed.
#' Negative approximations generate a warning and are not clamped. Prefer
#' \code{woolhouse = "none"} for exact payments under the chosen fractional
#' assumption. At \code{n = 0}, APV and the annuity value at start are zero,
#' while the zero-term pure endowment factor is one.
#'
#' @seealso \code{\link{insurance_x}}, \code{\link{premium_x}},
#'   \code{\link{reserve_x}}, \code{\link{t_px}}, \code{\link{t_Ex}}
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:65,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000)
#' )
#'
#' # Annual annuity-immediate
#' annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   timing = "immediate"
#' )
#'
#' # Annual annuity-due
#' annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   timing = "due"
#' )
#'
#' # Temporary annuity
#' annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 3,
#'   timing = "due"
#' )
#'
#' # Deferred annuity
#' annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   h = 2,
#'   timing = "due"
#' )
#'
#' # Tidy output
#' annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   n = 3,
#'   timing = "due",
#'   tidy = TRUE
#' )
#'
#' # Pipe workflow with a life contract
#' life_contract(lt = lt, lives = "single", x = 60, i = 0.06) |>
#'   annuity_x(n = 3, timing = "due")
#'
#' # Parallel single-life scenarios, suitable for mutate()
#' policies <- tibble::tibble(age = c(60, 61, 62), rate = c(0.03, 0.04, 0.05))
#' dplyr::mutate(policies,
#'   apv = annuity_x(lt, x = age, i = rate, n = 2, k = 12, timing = "due"))
#' annuity_x(lt, x = c(60, 61), i = 0.04, n = c(2, 3), tidy = TRUE)
#'
#' @export
annuity_x <- function(
    lt,
    x,
    i,
    i_type = "effective",
    m = 1L,
    n = NULL,
    h = 0L,
    k = 1L,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    frac = NULL,
    tidy = FALSE
) {
  if (missing(timing)) timing <- "immediate"
  if (missing(woolhouse)) woolhouse <- "none"

  # -------------------------------------------------------------------------
  # Pipe support: allow a tidyact_life_contract as first argument
  # -------------------------------------------------------------------------

  if (.as_life_contract(lt)) {
    contract <- lt

    if (!identical(contract$lives, "single")) {
      stop(
        "`annuity_x()` supports only single-life `life_contract()` objects.",
        call. = FALSE
      )
    }

    lt <- contract$lt

    if (missing(x) || is.null(x)) {
      x <- contract$x
    }

    if (missing(i) || is.null(i)) {
      i <- contract$i
    }

    if (missing(i_type) || is.null(i_type)) {
      i_type <- contract$i_type
    }

    if (missing(m) || is.null(m)) {
      m <- contract$m
    }
  }

  # -------------------------------------------------------------------------
  # Basic validation
  # -------------------------------------------------------------------------

  if (!missing(x) && !missing(i)) {
    scenarios <- .single_life_scenarios(annuity_x, lt,
      list(x = x, i = i, i_type = i_type, m = m, n = n, h = h, k = k,
           timing = timing, woolhouse = woolhouse, frac = frac), tidy)
    if (!is.null(scenarios)) return(scenarios)
  }
  timing <- match.arg(timing, c("immediate", "due"))
  woolhouse <- match.arg(woolhouse, c("none", "first", "second"))

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (!is.data.frame(lt)) {
    stop("`lt` must be a data frame, tibble, or single-life contract.", call. = FALSE)
  }

  if (!all(c("x", "lx") %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
  }

  if (missing(x) ||
      !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      abs(x - round(x)) > 1e-10 || abs(x) > .Machine$integer.max) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (missing(i) ||
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
      abs(h - round(h)) > 1e-10 || h > .Machine$integer.max) {
    stop("`h` must be a single nonnegative integer.", call. = FALSE)
  }

  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k < 1 ||
      abs(k - round(k)) > 1e-10) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  if (!is.null(n) &&
      (!is.numeric(n) ||
       length(n) != 1L ||
       is.na(n) ||
       n < 0 ||
       (!is.infinite(n) && !is.finite(n)))) {
    stop(
      "`n` must be `NULL`, `Inf`, or a single nonnegative finite value.",
      call. = FALSE
    )
  }

  if (is.null(frac)) {
    frac <- attr(lt, "frac")

    if (is.null(frac) || length(frac) != 1L || is.na(frac)) {
      frac <- "UDD"
    }
  }

  if (!is.character(frac) || length(frac) != 1L || is.na(frac)) {
    stop("`frac` must be a single character string.", call. = FALSE)
  }

  frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))

  x <- as.integer(round(x))
  m <- if (m <= .Machine$integer.max) as.integer(round(m)) else round(m)
  h <- as.integer(round(h))
  k <- if (k <= .Machine$integer.max) as.integer(round(k)) else round(k)

  if (!is.null(n) && !is.infinite(n)) {
    if (identical(woolhouse, "none")) {
      n_payments_check <- n * k

      if (!is.finite(n_payments_check) ||
          n_payments_check > .Machine$integer.max) {
        stop("The number of payments exceeds the supported integer range.", call. = FALSE)
      }
      if (abs(n_payments_check - round(n_payments_check)) > 1e-10) {
        stop(
          "For exact k-thly valuation, `n * k` must be an integer so that ",
          "the term contains a whole number of payments.",
          call. = FALSE
        )
      }

      n <- round(n_payments_check) / k
    } else {
      if (abs(n - round(n)) > 1e-10) {
        stop(
          "Woolhouse approximations require `n` to be an integer number of years.",
          call. = FALSE
        )
      }

      n <- if (n <= .Machine$integer.max) as.integer(round(n)) else round(n)
    }
  }

  # -------------------------------------------------------------------------
  # Interest conversion
  # -------------------------------------------------------------------------

  basis <- .interest_basis(i_type, i, m)
  i_effective <- basis$effective

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
    stop("Column `x` in `lt` must contain finite non-missing values.", call. = FALSE)
  }

  if (!nrow(lt) || any(abs(lt$x) > .Machine$integer.max) ||
      any(abs(lt$x - round(lt$x)) > 1e-10)) {
    stop("Column `x` in `lt` must contain integer ages.", call. = FALSE)
  }

  if (anyDuplicated(lt$x)) {
    stop("Life table ages in column `x` must be unique.", call. = FALSE)
  }

  if (any(is.na(lt$lx)) || any(!is.finite(lt$lx)) || any(lt$lx < 0)) {
    stop("Column `lx` in `lt` must contain finite nonnegative values.", call. = FALSE)
  }

  ages <- as.integer(round(lt$x))
  lx <- as.numeric(lt$lx)
  if (any(diff(lx) > 0)) {
    stop("Column `lx` must be nonincreasing.", call. = FALSE)
  }
  if (any(diff(ages) != 1)) {
    stop("Life table ages must form an annual age grid.", call. = FALSE)
  }
  omega <- max(ages)

  get_lx <- function(current_age) {
    idx <- match(current_age, ages)

    if (!is.na(idx)) {
      return(lx[[idx]])
    }

    if (current_age == as.double(omega) + 1) {
      return(0)
    }

    NA_real_
  }

  log_p_int <- function(current_age, tt) {
    if (tt == 0) {
      return(0)
    }

    l0 <- get_lx(current_age)
    l1 <- get_lx(as.double(current_age) + tt)

    if (is.na(l0) || is.na(l1) || l0 <= 0) {
      return(NA_real_)
    }

    .single_life_log_ratio(l1, l0)
  }

  log_p_frac <- function(current_age, u) {
    if (u < 0) {
      return(NA_real_)
    }

    if (u == 0) {
      return(0)
    }

    tt <- floor(u)
    s <- u - tt

    pt <- log_p_int(current_age, tt)

    if (is.na(pt)) {
      return(NA_real_)
    }

    if (s == 0) {
      return(pt)
    }
    if (pt == -Inf) return(-Inf)

    y_age <- current_age + tt

    ly <- get_lx(y_age)
    ly1 <- get_lx(y_age + 1L)

    if (is.na(ly) || is.na(ly1) || ly <= 0) {
      return(NA_real_)
    }

    p_y <- ly1 / ly
    q_y <- 1 - p_y

    log_p_s <- switch(
      frac,
      UDD = log1p(-s * q_y),
      CF = s * .single_life_log_ratio(ly1, ly),
      CML = s * .single_life_log_ratio(ly1, ly),
      Balducci = .single_life_log_ratio(ly1, ly) - log(p_y + s * q_y)
    )

    pt + log_p_s
  }

  # -------------------------------------------------------------------------
  # Deferment
  # -------------------------------------------------------------------------

  if (is.na(get_lx(x)) || get_lx(x) <= 0) {
    stop("The issue age `x` must have positive `lx` in the life table.", call. = FALSE)
  }
  log_survival_h <- log_p_int(x, h)
  log_deferment <- if (isTRUE(log_survival_h == -Inf)) -Inf else
    log_survival_h - basis$growth * h
  deferment_factor <- exp(log_deferment)

  if (is.na(deferment_factor)) {
    stop(
      "The deferred age `x + h` is outside the life table or `lx(x)` is zero.",
      call. = FALSE
    )
  }

  start_age <- as.double(x) + h
  if (abs(start_age) <= .Machine$integer.max) start_age <- as.integer(start_age)

  # -------------------------------------------------------------------------
  # Term
  # -------------------------------------------------------------------------

  max_years <- max(0, (as.double(omega) + 1) - start_age)
  if (max_years <= .Machine$integer.max) max_years <- as.integer(max_years)

  if (is.null(n) || is.infinite(n)) {
    n_used <- max_years
    n_display <- if (is.null(n)) Inf else n
  } else {
    n_used <- n
    n_display <- n

    if (n_used > max_years) {
      stop(
        "`n` exceeds the horizon allowed by the life table. ",
        "The table must support ages up to `x + h + n`.",
        call. = FALSE
      )
    }
  }

  if (n_used == 0L) {
    result <- 0

    if (!tidy) {
      return(result)
    }

    return(tibble::tibble(
      x = x,
      i = i,
      i_type = i_type,
      m = m,
      i_effective = i_effective,
      n = n_display,
      n_used = n_used,
      h = h,
      x_h = start_age,
      k = k,
      timing = timing,
      woolhouse = woolhouse,
      frac = frac,
      deferment_factor = deferment_factor,
      pure_endowment_factor = 1,
      annuity_value_at_start = 0,
      apv = result
    ))
  }

  if (n_used * k > .Machine$integer.max) {
    stop("The number of payments exceeds the supported integer range.", call. = FALSE)
  }
  log_endowment_survival <- log_p_frac(start_age, n_used)
  log_endowment <- if (isTRUE(log_endowment_survival == -Inf)) -Inf else
    log_endowment_survival - basis$growth * n_used
  pure_endowment_factor <- exp(log_endowment)

  if (is.na(pure_endowment_factor)) {
    pure_endowment_factor <- 0
  }

  # -------------------------------------------------------------------------
  # Annual and k-thly computation helpers
  # -------------------------------------------------------------------------

  annual_exact <- function(current_age, nn, tim, offset = 0) {
    if (tim == "due") {
      times <- 0:(nn - 1L)
    } else {
      times <- 1:nn
    }

    survival <- vapply(
      times,
      function(tt) log_p_int(current_age, offset + tt),
      numeric(1L)
    )

    if (anyNA(survival)) {
      stop(
        "The life table does not support the required ages for annual payments.",
        call. = FALSE
      )
    }

    log_values <- survival - basis$growth * (offset + times)
    log_values[survival == -Inf] <- -Inf
    sum(exp(log_values))
  }

  kthly_exact <- function(current_age, nn, kk, tim, offset = 0) {
    n_payments <- as.integer(round(kk * nn))

    if (tim == "due") {
      j <- 0:(n_payments - 1L)
    } else {
      j <- 1:n_payments
    }

    u <- j / kk

    survival <- vapply(
      u,
      function(uu) log_p_frac(current_age, offset + uu),
      numeric(1L)
    )

    if (anyNA(survival)) {
      stop(
        "The life table does not support the required ages for k-thly payments.",
        call. = FALSE
      )
    }

    log_values <- survival - basis$growth * (offset + u)
    log_values[survival == -Inf] <- -Inf
    sum(.cash_flow_value(rep(1 / kk, length(u)), log_values))
  }

  # -------------------------------------------------------------------------
  # Main APV computation
  # -------------------------------------------------------------------------

  if (log_survival_h == -Inf) {
    annuity_value_at_start <- NA_real_
    result <- 0
  } else if (k == 1L) {
    annuity_value_at_start <- annual_exact(
      current_age = start_age,
      nn = n_used,
      tim = timing
    )
    result <- annual_exact(x, n_used, timing, offset = h)
  } else if (woolhouse == "none") {
    annuity_value_at_start <- kthly_exact(
      current_age = start_age,
      nn = n_used,
      kk = k,
      tim = timing
    )
    result <- kthly_exact(x, n_used, k, timing, offset = h)
  } else {
    annual_due <- annual_exact(
      current_age = start_age,
      nn = n_used,
      tim = "due"
    )

    endowment_loss <- if (is.na(log_endowment)) 1 else -expm1(log_endowment)
    adj1 <- (k - 1) / (2 * k) * endowment_loss

    if (woolhouse == "first") {
      due_k <- annual_due - adj1
    } else {
      delta <- basis$growth

      ly <- get_lx(start_age)
      ly1 <- get_lx(as.double(start_age) + 1)

      if (is.na(ly) || is.na(ly1) || ly <= 0) {
        stop("Cannot compute the force approximation at the starting age.", call. = FALSE)
      }

      if (ly1 == 0) {
        stop("Second-order Woolhouse requires positive one-year survival; ",
             "use `woolhouse = 'none'` for exact fractional valuation.", call. = FALSE)
      }
      mu_y <- -.single_life_log_ratio(ly1, ly)

      lyn <- get_lx(as.double(start_age) + n_used)
      lyn1 <- get_lx(as.double(start_age) + n_used + 1)

      if (!is.na(lyn) && !is.na(lyn1) && lyn > 0) {
        if (lyn1 == 0 && pure_endowment_factor > 0) {
          stop("Second-order Woolhouse requires positive endpoint survival; ",
               "use `woolhouse = 'none'` for exact fractional valuation.", call. = FALSE)
        }
        mu_yn <- if (lyn1 > 0) -.single_life_log_ratio(lyn1, lyn) else 0
      } else {
        mu_yn <- 0
      }

      adj2 <- (k^2 - 1) /
        (12 * k^2) *
        (
          mu_y - mu_yn + (delta + mu_yn) * endowment_loss
        )

      due_k <- annual_due - adj1 - adj2
    }

    annuity_value_at_start <- if (timing == "due") {
      due_k
    } else {
      due_k - (1 / k) * endowment_loss
    }
    result <- .cash_flow_value(annuity_value_at_start, log_deferment)
    if (is.na(result)) {
      stop("Woolhouse approximation is undefined; use `woolhouse = 'none'`.", call. = FALSE)
    }
    if (result < 0) {
      warning("Woolhouse approximation is negative; use `woolhouse = 'none'` ",
              "for exact fractional valuation.", call. = FALSE)
    }
  }

  if (!tidy) {
    return(result)
  }

  tibble::tibble(
    x = x,
    i = i,
    i_type = i_type,
    m = m,
    i_effective = i_effective,
    n = n_display,
    n_used = n_used,
    h = h,
    x_h = start_age,
    k = k,
    timing = timing,
    woolhouse = woolhouse,
    frac = frac,
    deferment_factor = deferment_factor,
    pure_endowment_factor = pure_endowment_factor,
    annuity_value_at_start = annuity_value_at_start,
    apv = result
  )
}
