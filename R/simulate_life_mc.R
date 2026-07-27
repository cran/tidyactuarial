#' Monte Carlo simulation of a life annuity
#'
#' Simulates present values of a discrete single-life annuity from a life table.
#' Annual and subannual payments are supported. For subannual payments, a
#' complete future lifetime is generated inside the year of death under a
#' fractional-age assumption.
#'
#' @param lt A life table or a \code{tidyact_life_contract}. A life table must
#'   contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age. Optional for a life contract.
#' @param i Numeric scalar interest-rate input.
#' @param i_type Interest-rate convention.
#' @param m Positive integer nominal conversion frequency.
#' @param n Positive integer term in years, or \code{Inf} for whole life.
#' @param k Positive integer number of annuity payments per year.
#' @param timing Payment timing: \code{"due"} or \code{"immediate"}.
#' @param payment Nonnegative amount paid at each payment date. To simulate an
#'   annualized payment rate of 1 with \code{k} payments per year, use
#'   \code{payment = 1 / k}.
#' @param frac Fractional-age assumption used to simulate the complete lifetime
#'   within the year of death: \code{"udd"}, \code{"cml"}, or
#'   \code{"balducci"}.
#' @param method Simulation method. Currently only \code{"inverse"}.
#' @param n_sim Positive integer number of simulations.
#' @param seed Optional nonnegative integer seed.
#' @param output \code{"simulations"} or \code{"summary"}.
#' @param ... Transitional compatibility for \code{mortality_table},
#'   \code{age}, \code{rate}, \code{rate_type}, \code{term_years}, and
#'   \code{payments_per_year}.
#'
#' @details
#' The amount \code{payment} is a per-payment amount. If \eqn{N} payments are
#' made at times \eqn{t_j}, the simulated present value is
#' \deqn{
#' \sum_{j=1}^{N} \mathrm{payment}\,v^{t_j}.
#' }
#'
#' For \code{k > 1}, the curtate lifetime alone is insufficient because
#' payments may occur during the year of death. The function therefore
#' simulates a complete lifetime \eqn{T_x=K_x+S} under \code{frac}.
#'
#' @return A tibble with one row per simulation, or a summary from
#'   \code{\link{summary_mc}}. Simulation output includes
#'   \code{payment_per_payment}, \code{payment_annualized}, \code{Kx},
#'   \code{Tx}, \code{n_payments}, \code{present_value}, and
#'   \code{pv_annuity}.
#'
#' @family simulation
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 0)
#' )
#'
#' simulate_annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   n_sim = 25,
#'   seed = 123
#' )
#'
#' simulate_annuity_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   n = 5,
#'   k = 12,
#'   payment = 1 / 12,
#'   frac = "udd",
#'   n_sim = 25,
#'   seed = 123
#' )
#'
#' @export
simulate_annuity_x <- function(
    lt,
    x = NULL,
    i = NULL,
    i_type = NULL,
    m = NULL,
    n = Inf,
    k = 1L,
    timing = c("due", "immediate"),
    payment = 1,
    frac = c("udd", "cml", "balducci"),
    method = c("inverse"),
    n_sim = 10000L,
    seed = NULL,
    output = c("simulations", "summary"),
    ...
) {
  lt_missing <- missing(lt)
  x_missing <- missing(x)
  i_missing <- missing(i)
  i_type_missing <- missing(i_type)
  n_missing <- missing(n)
  k_missing <- missing(k)

  timing <- match.arg(timing)
  frac <- match.arg(frac)
  method <- match.arg(method)
  output <- match.arg(output)

  old <- .life_mc_collect_old_args(
    dots = list(...),
    allowed = c(
      "mortality_table",
      "age",
      "rate",
      "rate_type",
      "term_years",
      "payments_per_year"
    )
  )

  if ("mortality_table" %in% names(old)) {
    if (!lt_missing) {
      stop(
        "Provide only one of `lt` or deprecated `mortality_table`.",
        call. = FALSE
      )
    }

    lt <- old[["mortality_table"]]
    lt_missing <- FALSE
  }

  if ("age" %in% names(old)) {
    if (!x_missing) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }

    x <- old[["age"]]
  }

  if ("rate" %in% names(old)) {
    if (!i_missing) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- old[["rate"]]
  }

  if ("rate_type" %in% names(old)) {
    if (!i_type_missing) {
      stop(
        "Provide only one of `i_type` or deprecated `rate_type`.",
        call. = FALSE
      )
    }

    i_type <- old[["rate_type"]]
  }

  if ("term_years" %in% names(old)) {
    if (!n_missing) {
      stop(
        "Provide only one of `n` or deprecated `term_years`.",
        call. = FALSE
      )
    }

    n <- old[["term_years"]]
  }

  if ("payments_per_year" %in% names(old)) {
    if (!k_missing) {
      stop(
        "Provide only one of `k` or deprecated `payments_per_year`.",
        call. = FALSE
      )
    }

    k <- old[["payments_per_year"]]
  }

  if (lt_missing) {
    stop("`lt` must be provided.", call. = FALSE)
  }

  inputs <- .resolve_life_mc_inputs(
    lt = lt,
    x = x,
    i = i,
    i_type = i_type,
    m = m
  )

  lt <- inputs$lt
  x <- inputs$x
  i <- inputs$i
  i_type <- inputs$i_type
  m <- inputs$m

  .validate_mc_common(
    x = x,
    i = i,
    m = m,
    n_sim = n_sim,
    seed = seed,
    n = n
  )

  .mc_assert_positive_integer(k, "k")
  .mc_assert_numeric_scalar(payment, "payment", min = 0)

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  k <- as.integer(round(k))
  n_sim <- as.integer(round(n_sim))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  i_effective <- .mc_effective_rate(
    i = i,
    i_type = i_type,
    m = m
  )

  lifetime <- .simulate_lifetime_inverse_lx(
    lt = lt,
    x = x,
    n = n,
    n_sim = n_sim,
    seed = seed,
    frac = frac
  )

  n_payments <- .annuity_payment_count(
    Kx = lifetime$Kx,
    Tx = lifetime$Tx,
    died_within_horizon = lifetime$died_within_horizon,
    n = n,
    k = k,
    timing = timing
  )

  annuity_factor <- .annuity_factor_count(
    n_payments = n_payments,
    i_effective = i_effective,
    k = k,
    timing = timing
  )

  present_value <- payment * annuity_factor

  first_payment_time <- ifelse(
    n_payments > 0L,
    if (timing == "due") 0 else 1 / k,
    NA_real_
  )

  last_payment_time <- ifelse(
    n_payments > 0L,
    if (timing == "due") {
      (n_payments - 1L) / k
    } else {
      n_payments / k
    },
    NA_real_
  )

  out <- tibble::tibble(
    sim_id = seq_len(n_sim),
    simulation_id = seq_len(n_sim),
    contract_type = "annuity",
    x = x,
    age = x,
    n = n,
    term_years = n,
    k = k,
    payments_per_year = k,
    timing = timing,
    payment = payment,
    payment_per_payment = payment,
    payment_annualized = k * payment,
    frac = frac,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective,
    Kx = lifetime$Kx,
    curtate_lifetime = lifetime$Kx,
    Tx = lifetime$Tx,
    complete_lifetime = lifetime$Tx,
    death_age = lifetime$death_age,
    died_within_horizon = lifetime$died_within_horizon,
    n_payments = n_payments,
    first_payment_time = first_payment_time,
    last_payment_time = last_payment_time,
    present_value = present_value,
    pv_annuity = present_value
  )

  if (output == "summary") {
    return(
      summary_mc(
        out,
        value_col = "present_value"
      )
    )
  }

  out
}


#' Monte Carlo simulation of a life insurance
#'
#' Simulates present values of a discrete single-life insurance payable at the
#' end of the year of death, or at maturity for an endowment insurance.
#'
#' @param lt A life table or a \code{tidyact_life_contract}.
#' @param x Integer actuarial age.
#' @param i Numeric scalar interest-rate input.
#' @param i_type Interest-rate convention.
#' @param m Positive integer nominal conversion frequency.
#' @param type \code{"whole"}, \code{"term"}, or \code{"endowment"}.
#' @param n Positive integer term, or \code{Inf} for whole life.
#' @param benefit Nonnegative insurance benefit.
#' @param method Simulation method. Currently only \code{"inverse"}.
#' @param n_sim Positive integer number of simulations.
#' @param seed Optional nonnegative integer seed.
#' @param output \code{"simulations"} or \code{"summary"}.
#' @param ... Transitional compatibility for \code{mortality_table},
#'   \code{age}, \code{rate}, \code{rate_type}, \code{insurance_type}, and
#'   \code{term_years}.
#'
#' @return A tibble with one row per simulation, or a summary. Simulation
#'   output includes \code{benefit_indicator}, \code{benefit_time},
#'   \code{present_value}, and \code{pv_benefit}.
#'
#' @family simulation
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 0)
#' )
#'
#' simulate_insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   type = "term",
#'   n = 5,
#'   benefit = 100000,
#'   n_sim = 25,
#'   seed = 123
#' )
#'
#' @export
simulate_insurance_x <- function(
    lt,
    x = NULL,
    i = NULL,
    i_type = NULL,
    m = NULL,
    type = c("whole", "term", "endowment"),
    n = Inf,
    benefit = 1,
    method = c("inverse"),
    n_sim = 10000L,
    seed = NULL,
    output = c("simulations", "summary"),
    ...
) {
  lt_missing <- missing(lt)
  x_missing <- missing(x)
  i_missing <- missing(i)
  i_type_missing <- missing(i_type)
  type_missing <- missing(type)
  n_missing <- missing(n)

  method <- match.arg(method)
  output <- match.arg(output)

  old <- .life_mc_collect_old_args(
    dots = list(...),
    allowed = c(
      "mortality_table",
      "age",
      "rate",
      "rate_type",
      "insurance_type",
      "term_years"
    )
  )

  if ("mortality_table" %in% names(old)) {
    if (!lt_missing) {
      stop(
        "Provide only one of `lt` or deprecated `mortality_table`.",
        call. = FALSE
      )
    }

    lt <- old[["mortality_table"]]
    lt_missing <- FALSE
  }

  if ("age" %in% names(old)) {
    if (!x_missing) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }

    x <- old[["age"]]
  }

  if ("rate" %in% names(old)) {
    if (!i_missing) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- old[["rate"]]
  }

  if ("rate_type" %in% names(old)) {
    if (!i_type_missing) {
      stop(
        "Provide only one of `i_type` or deprecated `rate_type`.",
        call. = FALSE
      )
    }

    i_type <- old[["rate_type"]]
  }

  if ("insurance_type" %in% names(old)) {
    if (!type_missing) {
      stop(
        "Provide only one of `type` or deprecated `insurance_type`.",
        call. = FALSE
      )
    }

    type <- old[["insurance_type"]]
  }

  if ("term_years" %in% names(old)) {
    if (!n_missing) {
      stop(
        "Provide only one of `n` or deprecated `term_years`.",
        call. = FALSE
      )
    }

    n <- old[["term_years"]]
  }

  if (lt_missing) {
    stop("`lt` must be provided.", call. = FALSE)
  }

  type <- match.arg(
    type,
    choices = c("whole", "term", "endowment")
  )

  inputs <- .resolve_life_mc_inputs(
    lt = lt,
    x = x,
    i = i,
    i_type = i_type,
    m = m
  )

  lt <- inputs$lt
  x <- inputs$x
  i <- inputs$i
  i_type <- inputs$i_type
  m <- inputs$m

  .validate_mc_common(
    x = x,
    i = i,
    m = m,
    n_sim = n_sim,
    seed = seed,
    n = n
  )

  .mc_assert_numeric_scalar(
    benefit,
    "benefit",
    min = 0
  )

  if (type %in% c("term", "endowment") && is.infinite(n)) {
    stop(
      "`n` must be finite for term and endowment insurance.",
      call. = FALSE
    )
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  n_sim <- as.integer(round(n_sim))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  i_effective <- .mc_effective_rate(
    i = i,
    i_type = i_type,
    m = m
  )

  lifetime <- .simulate_lifetime_inverse_lx(
    lt = lt,
    x = x,
    n = n,
    n_sim = n_sim,
    seed = seed,
    frac = NULL
  )

  Kx <- lifetime$Kx
  v <- 1 / (1 + i_effective)

  if (type == "whole") {
    benefit_indicator <- as.numeric(
      lifetime$died_within_horizon
    )
    benefit_time <- ifelse(
      benefit_indicator == 1,
      Kx + 1L,
      NA_real_
    )
  } else if (type == "term") {
    benefit_indicator <- as.numeric(
      lifetime$died_within_horizon & Kx < n
    )
    benefit_time <- ifelse(
      benefit_indicator == 1,
      Kx + 1L,
      NA_real_
    )
  } else {
    death_before_term <-
      lifetime$died_within_horizon & Kx < n

    benefit_indicator <- rep(1, n_sim)
    benefit_time <- ifelse(
      death_before_term,
      Kx + 1L,
      n
    )
  }

  present_value <- ifelse(
    benefit_indicator == 1,
    benefit * v^benefit_time,
    0
  )

  out <- tibble::tibble(
    sim_id = seq_len(n_sim),
    simulation_id = seq_len(n_sim),
    contract_type = paste0(type, "_insurance"),
    x = x,
    age = x,
    n = n,
    term_years = n,
    type = type,
    insurance_type = type,
    benefit = benefit,
    i = i,
    rate = i,
    i_type = i_type,
    rate_type = i_type,
    m = m,
    i_effective = i_effective,
    Kx = lifetime$Kx,
    curtate_lifetime = lifetime$Kx,
    death_age = lifetime$death_age,
    died_within_horizon = lifetime$died_within_horizon,
    covered = as.logical(benefit_indicator),
    benefit_indicator = benefit_indicator,
    payment_time = benefit_time,
    benefit_time = benefit_time,
    present_value = present_value,
    pv_benefit = present_value
  )

  if (output == "summary") {
    return(
      summary_mc(
        out,
        value_col = "present_value"
      )
    )
  }

  out
}


#' Internal helper: collect deprecated Monte Carlo argument names
#'
#' @keywords internal
.life_mc_collect_old_args <- function(dots, allowed) {
  if (!is.list(dots)) {
    stop("`dots` must be a list.", call. = FALSE)
  }

  if (length(dots) == 0L) {
    return(dots)
  }

  dot_names <- names(dots)

  if (is.null(dot_names) ||
      anyNA(dot_names) ||
      any(!nzchar(dot_names))) {
    stop(
      "All arguments supplied through `...` must be named.",
      call. = FALSE
    )
  }

  bad <- setdiff(dot_names, allowed)

  if (length(bad) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  dots
}


#' Internal helper: resolve inputs for single-life Monte Carlo functions
#'
#' @keywords internal
.resolve_life_mc_inputs <- function(
    lt,
    x = NULL,
    i = NULL,
    i_type = NULL,
    m = NULL
) {
  `%||%` <- function(a, b) {
    if (!is.null(a)) a else b
  }

  if (inherits(lt, "tidyact_life_contract")) {
    contract <- lt

    if (!identical(contract$lives, "single")) {
      stop(
        "Only `lives = 'single'` is supported by these simulation ",
        "functions.",
        call. = FALSE
      )
    }

    lt <- contract$mortality_table %||% contract$lt

    if (is.null(x)) {
      x <- contract$x %||% contract$age
    }

    if (is.null(i)) {
      i <- contract$i %||% contract$rate
    }

    if (is.null(i_type)) {
      i_type <- contract$i_type %||% contract$rate_type
    }

    if (is.null(m)) {
      m <- contract$m
    }
  }

  if (is.null(i_type)) {
    i_type <- "effective"
  }

  if (is.null(m)) {
    m <- 1L
  }

  i_type <- .mc_normalize_i_type(i_type)

  list(
    lt = lt,
    x = x,
    i = i,
    i_type = i_type,
    m = m
  )
}


#' Internal helper: validate common Monte Carlo arguments
#'
#' @keywords internal
.validate_mc_common <- function(
    x,
    i,
    m,
    n_sim,
    seed,
    n
) {
  if (is.null(x) ||
      !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0 ||
      abs(x - round(x)) > 1e-10) {
    stop(
      "`x` must be a single nonnegative integer age.",
      call. = FALSE
    )
  }

  .mc_assert_numeric_scalar(i, "i")
  .mc_assert_positive_integer(m, "m")
  .mc_assert_positive_integer(n_sim, "n_sim")

  if (!is.null(seed)) {
    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed < 0 ||
        abs(seed - round(seed)) > 1e-10 ||
        seed > .Machine$integer.max) {
      stop(
        "`seed` must be `NULL` or a nonnegative integer within the ",
        "supported integer range.",
        call. = FALSE
      )
    }
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n <= 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) ||
        abs(n - round(n)) > 1e-10))) {
    stop(
      "`n` must be `Inf` or a single positive integer.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Internal helper: simulate future lifetimes from lx
#'
#' @param lt Life table with columns \code{x} and \code{lx}.
#' @param x Integer age.
#' @param n Positive integer horizon or \code{Inf}.
#' @param n_sim Positive integer number of simulations.
#' @param seed Optional seed.
#' @param frac Optional fractional-age assumption. If \code{NULL}, only the
#'   curtate lifetime is generated.
#'
#' @return A tibble with simulated curtate and complete future lifetimes.
#'
#' @keywords internal
.simulate_lifetime_inverse_lx <- function(
    lt,
    x,
    n,
    n_sim,
    seed = NULL,
    frac = NULL
) {
  if (!is.data.frame(lt)) {
    stop("`lt` must be a data frame or tibble.", call. = FALSE)
  }

  if (!all(c("x", "lx") %in% names(lt))) {
    stop(
      "`lt` must contain columns `x` and `lx`.",
      call. = FALSE
    )
  }

  if (!is.null(frac)) {
    frac <- match.arg(
      frac,
      c("udd", "cml", "balducci")
    )
  }

  lt <- lt[order(lt$x), , drop = FALSE]

  if (!is.numeric(lt$x) || !is.numeric(lt$lx)) {
    stop(
      "Columns `x` and `lx` must be numeric.",
      call. = FALSE
    )
  }

  if (anyNA(lt$x) ||
      any(!is.finite(lt$x)) ||
      any(abs(lt$x - round(lt$x)) > 1e-10)) {
    stop(
      "Column `x` must contain finite integer ages.",
      call. = FALSE
    )
  }

  if (anyNA(lt$lx) ||
      any(!is.finite(lt$lx)) ||
      any(lt$lx < 0)) {
    stop(
      "Column `lx` must contain finite nonnegative values.",
      call. = FALSE
    )
  }

  if (anyDuplicated(lt$x)) {
    stop(
      "Column `x` must not contain duplicated ages.",
      call. = FALSE
    )
  }

  if (any(diff(lt$lx) > 1e-10)) {
    stop(
      "Column `lx` must be nonincreasing with age.",
      call. = FALSE
    )
  }

  ages <- as.integer(round(lt$x))
  lx <- as.numeric(lt$lx)

  idx_x <- match(x, ages)

  if (is.na(idx_x)) {
    stop("`x` must be present in `lt$x`.", call. = FALSE)
  }

  if (lx[[idx_x]] <= 0) {
    stop(
      "`lx` at age `x` must be strictly positive.",
      call. = FALSE
    )
  }

  omega <- max(ages)

  if (is.infinite(n)) {
    horizon <- omega - x

    if (horizon <= 0) {
      stop(
        "The life table does not provide a future lifetime horizon.",
        call. = FALSE
      )
    }
  } else {
    horizon <- as.integer(round(n))
  }

  needed_ages <- x:(x + horizon)

  if (!all(needed_ages %in% ages)) {
    stop(
      "`lt` does not contain every integer age required for the ",
      "simulation horizon.",
      call. = FALSE
    )
  }

  lx_path <- lx[match(needed_ages, ages)]

  if (is.infinite(n) && tail(lx_path, 1L) > 0) {
    stop(
      "Whole-life simulation requires the life table to end with ",
      "`lx = 0` at or after age `x`.",
      call. = FALSE
    )
  }

  l0 <- lx_path[[1L]]
  death_counts <- lx_path[-length(lx_path)] -
    lx_path[-1L]

  if (any(death_counts < -1e-10)) {
    stop(
      "The life table implies negative death counts.",
      call. = FALSE
    )
  }

  death_counts <- pmax(death_counts, 0)
  survival_to_horizon <- tail(lx_path, 1L)

  if (is.infinite(n)) {
    probs <- death_counts / l0
  } else {
    probs <- c(
      death_counts / l0,
      survival_to_horizon / l0
    )
  }

  if (anyNA(probs) ||
      any(!is.finite(probs)) ||
      any(probs < 0) ||
      sum(probs) <= 0) {
    stop(
      "The simulated lifetime distribution has invalid probabilities.",
      call. = FALSE
    )
  }

  probability_error <- abs(sum(probs) - 1)

  if (probability_error > 1e-8) {
    stop(
      "The life-table probabilities do not sum to 1 over the requested ",
      "horizon.",
      call. = FALSE
    )
  }

  probs <- probs / sum(probs)

  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (!is.null(seed)) {
    if (had_seed) {
      old_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }

    on.exit(
      {
        if (had_seed) {
          assign(
            ".Random.seed",
            old_seed,
            envir = .GlobalEnv
          )
        } else if (exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )) {
          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )

    set.seed(as.integer(round(seed)))
  }

  sampled <- sample.int(
    n = length(probs),
    size = n_sim,
    replace = TRUE,
    prob = probs
  )

  if (is.infinite(n)) {
    Kx <- sampled - 1L
    died_within_horizon <- rep(TRUE, n_sim)
  } else {
    survival_category <- length(probs)
    died_within_horizon <-
      sampled != survival_category

    Kx <- ifelse(
      died_within_horizon,
      sampled - 1L,
      horizon
    )
  }

  Kx <- as.integer(Kx)
  Tx <- rep(NA_real_, n_sim)

  if (!is.null(frac) && any(died_within_horizon)) {
    death_rows <- which(died_within_horizon)
    death_ages <- x + Kx[death_rows]

    l_start <- lx[match(death_ages, ages)]
    l_end <- lx[match(death_ages + 1L, ages)]

    p_year <- l_end / l_start
    q_year <- 1 - p_year

    u <- stats::runif(length(death_rows))

    s <- if (identical(frac, "udd")) {
      u
    } else if (identical(frac, "cml")) {
      out <- numeric(length(u))
      degenerate <- p_year <= 0

      out[degenerate] <- 0
      regular <- !degenerate

      out[regular] <-
        log1p(-u[regular] * q_year[regular]) /
        log(p_year[regular])

      out
    } else {
      out <- numeric(length(u))
      degenerate <- p_year <= 0

      out[degenerate] <- 0
      regular <- !degenerate

      out[regular] <-
        u[regular] * p_year[regular] /
        (1 - u[regular] * q_year[regular])

      out
    }

    s <- pmin(pmax(s, 0), 1 - .Machine$double.eps)
    Tx[death_rows] <- Kx[death_rows] + s
  }

  tibble::tibble(
    Kx = Kx,
    curtate_lifetime = Kx,
    Tx = Tx,
    complete_lifetime = Tx,
    death_age = ifelse(
      died_within_horizon,
      x + Kx + 1L,
      NA_real_
    ),
    died_within_horizon = died_within_horizon
  )
}


#' Internal helper: count annuity payments
#'
#' @keywords internal
.annuity_payment_count <- function(
    Kx,
    died_within_horizon,
    n,
    k,
    timing,
    Tx = NULL
) {
  timing <- match.arg(
    timing,
    c("due", "immediate")
  )

  .mc_assert_positive_integer(k, "k")
  k <- as.integer(round(k))

  if (length(Kx) != length(died_within_horizon)) {
    stop(
      "`Kx` and `died_within_horizon` must have the same length.",
      call. = FALSE
    )
  }

  if (k == 1L) {
    years_lived <- if (is.infinite(n)) {
      Kx
    } else {
      pmin(Kx, n)
    }

    out <- if (timing == "due") {
      years_lived + 1L
    } else {
      years_lived
    }

    if (!is.infinite(n)) {
      out <- pmin(out, as.integer(n))
    }

    return(as.integer(pmax(out, 0L)))
  }

  if (is.null(Tx) ||
      !is.numeric(Tx) ||
      length(Tx) != length(Kx)) {
    stop(
      "`Tx` must be supplied with the same length as `Kx` when `k > 1`.",
      call. = FALSE
    )
  }

  vapply(
    seq_along(Kx),
    function(idx) {
      if (isTRUE(died_within_horizon[[idx]])) {
        T_schedule <- Tx[[idx]]

        if (is.na(T_schedule) || !is.finite(T_schedule)) {
          stop(
            "A complete lifetime is required for every simulated death ",
            "when `k > 1`.",
            call. = FALSE
          )
        }
      } else {
        if (is.infinite(n)) {
          stop(
            "A whole-life simulation cannot contain a censored lifetime.",
            call. = FALSE
          )
        }

        T_schedule <- n + 1 / k
      }

      times <- .mc_annuity_payment_times_fractional(
        T = T_schedule,
        type = if (is.infinite(n)) "whole" else "temporary",
        n = if (is.infinite(n)) NULL else n,
        h = 0,
        n_guar = NULL,
        timing = timing,
        payment_interval = 1 / k
      )

      length(times)
    },
    integer(1)
  )
}


#' Internal helper: convert payment counts into annuity factors
#'
#' @keywords internal
.annuity_factor_count <- function(
    n_payments,
    i_effective,
    k,
    timing
) {
  timing <- match.arg(
    timing,
    c("due", "immediate")
  )

  .mc_assert_positive_integer(k, "k")
  k <- as.integer(round(k))

  if (!is.numeric(n_payments) ||
      anyNA(n_payments) ||
      any(!is.finite(n_payments)) ||
      any(n_payments < 0) ||
      any(abs(n_payments - round(n_payments)) > 1e-10)) {
    stop(
      "`n_payments` must contain nonnegative integers.",
      call. = FALSE
    )
  }

  i_period <- (1 + i_effective)^(1 / k) - 1
  v_period <- 1 / (1 + i_period)

  out <- numeric(length(n_payments))
  positive <- n_payments > 0L

  if (!any(positive)) {
    return(out)
  }

  if (abs(i_period) < 1e-12) {
    out[positive] <- n_payments[positive]
    return(out)
  }

  due_factor <-
    (1 - v_period^n_payments[positive]) /
    (1 - v_period)

  if (timing == "due") {
    out[positive] <- due_factor
  } else {
    out[positive] <- v_period * due_factor
  }

  out
}
