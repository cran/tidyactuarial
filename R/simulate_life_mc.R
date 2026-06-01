#' Monte Carlo simulation of a life annuity
#'
#' Simulates the present value of a discrete single-life annuity using a life
#' table and annual curtate future lifetimes, with compact actuarial notation.
#'
#' The function supports direct use with \code{lt}, \code{x}, and \code{i}, and
#' pipe workflows with \code{\link{life_contract}}.
#'
#' @param lt A life table or a \code{tidyact_life_contract} object. A life table
#'   must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age. Optional when \code{lt} is a
#'   \code{tidyact_life_contract}.
#' @param i Numeric scalar. Annual interest-rate input. Optional when
#'   \code{lt} is a \code{tidyact_life_contract}.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param n Term in years. Use \code{Inf} for whole life.
#' @param k Positive integer. Number of annuity payments per year.
#' @param timing Character string. Either \code{"due"} or \code{"immediate"}.
#' @param payment Numeric scalar. Level payment per payment period.
#' @param method Simulation method. Currently \code{"inverse"} is supported.
#' @param n_sim Positive integer. Number of simulations.
#' @param seed Optional seed for reproducibility.
#' @param output Character string. Use \code{"simulations"} for the simulated
#'   data or \code{"summary"} for summary statistics using
#'   \code{\link{summary_mc}}.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{term_years}, and \code{payments_per_year}.
#'
#' @return A tibble. If \code{output = "simulations"}, the tibble contains one
#'   row per simulation. If \code{output = "summary"}, the tibble is the output
#'   of \code{\link{summary_mc}} applied to the simulated present values.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age,
#' \code{i} is the interest-rate input, \code{i_type} is the interest-rate type,
#' \code{m} is the conversion frequency for nominal rates, \code{n} is the term,
#' and \code{k} is the payment frequency.
#'
#' The simulation is based on the curtate future lifetime \eqn{K_x}. For annual
#' payments, an annuity-due pays while the life is alive at the beginning of the
#' year, and an annuity-immediate pays at the end of completed years.
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
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   n = 5,
#'   k = 12,
#'   n_sim = 25,
#'   seed = 123,
#'   output = "summary"
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
    method = c("inverse"),
    n_sim = 10000L,
    seed = NULL,
    output = c("simulations", "summary"),
    ...
) {
  timing <- match.arg(timing)
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

  if (!is.null(old$mortality_table)) {
    if (!missing(lt)) {
      stop("Provide only one of `lt` or deprecated `mortality_table`.", call. = FALSE)
    }
    lt <- old$mortality_table
  }

  if (!is.null(old$age)) {
    if (!is.null(x)) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }
    x <- old$age
  }

  if (!is.null(old$rate)) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }
    i <- old$rate
  }

  if (!is.null(old$rate_type)) {
    if (!is.null(i_type)) {
      stop("Provide only one of `i_type` or deprecated `rate_type`.", call. = FALSE)
    }
    i_type <- old$rate_type
  }

  if (!is.null(old$term_years)) {
    if (!is.infinite(n)) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }
    n <- old$term_years
  }

  if (!is.null(old$payments_per_year)) {
    if (!identical(k, 1L) && !identical(k, 1)) {
      stop("Provide only one of `k` or deprecated `payments_per_year`.", call. = FALSE)
    }
    k <- old$payments_per_year
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

  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k < 1 ||
      abs(k - round(k)) > 1e-10) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(payment) ||
      length(payment) != 1L ||
      is.na(payment) ||
      !is.finite(payment)) {
    stop("`payment` must be a single finite numeric value.", call. = FALSE)
  }

  x <- as.integer(round(x))
  m <- as.integer(round(m))
  k <- as.integer(round(k))
  n_sim <- as.integer(round(n_sim))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (!is.finite(i_effective) || i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  lifetime <- .simulate_lifetime_inverse_lx(
    lt = lt,
    x = x,
    n = n,
    n_sim = n_sim,
    seed = seed
  )

  Kx <- lifetime$Kx

  n_payments <- .annuity_payment_count(
    Kx = Kx,
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
    n_payments = n_payments,
    present_value = present_value
  )

  if (output == "summary") {
    return(summary_mc(out, value_col = "present_value"))
  }

  out
}


#' Monte Carlo simulation of a life insurance
#'
#' Simulates the present value of a discrete single-life insurance using a life
#' table and annual curtate future lifetimes, with compact actuarial notation.
#'
#' The function supports direct use with \code{lt}, \code{x}, and \code{i}, and
#' pipe workflows with \code{\link{life_contract}}.
#'
#' @param lt A life table or a \code{tidyact_life_contract} object. A life table
#'   must contain columns \code{x} and \code{lx}.
#' @param x Integer actuarial age. Optional when \code{lt} is a
#'   \code{tidyact_life_contract}.
#' @param i Numeric scalar. Annual interest-rate input. Optional when
#'   \code{lt} is a \code{tidyact_life_contract}.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param type Character string. One of \code{"whole"}, \code{"term"}, or
#'   \code{"endowment"}.
#' @param n Term in years. Required as finite for term and endowment insurance.
#'   Use \code{Inf} for whole life.
#' @param benefit Numeric scalar. Insurance benefit.
#' @param method Simulation method. Currently \code{"inverse"} is supported.
#' @param n_sim Positive integer. Number of simulations.
#' @param seed Optional seed for reproducibility.
#' @param output Character string. Use \code{"simulations"} for the simulated
#'   data or \code{"summary"} for summary statistics using
#'   \code{\link{summary_mc}}.
#' @param ... Transitional compatibility for older calls using
#'   \code{mortality_table}, \code{age}, \code{rate}, \code{rate_type},
#'   \code{insurance_type}, and \code{term_years}.
#'
#' @return A tibble. If \code{output = "simulations"}, the tibble contains one
#'   row per simulation. If \code{output = "summary"}, the tibble is the output
#'   of \code{\link{summary_mc}} applied to the simulated present values.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age,
#' \code{i} is the interest-rate input, \code{i_type} is the interest-rate type,
#' \code{m} is the conversion frequency for nominal rates, and \code{n} is the
#' insurance term.
#'
#' The benefit is paid at the end of the year of death for whole-life and term
#' insurance. For endowment insurance, the benefit is paid at death within the
#' term or at maturity if the life survives.
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
#'   type = "whole",
#'   n_sim = 25,
#'   seed = 123
#' )
#'
#' simulate_insurance_x(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   type = "term",
#'   n = 5,
#'   benefit = 100000,
#'   n_sim = 25,
#'   seed = 123,
#'   output = "summary"
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
  type <- match.arg(type)
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

  if (!is.null(old$mortality_table)) {
    if (!missing(lt)) {
      stop("Provide only one of `lt` or deprecated `mortality_table`.", call. = FALSE)
    }
    lt <- old$mortality_table
  }

  if (!is.null(old$age)) {
    if (!is.null(x)) {
      stop("Provide only one of `x` or deprecated `age`.", call. = FALSE)
    }
    x <- old$age
  }

  if (!is.null(old$rate)) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }
    i <- old$rate
  }

  if (!is.null(old$rate_type)) {
    if (!is.null(i_type)) {
      stop("Provide only one of `i_type` or deprecated `rate_type`.", call. = FALSE)
    }
    i_type <- old$rate_type
  }

  if (!is.null(old$insurance_type)) {
    type <- match.arg(old$insurance_type, c("whole", "term", "endowment"))
  }

  if (!is.null(old$term_years)) {
    if (!is.infinite(n)) {
      stop("Provide only one of `n` or deprecated `term_years`.", call. = FALSE)
    }
    n <- old$term_years
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
  m <- as.integer(round(m))
  n_sim <- as.integer(round(n_sim))

  if (!is.infinite(n)) {
    n <- as.integer(round(n))
  }

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (!is.finite(i_effective) || i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  lifetime <- .simulate_lifetime_inverse_lx(
    lt = lt,
    x = x,
    n = n,
    n_sim = n_sim,
    seed = seed
  )

  Kx <- lifetime$Kx
  v <- 1 / (1 + i_effective)

  if (type == "whole") {
    covered <- lifetime$died_within_horizon
    payment_time <- Kx + 1L
    present_value <- ifelse(covered, benefit * v^payment_time, 0)
  } else if (type == "term") {
    covered <- lifetime$died_within_horizon & Kx < n
    payment_time <- Kx + 1L
    present_value <- ifelse(covered, benefit * v^payment_time, 0)
  } else {
    death_covered <- lifetime$died_within_horizon & Kx < n
    death_payment_time <- Kx + 1L
    maturity_payment_time <- n

    present_value <- ifelse(
      death_covered,
      benefit * v^death_payment_time,
      benefit * v^maturity_payment_time
    )

    covered <- rep(TRUE, n_sim)
    payment_time <- ifelse(death_covered, death_payment_time, maturity_payment_time)
  }

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
    covered = covered,
    payment_time = payment_time,
    present_value = present_value
  )

  if (output == "summary") {
    return(summary_mc(out, value_col = "present_value"))
  }

  out
}


#' Internal helper: collect old Monte Carlo argument names
#'
#' @param dots List of arguments from `...`.
#' @param allowed Character vector with allowed deprecated names.
#'
#' @return A list.
#'
#' @keywords internal
.life_mc_collect_old_args <- function(dots, allowed) {
  bad <- setdiff(names(dots), allowed)

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
#' @param lt A life table or tidyactuarial life contract.
#' @param x Integer age.
#' @param i Interest-rate input.
#' @param i_type Interest-rate type.
#' @param m Nominal conversion frequency.
#'
#' @return A list with resolved arguments.
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
        "Only `lives = 'single'` is currently supported by individual simulation functions.",
        call. = FALSE
      )
    }

    lt <- contract$mortality_table

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
#' @param x Integer age.
#' @param i Interest-rate input.
#' @param m Nominal conversion frequency.
#' @param n_sim Number of simulations.
#' @param seed Optional seed.
#' @param n Term in years.
#'
#' @return Invisibly returns `TRUE`.
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
      abs(x - round(x)) > 1e-10) {
    stop("`x` must be a single integer age.", call. = FALSE)
  }

  if (is.null(i) ||
      !is.numeric(i) ||
      length(i) != 1L ||
      is.na(i) ||
      !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      is.na(m) ||
      !is.finite(m) ||
      m < 1 ||
      abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(n_sim) ||
      length(n_sim) != 1L ||
      is.na(n_sim) ||
      !is.finite(n_sim) ||
      n_sim < 1 ||
      abs(n_sim - round(n_sim)) > 1e-10) {
    stop("`n_sim` must be a single positive integer.", call. = FALSE)
  }

  if (!is.null(seed) &&
      (!is.numeric(seed) ||
       length(seed) != 1L ||
       is.na(seed) ||
       !is.finite(seed))) {
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n <= 0 ||
      (!is.infinite(n) &&
       (!is.finite(n) ||
        abs(n - round(n)) > 1e-10))) {
    stop("`n` must be `Inf` or a single positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Internal helper: simulate curtate future lifetimes from lx
#'
#' @param lt Life table with columns `x` and `lx`.
#' @param x Integer age.
#' @param n Term in years. Use `Inf` for whole life.
#' @param n_sim Number of simulations.
#' @param seed Optional seed.
#'
#' @return A tibble with simulated curtate lifetimes.
#'
#' @keywords internal
.simulate_lifetime_inverse_lx <- function(
    lt,
    x,
    n,
    n_sim,
    seed = NULL
) {
  if (!is.data.frame(lt)) {
    stop("`lt` must be a data.frame or tibble.", call. = FALSE)
  }

  if (!all(c("x", "lx") %in% names(lt))) {
    stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
  }

  lt <- lt[order(lt$x), , drop = FALSE]

  if (!is.numeric(lt$x) || !is.numeric(lt$lx)) {
    stop("Columns `x` and `lx` must be numeric.", call. = FALSE)
  }

  if (any(is.na(lt$x)) || any(!is.finite(lt$x)) ||
      any(abs(lt$x - round(lt$x)) > 1e-10)) {
    stop("Column `x` must contain finite integer ages.", call. = FALSE)
  }

  if (any(is.na(lt$lx)) || any(!is.finite(lt$lx)) || any(lt$lx < 0)) {
    stop("Column `lx` must contain finite nonnegative values.", call. = FALSE)
  }

  if (anyDuplicated(lt$x)) {
    stop("Column `x` must not contain duplicated ages.", call. = FALSE)
  }

  ages <- as.integer(round(lt$x))
  lx <- as.numeric(lt$lx)

  idx_x <- match(x, ages)

  if (is.na(idx_x)) {
    stop("`x` must be present in `lt$x`.", call. = FALSE)
  }

  if (lx[[idx_x]] <= 0) {
    stop("`lx` at age `x` must be strictly positive.", call. = FALSE)
  }

  omega <- max(ages)

  if (is.infinite(n)) {
    horizon <- omega - x

    if (horizon <= 0) {
      stop("The life table does not provide a future lifetime horizon.", call. = FALSE)
    }
  } else {
    horizon <- as.integer(round(n))
  }

  needed_ages <- x:(x + horizon)

  if (!all(needed_ages %in% ages)) {
    stop(
      "`lt` does not contain enough integer ages for the requested simulation horizon.",
      call. = FALSE
    )
  }

  lx_path <- lx[match(needed_ages, ages)]

  if (is.infinite(n) && tail(lx_path, 1) > 0) {
    stop(
      "Whole-life simulation requires the life table to end with `lx = 0` at or after `x`.",
      call. = FALSE
    )
  }

  l0 <- lx_path[[1]]

  death_counts <- pmax(lx_path[-length(lx_path)] - lx_path[-1], 0)
  survival_to_horizon <- tail(lx_path, 1)

  if (is.infinite(n)) {
    probs <- death_counts / l0
    categories <- seq_along(death_counts)
  } else {
    probs <- c(death_counts / l0, survival_to_horizon / l0)
    categories <- seq_along(probs)
  }

  if (sum(probs) <= 0 || any(is.na(probs))) {
    stop("The simulated lifetime distribution has invalid probabilities.", call. = FALSE)
  }

  if (abs(sum(probs) - 1) > 1e-8) {
    probs <- probs / sum(probs)
  }

  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      stats::runif(1)
    }

    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    on.exit(
      assign(".Random.seed", old_seed, envir = .GlobalEnv),
      add = TRUE
    )

    set.seed(seed)
  }

  sampled <- sample(
    x = categories,
    size = n_sim,
    replace = TRUE,
    prob = probs
  )

  if (is.infinite(n)) {
    Kx <- sampled - 1L
    died_within_horizon <- rep(TRUE, n_sim)
  } else {
    survival_category <- length(probs)
    died_within_horizon <- sampled != survival_category
    Kx <- ifelse(died_within_horizon, sampled - 1L, horizon)
  }

  tibble::tibble(
    Kx = as.integer(Kx),
    curtate_lifetime = as.integer(Kx),
    death_age = ifelse(
      died_within_horizon,
      x + as.integer(Kx) + 1L,
      NA_real_
    ),
    died_within_horizon = died_within_horizon
  )
}


#' Internal helper: count annuity payments from simulated curtate lifetimes
#'
#' @param Kx Simulated curtate future lifetimes.
#' @param died_within_horizon Logical vector.
#' @param n Term in years.
#' @param k Payment frequency.
#' @param timing Payment timing.
#'
#' @return Integer vector with number of payments.
#'
#' @keywords internal
.annuity_payment_count <- function(
    Kx,
    died_within_horizon,
    n,
    k,
    timing
) {
  kk <- as.integer(k)

  if (is.infinite(n)) {
    years_lived <- Kx
    cap_payments <- Inf
  } else {
    years_lived <- pmin(Kx, n)
    cap_payments <- as.integer(round(n * kk))
  }

  if (timing == "due") {
    out <- pmax(years_lived * kk + 1L, 0L)
  } else {
    out <- pmax(years_lived * kk, 0L)
  }

  # For finite temporary annuities, the maximum number of payments is n * k.
  if (is.finite(cap_payments)) {
    out <- pmin(out, cap_payments)
  }

  as.integer(out)
}


#' Internal helper: convert number of payments into annuity factors
#'
#' @param n_payments Integer vector with number of payments.
#' @param i_effective Annual effective interest rate.
#' @param k Payment frequency.
#' @param timing Payment timing.
#'
#' @return Numeric vector of annuity factors.
#'
#' @keywords internal
.annuity_factor_count <- function(
    n_payments,
    i_effective,
    k,
    timing
) {
  i_period <- (1 + i_effective)^(1 / k) - 1
  v_period <- 1 / (1 + i_period)

  out <- numeric(length(n_payments))
  positive <- n_payments > 0L

  if (!any(positive)) {
    return(out)
  }

  if (abs(i_period) < 1e-12) {
    out[positive] <- n_payments[positive] / k
    return(out)
  }

  due_factor <- (1 / k) *
    (1 - v_period^n_payments[positive]) /
    (1 - v_period)

  if (timing == "due") {
    out[positive] <- due_factor
  } else {
    out[positive] <- v_period * due_factor
  }

  out
}
