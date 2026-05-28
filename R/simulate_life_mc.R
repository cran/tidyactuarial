#' Monte Carlo simulation of a life annuity
#'
#' Simulates the present value of a discrete life annuity using a mortality
#' table and annual curtate future lifetimes.
#'
#' The function supports direct use with `mortality_table`, `age`, and `rate`,
#' and pipe workflows with [life_contract()].
#'
#' @param mortality_table A life table or a `tidyact_life_contract` object.
#'   A life table must contain columns `x` and `lx`.
#' @param age Integer actuarial age. Optional when `mortality_table` is a
#'   `tidyact_life_contract`.
#' @param rate Numeric scalar. Annual interest-rate input. Optional when
#'   `mortality_table` is a `tidyact_life_contract`.
#' @param rate_type Character string indicating the rate type.
#' @param m Positive integer. Compounding frequency for nominal rates.
#' @param term_years Term in years. Use `Inf` for whole life.
#' @param payments_per_year Positive integer. Number of annuity payments per year.
#' @param timing `"due"` or `"immediate"`.
#' @param payment Numeric scalar. Level payment per payment period.
#' @param method Simulation method. Currently `"inverse"` is supported.
#' @param n_sim Positive integer. Number of simulations.
#' @param seed Optional seed for reproducibility.
#' @param output `"simulations"` for the simulated data or `"summary"` for
#'   summary statistics using [summary_mc()].
#'
#' @return A tibble.
#'
#' @family simulation
#' @family life-contingencies
#'
#' @export
simulate_annuity_x <- function(
    mortality_table,
    age = NULL,
    rate = NULL,
    rate_type = NULL,
    m = NULL,
    term_years = Inf,
    payments_per_year = 1L,
    timing = c("due", "immediate"),
    payment = 1,
    method = c("inverse"),
    n_sim = 10000L,
    seed = NULL,
    output = c("simulations", "summary")
) {
  timing <- match.arg(timing)
  method <- match.arg(method)
  output <- match.arg(output)

  inputs <- .resolve_life_mc_inputs(
    mortality_table = mortality_table,
    age = age,
    rate = rate,
    rate_type = rate_type,
    m = m
  )

  mortality_table <- inputs$mortality_table
  age <- inputs$age
  rate <- inputs$rate
  rate_type <- inputs$rate_type
  m <- inputs$m

  .validate_mc_common(
    age = age,
    rate = rate,
    m = m,
    n_sim = n_sim,
    seed = seed,
    term_years = term_years
  )

  if (!is.numeric(payments_per_year) ||
      length(payments_per_year) != 1L ||
      is.na(payments_per_year) ||
      !is.finite(payments_per_year) ||
      payments_per_year < 1 ||
      abs(payments_per_year - round(payments_per_year)) > 1e-10) {
    stop("`payments_per_year` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(payment) ||
      length(payment) != 1L ||
      is.na(payment) ||
      !is.finite(payment)) {
    stop("`payment` must be a single finite numeric value.", call. = FALSE)
  }

  payments_per_year <- as.integer(round(payments_per_year))

  i_effective <- standardize_interest(
    type = rate_type,
    rate = rate,
    m = m
  )

  lifetime <- .simulate_lifetime_inverse_lx(
    mortality_table = mortality_table,
    age = age,
    term_years = term_years,
    n_sim = n_sim,
    seed = seed
  )

  k <- lifetime$curtate_lifetime

  n_payments <- .annuity_payment_count(
    curtate_lifetime = k,
    died_within_horizon = lifetime$died_within_horizon,
    term_years = term_years,
    payments_per_year = payments_per_year,
    timing = timing
  )

  annuity_factor <- .annuity_factor_count(
    n_payments = n_payments,
    i_effective = i_effective,
    payments_per_year = payments_per_year,
    timing = timing
  )

  present_value <- payment * annuity_factor

  out <- tibble::tibble(
    simulation_id = seq_len(n_sim),
    contract_type = "annuity",
    age = age,
    term_years = term_years,
    payments_per_year = payments_per_year,
    timing = timing,
    payment = payment,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    curtate_lifetime = lifetime$curtate_lifetime,
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
#' Simulates the present value of a discrete life insurance using a mortality
#' table and annual curtate future lifetimes.
#'
#' The function supports direct use with `mortality_table`, `age`, and `rate`,
#' and pipe workflows with [life_contract()].
#'
#' @param mortality_table A life table or a `tidyact_life_contract` object.
#' @param age Integer actuarial age. Optional when `mortality_table` is a
#'   `tidyact_life_contract`.
#' @param rate Numeric scalar. Annual interest-rate input. Optional when
#'   `mortality_table` is a `tidyact_life_contract`.
#' @param rate_type Character string indicating the rate type.
#' @param m Positive integer. Compounding frequency for nominal rates.
#' @param insurance_type `"whole"`, `"term"`, or `"endowment"`.
#' @param term_years Term in years. Required for term and endowment insurance.
#'   Use `Inf` for whole life.
#' @param benefit Numeric scalar. Insurance benefit.
#' @param method Simulation method. Currently `"inverse"` is supported.
#' @param n_sim Positive integer. Number of simulations.
#' @param seed Optional seed for reproducibility.
#' @param output `"simulations"` for the simulated data or `"summary"` for
#'   summary statistics using [summary_mc()].
#'
#' @return A tibble.
#'
#' @family simulation
#' @family life-contingencies
#'
#' @export
simulate_insurance_x <- function(
    mortality_table,
    age = NULL,
    rate = NULL,
    rate_type = NULL,
    m = NULL,
    insurance_type = c("whole", "term", "endowment"),
    term_years = Inf,
    benefit = 1,
    method = c("inverse"),
    n_sim = 10000L,
    seed = NULL,
    output = c("simulations", "summary")
) {
  insurance_type <- match.arg(insurance_type)
  method <- match.arg(method)
  output <- match.arg(output)

  inputs <- .resolve_life_mc_inputs(
    mortality_table = mortality_table,
    age = age,
    rate = rate,
    rate_type = rate_type,
    m = m
  )

  mortality_table <- inputs$mortality_table
  age <- inputs$age
  rate <- inputs$rate
  rate_type <- inputs$rate_type
  m <- inputs$m

  .validate_mc_common(
    age = age,
    rate = rate,
    m = m,
    n_sim = n_sim,
    seed = seed,
    term_years = term_years
  )

  if (!is.numeric(benefit) ||
      length(benefit) != 1L ||
      is.na(benefit) ||
      !is.finite(benefit) ||
      benefit < 0) {
    stop("`benefit` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (insurance_type %in% c("term", "endowment") && is.infinite(term_years)) {
    stop(
      "`term_years` must be finite for term and endowment insurance.",
      call. = FALSE
    )
  }

  i_effective <- standardize_interest(
    type = rate_type,
    rate = rate,
    m = m
  )

  lifetime <- .simulate_lifetime_inverse_lx(
    mortality_table = mortality_table,
    age = age,
    term_years = term_years,
    n_sim = n_sim,
    seed = seed
  )

  k <- lifetime$curtate_lifetime
  v <- 1 / (1 + i_effective)

  if (insurance_type == "whole") {
    covered <- lifetime$died_within_horizon
    payment_time <- k + 1L
    present_value <- ifelse(covered, benefit * v^payment_time, 0)
  } else if (insurance_type == "term") {
    covered <- lifetime$died_within_horizon & k < term_years
    payment_time <- k + 1L
    present_value <- ifelse(covered, benefit * v^payment_time, 0)
  } else {
    death_covered <- lifetime$died_within_horizon & k < term_years
    death_payment_time <- k + 1L
    maturity_payment_time <- term_years

    present_value <- ifelse(
      death_covered,
      benefit * v^death_payment_time,
      benefit * v^maturity_payment_time
    )

    covered <- rep(TRUE, n_sim)
    payment_time <- ifelse(death_covered, death_payment_time, maturity_payment_time)
  }

  out <- tibble::tibble(
    simulation_id = seq_len(n_sim),
    contract_type = paste0(insurance_type, "_insurance"),
    age = age,
    term_years = term_years,
    benefit = benefit,
    rate = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    curtate_lifetime = lifetime$curtate_lifetime,
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


.resolve_life_mc_inputs <- function(
    mortality_table,
    age = NULL,
    rate = NULL,
    rate_type = NULL,
    m = NULL
) {
  if (inherits(mortality_table, "tidyact_life_contract")) {
    contract <- mortality_table

    if (!identical(contract$lives, "single")) {
      stop(
        "Only `lives = 'single'` is currently supported by individual simulation functions.",
        call. = FALSE
      )
    }

    mortality_table <- contract$mortality_table

    if (is.null(age)) {
      age <- contract$age
    }

    if (is.null(rate)) {
      rate <- contract$rate
    }

    if (is.null(rate_type)) {
      rate_type <- contract$rate_type
    }

    if (is.null(m)) {
      m <- contract$m
    }
  }

  if (is.null(rate_type)) {
    rate_type <- "effective"
  }

  if (is.null(m)) {
    m <- 1L
  }

  list(
    mortality_table = mortality_table,
    age = age,
    rate = rate,
    rate_type = rate_type,
    m = m
  )
}


.validate_mc_common <- function(
    age,
    rate,
    m,
    n_sim,
    seed,
    term_years
) {
  if (is.null(age) ||
      !is.numeric(age) ||
      length(age) != 1L ||
      is.na(age) ||
      !is.finite(age) ||
      abs(age - round(age)) > 1e-10) {
    stop("`age` must be a single integer age.", call. = FALSE)
  }

  if (is.null(rate) ||
      !is.numeric(rate) ||
      length(rate) != 1L ||
      is.na(rate) ||
      !is.finite(rate)) {
    stop("`rate` must be a single finite numeric value.", call. = FALSE)
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

  if (!is.numeric(term_years) ||
      length(term_years) != 1L ||
      is.na(term_years) ||
      term_years <= 0 ||
      (!is.infinite(term_years) &&
       (!is.finite(term_years) ||
        abs(term_years - round(term_years)) > 1e-10))) {
    stop("`term_years` must be `Inf` or a single positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}


.simulate_lifetime_inverse_lx <- function(
    mortality_table,
    age,
    term_years,
    n_sim,
    seed = NULL
) {
  if (!is.data.frame(mortality_table)) {
    stop("`mortality_table` must be a data.frame or tibble.", call. = FALSE)
  }

  if (!all(c("x", "lx") %in% names(mortality_table))) {
    stop("`mortality_table` must contain columns `x` and `lx`.", call. = FALSE)
  }

  lt <- mortality_table[order(mortality_table$x), , drop = FALSE]

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

  idx_age <- match(age, ages)

  if (is.na(idx_age)) {
    stop("`age` must be present in `mortality_table$x`.", call. = FALSE)
  }

  if (lx[[idx_age]] <= 0) {
    stop("`lx` at `age` must be strictly positive.", call. = FALSE)
  }

  max_age <- max(ages)

  if (is.infinite(term_years)) {
    horizon <- max_age - age

    if (horizon <= 0) {
      stop("The mortality table does not provide a future lifetime horizon.", call. = FALSE)
    }
  } else {
    horizon <- as.integer(round(term_years))
  }

  needed_ages <- age:(age + horizon)

  if (!all(needed_ages %in% ages)) {
    stop(
      "`mortality_table` does not contain enough integer ages for the requested simulation horizon.",
      call. = FALSE
    )
  }

  lx_path <- lx[match(needed_ages, ages)]

  if (is.infinite(term_years) && tail(lx_path, 1) > 0) {
    stop(
      "Whole-life simulation requires the mortality table to end with `lx = 0` at or after `age`.",
      call. = FALSE
    )
  }

  l0 <- lx_path[[1]]

  death_counts <- pmax(lx_path[-length(lx_path)] - lx_path[-1], 0)
  survival_to_horizon <- tail(lx_path, 1)

  if (is.infinite(term_years)) {
    probs <- death_counts / l0
    categories <- seq_along(death_counts)
  } else {
    probs <- c(death_counts / l0, survival_to_horizon / l0)
    categories <- seq_along(probs)
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

  if (is.infinite(term_years)) {
    curtate_lifetime <- sampled - 1L
    died_within_horizon <- rep(TRUE, n_sim)
  } else {
    survival_category <- length(probs)
    died_within_horizon <- sampled != survival_category
    curtate_lifetime <- ifelse(died_within_horizon, sampled - 1L, horizon)
  }

  tibble::tibble(
    curtate_lifetime = as.integer(curtate_lifetime),
    death_age = ifelse(
      died_within_horizon,
      age + as.integer(curtate_lifetime) + 1L,
      NA_real_
    ),
    died_within_horizon = died_within_horizon
  )
}


.annuity_payment_count <- function(
    curtate_lifetime,
    died_within_horizon,
    term_years,
    payments_per_year,
    timing
) {
  k_month <- as.integer(payments_per_year)

  if (is.infinite(term_years)) {
    years_lived <- curtate_lifetime
  } else {
    years_lived <- pmin(curtate_lifetime, term_years)
  }

  if (timing == "due") {
    pmax(years_lived * k_month + 1L, 0L)
  } else {
    pmax(years_lived * k_month, 0L)
  }
}


.annuity_factor_count <- function(
    n_payments,
    i_effective,
    payments_per_year,
    timing
) {
  i_period <- (1 + i_effective)^(1 / payments_per_year) - 1
  v_period <- 1 / (1 + i_period)

  out <- numeric(length(n_payments))

  positive <- n_payments > 0L

  if (!any(positive)) {
    return(out)
  }

  if (abs(i_period) < 1e-12) {
    out[positive] <- n_payments[positive] / payments_per_year
    return(out)
  }

  due_factor <- (1 / payments_per_year) *
    (1 - v_period^n_payments[positive]) /
    (1 - v_period)

  if (timing == "due") {
    out[positive] <- due_factor
  } else {
    out[positive] <- v_period * due_factor
  }

  out
}
