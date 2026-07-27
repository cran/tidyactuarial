#' Add an insurance benefit specification to a life contract
#'
#' Adds the benefit side of a single-life insurance contract to an existing
#' [life_contract()] object. This function stores the specification only; it
#' does not calculate an actuarial present value.
#'
#' @param contract A single-life `tidyact_life_contract` object.
#' @param type Character string. Insurance type: `"whole"`, `"term"`,
#'   `"endowment"`, or `"variable_k"`.
#' @param benefit Benefit amount. For standard products, a single nonnegative
#'   numeric value. For `type = "variable_k"`, a numeric vector or a function
#'   of time may be supplied.
#' @param n Insurance term in years. Use `Inf` for whole-life insurance.
#'   A finite positive value is required for term and endowment insurance.
#' @param h Nonnegative integer deferment period in years.
#' @param k Positive integer. Frequency used for benefits evaluated within
#'   each year when `type = "variable_k"`.
#' @param frac Fractional-age assumption. One of `"UDD"`, `"CF"`, `"CML"`,
#'   or `"Balducci"`.
#'
#' @return The original contract with an `insurance` component.
#'
#' @examples
#' lt <- data.frame(
#'   x = 40:90,
#'   lx = round(100000 * exp(-0.018 * (0:50)^1.35))
#' )
#'
#' life_contract(
#'   lt = lt,
#'   lives = "single",
#'   x = 40,
#'   i = 0.05
#' ) |>
#'   add_insurance(
#'     type = "term",
#'     benefit = 100000,
#'     n = 20
#'   )
#'
#' @family life-contingencies
#'
#' @export
add_insurance <- function(
    contract,
    type = c("whole", "term", "endowment", "variable_k"),
    benefit = 1,
    n = Inf,
    h = 0L,
    k = 1L,
    frac = c("UDD", "CF", "CML", "Balducci")
) {
  if (!inherits(contract, "tidyact_life_contract")) {
    stop("`contract` must be a `tidyact_life_contract` object.", call. = FALSE)
  }

  if (!identical(contract$lives, "single")) {
    stop(
      "`add_insurance()` currently supports only single-life contracts.",
      call. = FALSE
    )
  }

  type <- match.arg(type)
  frac <- match.arg(frac)

  if (identical(frac, "CML")) {
    frac <- "CF"
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n <= 0 ||
      (!is.infinite(n) && !is.finite(n))) {
    stop("`n` must be `Inf` or a single positive finite value.", call. = FALSE)
  }

  if (type %in% c("term", "endowment") && is.infinite(n)) {
    stop("`n` must be finite for term and endowment insurance.", call. = FALSE)
  }

  if (!is.numeric(h) ||
      length(h) != 1L ||
      is.na(h) ||
      !is.finite(h) ||
      h < 0 ||
      abs(h - round(h)) > 1e-10) {
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

  if (type == "variable_k") {
    if (!(is.function(benefit) ||
          (is.numeric(benefit) &&
           length(benefit) >= 1L &&
           all(is.finite(benefit)) &&
           all(benefit >= 0)))) {
      stop(
        "For `type = 'variable_k'`, `benefit` must be a nonnegative numeric ",
        "vector or a function of time.",
        call. = FALSE
      )
    }

    if (is.function(benefit) && is.infinite(n)) {
      stop(
        "For a functional variable benefit, `n` must be finite.",
        call. = FALSE
      )
    }
  } else {
    if (!is.numeric(benefit) ||
        length(benefit) != 1L ||
        is.na(benefit) ||
        !is.finite(benefit) ||
        benefit < 0) {
      stop(
        "For standard products, `benefit` must be a single finite ",
        "nonnegative number.",
        call. = FALSE
      )
    }
  }

  contract$insurance <- list(
    type = type,
    benefit = benefit,
    n = n,
    h = as.integer(round(h)),
    k = as.integer(round(k)),
    frac = frac
  )

  contract
}


#' Add a contingent premium-payment schedule to a life contract
#'
#' Adds the premium-payment side of a single-life insurance contract. The
#' schedule defines when premiums are payable while the insured is alive.
#' This function stores the specification only; [premium_x()] later applies
#' the equivalence principle.
#'
#' @param contract A single-life `tidyact_life_contract` object.
#' @param k Positive integer. Number of premium payments per year.
#' @param timing Timing of premium payments: `"due"` or `"immediate"`.
#' @param premium_start Start of premium payments: `"issue"` or `"deferred"`.
#' @param n_prem Premium-paying term in years. Use `NULL` to infer the term
#'   from the insurance component. A finite fractional value is allowed when
#'   `n_prem * k` is an integer.
#' @param woolhouse Woolhouse approximation order: `"none"`, `"first"`, or
#'   `"second"`.
#' @param frac Fractional-age assumption for exact k-thly valuation. One of
#'   `"UDD"`, `"CF"`, `"CML"`, or `"Balducci"`. If `NULL`, the insurance
#'   component's assumption is used when available; otherwise `"UDD"`.
#'
#' @return The original contract with a `premium_schedule` component.
#'
#' @examples
#' lt <- data.frame(
#'   x = 40:90,
#'   lx = round(100000 * exp(-0.018 * (0:50)^1.35))
#' )
#'
#' life_contract(
#'   lt = lt,
#'   lives = "single",
#'   x = 40,
#'   i = 0.05
#' ) |>
#'   add_insurance(
#'     type = "term",
#'     benefit = 100000,
#'     n = 20
#'   ) |>
#'   add_premium_schedule(
#'     k = 12,
#'     n_prem = 10,
#'     timing = "due"
#'   )
#'
#' @family life-contingencies
#'
#' @export
add_premium_schedule <- function(
    contract,
    k = 1L,
    timing = c("due", "immediate"),
    premium_start = c("issue", "deferred"),
    n_prem = NULL,
    woolhouse = c("none", "first", "second"),
    frac = NULL
) {
  if (!inherits(contract, "tidyact_life_contract")) {
    stop("`contract` must be a `tidyact_life_contract` object.", call. = FALSE)
  }

  if (!identical(contract$lives, "single")) {
    stop(
      "`add_premium_schedule()` currently supports only single-life contracts.",
      call. = FALSE
    )
  }

  timing <- match.arg(timing)
  premium_start <- match.arg(premium_start)
  woolhouse <- match.arg(woolhouse)

  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k < 1 ||
      abs(k - round(k)) > 1e-10) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  k <- as.integer(round(k))

  if (!is.null(n_prem)) {
    if (!is.numeric(n_prem) ||
        length(n_prem) != 1L ||
        is.na(n_prem) ||
        n_prem <= 0 ||
        (!is.infinite(n_prem) && !is.finite(n_prem))) {
      stop(
        "`n_prem` must be `NULL`, `Inf`, or a single positive finite value.",
        call. = FALSE
      )
    }

    if (!is.infinite(n_prem)) {
      n_payments <- n_prem * k

      if (abs(n_payments - round(n_payments)) > 1e-10) {
        stop(
          "`n_prem * k` must be an integer so that the schedule contains ",
          "a whole number of premium payments.",
          call. = FALSE
        )
      }

      n_prem <- round(n_payments) / k
    }
  }

  if (is.null(frac)) {
    frac <- if (!is.null(contract$insurance$frac)) {
      contract$insurance$frac
    } else {
      "UDD"
    }
  }

  if (!is.character(frac) || length(frac) != 1L || is.na(frac)) {
    stop("`frac` must be `NULL` or a single character string.", call. = FALSE)
  }

  frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))

  if (identical(frac, "CML")) {
    frac <- "CF"
  }

  contract$premium_schedule <- list(
    k = k,
    timing = timing,
    premium_start = premium_start,
    n_prem = n_prem,
    woolhouse = woolhouse,
    frac = frac
  )

  contract
}
