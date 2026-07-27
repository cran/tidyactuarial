#' Add an insurance benefit specification to a life contract
#'
#' Adds the benefit side of a single-life or two-life insurance contract to an
#' existing [life_contract()] object. This function stores the specification
#' only; it does not calculate an actuarial present value.
#'
#' For two-life contracts, the insurance status is inferred from
#' `contract$lives`: `"joint"` becomes joint-life status and
#' `"last_survivor"` becomes last-survivor status.
#'
#' @param contract A `tidyact_life_contract` object.
#' @param type Insurance type. Single-life contracts support `"whole"`,
#'   `"term"`, `"endowment"`, and `"variable_k"`. Two-life contracts support
#'   `"whole"`, `"term"`, `"endowment"`, and `"pure_endowment"`.
#' @param benefit Benefit amount. For standard products, a single nonnegative
#'   numeric value. For single-life `type = "variable_k"`, a numeric vector or
#'   a function of time may be supplied.
#' @param n Insurance term in years. Use `Inf` for whole-life insurance.
#' @param h Nonnegative integer deferment period in years.
#' @param k Positive integer. Benefit frequency for single-life
#'   `type = "variable_k"`.
#' @param frac Fractional-age assumption: `"UDD"`, `"CF"`, `"CML"`, or
#'   `"Balducci"`.
#'
#' @return The original contract with an `insurance` component.
#'
#' @family life-contingencies
#'
#' @export
add_insurance <- function(
    contract,
    type = c(
      "whole",
      "term",
      "endowment",
      "pure_endowment",
      "variable_k"
    ),
    benefit = 1,
    n = Inf,
    h = 0L,
    k = 1L,
    frac = c("UDD", "CF", "CML", "Balducci")
) {
  if (!inherits(contract, "tidyact_life_contract")) {
    stop(
      "`contract` must be a `tidyact_life_contract` object.",
      call. = FALSE
    )
  }

  type <- match.arg(type)
  frac <- match.arg(frac)

  if (identical(frac, "CML")) {
    frac <- "CF"
  }

  is_single <- identical(contract$lives, "single")
  is_two_life <- contract$lives %in% c("joint", "last_survivor")

  if (!is_single && !is_two_life) {
    stop("Unsupported life-contract status.", call. = FALSE)
  }

  allowed_types <- if (is_single) {
    c("whole", "term", "endowment", "variable_k")
  } else {
    c("whole", "term", "endowment", "pure_endowment")
  }

  if (!type %in% allowed_types) {
    stop(
      "`type = '", type, "'` is not supported for a ",
      if (is_single) "single-life" else "two-life",
      " contract.",
      call. = FALSE
    )
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      n <= 0 ||
      (!is.infinite(n) && !is.finite(n))) {
    stop(
      "`n` must be `Inf` or a single positive finite value.",
      call. = FALSE
    )
  }

  if (type %in% c("term", "endowment", "pure_endowment") &&
      is.infinite(n)) {
    stop(
      "`n` must be finite for term, endowment, and pure endowment insurance.",
      call. = FALSE
    )
  }

  if (type != "variable_k" &&
      !is.infinite(n) &&
      abs(n - round(n)) > 1e-10) {
    stop(
      "Standard insurance terms must be an integer number of years.",
      call. = FALSE
    )
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

  if (identical(type, "variable_k")) {
    if (!(is.function(benefit) ||
          (is.numeric(benefit) &&
           length(benefit) >= 1L &&
           all(is.finite(benefit)) &&
           all(benefit >= 0)))) {
      stop(
        "For `type = 'variable_k'`, `benefit` must be a nonnegative ",
        "numeric vector or a function of time.",
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

  status <- if (is_single) {
    NULL
  } else if (identical(contract$lives, "joint")) {
    "joint"
  } else {
    "last"
  }

  contract$insurance <- list(
    type = type,
    benefit = benefit,
    n = if (!is.infinite(n) && type != "variable_k") {
      as.integer(round(n))
    } else {
      n
    },
    h = as.integer(round(h)),
    k = as.integer(round(k)),
    frac = frac,
    status = status
  )

  contract
}


#' Add a contingent premium-payment schedule to a life contract
#'
#' Adds the premium-payment side of a single-life or two-life insurance
#' contract. The schedule defines when premiums are payable while the selected
#' life status is in force. This function stores the specification only.
#'
#' @param contract A `tidyact_life_contract` object.
#' @param k Positive integer. Number of premium payments per year.
#' @param timing Timing of premium payments: `"due"` or `"immediate"`.
#' @param premium_start Start of premium payments: `"issue"` or `"deferred"`.
#' @param n_prem Premium-paying term in years. Use `NULL` to infer the term
#'   from the insurance component. A fractional value is allowed when
#'   `n_prem * k` is an integer.
#' @param woolhouse Woolhouse approximation order: `"none"`, `"first"`, or
#'   `"second"`.
#' @param frac Fractional-age assumption. If `NULL`, the insurance component's
#'   assumption is used when available; otherwise `"UDD"`.
#'
#' @return The original contract with a `premium_schedule` component.
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
    stop(
      "`contract` must be a `tidyact_life_contract` object.",
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

  if (!is.character(frac) ||
      length(frac) != 1L ||
      is.na(frac)) {
    stop(
      "`frac` must be `NULL` or a single character string.",
      call. = FALSE
    )
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
