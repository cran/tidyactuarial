#' Simulate future lifetimes from a life table
#'
#' Simulates curtate and, optionally, complete future lifetimes from a life
#' table containing one-year death probabilities.
#'
#' @param lt A data frame or tibble containing the life table.
#' @param x Nonnegative integer initial actuarial age. The age must appear
#'   explicitly in the column selected by \code{x_col}.
#' @param n_sim Positive integer number of simulations.
#' @param x_col Character scalar naming the age column.
#' @param qx_col Character scalar naming the one-year death-probability column.
#' @param method Simulation method for the curtate future lifetime:
#'   \code{"inverse"}, \code{"multinomial"}, or \code{"antithetic"}.
#' @param frac Fractional-age assumption used inside the year of death:
#'   \code{"udd"}, \code{"cml"}, \code{"balducci"}, or \code{"none"}.
#'   The legacy value \code{"constant_force"} is accepted as an alias for
#'   \code{"cml"}.
#' @param seed Optional nonnegative integer seed. When supplied, the caller's
#'   random-number state is restored after simulation.
#' @param include_distribution Logical scalar. If \code{TRUE}, the probability
#'   distribution used for simulation is attached as a list-column.
#' @param truncation Treatment of a life table whose available death
#'   probabilities do not exhaust the lifetime distribution. With
#'   \code{"conditional"}, the historical behavior is retained and the
#'   probabilities are normalized conditional on death within the available
#'   ages, with a warning. With \code{"error"}, the simulation stops.
#' @param tol Nonnegative numeric tolerance used for age-grid and probability
#'   checks.
#'
#' @details
#' For a life aged \eqn{x}, the curtate future lifetime has probability mass
#' function
#' \deqn{
#' P(K_x=k)={}_kp_xq_{x+k}.
#' }
#' The selected life-table ages must therefore start at \eqn{x} and proceed in
#' consecutive one-year steps. Missing or duplicated ages are rejected because
#' the survival recursion would otherwise no longer represent
#' \eqn{{}_kp_x}.
#'
#' Under UDD, the fractional part of \eqn{T_x} is uniform conditional on the
#' year of death. Under CML, it follows the conditional constant-force
#' distribution. Under Balducci, it follows the corresponding conditional
#' Balducci distribution.
#'
#' @return A tibble with one row per simulation and standardized columns
#' including \code{sim_id}, \code{Kx}, \code{curtate_lifetime}, \code{Tx},
#' \code{complete_lifetime}, \code{death_age}, \code{method}, and \code{frac}.
#' The logical column \code{distribution_conditioned} indicates whether a
#' truncated table was normalized conditionally.
#'
#' @seealso \code{\link{simulate_lifetimes}}, \code{\link{mc_insurance}},
#'   \code{\link{mc_annuity}}, \code{\link{mc_premium}},
#'   \code{\link{mc_loss}}, \code{\link{mc_reserve}}
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). \emph{Actuarial Mathematics}. Second Edition.
#' Society of Actuaries.
#'
#' @family monte-carlo
#'
#' @examples
#' lt <- tibble::tibble(
#'   x = 40:43,
#'   qx = c(0.10, 0.20, 0.30, 1)
#' )
#'
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     method = "inverse",
#'     frac = "udd",
#'     seed = 123
#'   )
#'
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     method = "antithetic",
#'     frac = "cml",
#'     seed = 123
#'   )
#'
#' @export
simulate_lifetime <- function(
    lt,
    x,
    n_sim = 10000L,
    x_col = "x",
    qx_col = "qx",
    method = c("inverse", "multinomial", "antithetic"),
    frac = c(
      "udd",
      "cml",
      "balducci",
      "none",
      "constant_force"
    ),
    seed = NULL,
    include_distribution = FALSE,
    truncation = c("conditional", "error"),
    tol = 1e-10
) {
  method <- match.arg(method)
  frac <- match.arg(frac)
  truncation <- match.arg(truncation)

  if (identical(frac, "constant_force")) {
    frac <- "cml"
  }

  if (!is.data.frame(lt)) {
    stop("`lt` must be a data frame or tibble.", call. = FALSE)
  }

  if (nrow(lt) == 0L) {
    stop("`lt` must contain at least one life-table row.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(
    x,
    "x",
    min = 0
  )

  if (abs(x - round(x)) > tol) {
    stop("`x` must be a nonnegative integer age.", call. = FALSE)
  }

  x <- as.integer(round(x))

  .mc_assert_positive_integer(n_sim, "n_sim")

  if (n_sim > .Machine$integer.max) {
    stop(
      "`n_sim` exceeds the supported integer range.",
      call. = FALSE
    )
  }

  n_sim <- as.integer(round(n_sim))

  .mc_assert_character_scalar(x_col, "x_col")
  .mc_assert_character_scalar(qx_col, "qx_col")

  if (identical(x_col, qx_col)) {
    stop("`x_col` and `qx_col` must identify different columns.", call. = FALSE)
  }

  .mc_assert_numeric_column(
    lt,
    x_col,
    "x_col"
  )

  .mc_assert_numeric_column(
    lt,
    qx_col,
    "qx_col"
  )

  if (!is.logical(include_distribution) ||
      length(include_distribution) != 1L ||
      is.na(include_distribution)) {
    stop(
      "`include_distribution` must be a logical scalar.",
      call. = FALSE
    )
  }

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol < 0) {
    stop(
      "`tol` must be a single nonnegative finite number.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed < 0 ||
        abs(seed - round(seed)) > tol ||
        seed > .Machine$integer.max) {
      stop(
        "`seed` must be `NULL` or a nonnegative integer within the ",
        "supported integer range.",
        call. = FALSE
      )
    }

    seed <- as.integer(round(seed))
  }

  age_vec <- lt[[x_col]]
  qx_vec <- lt[[qx_col]]

  if (anyNA(age_vec) ||
      any(!is.finite(age_vec)) ||
      any(age_vec < 0) ||
      any(abs(age_vec - round(age_vec)) > tol)) {
    stop(
      "The column selected by `x_col` must contain finite nonnegative ",
      "integer ages.",
      call. = FALSE
    )
  }

  if (anyDuplicated(age_vec)) {
    stop(
      "The column selected by `x_col` must not contain duplicated ages.",
      call. = FALSE
    )
  }

  if (anyNA(qx_vec) ||
      any(!is.finite(qx_vec)) ||
      any(qx_vec < 0 | qx_vec > 1)) {
    stop(
      "The column selected by `qx_col` must contain finite probabilities ",
      "between 0 and 1.",
      call. = FALSE
    )
  }

  age_vec <- as.integer(round(age_vec))
  order_age <- order(age_vec)
  age_vec <- age_vec[order_age]
  qx_vec <- qx_vec[order_age]

  idx_x <- match(x, age_vec)

  if (is.na(idx_x)) {
    stop(
      "`x` must appear explicitly in the age column selected by `x_col`.",
      call. = FALSE
    )
  }

  attained_age <- age_vec[idx_x:length(age_vec)]
  qx_selected <- qx_vec[idx_x:length(qx_vec)]

  expected_ages <- seq.int(
    from = x,
    length.out = length(attained_age)
  )

  if (!identical(attained_age, expected_ages)) {
    stop(
      "The selected life-table ages must be consecutive one-year ages ",
      "starting at `x`.",
      call. = FALSE
    )
  }

  px_selected <- 1 - qx_selected

  survival_to_k <- c(
    1,
    if (length(px_selected) > 1L) {
      cumprod(px_selected[-length(px_selected)])
    } else {
      numeric(0)
    }
  )

  prob_unconditional <- survival_to_k * qx_selected
  total_prob <- sum(prob_unconditional)

  if (!is.finite(total_prob) || total_prob <= 0) {
    stop(
      "The simulated lifetime distribution has zero or invalid total ",
      "probability.",
      call. = FALSE
    )
  }

  conditioned <- abs(total_prob - 1) > tol

  if (conditioned && identical(truncation, "error")) {
    stop(
      "The available life table does not exhaust the lifetime ",
      "distribution. Extend the table to a terminal `qx = 1`, or use ",
      "`truncation = \"conditional\"` explicitly.",
      call. = FALSE
    )
  }

  if (conditioned) {
    warning(
      "The life table is truncated. Simulation is conditional on death ",
      "within the available ages.",
      call. = FALSE
    )
  }

  prob <- prob_unconditional / total_prob

  # Protect the inverse-transform grid against platform-specific floating-point
  # accumulation. On some systems, cumsum(prob) can end slightly above 1
  # (for example, 1 + 7e-16). Replacing only the last value by 1 would then
  # make the CDF decrease at its endpoint and `findInterval()` would fail.
  cdf <- cumsum(prob)
  cdf <- pmin(cdf, 1)
  cdf <- cummax(cdf)
  cdf[[length(cdf)]] <- 1

  dist <- tibble::tibble(
    attained_age = attained_age,
    k = seq_along(attained_age) - 1L,
    qx = qx_selected,
    px = px_selected,
    survival_to_k = survival_to_k,
    prob_unconditional = prob_unconditional,
    prob = prob,
    cdf = cdf
  )

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

    set.seed(seed)
  }

  antithetic_uniforms <- function(size) {
    n_base <- ceiling(size / 2)
    u_base <- stats::runif(n_base)

    as.vector(
      rbind(
        u_base,
        1 - u_base
      )
    )[seq_len(size)]
  }

  inverse_index <- function(u) {
    index <- findInterval(
      x = u,
      vec = c(0, cdf),
      rightmost.closed = TRUE
    )

    pmin(
      pmax(index, 1L),
      length(cdf)
    )
  }

  if (identical(method, "multinomial")) {
    index <- sample.int(
      n = nrow(dist),
      size = n_sim,
      replace = TRUE,
      prob = prob
    )
  } else {
    u <- if (identical(method, "antithetic")) {
      antithetic_uniforms(n_sim)
    } else {
      stats::runif(n_sim)
    }

    index <- inverse_index(u)
  }

  Kx <- dist$k[index]
  q_death <- dist$qx[index]
  px_death <- 1 - q_death

  Tx <- rep(NA_real_, n_sim)

  if (!identical(frac, "none")) {
    u_frac <- if (identical(method, "antithetic")) {
      antithetic_uniforms(n_sim)
    } else {
      stats::runif(n_sim)
    }

    fractional_year <- if (identical(frac, "udd")) {
      u_frac
    } else if (identical(frac, "cml")) {
      out <- numeric(n_sim)
      degenerate <- px_death <= tol
      out[degenerate] <- 0

      regular <- !degenerate

      out[regular] <-
        log1p(-u_frac[regular] * q_death[regular]) /
        log(px_death[regular])

      out
    } else {
      out <- numeric(n_sim)
      degenerate <- px_death <= tol
      out[degenerate] <- 0

      regular <- !degenerate

      out[regular] <-
        u_frac[regular] * px_death[regular] /
        (1 - u_frac[regular] * q_death[regular])

      out
    }

    fractional_year <- pmin(
      pmax(fractional_year, 0),
      1 - .Machine$double.eps
    )

    Tx <- Kx + fractional_year
  }

  out <- tibble::tibble(
    sim_id = seq_len(n_sim),
    simulation_id = seq_len(n_sim),
    x = x,
    age = x,
    method = method,
    frac = frac,
    fractional = frac,
    Kx = as.integer(Kx),
    curtate_lifetime = as.integer(Kx),
    Tx = Tx,
    complete_lifetime = Tx,
    death_age = x + as.integer(Kx) + 1L,
    qx_at_death_year = q_death,
    distribution_conditioned = conditioned
  )

  if (include_distribution) {
    out <- out |>
      dplyr::mutate(
        distribution = rep(
          list(dist),
          n_sim
        )
      )
  }

  out
}
