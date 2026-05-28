#' Simulate future lifetimes from a life table
#'
#' Simulates curtate and, optionally, complete future lifetimes from a life
#' table containing one-year death probabilities.
#'
#' This function is designed as the simulation engine for Monte Carlo life
#' contingency calculations in `tidyactuarial`. It generates simulated values
#' of the curtate future lifetime \eqn{K_x} and a simulated complete future
#' lifetime \eqn{T_x}. The output is a tidy tibble that can be used naturally
#' with pipes and downstream functions such as [mc_insurance()],
#' [mc_annuity()], [mc_premium()], and [mc_loss()].
#'
#' For a life aged \eqn{x}, the curtate future lifetime \eqn{K_x} follows
#'
#' \deqn{
#'   P(K_x = k) = {}_k p_x q_{x+k},
#' }
#'
#' where \eqn{{}_k p_x} is the probability of surviving \eqn{k} complete years
#' from age \eqn{x}, and \eqn{q_{x+k}} is the one-year probability of death
#' at attained age \eqn{x+k}.
#'
#' @param data A data frame or tibble containing the life table.
#' @param age Numeric scalar. Initial age \eqn{x} of the individual.
#' @param n_sim Positive integer. Number of Monte Carlo simulations.
#'   Default is `10000`.
#' @param age_col Character string. Name of the age column in `data`.
#'   Default is `"age"`.
#' @param qx_col Character string. Name of the one-year death probability
#'   column in `data`. Default is `"qx"`.
#' @param method Character string specifying the simulation method for
#'   \eqn{K_x}. Available options are `"inverse"`, `"multinomial"`, and
#'   `"antithetic"`. Default is `"inverse"`.
#' @param fractional Character string specifying how the fractional part of
#'   the complete future lifetime is generated within the year of death.
#'   Available options are `"udd"`, `"constant_force"`, and `"none"`.
#'   Default is `"udd"`.
#' @param seed Optional integer seed for reproducibility. Default is `NULL`.
#' @param include_distribution Logical. If `TRUE`, the probability mass
#'   function used to simulate \eqn{K_x} is attached as a list-column named
#'   `distribution`. Default is `FALSE`.
#'
#' @details
#' The function first constructs the conditional distribution of \eqn{K_x}
#' from the selected age onward. Then it generates simulated values according
#' to the selected method.
#'
#' The available simulation methods are:
#'
#' * `"inverse"`: inverse transform simulation using the cumulative
#'   distribution of \eqn{K_x}.
#' * `"multinomial"`: direct sampling from the probability mass function of
#'   \eqn{K_x}.
#' * `"antithetic"`: inverse transform simulation using antithetic uniforms
#'   \eqn{U} and \eqn{1-U}, a basic variance reduction technique.
#'
#' The `fractional` argument controls the simulated complete future lifetime
#' \eqn{T_x}. If `fractional = "udd"`, a uniform fractional lifetime is added
#' to \eqn{K_x}. If `fractional = "constant_force"`, the fractional lifetime
#' within the year of death is generated under a constant force of mortality
#' assumption conditional on death during that year. If `fractional = "none"`,
#' the complete lifetime is returned as `NA`.
#'
#' The probabilities are normalized internally to handle finite life tables.
#' If the life table is truncated, the simulated distribution is conditional
#' on death occurring within the available range of ages. In practical work,
#' it is recommended that the last available death probability be equal to 1.
#'
#' @return A tibble with one row per simulation and the following columns:
#'
#' * `sim_id`: simulation identifier.
#' * `age`: initial age.
#' * `method`: simulation method used.
#' * `fractional`: fractional age assumption used for \eqn{T_x}.
#' * `Kx`: simulated curtate future lifetime.
#' * `Tx`: simulated complete future lifetime. If `fractional = "none"`,
#'   this column contains `NA`.
#'
#' If `include_distribution = TRUE`, the output also includes a list-column
#' named `distribution` containing the probability mass function used for the
#' simulation.
#'
#' @seealso
#' [mc_insurance()], [mc_annuity()], [mc_premium()], [mc_loss()]
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). *Actuarial Mathematics*. Second Edition.
#' Society of Actuaries.
#'
#' @examples
#' life_table <- tibble::tibble(
#'   age = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' # Basic simulation using inverse transform sampling
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     method = "inverse",
#'     seed = 123
#'   )
#'
#' # Antithetic simulation
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     method = "antithetic",
#'     seed = 123
#'   )
#'
#' # Returning the distribution used for simulation
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 100,
#'     include_distribution = TRUE,
#'     seed = 123
#'   )
#'
#' # Full pipeline with a life insurance present value
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     benefit = 1
#'   )
#'
#' # Complete future lifetime for benefits payable at the moment of death
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     payment_timing = "moment_of_death",
#'     tx_col = "Tx"
#'   )
#'
#' @export
simulate_lifetime <- function(data,
                              age,
                              n_sim = 10000,
                              age_col = "age",
                              qx_col = "qx",
                              method = c("inverse", "multinomial", "antithetic"),
                              fractional = c("udd", "constant_force", "none"),
                              seed = NULL,
                              include_distribution = FALSE) {
  method <- match.arg(method)
  fractional <- match.arg(fractional)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(age, "age")
  .mc_assert_positive_integer(n_sim, "n_sim")
  .mc_assert_numeric_column(data, age_col, "age_col")
  .mc_assert_numeric_column(data, qx_col, "qx_col")

  if (!is.logical(include_distribution) || length(include_distribution) != 1) {
    stop("`include_distribution` must be a logical scalar.", call. = FALSE)
  }

  if (!is.null(seed)) {
    .mc_assert_positive_integer(seed, "seed")
    set.seed(seed)
  }

  age_vec <- data[[age_col]]
  qx_vec <- data[[qx_col]]

  keep <- age_vec >= age

  if (!any(keep)) {
    stop(
      "No ages greater than or equal to `age` were found in `data`.",
      call. = FALSE
    )
  }

  dist <- data.frame(
    attained_age = age_vec[keep],
    qx = qx_vec[keep]
  )

  dist <- dist[order(dist$attained_age), , drop = FALSE]

  if (any(is.na(dist$qx))) {
    stop("The selected `qx_col` contains missing values.", call. = FALSE)
  }

  if (any(dist$qx < 0 | dist$qx > 1)) {
    stop("All death probabilities in `qx_col` must be between 0 and 1.",
         call. = FALSE)
  }

  dist$k <- dist$attained_age - age
  dist$px <- 1 - dist$qx

  survival_to_k <- c(1, cumprod(dist$px[-length(dist$px)]))
  dist$survival_to_k <- survival_to_k
  dist$prob <- dist$survival_to_k * dist$qx

  total_prob <- sum(dist$prob)

  if (total_prob <= 0 || is.na(total_prob)) {
    stop(
      "The simulated lifetime distribution has zero total probability.",
      call. = FALSE
    )
  }

  if (abs(total_prob - 1) > 1e-8) {
    warning(
      "The life table appears to be finite or truncated. ",
      "The simulated probabilities were normalized over the available ages.",
      call. = FALSE
    )
  }

  dist$prob <- dist$prob / total_prob
  dist$cdf <- cumsum(dist$prob)
  dist$cdf[length(dist$cdf)] <- 1

  if (method == "multinomial") {
    Kx <- sample(
      x = dist$k,
      size = n_sim,
      replace = TRUE,
      prob = dist$prob
    )
  }

  if (method == "inverse") {
    u <- stats::runif(n_sim)

    index <- findInterval(
      x = u,
      vec = c(0, dist$cdf),
      rightmost.closed = TRUE
    )

    index[index < 1] <- 1
    index[index > nrow(dist)] <- nrow(dist)

    Kx <- dist$k[index]
  }

  if (method == "antithetic") {
    n_half <- ceiling(n_sim / 2)

    u <- stats::runif(n_half)
    u <- c(u, 1 - u)[seq_len(n_sim)]

    index <- findInterval(
      x = u,
      vec = c(0, dist$cdf),
      rightmost.closed = TRUE
    )

    index[index < 1] <- 1
    index[index > nrow(dist)] <- nrow(dist)

    Kx <- dist$k[index]
  }

  Tx <- rep(NA_real_, n_sim)

  if (fractional == "udd") {
    Tx <- Kx + stats::runif(n_sim)
  }

  if (fractional == "constant_force") {
    qx_by_k <- stats::setNames(dist$qx, dist$k)
    q_death <- unname(qx_by_k[as.character(Kx)])

    u_frac <- stats::runif(n_sim)

    Tx <- Kx + ifelse(
      q_death <= 0,
      0,
      log(1 - u_frac * q_death) / log(1 - q_death)
    )
  }

  out <- tibble::tibble(
    sim_id = seq_len(n_sim),
    age = age,
    method = method,
    fractional = fractional,
    Kx = Kx,
    Tx = Tx
  )

  if (include_distribution) {
    dist_tbl <- tibble::tibble(
      attained_age = dist$attained_age,
      k = dist$k,
      qx = dist$qx,
      px = dist$px,
      survival_to_k = dist$survival_to_k,
      prob = dist$prob,
      cdf = dist$cdf
    )

    out <- dplyr::mutate(
      out,
      distribution = list(dist_tbl)
    )
  }

  out
}
