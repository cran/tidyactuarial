#' Simulate future lifetimes from a life table
#'
#' Simulates curtate and, optionally, complete future lifetimes from a life table
#' containing one-year death probabilities, using compact actuarial notation.
#'
#' This function is designed as the simulation engine for Monte Carlo life
#' contingency calculations in \code{tidyactuarial}. It generates simulated
#' values of the curtate future lifetime \eqn{K_x} and a simulated complete
#' future lifetime \eqn{T_x}. The output is a tidy tibble that can be used
#' naturally with pipes and downstream Monte Carlo functions.
#'
#' For a life aged \eqn{x}, the curtate future lifetime \eqn{K_x} follows
#'
#' \deqn{
#'   P(K_x = k) = {}_k p_x q_{x+k},
#' }
#'
#' where \eqn{{}_k p_x} is the probability of surviving \eqn{k} complete years
#' from age \eqn{x}, and \eqn{q_{x+k}} is the one-year probability of death at
#' attained age \eqn{x+k}.
#'
#' @param lt A data frame or tibble containing the life table.
#' @param x Numeric scalar. Initial actuarial age of the individual.
#' @param n_sim Positive integer. Number of Monte Carlo simulations. Default is
#'   \code{10000}.
#' @param x_col Character string. Name of the age column in \code{lt}. Default
#'   is \code{"x"}.
#' @param qx_col Character string. Name of the one-year death probability
#'   column in \code{lt}. Default is \code{"qx"}.
#' @param method Character string specifying the simulation method for
#'   \eqn{K_x}. Available options are \code{"inverse"}, \code{"multinomial"},
#'   and \code{"antithetic"}.
#' @param frac Character string specifying how the fractional part of the
#'   complete future lifetime is generated within the year of death. Available
#'   options are \code{"udd"}, \code{"constant_force"}, and \code{"none"}.
#' @param seed Optional integer seed for reproducibility. Default is
#'   \code{NULL}.
#' @param include_distribution Logical. If \code{TRUE}, the probability mass
#'   function used to simulate \eqn{K_x} is attached as a list-column named
#'   \code{distribution}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} denotes the life table, \code{x} denotes the
#' actuarial age, and \code{frac} denotes the fractional-age simulation
#' assumption.
#'
#' The function first constructs the conditional distribution of \eqn{K_x} from
#' the selected age onward. Then it generates simulated values according to the
#' selected method.
#'
#' The available simulation methods are:
#' \itemize{
#'   \item \code{"inverse"}: inverse transform simulation using the cumulative
#'   distribution of \eqn{K_x}.
#'   \item \code{"multinomial"}: direct sampling from the probability mass
#'   function of \eqn{K_x}.
#'   \item \code{"antithetic"}: inverse transform simulation using antithetic
#'   uniforms \eqn{U} and \eqn{1-U}.
#' }
#'
#' The \code{frac} argument controls the simulated complete future lifetime
#' \eqn{T_x}. If \code{frac = "udd"}, a uniform fractional lifetime is added to
#' \eqn{K_x}. If \code{frac = "constant_force"}, the fractional lifetime within
#' the year of death is generated under a constant force of mortality
#' assumption conditional on death during that year. If \code{frac = "none"},
#' the complete lifetime is returned as \code{NA}.
#'
#' The probabilities are normalized internally to handle finite life tables. If
#' the life table is truncated, the simulated distribution is conditional on
#' death occurring within the available range of ages. In practical work, it is
#' recommended that the last available death probability be equal to 1.
#'
#' @return A tibble with one row per simulation and columns:
#' \describe{
#'   \item{sim_id}{Simulation identifier.}
#'   \item{x}{Initial actuarial age.}
#'   \item{age}{Compatibility column equal to \code{x}.}
#'   \item{method}{Simulation method used.}
#'   \item{frac}{Fractional-age simulation assumption.}
#'   \item{fractional}{Compatibility column equal to \code{frac}.}
#'   \item{Kx}{Simulated curtate future lifetime.}
#'   \item{Tx}{Simulated complete future lifetime. If \code{frac = "none"},
#'   this column contains \code{NA}.}
#' }
#'
#' If \code{include_distribution = TRUE}, the output also includes a list-column
#' named \code{distribution} containing the probability mass function used for
#' the simulation.
#'
#' @seealso \code{\link{mc_insurance}}, \code{\link{mc_annuity}},
#'   \code{\link{mc_premium}}, \code{\link{mc_loss}}
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
#'   x = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' # Basic simulation using inverse transform sampling
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     method = "inverse",
#'     seed = 123
#'   )
#'
#' # Antithetic simulation
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     method = "antithetic",
#'     seed = 123
#'   )
#'
#' # Returning the distribution used for simulation
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     include_distribution = TRUE,
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
    frac = c("udd", "constant_force", "none"),
    seed = NULL,
    include_distribution = FALSE
) {
  method <- match.arg(method)
  frac <- match.arg(frac)

  if (!is.data.frame(lt)) {
    stop("`lt` must be a data frame or tibble.", call. = FALSE)
  }

  if (missing(x) ||
      !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x)) {
    stop("`x` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(n_sim) ||
      length(n_sim) != 1L ||
      is.na(n_sim) ||
      !is.finite(n_sim) ||
      n_sim <= 0 ||
      abs(n_sim - round(n_sim)) > 1e-10) {
    stop("`n_sim` must be a single positive integer.", call. = FALSE)
  }

  n_sim <- as.integer(round(n_sim))

  if (!is.character(x_col) || length(x_col) != 1L || is.na(x_col)) {
    stop("`x_col` must be a single character string.", call. = FALSE)
  }

  if (!is.character(qx_col) || length(qx_col) != 1L || is.na(qx_col)) {
    stop("`qx_col` must be a single character string.", call. = FALSE)
  }

  if (!x_col %in% names(lt)) {
    stop("`x_col` must identify a column in `lt`.", call. = FALSE)
  }

  if (!qx_col %in% names(lt)) {
    stop("`qx_col` must identify a column in `lt`.", call. = FALSE)
  }

  if (!is.numeric(lt[[x_col]])) {
    stop("The column identified by `x_col` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(lt[[qx_col]])) {
    stop("The column identified by `qx_col` must be numeric.", call. = FALSE)
  }

  if (!is.logical(include_distribution) ||
      length(include_distribution) != 1L ||
      is.na(include_distribution)) {
    stop("`include_distribution` must be a logical scalar.", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed <= 0 ||
        abs(seed - round(seed)) > 1e-10) {
      stop("`seed` must be a single positive integer or NULL.", call. = FALSE)
    }

    set.seed(as.integer(round(seed)))
  }

  age_vec <- lt[[x_col]]
  qx_vec <- lt[[qx_col]]

  keep <- age_vec >= x

  if (!any(keep)) {
    stop(
      "No ages greater than or equal to `x` were found in `lt`.",
      call. = FALSE
    )
  }

  dist <- data.frame(
    attained_age = age_vec[keep],
    qx = qx_vec[keep]
  )

  dist <- dist[order(dist$attained_age), , drop = FALSE]

  if (any(is.na(dist$attained_age)) || any(!is.finite(dist$attained_age))) {
    stop("The selected `x_col` contains non-finite or missing ages.", call. = FALSE)
  }

  if (any(is.na(dist$qx))) {
    stop("The selected `qx_col` contains missing values.", call. = FALSE)
  }

  if (any(!is.finite(dist$qx)) || any(dist$qx < 0 | dist$qx > 1)) {
    stop("All death probabilities in `qx_col` must be finite and between 0 and 1.",
         call. = FALSE)
  }

  dist$k <- dist$attained_age - x

  if (any(dist$k < 0) || any(abs(dist$k - round(dist$k)) > 1e-10)) {
    stop("The selected ages must define nonnegative integer curtate lifetimes.",
         call. = FALSE)
  }

  dist$k <- as.integer(round(dist$k))
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

    index[index < 1L] <- 1L
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

    index[index < 1L] <- 1L
    index[index > nrow(dist)] <- nrow(dist)

    Kx <- dist$k[index]
  }

  Tx <- rep(NA_real_, n_sim)

  if (frac == "udd") {
    Tx <- Kx + stats::runif(n_sim)
  }

  if (frac == "constant_force") {
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
    x = x,
    age = x,
    method = method,
    frac = frac,
    fractional = frac,
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
