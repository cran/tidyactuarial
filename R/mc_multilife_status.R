#' Construct multiple-life status lifetimes from simulated lives
#'
#' Constructs simulated multiple-life status lifetimes from simulated individual
#' lifetimes in long format.
#'
#' This function is designed to be used after [simulate_lifetimes()]. It takes
#' one row per simulation and per life, and returns one row per simulation with
#' simulated curtate and complete lifetimes for the selected multiple-life
#' status.
#'
#' The key output columns are `K_status` and `T_status`. They are intentionally
#' named so they can be used directly in downstream functions such as
#' [mc_insurance()] and [mc_annuity()] through their `k_col` and `tx_col`
#' arguments.
#'
#' @param data A data frame or tibble containing simulated multiple-life
#'   lifetimes, typically returned by [simulate_lifetimes()].
#' @param status Character string specifying the multiple-life status.
#'   Available options are `"joint_life"`, `"last_survivor"`,
#'   `"first_death"`, `"last_death"`, `"kth_death"`, and
#'   `"at_least_k_alive"`.
#' @param k Optional positive integer. Required when `status = "kth_death"` or
#'   `status = "at_least_k_alive"`.
#' @param sim_col Character string. Name of the simulation identifier column.
#'   Default is `"sim_id"`.
#' @param life_col Character string. Name of the life identifier column.
#'   Default is `"life_id"`.
#' @param k_col Character string. Name of the column containing simulated
#'   curtate future lifetimes. Default is `"Kx"`.
#' @param tx_col Character string. Name of the column containing simulated
#'   complete future lifetimes. Default is `"Tx"`.
#' @param k_status_col Character string. Name of the output column containing
#'   the curtate lifetime of the multiple-life status. Default is `"K_status"`.
#' @param t_status_col Character string. Name of the output column containing
#'   the complete lifetime of the multiple-life status. Default is `"T_status"`.
#' @param keep_lifetimes Logical. If `TRUE`, the output includes a list-column
#'   named `lifetimes` containing the simulated individual lifetimes used in
#'   each simulation. Default is `FALSE`.
#'
#' @details
#' Suppose that, for a given simulation, the complete future lifetimes of
#' \eqn{r} lives are
#'
#' \deqn{
#'   T_1, T_2, \ldots, T_r.
#' }
#'
#' Let
#'
#' \deqn{
#'   T_{(1)} \le T_{(2)} \le \cdots \le T_{(r)}
#' }
#'
#' denote the ordered lifetimes. The selected multiple-life status is converted
#' into an order statistic:
#'
#' * `"joint_life"`: fails at the first death, \eqn{T_{(1)}}.
#' * `"first_death"`: time until the first death, \eqn{T_{(1)}}.
#' * `"last_survivor"`: fails at the last death, \eqn{T_{(r)}}.
#' * `"last_death"`: time until the last death, \eqn{T_{(r)}}.
#' * `"kth_death"`: time until the `k`-th death, \eqn{T_{(k)}}.
#' * `"at_least_k_alive"`: status remains active while at least `k` lives are
#'   alive, and fails at the \eqn{(r-k+1)}-th death.
#'
#' The same order-statistic logic is applied to the curtate future lifetimes
#' \eqn{K_1, K_2, \ldots, K_r}, producing `K_status`.
#'
#' This construction allows multiple-life quantities to be evaluated using the
#' single-life Monte Carlo functions. For example, a joint-life annuity can be
#' computed by passing `k_col = "K_status"` to [mc_annuity()], and a first-death
#' insurance can be computed by passing `k_col = "K_status"` and
#' `tx_col = "T_status"` to [mc_insurance()].
#'
#' The function assumes the individual lifetimes have already been simulated.
#' If the input comes from [simulate_lifetimes()], the default model is
#' independence across lives.
#'
#' @return A tibble with one row per simulation and the following columns:
#'
#' * `sim_id`: simulation identifier, or the column named by `sim_col`.
#' * `status`: multiple-life status.
#' * `n_lives`: number of lives in the simulation.
#' * `k`: value of `k`, when applicable.
#' * `status_rank`: order statistic rank used to define the status.
#' * `K_status`: simulated curtate lifetime of the status, or another name
#'   supplied through `k_status_col`.
#' * `T_status`: simulated complete lifetime of the status, or another name
#'   supplied through `t_status_col`. If complete lifetimes are unavailable,
#'   this column contains `NA`.
#' * `event_life_id`: identifier of the life or lives associated with the
#'   status event time.
#' * `n_event_lives`: number of lives tied at the status event time.
#'
#' If `keep_lifetimes = TRUE`, a list-column named `lifetimes` is also returned.
#'
#' @seealso
#' [simulate_lifetimes()], [simulate_lifetime()], [mc_insurance()],
#' [mc_annuity()], [summary_mc()]
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). *Actuarial Mathematics*. Second Edition.
#' Society of Actuaries.
#'
#' @examples
#' life_table <- tibble::tibble(
#'   age = 60:90,
#'   qx = seq(0.01, 1, length.out = 31)
#' )
#'
#' # Simulated lifetimes for two independent lives
#' sim_two <- life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58),
#'     n_sim = 25,
#'     seed = 123
#'   )
#'
#' # Joint-life status: time until the first death
#' mc_multilife_status(
#'   sim_two,
#'   status = "joint_life"
#' )
#'
#' # Last-survivor status: time until the last death
#' mc_multilife_status(
#'   sim_two,
#'   status = "last_survivor"
#' )
#'
#' # Simulated lifetimes for three independent lives
#' sim_three <- life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58, 55),
#'     life_id = c("x", "y", "z"),
#'     n_sim = 25,
#'     seed = 123
#'   )
#'
#' # Time until the second death among three lives
#' mc_multilife_status(
#'   sim_three,
#'   status = "kth_death",
#'   k = 2
#' )
#'
#' # Status active while at least two lives are alive
#' mc_multilife_status(
#'   sim_three,
#'   status = "at_least_k_alive",
#'   k = 2
#' )
#'
#' @export
mc_multilife_status <- function(data,
                                status = c(
                                  "joint_life",
                                  "last_survivor",
                                  "first_death",
                                  "last_death",
                                  "kth_death",
                                  "at_least_k_alive"
                                ),
                                k = NULL,
                                sim_col = "sim_id",
                                life_col = "life_id",
                                k_col = "Kx",
                                tx_col = "Tx",
                                k_status_col = "K_status",
                                t_status_col = "T_status",
                                keep_lifetimes = FALSE) {
  status <- match.arg(status)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_column(data, sim_col, "sim_col")
  .mc_assert_column(data, life_col, "life_col")
  .mc_assert_numeric_column(data, k_col, "k_col")
  .mc_assert_character_scalar(k_status_col, "k_status_col")
  .mc_assert_character_scalar(t_status_col, "t_status_col")

  has_tx <- tx_col %in% names(data) &&
    is.numeric(data[[tx_col]]) &&
    !all(is.na(data[[tx_col]]))

  if (!is.logical(keep_lifetimes) || length(keep_lifetimes) != 1 ||
      is.na(keep_lifetimes)) {
    stop("`keep_lifetimes` must be a logical scalar.", call. = FALSE)
  }

  if (status %in% c("kth_death", "at_least_k_alive")) {
    .mc_assert_positive_integer(k, "k")
  }

  summarise_status_one <- function(df) {
    k_values <- df[[k_col]]
    life_ids <- as.character(df[[life_col]])

    if (anyNA(k_values)) {
      stop(
        "The selected curtate lifetime column contains missing values within a simulation.",
        call. = FALSE
      )
    }

    n_lives <- length(k_values)

    if (n_lives == 0) {
      stop("Each simulation must contain at least one life.", call. = FALSE)
    }

    status_rank <- .mc_multilife_status_rank(
      status = status,
      k = k,
      n_lives = n_lives
    )

    k_used <- if (status %in% c("kth_death", "at_least_k_alive")) {
      as.integer(k)
    } else {
      NA_integer_
    }

    k_status <- sort(k_values)[[status_rank]]

    if (has_tx) {
      t_values <- df[[tx_col]]

      if (anyNA(t_values)) {
        stop(
          "The selected complete lifetime column contains missing values within a simulation.",
          call. = FALSE
        )
      }

      t_status <- sort(t_values)[[status_rank]]

      tied <- abs(t_values - t_status) <= sqrt(.Machine$double.eps)
    } else {
      t_status <- NA_real_
      tied <- abs(k_values - k_status) <= sqrt(.Machine$double.eps)
    }

    event_life_id <- paste(life_ids[tied], collapse = ",")
    n_event_lives <- sum(tied)

    out <- tibble::tibble(
      status = status,
      n_lives = n_lives,
      k = k_used,
      status_rank = status_rank,
      "{k_status_col}" := k_status,
      "{t_status_col}" := t_status,
      event_life_id = event_life_id,
      n_event_lives = n_event_lives
    )

    if (keep_lifetimes) {
      lifetimes_tbl <- tibble::tibble(
        life_id = life_ids,
        Kx = k_values,
        Tx = if (has_tx) df[[tx_col]] else NA_real_
      )

      out <- dplyr::mutate(out, lifetimes = list(lifetimes_tbl))
    }

    out
  }

  data |>
    dplyr::ungroup() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sim_col))) |>
    dplyr::group_modify(~ summarise_status_one(.x)) |>
    dplyr::ungroup()
}


#' Internal helper: determine order statistic rank for a multiple-life status
#'
#' @param status Character string with the multiple-life status.
#' @param k Optional positive integer.
#' @param n_lives Positive integer. Number of lives.
#'
#' @return Positive integer giving the order statistic rank.
#'
#' @keywords internal
.mc_multilife_status_rank <- function(status, k, n_lives) {
  if (status %in% c("joint_life", "first_death")) {
    return(1L)
  }

  if (status %in% c("last_survivor", "last_death")) {
    return(as.integer(n_lives))
  }

  if (status == "kth_death") {
    if (k > n_lives) {
      stop(
        "`k` must be less than or equal to the number of lives in each simulation.",
        call. = FALSE
      )
    }

    return(as.integer(k))
  }

  if (status == "at_least_k_alive") {
    if (k > n_lives) {
      stop(
        "`k` must be less than or equal to the number of lives in each simulation.",
        call. = FALSE
      )
    }

    return(as.integer(n_lives - k + 1))
  }

  stop("Unknown multiple-life status.", call. = FALSE)
}

