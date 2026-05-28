#' Simulate multiple independent future lifetimes
#'
#' Simulates curtate and complete future lifetimes for several lives under an
#' independence assumption.
#'
#' This function extends [simulate_lifetime()] to multiple lives. It is designed
#' as the simulation engine for Monte Carlo multiple-life actuarial calculations.
#' Each life is simulated independently, and the output is returned in tidy long
#' format with one row per simulation and per life.
#'
#' @param data A life table data frame or a list of life table data frames.
#'   If a single data frame is supplied, the same life table is used for all
#'   lives. If a list is supplied, it must have length 1 or the same length as
#'   `ages`.
#' @param ages Numeric vector. Initial ages of the lives to simulate.
#' @param n_sim Positive integer. Number of Monte Carlo simulations per life.
#'   Default is `10000`.
#' @param life_id Optional character vector identifying each life. If `NULL`,
#'   IDs are created automatically as `"life_1"`, `"life_2"`, and so on. If
#'   `data` is a named list and `life_id = NULL`, the list names are used when
#'   possible.
#' @param age_col Character string. Name of the age column in the life table or
#'   life tables. Default is `"age"`.
#' @param qx_col Character string. Name of the one-year death probability
#'   column in the life table or life tables. Default is `"qx"`.
#' @param method Character string specifying the simulation method for each
#'   curtate future lifetime. Available options are `"inverse"`,
#'   `"multinomial"`, and `"antithetic"`. Default is `"inverse"`.
#' @param fractional Character string specifying how the fractional part of the
#'   complete future lifetime is generated within the year of death. Available
#'   options are `"udd"`, `"constant_force"`, and `"none"`. Default is `"udd"`.
#' @param seed Optional integer seed for reproducibility. Default is `NULL`.
#'
#' @details
#' For lives aged \eqn{x_1, x_2, \ldots, x_r}, the function simulates
#'
#' \deqn{
#'   K_{x_1}, K_{x_2}, \ldots, K_{x_r}
#' }
#'
#' independently. The same applies to the complete future lifetimes
#' \eqn{T_{x_1}, T_{x_2}, \ldots, T_{x_r}} when `fractional` is not `"none"`.
#'
#' The independence assumption means that no common shock, copula, frailty, or
#' other dependence structure is imposed. This is a natural first model for
#' multiple-life Monte Carlo calculations and allows direct construction of
#' joint-life, last-survivor, first-death, last-death, and k-th death quantities
#' through downstream functions.
#'
#' The output is intentionally long:
#'
#' ```
#' sim_id | life_id | life_index | age | Kx | Tx
#' ```
#'
#' This format works naturally with [dplyr::group_by()],
#' [dplyr::summarise()], [tidyr::pivot_wider()], and downstream functions such
#' as `mc_multilife_status()`.
#'
#' If different mortality tables are needed for different lives, pass `data` as
#' a list of life tables. For example, one table may be used for a male life and
#' another for a female life.
#'
#' @return A tibble with one row per simulation and per life. It contains:
#'
#' * `sim_id`: simulation identifier.
#' * `life_id`: life identifier.
#' * `life_index`: position of the life in `ages`.
#' * `age`: initial age of the life.
#' * `method`: simulation method used.
#' * `fractional`: fractional age assumption.
#' * `Kx`: simulated curtate future lifetime.
#' * `Tx`: simulated complete future lifetime. If `fractional = "none"`,
#'   this column contains `NA`.
#'
#' @seealso
#' [simulate_lifetime()], [summary_mc()]
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). *Actuarial Mathematics*. Second Edition.
#' Society of Actuaries.
#'
#' @examples
#' life_table <- tibble::tibble(
#'   age = 40:110,
#'   qx = seq(0.002, 1, length.out = 71)
#' )
#'
#' # Two independent lives using the same life table
#' life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58),
#'     n_sim = 1000,
#'     seed = 123
#'   )
#'
#' # Three independent lives
#' life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58, 55),
#'     life_id = c("x", "y", "z"),
#'     n_sim = 1000,
#'     seed = 123
#'   )
#'
#' # Different life tables for different lives
#' male_table <- tibble::tibble(
#'   age = 40:110,
#'   qx = seq(0.003, 1, length.out = 71)
#' )
#'
#' female_table <- tibble::tibble(
#'   age = 40:110,
#'   qx = seq(0.002, 1, length.out = 71)
#' )
#'
#' list(male = male_table, female = female_table) |>
#'   simulate_lifetimes(
#'     ages = c(60, 58),
#'     n_sim = 1000,
#'     seed = 123
#'   )
#'
#' @export
simulate_lifetimes <- function(data,
                               ages,
                               n_sim = 10000,
                               life_id = NULL,
                               age_col = "age",
                               qx_col = "qx",
                               method = c("inverse", "multinomial", "antithetic"),
                               fractional = c("udd", "constant_force", "none"),
                               seed = NULL) {
  method <- match.arg(method)
  fractional <- match.arg(fractional)

  if (!is.numeric(ages) || length(ages) == 0 ||
      anyNA(ages) || any(!is.finite(ages))) {
    stop("`ages` must be a numeric vector without missing values.", call. = FALSE)
  }

  if (any(ages < 0)) {
    stop("All values in `ages` must be non-negative.", call. = FALSE)
  }

  .mc_assert_positive_integer(n_sim, "n_sim")
  .mc_assert_character_scalar(age_col, "age_col")
  .mc_assert_character_scalar(qx_col, "qx_col")

  n_lives <- length(ages)

  if (is.data.frame(data)) {
    life_tables <- rep(list(data), n_lives)
  } else if (is.list(data)) {
    if (!length(data) %in% c(1, n_lives)) {
      stop(
        "`data` must be a data frame or a list of length 1 or length equal ",
        "to `length(ages)`.",
        call. = FALSE
      )
    }

    if (!all(vapply(data, is.data.frame, logical(1)))) {
      stop("Every element of `data` must be a data frame or tibble.",
           call. = FALSE)
    }

    life_tables <- if (length(data) == 1) {
      rep(data, n_lives)
    } else {
      data
    }
  } else {
    stop(
      "`data` must be a life table data frame or a list of life table data frames.",
      call. = FALSE
    )
  }

  if (is.null(life_id)) {
    data_names <- names(data)

    if (is.list(data) && length(data) == n_lives &&
        !is.null(data_names) && all(nzchar(data_names))) {
      life_id <- data_names
    } else {
      life_id <- paste0("life_", seq_len(n_lives))
    }
  }

  if (!is.character(life_id) || length(life_id) != n_lives || anyNA(life_id)) {
    stop(
      "`life_id` must be NULL or a character vector with length equal to ",
      "`length(ages)`.",
      call. = FALSE
    )
  }

  if (anyDuplicated(life_id)) {
    stop("`life_id` values must be unique.", call. = FALSE)
  }

  if (!is.null(seed)) {
    .mc_assert_positive_integer(seed, "seed")
    set.seed(seed)
  }

  simulated_list <- lapply(
    seq_len(n_lives),
    function(i) {
      simulate_lifetime(
        data = life_tables[[i]],
        age = ages[[i]],
        n_sim = n_sim,
        age_col = age_col,
        qx_col = qx_col,
        method = method,
        fractional = fractional,
        seed = NULL,
        include_distribution = FALSE
      ) |>
        dplyr::mutate(
          life_id = life_id[[i]],
          life_index = i,
          .before = "age"
        ) |>
        dplyr::select(
          .data$sim_id,
          .data$life_id,
          .data$life_index,
          dplyr::everything()
        )
    }
  )

  dplyr::bind_rows(simulated_list)
}
