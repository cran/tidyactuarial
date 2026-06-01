#' Compute multiple-life simulated status variables
#'
#' @description
#' Combines simulated future lifetimes from multiple lives into a single
#' multiple-life status. For a joint-life status, the status fails at the first
#' death. For a last-survivor status, the status fails at the last death.
#'
#' @param data A data frame or tibble produced by `simulate_lifetimes()`.
#' @param status Multiple-life status. One of `"joint"`, `"joint_life"`,
#'   `"last"`, `"last_survivor"`, or `"last_survivor_life"`.
#' @param col_sim Name of the simulation identifier column.
#' @param col_life Name of the life identifier column.
#' @param col_K Name of the curtate future lifetime column.
#' @param col_T Name of the complete future lifetime column.
#'
#' @return A tibble with one row per simulation and the columns `K_status`
#'   and `T_status`.
#'
#' @examples
#' lt <- tibble::tibble(
#'   x = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' lt |>
#'   simulate_lifetimes(
#'     x = c(60, 58),
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint")
#'
#' @export
mc_multilife_status <- function(data,
                                status = c("joint", "last_survivor"),
                                col_sim = "sim_id",
                                col_life = "life_id",
                                col_K = "Kx",
                                col_T = "Tx") {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  if (length(status) != 1L) {
    status <- status[1L]
  }

  status <- match.arg(
    status,
    choices = c(
      "joint",
      "joint_life",
      "last",
      "last_survivor",
      "last_survivor_life"
    )
  )

  required_cols <- c(col_sim, col_K, col_T)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0L) {
    stop(
      paste0(
        "`",
        paste(missing_cols, collapse = "`, `"),
        "` must identify a column in `data`."
      ),
      call. = FALSE
    )
  }

  if (!is.null(col_life) && !(col_life %in% names(data))) {
    col_life <- NULL
  }

  sim_values <- data[[col_sim]]
  K_values <- data[[col_K]]
  T_values <- data[[col_T]]

  if (any(is.na(sim_values))) {
    stop("`col_sim` must not contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(K_values))) {
    stop("`col_K` must contain only finite values.", call. = FALSE)
  }

  if (any(!is.finite(T_values))) {
    stop("`col_T` must contain only finite values.", call. = FALSE)
  }

  sim_levels <- unique(sim_values)

  compute_one_sim <- function(s) {
    idx <- which(sim_values == s)

    K_i <- K_values[idx]
    T_i <- T_values[idx]

    if (status %in% c("joint", "joint_life")) {
      K_status <- min(K_i)
      T_status <- min(T_i)
    } else {
      K_status <- max(K_i)
      T_status <- max(T_i)
    }

    tibble::tibble(
      sim_id = s,
      K_status = as.integer(K_status),
      T_status = as.numeric(T_status),
      n_lives = length(idx),
      status = status
    )
  }

  out <- lapply(sim_levels, compute_one_sim)
  out <- dplyr::bind_rows(out)

  # Keep a generic simulation alias for compatibility with functions that may
  # use either `sim_id` or `sim`.
  out$sim <- out$sim_id

  out
}
