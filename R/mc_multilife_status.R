#' Compute multiple-life simulated status variables
#'
#' Combines simulated future lifetimes for several lives into one simulated
#' multiple-life status. A joint-life status terminates at the first death,
#' whereas a last-survivor status terminates at the last death.
#'
#' @param data A data frame or tibble, typically returned by
#'   \code{\link{simulate_lifetimes}}.
#' @param status Multiple-life status. Canonical values are \code{"joint"} and
#'   \code{"last_survivor"}. The aliases \code{"joint_life"},
#'   \code{"last"}, and \code{"last_survivor_life"} are accepted.
#' @param col_sim Character scalar naming the simulation identifier column.
#' @param col_life Character scalar naming the life identifier column. The
#'   default is \code{"life_id"}. Supply \code{NULL} only when no life
#'   identifier exists; in that case the function can validate equal row counts
#'   across simulations but cannot detect duplicated or substituted lives.
#' @param col_K Character scalar naming the curtate future lifetime column.
#' @param col_T Character scalar naming the complete future lifetime column.
#' @param tol Nonnegative numeric tolerance used when checking that each
#'   complete future lifetime is compatible with its curtate lifetime.
#'
#' @details
#' For simulated complete future lifetimes \eqn{T_1,\ldots,T_r}, the
#' multiple-life complete future lifetime is
#' \deqn{
#' T_{\mathrm{joint}}=\min(T_1,\ldots,T_r)
#' }
#' for the joint-life status and
#' \deqn{
#' T_{\mathrm{last}}=\max(T_1,\ldots,T_r)
#' }
#' for the last-survivor status.
#'
#' The corresponding curtate future lifetime is obtained with the same minimum
#' or maximum operation on \eqn{K_1,\ldots,K_r}. Because the function validates
#' \eqn{K_j \le T_j < K_j+1}, the returned values satisfy the same
#' curtate-complete relationship.
#'
#' When \code{col_life} is available, every simulation must contain exactly one
#' row for every life and the same set of lives must appear in every
#' simulation. This prevents a missing or duplicated life from being silently
#' treated as a different multiple-life contract.
#'
#' @return A tibble with one row per simulation and:
#' \describe{
#'   \item{sim_id}{Standard simulation identifier.}
#'   \item{sim}{Compatibility alias for the simulation identifier.}
#'   \item{K_status}{Curtate future lifetime of the status.}
#'   \item{T_status}{Complete future lifetime of the status.}
#'   \item{n_lives}{Number of lives represented in each simulation.}
#'   \item{status}{Canonical status: \code{"joint"} or
#'   \code{"last_survivor"}.}
#' }
#'
#' If \code{col_sim} is neither \code{"sim_id"} nor \code{"sim"}, that
#' original identifier column is also retained.
#'
#' @seealso \code{\link{simulate_lifetimes}}, \code{\link{mc_insurance}},
#'   \code{\link{mc_annuity}}, \code{\link{mc_reserve}}
#'
#' @family monte-carlo
#'
#' @examples
#' simulated <- tibble::tibble(
#'   sim_id = rep(1:2, each = 2),
#'   life_id = rep(c("x", "y"), times = 2),
#'   Kx = c(2, 5, 4, 1),
#'   Tx = c(2.4, 5.2, 4.5, 1.3)
#' )
#'
#' simulated |>
#'   mc_multilife_status(status = "joint")
#'
#' simulated |>
#'   mc_multilife_status(status = "last_survivor")
#'
#' @export
mc_multilife_status <- function(
    data,
    status = c("joint", "last_survivor"),
    col_sim = "sim_id",
    col_life = "life_id",
    col_K = "Kx",
    col_T = "Tx",
    tol = 1e-10
) {
  status_missing <- missing(status)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  if (nrow(data) == 0L) {
    stop("`data` must contain at least one simulated contract.", call. = FALSE)
  }

  assert_column_name <- function(value, name, allow_null = FALSE) {
    if (allow_null && is.null(value)) {
      return(invisible(NULL))
    }

    if (!is.character(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !nzchar(value)) {
      stop(
        "`", name, "` must be ",
        if (allow_null) "`NULL` or " else "",
        "a single nonmissing column name.",
        call. = FALSE
      )
    }

    invisible(NULL)
  }

  assert_column_name(col_sim, "col_sim")
  assert_column_name(col_life, "col_life", allow_null = TRUE)
  assert_column_name(col_K, "col_K")
  assert_column_name(col_T, "col_T")

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  if (status_missing) {
    status <- "joint"
  } else {
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
  }

  status <- if (status %in% c("joint", "joint_life")) {
    "joint"
  } else {
    "last_survivor"
  }

  required_cols <- c(col_sim, col_K, col_T)

  if (!is.null(col_life)) {
    required_cols <- c(required_cols, col_life)
  }

  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(sprintf("`%s`", missing_cols), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  sim_values <- data[[col_sim]]
  K_values <- data[[col_K]]
  T_values <- data[[col_T]]

  if (is.list(sim_values) && !is.data.frame(sim_values)) {
    stop("`col_sim` must identify an atomic simulation identifier.", call. = FALSE)
  }

  if (anyNA(sim_values)) {
    stop("`col_sim` must not contain missing values.", call. = FALSE)
  }

  if (!is.numeric(K_values)) {
    stop("`col_K` must identify a numeric column.", call. = FALSE)
  }

  if (!is.numeric(T_values)) {
    stop("`col_T` must identify a numeric column.", call. = FALSE)
  }

  if (anyNA(K_values) || any(!is.finite(K_values))) {
    stop("`col_K` must contain only finite nonmissing values.", call. = FALSE)
  }

  if (anyNA(T_values) || any(!is.finite(T_values))) {
    stop("`col_T` must contain only finite nonmissing values.", call. = FALSE)
  }

  if (any(K_values < 0)) {
    stop("`col_K` must contain nonnegative curtate lifetimes.", call. = FALSE)
  }

  if (any(T_values < 0)) {
    stop("`col_T` must contain nonnegative complete lifetimes.", call. = FALSE)
  }

  if (any(abs(K_values - round(K_values)) > tol)) {
    stop("`col_K` must contain whole-number curtate lifetimes.", call. = FALSE)
  }

  K_values <- round(K_values)

  incompatible_lifetime <- T_values < K_values - tol |
    T_values >= K_values + 1 + tol

  if (any(incompatible_lifetime)) {
    stop(
      "Each complete lifetime must satisfy `K <= T < K + 1` within ",
      "the specified tolerance.",
      call. = FALSE
    )
  }

  sim_levels <- unique(sim_values)
  group_id <- match(sim_values, sim_levels)
  rows_by_sim <- split(seq_len(nrow(data)), group_id)

  n_lives_by_sim <- vapply(rows_by_sim, length, integer(1))

  if (any(n_lives_by_sim < 2L)) {
    stop(
      "Every simulation must contain at least two lives.",
      call. = FALSE
    )
  }

  if (!is.null(col_life)) {
    life_values <- data[[col_life]]

    if (is.list(life_values) && !is.data.frame(life_values)) {
      stop("`col_life` must identify an atomic life identifier.", call. = FALSE)
    }

    if (anyNA(life_values)) {
      stop("`col_life` must not contain missing values.", call. = FALSE)
    }

    reference_lives <- NULL

    for (rows in rows_by_sim) {
      lives_i <- as.character(life_values[rows])

      if (anyDuplicated(lives_i)) {
        stop(
          "Each life identifier must appear exactly once within every ",
          "simulation.",
          call. = FALSE
        )
      }

      lives_i <- sort(lives_i)

      if (is.null(reference_lives)) {
        reference_lives <- lives_i
      } else if (!identical(lives_i, reference_lives)) {
        stop(
          "Every simulation must contain the same set of life identifiers.",
          call. = FALSE
        )
      }
    }
  } else if (length(unique(n_lives_by_sim)) != 1L) {
    stop(
      "Without `col_life`, every simulation must contain the same number ",
      "of life rows.",
      call. = FALSE
    )
  }

  K_status <- vapply(
    rows_by_sim,
    function(rows) {
      if (identical(status, "joint")) {
        min(K_values[rows])
      } else {
        max(K_values[rows])
      }
    },
    numeric(1)
  )

  T_status <- vapply(
    rows_by_sim,
    function(rows) {
      if (identical(status, "joint")) {
        min(T_values[rows])
      } else {
        max(T_values[rows])
      }
    },
    numeric(1)
  )

  out <- tibble::tibble(
    "{col_sim}" := sim_levels,
    K_status = as.integer(K_status),
    T_status = as.numeric(T_status),
    n_lives = as.integer(n_lives_by_sim),
    status = status
  )

  if (!"sim_id" %in% names(out)) {
    out <- out |>
      dplyr::mutate(
        sim_id = .data[[col_sim]],
        .before = 1
      )
  }

  if (!"sim" %in% names(out)) {
    out <- out |>
      dplyr::mutate(
        sim = .data[[col_sim]],
        .after = "sim_id"
      )
  }

  out
}
