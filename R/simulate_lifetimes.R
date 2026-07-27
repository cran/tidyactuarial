#' Simulate future lifetimes for multiple lives
#'
#' Simulates independent curtate and complete future lifetimes for several
#' initial ages from one life table. The function uses
#' \code{\link{simulate_lifetime}} as its single-life simulation engine so that
#' age validation, truncation treatment, fractional-age assumptions, and
#' random-number conventions remain consistent across the Monte Carlo module.
#'
#' @param data A data frame or tibble containing a life table.
#' @param x Nonempty numeric vector of nonnegative integer initial ages. Each
#'   position represents one life; repeated ages are allowed.
#' @param n_sim Positive integer number of simulations per life.
#' @param frac Fractional-age assumption. Canonical values are \code{"udd"},
#'   \code{"cml"}, and \code{"balducci"}. The historical values
#'   \code{"constant"}, \code{"cfm"}, and \code{"constant_force"} are accepted
#'   as aliases for \code{"cml"}.
#' @param method Simulation method passed to \code{\link{simulate_lifetime}}:
#'   \code{"inverse"}, \code{"multinomial"}, or \code{"antithetic"}.
#' @param seed Optional nonnegative integer seed. When supplied, one random
#'   stream is initialized for the complete multiple-life simulation and the
#'   caller's previous random-number state is restored afterward.
#' @param truncation Treatment of a mortality table that does not exhaust the
#'   lifetime distribution: \code{"conditional"} retains the historical
#'   conditional simulation with one warning, while \code{"error"} stops.
#' @param tol Nonnegative numeric tolerance used for age and mortality checks.
#'
#' @details
#' The life table must contain an age column named \code{x}, \code{age}, or
#' \code{Age}, and at least one mortality basis among \code{qx}, \code{px}, or
#' \code{lx}. The precedence is \code{qx}, then \code{px}, then \code{lx}.
#'
#' Unlike the previous implementation, invalid mortality values are not silently
#' replaced by zero and the last death probability is not forcibly overwritten.
#' Such repairs can materially change the simulated lifetime distribution.
#'
#' If \code{lx} is supplied, one-year death probabilities are derived as
#' \deqn{
#' q_y = 1 - \frac{l_{y+1}}{l_y}
#' }
#' for ages with positive exposure. A terminal zero in \code{lx} therefore
#' generates the appropriate terminal \code{qx = 1}. If the final available
#' \code{lx} remains positive, the resulting distribution is recognized as
#' truncated and handled according to \code{truncation}.
#'
#' Conditional on the supplied life table, the simulated lives are independent.
#' A single seed is set once for the entire call; the seed is not reset for each
#' life, which avoids introducing artificial perfect dependence between lives
#' of the same age.
#'
#' @return A tibble with one row per simulation and life. Standard columns
#' include \code{sim_id}, \code{simulation_id}, \code{life_id}, \code{x},
#' \code{Kx}, \code{curtate_lifetime}, \code{Tx},
#' \code{complete_lifetime}, \code{death_age}, \code{method}, and \code{frac}.
#' Compatibility aliases \code{sim}, \code{life}, \code{K}, and \code{T} are
#' retained.
#'
#' @seealso \code{\link{simulate_lifetime}},
#'   \code{\link{mc_multilife_status}}, \code{\link{mc_insurance}},
#'   \code{\link{mc_annuity}}, \code{\link{mc_reserve}}
#'
#' @family monte-carlo
#'
#' @examples
#' lt <- tibble::tibble(
#'   x = 40:43,
#'   qx = c(0.10, 0.20, 0.30, 1)
#' )
#'
#' simulate_lifetimes(
#'   data = lt,
#'   x = c(40, 41),
#'   n_sim = 25,
#'   frac = "udd",
#'   seed = 123
#' )
#'
#' @export
simulate_lifetimes <- function(
    data,
    x,
    n_sim = 1000L,
    frac = c(
      "udd",
      "cml",
      "balducci",
      "constant",
      "cfm",
      "constant_force"
    ),
    method = c(
      "inverse",
      "multinomial",
      "antithetic"
    ),
    seed = NULL,
    truncation = c("conditional", "error"),
    tol = 1e-10
) {
  frac <- match.arg(frac)
  method <- match.arg(method)
  truncation <- match.arg(truncation)

  if (frac %in% c(
    "constant",
    "cfm",
    "constant_force"
  )) {
    frac <- "cml"
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  if (nrow(data) == 0L) {
    stop("`data` must contain at least one life-table row.", call. = FALSE)
  }

  if (!is.numeric(x) ||
      length(x) == 0L ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x < 0) ||
      any(abs(x - round(x)) > tol)) {
    stop(
      "`x` must be a nonempty vector of nonnegative integer ages.",
      call. = FALSE
    )
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

  age_candidates <- intersect(
    c("x", "age", "Age"),
    names(data)
  )

  if (length(age_candidates) == 0L) {
    stop(
      "`data` must contain an age column named `x`, `age`, or `Age`.",
      call. = FALSE
    )
  }

  age_col <- age_candidates[[1L]]
  ages <- data[[age_col]]

  if (!is.numeric(ages) ||
      anyNA(ages) ||
      any(!is.finite(ages)) ||
      any(ages < 0) ||
      any(abs(ages - round(ages)) > tol)) {
    stop(
      "The age column must contain finite nonnegative integer ages.",
      call. = FALSE
    )
  }

  ages <- as.integer(round(ages))

  if (anyDuplicated(ages)) {
    stop(
      "The age column must not contain duplicated ages.",
      call. = FALSE
    )
  }

  ord <- order(ages)
  ages <- ages[ord]
  data <- data[ord, , drop = FALSE]

  mortality_basis <- NULL

  if ("qx" %in% names(data)) {
    qx <- data[["qx"]]

    if (!is.numeric(qx) ||
        anyNA(qx) ||
        any(!is.finite(qx)) ||
        any(qx < -tol | qx > 1 + tol)) {
      stop(
        "Column `qx` must contain finite probabilities between 0 and 1.",
        call. = FALSE
      )
    }

    qx <- pmin(pmax(qx, 0), 1)

    mortality_basis <- tibble::tibble(
      x = ages,
      qx = qx
    )
  } else if ("px" %in% names(data)) {
    px <- data[["px"]]

    if (!is.numeric(px) ||
        anyNA(px) ||
        any(!is.finite(px)) ||
        any(px < -tol | px > 1 + tol)) {
      stop(
        "Column `px` must contain finite probabilities between 0 and 1.",
        call. = FALSE
      )
    }

    px <- pmin(pmax(px, 0), 1)

    mortality_basis <- tibble::tibble(
      x = ages,
      qx = 1 - px
    )
  } else if ("lx" %in% names(data)) {
    lx <- data[["lx"]]

    if (!is.numeric(lx) ||
        anyNA(lx) ||
        any(!is.finite(lx)) ||
        any(lx < 0)) {
      stop(
        "Column `lx` must contain finite nonnegative values.",
        call. = FALSE
      )
    }

    if (length(lx) < 2L) {
      stop(
        "At least two consecutive `lx` values are required.",
        call. = FALSE
      )
    }

    if (any(diff(lx) > tol)) {
      stop(
        "Column `lx` must be nonincreasing with age.",
        call. = FALSE
      )
    }

    positive_rows <- which(lx > tol)

    if (length(positive_rows) == 0L) {
      stop(
        "Column `lx` must contain positive exposure at some age.",
        call. = FALSE
      )
    }

    last_positive <- max(positive_rows)

    if (last_positive == length(lx)) {
      derivation_rows <- seq_len(length(lx) - 1L)
    } else {
      derivation_rows <- seq_len(last_positive)
    }

    denominator <- lx[derivation_rows]
    numerator <- lx[derivation_rows + 1L]

    qx <- 1 - numerator / denominator

    if (any(!is.finite(qx)) ||
        any(qx < -tol | qx > 1 + tol)) {
      stop(
        "Column `lx` implies invalid one-year death probabilities.",
        call. = FALSE
      )
    }

    qx <- pmin(pmax(qx, 0), 1)

    mortality_basis <- tibble::tibble(
      x = ages[derivation_rows],
      qx = qx
    )
  } else {
    stop(
      "`data` must contain one of the mortality columns `qx`, `px`, or `lx`.",
      call. = FALSE
    )
  }

  missing_start_ages <- setdiff(
    unique(x),
    mortality_basis$x
  )

  if (length(missing_start_ages) > 0L) {
    stop(
      "Every initial age in `x` must appear explicitly in the usable ",
      "mortality basis. Missing age(s): ",
      paste(missing_start_ages, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

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

  conditioned_any <- FALSE

  simulate_one_life <- function(age_start, life_id) {
    result <- withCallingHandlers(
      simulate_lifetime(
        lt = mortality_basis,
        x = age_start,
        n_sim = n_sim,
        x_col = "x",
        qx_col = "qx",
        method = method,
        frac = frac,
        seed = NULL,
        include_distribution = FALSE,
        truncation = truncation,
        tol = tol
      ),
      warning = function(w) {
        if (grepl(
          "truncated",
          conditionMessage(w),
          ignore.case = TRUE
        )) {
          invokeRestart("muffleWarning")
        }
      }
    )

    conditioned_any <<- conditioned_any ||
      any(result$distribution_conditioned)

    result |>
      dplyr::mutate(
        life_id = life_id,
        life = life_id,
        sim = .data$sim_id,
        K = .data$Kx,
        T = .data$Tx
      ) |>
      dplyr::select(
        dplyr::all_of(
          c(
            "sim_id",
            "simulation_id",
            "life_id",
            "life",
            "x",
            "age",
            "method",
            "frac",
            "fractional",
            "Kx",
            "curtate_lifetime",
            "Tx",
            "complete_lifetime",
            "death_age",
            "qx_at_death_year",
            "distribution_conditioned",
            "sim",
            "K",
            "T"
          )
        )
      )
  }

  out <- lapply(
    seq_along(x),
    function(j) {
      simulate_one_life(
        age_start = x[[j]],
        life_id = j
      )
    }
  ) |>
    dplyr::bind_rows()

  if (conditioned_any &&
      identical(truncation, "conditional")) {
    warning(
      "The life table is truncated for at least one initial age. ",
      "Simulation is conditional on death within the available ages.",
      call. = FALSE
    )
  }

  out
}
