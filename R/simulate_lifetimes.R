#' Simulate future lifetimes from a life table
#'
#' @description
#' Simulates curtate and complete future lifetimes from a life table.
#'
#' This version is intentionally simple and defensive. It avoids using
#' `findInterval()` on a cumulative distribution that may contain missing values,
#' while preserving the column names expected by the Monte Carlo workflow:
#' `sim_id`, `life_id`, `Kx`, and `Tx`.
#'
#' @param data A data frame or tibble containing a life table.
#' @param x Numeric vector of initial ages.
#' @param n_sim Number of simulations per life.
#' @param frac Fractional age assumption. One of `"udd"`, `"constant"`,
#'   `"cfm"`, or `"balducci"`.
#' @param seed Optional random seed.
#'
#' @return A tibble with simulated curtate and complete future lifetimes.
#'
#' @details
#' The input life table must contain an age column named `x`, `age`, or `Age`,
#' and at least one mortality column among `qx`, `px`, or `lx`.
#'
#' If `qx` is not available, it is obtained from `px` as `1 - px`, or from
#' consecutive `lx` values as `1 - lx[x + 1] / lx[x]`.
#'
#' Invalid mortality values are handled defensively. Non-finite values, negative
#' values, and values greater than one are removed from the mortality basis.
#' Missing values are then treated as zero, and the final available age is forced
#' to have `qx = 1` so that the lifetime distribution is closed.
#'
#' @examples
#' lt <- tibble::tibble(
#'   x = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' simulate_lifetimes(
#'   data = lt,
#'   x = c(60, 58),
#'   n_sim = 25,
#'   frac = "udd",
#'   seed = 123
#' )
#'
#' @export
simulate_lifetimes <- function(data,
                               x,
                               n_sim = 1000,
                               frac = "udd",
                               seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  if (!is.numeric(x) || length(x) == 0L) {
    stop("`x` must be a non-empty numeric vector of ages.", call. = FALSE)
  }

  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) ||
      n_sim <= 0) {
    stop("`n_sim` must be a positive integer.", call. = FALSE)
  }

  n_sim <- as.integer(n_sim)

  frac <- match.arg(
    frac,
    choices = c("udd", "constant", "cfm", "balducci")
  )

  age_col <- intersect(c("x", "age", "Age"), names(data))[1]

  if (is.na(age_col)) {
    stop("`data` must contain an age column named `x`, `age`, or `Age`.",
         call. = FALSE)
  }

  ages <- as.numeric(data[[age_col]])

  if (any(!is.finite(ages))) {
    stop("The age column must contain only finite numeric values.",
         call. = FALSE)
  }

  ord <- order(ages)
  data <- data[ord, , drop = FALSE]
  ages <- ages[ord]

  if (any(duplicated(ages))) {
    stop("The age column must not contain duplicated ages.", call. = FALSE)
  }

  if ("qx" %in% names(data)) {
    qx <- as.numeric(data[["qx"]])
  } else if ("px" %in% names(data)) {
    qx <- 1 - as.numeric(data[["px"]])
  } else if ("lx" %in% names(data)) {
    lx <- as.numeric(data[["lx"]])
    qx <- rep(NA_real_, length(lx))

    if (length(lx) >= 2L) {
      denom <- lx[-length(lx)]
      numer <- lx[-1]

      valid <- is.finite(denom) & denom > 0 & is.finite(numer)
      qx[-length(lx)][valid] <- 1 - numer[valid] / denom[valid]
    }

    qx[length(qx)] <- 1
  } else {
    stop("`data` must contain one of the columns `qx`, `px`, or `lx`.",
         call. = FALSE)
  }

  qx <- as.numeric(qx)

  if (length(qx) != length(ages)) {
    stop("The mortality column must have the same length as the age column.",
         call. = FALSE)
  }

  qx[!is.finite(qx)] <- NA_real_
  qx[qx < 0] <- NA_real_
  qx[qx > 1] <- NA_real_

  if (all(is.na(qx))) {
    stop("A valid mortality column could not be constructed.", call. = FALSE)
  }

  qx[is.na(qx)] <- 0
  qx[length(qx)] <- 1

  simulate_one_life <- function(age_start, life_id) {
    if (!is.finite(age_start)) {
      stop("All initial ages in `x` must be finite.", call. = FALSE)
    }

    start_pos <- which(ages >= age_start)[1]

    if (is.na(start_pos)) {
      stop("An initial age in `x` is outside the life table range.",
           call. = FALSE)
    }

    q <- qx[start_pos:length(qx)]

    q[!is.finite(q)] <- 0
    q <- pmin(pmax(q, 0), 1)
    q[length(q)] <- 1

    surv_before <- c(1, cumprod(1 - q[-length(q)]))
    prob_death <- surv_before * q

    prob_death[!is.finite(prob_death)] <- 0
    prob_death <- pmax(prob_death, 0)

    total_prob <- sum(prob_death)

    if (!is.finite(total_prob) || total_prob <= 0) {
      stop("A valid lifetime distribution could not be constructed.",
           call. = FALSE)
    }

    prob_death <- prob_death / total_prob

    Kx <- sample.int(
      n = length(prob_death),
      size = n_sim,
      replace = TRUE,
      prob = prob_death
    ) - 1L

    q_death <- q[Kx + 1L]
    u <- stats::runif(n_sim)

    frac_time <- u

    if (frac %in% c("constant", "cfm")) {
      valid <- q_death > 0 & q_death < 1
      p_death <- 1 - q_death

      frac_time[valid] <- log(1 - u[valid] * q_death[valid]) /
        log(p_death[valid])

      frac_time[!is.finite(frac_time)] <- u[!is.finite(frac_time)]
    }

    if (frac == "balducci") {
      valid <- q_death > 0 & q_death < 1

      frac_time[valid] <- u[valid] * (1 - q_death[valid]) /
        (1 - u[valid] * q_death[valid])

      frac_time[!is.finite(frac_time)] <- u[!is.finite(frac_time)]
    }

    frac_time <- pmin(pmax(frac_time, 0), 1)

    Tx <- Kx + frac_time

    tibble::tibble(
      sim_id = seq_len(n_sim),
      life_id = life_id,
      x = age_start,
      age = age_start,
      method = "inverse",
      frac = frac,
      fractional = frac,
      Kx = as.integer(Kx),
      Tx = as.numeric(Tx),
      K = as.integer(Kx),
      T = as.numeric(Tx),
      sim = seq_len(n_sim),
      life = life_id
    )
  }

  out <- lapply(seq_along(x), function(j) {
    simulate_one_life(age_start = x[j], life_id = j)
  })

  dplyr::bind_rows(out)
}

