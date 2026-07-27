#' Accumulation factor under conventional or time-varying interest
#'
#' Computes the accumulation factor from time \code{s} to time \code{t} using
#' exactly one of three interest specifications:
#' \itemize{
#'   \item a conventional interest rate \code{i};
#'   \item an accumulation function \code{a(t)};
#'   \item a time-varying force of interest \code{delta(t)}.
#' }
#'
#' The accumulation factor is
#' \deqn{A(s,t) = (1+i)^{t-s}}
#' for a constant annual effective rate,
#' \deqn{A(s,t) = \frac{a(t)}{a(s)}}
#' for an accumulation function, and
#' \deqn{A(s,t) = \exp\left(\int_s^t \delta(u)\,du\right)}
#' for a time-varying force of interest.
#'
#' @param s Numeric vector of starting times. Defaults to 0.
#' @param t Numeric vector of ending times.
#' @param i Optional numeric vector of conventional interest-rate values.
#' @param i_type Character vector indicating the conventional interest-rate
#'   type. Allowed values are \code{"effective"},
#'   \code{"nominal_interest"}, \code{"nominal_discount"}, and
#'   \code{"force"}. Used only when \code{i} is supplied.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Used only when \code{i} is supplied.
#' @param a Optional accumulation function of one numeric time argument. The
#'   function must return one positive finite numeric value for each supplied
#'   time.
#' @param delta Optional force-of-interest function of one numeric time
#'   argument. Numerical integration is performed separately over each
#'   interval \eqn{[s,t]}.
#' @param subdivisions Positive integer giving the maximum number of
#'   subintervals used by \code{stats::integrate()} when \code{delta} is
#'   supplied.
#' @param rel.tol Positive relative tolerance passed to
#'   \code{stats::integrate()} when \code{delta} is supplied.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric vector. If
#'   \code{TRUE}, returns a tibble with the inputs and intermediate quantities.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of accumulation factors.
#'
#' If \code{tidy = TRUE}, a tibble containing the interval, model, relevant
#' intermediate quantities, and accumulation factor.
#'
#' @details
#' Exactly one of \code{i}, \code{a}, or \code{delta} must be supplied.
#'
#' Numeric arguments must have length 1 or a common length. Scalars are
#' recycled over the remaining scenarios, making the function suitable for
#' use inside \code{dplyr::mutate()} pipelines.
#'
#' Times must satisfy \eqn{0 \le s \le t}. Missing times propagate as missing
#' accumulation factors.
#'
#' @seealso \code{\link{standardize_interest}},
#'   \code{\link{present_value}}, \code{\link{future_value}}
#'
#' @family interest
#' @family time-value
#'
#' @examples
#' # Constant annual effective interest
#' accumulation_factor(t = 5, i = 0.07)
#'
#' # One rate recycled over several times
#' accumulation_factor(
#'   t = c(1, 2, 5),
#'   i = 0.07
#' )
#'
#' # Accumulation function
#' a_fun <- function(t) {
#'   exp(0.03 * t + 0.002 * t^2)
#' }
#'
#' accumulation_factor(
#'   s = 2,
#'   t = 5,
#'   a = a_fun
#' )
#'
#' # Time-varying force of interest
#' delta_fun <- function(t) {
#'   0.03 + 0.004 * t
#' }
#'
#' accumulation_factor(
#'   s = 2,
#'   t = 5,
#'   delta = delta_fun
#' )
#'
#' # Pipe-friendly use with a small tibble
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   scenarios <- tibble::tibble(
#'     i = 0.07,
#'     t = c(1, 2, 3, 5)
#'   )
#'
#'   scenarios |>
#'     dplyr::mutate(
#'       accumulation = accumulation_factor(
#'         t = t,
#'         i = i
#'       )
#'     )
#' }
#'
#' @export
accumulation_factor <- function(
    s = 0,
    t,
    i = NULL,
    i_type = "effective",
    m = 1,
    a = NULL,
    delta = NULL,
    subdivisions = 100L,
    rel.tol = 1e-8,
    tidy = FALSE
) {
  if (missing(t)) {
    stop("`t` must be provided.", call. = FALSE)
  }

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (!is.numeric(s)) {
    stop("`s` must be a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(t)) {
    stop("`t` must be a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(subdivisions) ||
      length(subdivisions) != 1L ||
      is.na(subdivisions) ||
      !is.finite(subdivisions) ||
      subdivisions < 1 ||
      subdivisions != floor(subdivisions)) {
    stop("`subdivisions` must be a positive integer scalar.", call. = FALSE)
  }

  if (!is.numeric(rel.tol) ||
      length(rel.tol) != 1L ||
      is.na(rel.tol) ||
      !is.finite(rel.tol) ||
      rel.tol <= 0) {
    stop("`rel.tol` must be a positive finite numeric scalar.", call. = FALSE)
  }

  model_supplied <- c(
    interest_rate = !is.null(i),
    accumulation_function = !is.null(a),
    force_function = !is.null(delta)
  )

  if (sum(model_supplied) != 1L) {
    stop(
      "Supply exactly one of `i`, `a`, or `delta`.",
      call. = FALSE
    )
  }

  model <- names(model_supplied)[model_supplied]

  if (identical(model, "interest_rate")) {
    if (!is.numeric(i)) {
      stop("`i` must be a numeric vector.", call. = FALSE)
    }

    if (!is.character(i_type)) {
      stop("`i_type` must be a character vector.", call. = FALSE)
    }

    if (!is.numeric(m)) {
      stop("`m` must be a numeric vector.", call. = FALSE)
    }

    size <- max(
      length(s),
      length(t),
      length(i),
      length(i_type),
      length(m),
      1L
    )

    valid_size <- function(x) {
      length(x) %in% c(1L, size)
    }

    if (!valid_size(s) ||
        !valid_size(t) ||
        !valid_size(i) ||
        !valid_size(i_type) ||
        !valid_size(m)) {
      stop(
        "`s`, `t`, `i`, `i_type`, and `m` must have length 1 or a common length.",
        call. = FALSE
      )
    }

    s <- rep_len(s, size)
    t <- rep_len(t, size)
    i <- rep_len(i, size)
    i_type <- rep_len(i_type, size)
    m <- rep_len(m, size)
  } else {
    fn <- if (identical(model, "accumulation_function")) a else delta

    if (!is.function(fn)) {
      argument_name <- if (identical(model, "accumulation_function")) {
        "a"
      } else {
        "delta"
      }

      stop(
        sprintf("`%s` must be a function.", argument_name),
        call. = FALSE
      )
    }

    size <- max(length(s), length(t), 1L)

    valid_size <- function(x) {
      length(x) %in% c(1L, size)
    }

    if (!valid_size(s) || !valid_size(t)) {
      stop(
        "`s` and `t` must have length 1 or a common length.",
        call. = FALSE
      )
    }

    s <- rep_len(s, size)
    t <- rep_len(t, size)
    i <- rep(NA_real_, size)
    i_type <- rep(NA_character_, size)
    m <- rep(NA_real_, size)
  }

  bad_s <- !is.na(s) & (!is.finite(s) | s < 0)
  if (any(bad_s)) {
    stop(
      "`s` must contain only finite values greater than or equal to 0, or NA.",
      call. = FALSE
    )
  }

  bad_t <- !is.na(t) & (!is.finite(t) | t < 0)
  if (any(bad_t)) {
    stop(
      "`t` must contain only finite values greater than or equal to 0, or NA.",
      call. = FALSE
    )
  }

  bad_order <- !is.na(s) & !is.na(t) & t < s
  if (any(bad_order)) {
    stop("Each ending time `t` must be greater than or equal to `s`.", call. = FALSE)
  }

  elapsed_time <- t - s
  i_effective <- rep(NA_real_, size)
  a_s <- rep(NA_real_, size)
  a_t <- rep(NA_real_, size)
  integral_delta <- rep(NA_real_, size)
  accumulation_out <- rep(NA_real_, size)

  if (identical(model, "interest_rate")) {
    i_effective <- standardize_interest(
      i_type = i_type,
      i = i,
      m = m
    )

    ok <- !is.na(elapsed_time) & !is.na(i_effective)

    accumulation_out[ok] <-
      (1 + i_effective[ok])^elapsed_time[ok]
  }

  if (identical(model, "accumulation_function")) {
    a_s <- .evaluate_scalar_time_function(
      fn = a,
      time = s,
      argument = "a"
    )

    a_t <- .evaluate_scalar_time_function(
      fn = a,
      time = t,
      argument = "a"
    )

    bad_a <- (!is.na(a_s) & a_s <= 0) |
      (!is.na(a_t) & a_t <= 0)

    if (any(bad_a)) {
      stop(
        "`a` must return positive finite numeric values.",
        call. = FALSE
      )
    }

    ok <- !is.na(a_s) & !is.na(a_t)
    accumulation_out[ok] <- a_t[ok] / a_s[ok]
  }

  if (identical(model, "force_function")) {
    ok <- !is.na(s) & !is.na(t)

    integral_delta[ok] <- vapply(
      which(ok),
      function(index) {
        .integrate_force_interval(
          delta = delta,
          lower = s[[index]],
          upper = t[[index]],
          subdivisions = as.integer(subdivisions),
          rel.tol = rel.tol
        )
      },
      numeric(1)
    )

    accumulation_out[ok] <- exp(integral_delta[ok])
  }

  bad_result <- !is.na(accumulation_out) &
    (!is.finite(accumulation_out) | accumulation_out <= 0)

  if (any(bad_result)) {
    stop(
      "The computed accumulation factors must be positive and finite.",
      call. = FALSE
    )
  }

  if (!tidy) {
    return(accumulation_out)
  }

  tibble::tibble(
    s = s,
    t = t,
    elapsed_time = elapsed_time,
    model = rep(model, size),
    i_input = i,
    i_type = i_type,
    m = m,
    i_effective = i_effective,
    a_s = a_s,
    a_t = a_t,
    integral_delta = integral_delta,
    accumulation_factor = accumulation_out
  )
}

.evaluate_scalar_time_function <- function(
    fn,
    time,
    argument
) {
  out <- rep(NA_real_, length(time))
  ok <- !is.na(time)

  if (!any(ok)) {
    return(out)
  }

  out[ok] <- vapply(
    time[ok],
    function(value) {
      result <- fn(value)

      if (!is.numeric(result) ||
          length(result) != 1L ||
          is.na(result) ||
          !is.finite(result)) {
        stop(
          sprintf(
            "`%s` must return one finite numeric value for each supplied time.",
            argument
          ),
          call. = FALSE
        )
      }

      as.numeric(result)
    },
    numeric(1)
  )

  out
}

.integrate_force_interval <- function(
    delta,
    lower,
    upper,
    subdivisions,
    rel.tol
) {
  if (identical(lower, upper)) {
    return(0)
  }

  integrand <- function(time) {
    vapply(
      time,
      function(value) {
        result <- delta(value)

        if (!is.numeric(result) ||
            length(result) != 1L ||
            is.na(result) ||
            !is.finite(result)) {
          stop(
            "`delta` must return one finite numeric value for each supplied time.",
            call. = FALSE
          )
        }

        as.numeric(result)
      },
      numeric(1)
    )
  }

  integral <- tryCatch(
    stats::integrate(
      f = integrand,
      lower = lower,
      upper = upper,
      subdivisions = subdivisions,
      rel.tol = rel.tol,
      stop.on.error = TRUE
    ),
    error = function(error) {
      stop(
        "Could not integrate `delta` over [",
        lower,
        ", ",
        upper,
        "]: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  integral$value
}
