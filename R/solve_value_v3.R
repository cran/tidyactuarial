#' Solve a scalar equation of value
#'
#' Solves for one unknown scalar argument of an existing actuarial or financial
#' function. The remaining arguments are supplied through `args`, and the
#' unknown is chosen with `solve_for`.
#'
#' The equation solved is
#' \deqn{g(x) = f(\ldots, x, \ldots) - \mathrm{target} = 0,}
#' where `f` is the function supplied in `fn` and `x` is the argument named in
#' `solve_for`.
#'
#' This function is intended to reuse the stable calculation functions already
#' available in tidyactuarial. It does not reproduce their actuarial formulas.
#'
#' @param fn A function, or the name of a function, whose scalar result is to be
#'   matched to `target`.
#' @param solve_for Character scalar. Name of the scalar argument to solve for.
#' @param target Finite numeric scalar. Desired value of the selected result.
#' @param args Named list containing all known arguments passed to `fn`. Do not
#'   supply a non-`NULL` value for the argument named in `solve_for`.
#' @param result Optional result extractor. Use `NULL` when `fn` returns a
#'   numeric scalar; a character scalar naming a column or list element; or a
#'   function that extracts a numeric scalar from the object returned by `fn`.
#' @param interval Optional finite numeric vector of length 2. It defines the
#'   admissible search interval. When omitted, a conservative interval is
#'   inferred from `solve_for`, `target`, and the numeric values in `args`.
#' @param start Optional finite numeric scalar used as the initial value for
#'   Newton--Raphson. It is also used as a fallback starting point when
#'   `method = "auto"` cannot bracket a root.
#' @param derivative Optional function of one argument returning the derivative
#'   of the equation residual with respect to the unknown. When omitted under
#'   Newton--Raphson, a central finite-difference derivative is used.
#' @param method Character scalar. One of `"auto"`, `"uniroot"`, or
#'   `"newton"`. The default first searches for bracketed roots and uses
#'   Newton--Raphson only when no bracket is found and `start` is available.
#' @param multiple Character scalar controlling multiple roots found inside the
#'   interval: `"error"`, `"all"`, or `"first"`.
#' @param tol Positive numeric scalar. Convergence and residual tolerance.
#' @param maxiter Positive integer. Maximum number of iterations.
#' @param scan_points Positive integer. Number of points used to detect sign
#'   changes and possible multiple roots.
#' @param max_expand Nonnegative integer. Maximum number of automatic interval
#'   expansions when an interval was inferred and a bracket is not found
#'   initially. An interval supplied explicitly by the user is never expanded.
#' @param output Character scalar. One of `"value"`, `"summary"`, or
#'   `"audit"`.
#'
#' @return
#' Depending on `output`:
#' \itemize{
#'   \item `"value"`: numeric scalar or vector containing the solution(s).
#'   \item `"summary"`: compact tibble with the unknown, solution, target,
#'     achieved value, residual, and method.
#'   \item `"audit"`: detailed tibble with convergence diagnostics and search
#'     information.
#' }
#'
#' @details
#' `solve_value()` is restricted to one scalar unknown. It is appropriate for
#' quantities such as an effective interest rate, a payment, a principal, a
#' present value, or another continuous scalar parameter.
#'
#' Bracketed root finding is preferred because it is generally more stable than
#' an unconstrained Newton iteration. The interval is scanned before solving;
#' therefore, possible multiple roots can be detected. This is important for
#' non-monotone equations of value.
#'
#' Discrete unknowns, such as an integer number of payments, should not be
#' solved by pretending they are continuous. A dedicated discrete solver or a
#' model-specific wrapper should be used for those cases.
#'
#' @seealso [stats::uniroot()], [irr_flow()], [irr_flow_multi()],
#'   [bond_ytm()], [present_value()], [a_angle()]
#'
#' @family time-value
#'
#' @examples
#' # Solve for the annual effective rate
#' solve_value(
#'   fn = present_value,
#'   solve_for = "i",
#'   target = 8000,
#'   args = list(
#'     C = 10000,
#'     t = 5,
#'     i_type = "effective"
#'   ),
#'   interval = c(0, 0.20)
#' )
#'
#' # Solve for a level annuity payment. Because a_angle(tidy = TRUE)
#' # returns a tibble, select the present_value column explicitly.
#' solve_value(
#'   fn = a_angle,
#'   solve_for = "payment",
#'   target = 100000,
#'   args = list(
#'     n = 10,
#'     i = 0.05,
#'     timing = "immediate",
#'     tidy = TRUE
#'   ),
#'   result = "present_value",
#'   interval = c(0, 20000),
#'   output = "summary"
#' )
#'
#' # Newton--Raphson for a simple custom equation
#' solve_value(
#'   fn = function(x) x^2,
#'   solve_for = "x",
#'   target = 2,
#'   start = 1,
#'   interval = c(0, 2),
#'   method = "newton"
#' )
#'
#' @export
solve_value <- function(
    fn,
    solve_for,
    target,
    args = list(),
    result = NULL,
    interval = NULL,
    start = NULL,
    derivative = NULL,
    method = c("auto", "uniroot", "newton"),
    multiple = c("error", "all", "first"),
    tol = 1e-10,
    maxiter = 1000L,
    scan_points = 401L,
    max_expand = 6L,
    output = c("value", "summary", "audit")
) {
  method <- match.arg(method)
  multiple <- match.arg(multiple)
  output <- match.arg(output)

  fn_expr <- substitute(fn)
  fn_name <- if (is.character(fn) && length(fn) == 1L) {
    fn
  } else {
    paste(deparse(fn_expr), collapse = "")
  }
  fn <- match.fun(fn)

  if (!is.character(solve_for) || length(solve_for) != 1L ||
      is.na(solve_for) || !nzchar(solve_for)) {
    stop("`solve_for` must be a single nonempty character string.", call. = FALSE)
  }

  if (!is.numeric(target) || length(target) != 1L || is.na(target) ||
      !is.finite(target)) {
    stop("`target` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.list(args)) {
    stop("`args` must be a named list.", call. = FALSE)
  }
  if (length(args) > 0L && (is.null(names(args)) || any(!nzchar(names(args))))) {
    stop("Every element of `args` must be named.", call. = FALSE)
  }
  if (anyDuplicated(names(args))) {
    stop("Names in `args` must be unique.", call. = FALSE)
  }
  if (solve_for %in% names(args) && !is.null(args[[solve_for]])) {
    stop(
      "Do not provide a known value for `", solve_for,
      "` in `args`; it is the unknown selected by `solve_for`.",
      call. = FALSE
    )
  }
  args[[solve_for]] <- NULL

  fn_formals <- names(formals(fn))
  if (!is.null(fn_formals) &&
      !(solve_for %in% fn_formals) &&
      !("..." %in% fn_formals)) {
    stop(
      "`solve_for = \"", solve_for,
      "\"` is not a formal argument of `fn`.",
      call. = FALSE
    )
  }

  if (!is.null(result) &&
      !is.function(result) &&
      !(is.character(result) && length(result) == 1L &&
        !is.na(result) && nzchar(result))) {
    stop(
      "`result` must be NULL, a single column/element name, or an extractor function.",
      call. = FALSE
    )
  }

  if (!is.null(start) &&
      (!is.numeric(start) || length(start) != 1L || is.na(start) ||
       !is.finite(start))) {
    stop("`start` must be NULL or a single finite numeric value.", call. = FALSE)
  }

  if (!is.null(derivative) && !is.function(derivative)) {
    stop("`derivative` must be NULL or a function.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a single finite positive number.", call. = FALSE)
  }

  .validate_positive_integer <- function(x, name, allow_zero = FALSE) {
    lower <- if (allow_zero) 0 else 1
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < lower || x != floor(x)) {
      qualifier <- if (allow_zero) "nonnegative" else "positive"
      stop("`", name, "` must be a single ", qualifier, " integer.", call. = FALSE)
    }
    as.integer(x)
  }

  maxiter <- .validate_positive_integer(maxiter, "maxiter")
  scan_points <- .validate_positive_integer(scan_points, "scan_points")
  max_expand <- .validate_positive_integer(max_expand, "max_expand", allow_zero = TRUE)

  if (scan_points < 3L) {
    stop("`scan_points` must be at least 3.", call. = FALSE)
  }

  .numeric_scale <- function(args, target) {
    values <- unlist(
      lapply(args, function(x) {
        if (is.numeric(x) && length(x) == 1L && is.finite(x)) x else NULL
      }),
      use.names = FALSE
    )
    max(c(1, abs(target), abs(values)), na.rm = TRUE)
  }

  .default_interval <- function(name, args, target) {
    name_lower <- tolower(name)
    scale <- .numeric_scale(args, target)

    rate_names <- c(
      "i", "rate", "interest", "interest_rate", "yield", "y",
      "discount_rate", "growth_rate", "irr"
    )
    time_names <- c("t", "time", "n", "term", "duration")
    positive_names <- c(
      "c", "cf", "payment", "p", "principal", "price", "face",
      "benefit", "pv", "fv", "present_value", "future_value", "r"
    )

    if (name_lower %in% rate_names) {
      return(c(-0.99, 1))
    }
    if (name_lower %in% time_names) {
      return(c(0, 100))
    }
    if (name_lower %in% positive_names) {
      return(c(0, 10 * scale))
    }

    c(-10 * scale, 10 * scale)
  }

  if (is.null(interval)) {
    interval <- .default_interval(solve_for, args, target)
    interval_inferred <- TRUE
  } else {
    interval_inferred <- FALSE
  }

  if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval) ||
      any(!is.finite(interval)) || interval[1L] == interval[2L]) {
    stop("`interval` must contain two distinct finite numeric values.", call. = FALSE)
  }
  interval <- sort(as.numeric(interval))

  if (!is.null(start) && (start < interval[1L] || start > interval[2L])) {
    stop("`start` must lie inside `interval`.", call. = FALSE)
  }

  .result_abort <- function(...) {
    rlang::abort(
      message = paste0(...),
      class = "tidyactuarial_result_error"
    )
  }

  .extract_scalar <- function(value) {
    extracted <- if (is.function(result)) {
      result(value)
    } else if (is.character(result)) {
      if (is.data.frame(value)) {
        if (!(result %in% names(value))) {
          .result_abort(
            "The result column `", result,
            "` was not returned by `fn`."
          )
        }
        if (nrow(value) != 1L) {
          .result_abort(
            "When `result` names a data-frame column, ",
            "`fn` must return exactly one row."
          )
        }
        value[[result]][1L]
      } else if (is.list(value)) {
        if (is.null(value[[result]])) {
          .result_abort(
            "The result element `", result,
            "` was not returned by `fn`."
          )
        }
        value[[result]]
      } else {
        .result_abort(
          "A character `result` can only extract from ",
          "a data frame or list."
        )
      }
    } else {
      value
    }

    if (!is.numeric(extracted) || length(extracted) != 1L ||
        is.na(extracted) || !is.finite(extracted)) {
      .result_abort(
        "The selected result of `fn` must be one finite numeric scalar. ",
        "Use `result` to select the intended quantity."
      )
    }

    as.numeric(extracted)
  }

  .evaluate_value <- function(x) {
    call_args <- args
    call_args[[solve_for]] <- x
    value <- do.call(fn, call_args)
    .extract_scalar(value)
  }

  .residual <- function(x) {
    .evaluate_value(x) - target
  }

  .safe_residual <- function(x) {
    tryCatch(
      .residual(x),
      error = function(e) {
        if (inherits(e, "tidyactuarial_result_error")) {
          stop(e)
        }

        NA_real_
      }
    )
  }

  .unique_numeric <- function(x, tolerance) {
    x <- sort(x[is.finite(x)])
    if (length(x) <= 1L) return(x)

    keep <- c(TRUE, diff(x) > tolerance * pmax(1, abs(x[-length(x)])))
    x[keep]
  }

  .scan_interval <- function(bounds) {
    grid <- seq(bounds[1L], bounds[2L], length.out = scan_points)
    values <- vapply(grid, .safe_residual, numeric(1))
    finite <- is.finite(values)

    exact_idx <- which(finite & abs(values) <= tol)
    exact_roots <- grid[exact_idx]

    brackets <- list()
    if (length(grid) >= 2L) {
      for (j in seq_len(length(grid) - 1L)) {
        if (!finite[j] || !finite[j + 1L]) next
        if (abs(values[j]) <= tol || abs(values[j + 1L]) <= tol) next
        if (sign(values[j]) != sign(values[j + 1L])) {
          brackets[[length(brackets) + 1L]] <- c(grid[j], grid[j + 1L])
        }
      }
    }

    list(
      grid = grid,
      values = values,
      exact_roots = exact_roots,
      brackets = brackets
    )
  }

  .expand_bounds <- function(bounds) {
    name_lower <- tolower(solve_for)
    rate_names <- c(
      "i", "rate", "interest", "interest_rate", "yield", "y",
      "discount_rate", "growth_rate", "irr"
    )
    positive_names <- c(
      "t", "time", "n", "term", "duration", "c", "cf", "payment",
      "p", "principal", "price", "face", "benefit", "pv", "fv",
      "present_value", "future_value", "r"
    )

    if (name_lower %in% rate_names) {
      return(c(max(-0.999999, bounds[1L]), max(2 * bounds[2L] + 0.01, 0.01)))
    }
    if (name_lower %in% positive_names && bounds[1L] >= 0) {
      return(c(0, max(10 * bounds[2L], 1)))
    }

    center <- mean(bounds)
    half_width <- diff(bounds) / 2
    c(center - 2 * half_width, center + 2 * half_width)
  }

  bounds_used <- interval
  scan <- .scan_interval(bounds_used)
  expansions <- 0L

  if (method %in% c("auto", "uniroot") &&
      interval_inferred &&
      length(scan$exact_roots) == 0L &&
      length(scan$brackets) == 0L &&
      max_expand > 0L) {
    for (expansion_idx in seq_len(max_expand)) {
      new_bounds <- .expand_bounds(bounds_used)
      if (identical(new_bounds, bounds_used)) break
      bounds_used <- new_bounds
      scan <- .scan_interval(bounds_used)
      expansions <- expansion_idx
      if (length(scan$exact_roots) > 0L || length(scan$brackets) > 0L) break
    }
  }

  .newton_solve <- function(x0) {
    x <- x0
    left <- bounds_used[1L]
    right <- bounds_used[2L]

    for (iteration in seq_len(maxiter)) {
      fx <- .residual(x)
      if (abs(fx) <= tol) {
        return(list(root = x, residual = fx, iter = iteration - 1L, converged = TRUE))
      }

      dfx <- if (is.function(derivative)) {
        derivative(x)
      } else {
        h <- sqrt(.Machine$double.eps) * max(1, abs(x))
        f_plus <- .safe_residual(x + h)
        f_minus <- .safe_residual(x - h)
        if (!is.finite(f_plus) || !is.finite(f_minus)) NA_real_ else
          (f_plus - f_minus) / (2 * h)
      }

      if (!is.numeric(dfx) || length(dfx) != 1L || !is.finite(dfx) ||
          abs(dfx) <= sqrt(.Machine$double.eps)) {
        return(list(root = x, residual = fx, iter = iteration, converged = FALSE))
      }

      step <- fx / dfx
      candidate <- x - step

      if (!is.finite(candidate)) {
        candidate <- (left + right) / 2
      } else if (candidate < left || candidate > right) {
        candidate <- min(max(candidate, left), right)
      }

      candidate_residual <- .safe_residual(candidate)
      backtrack <- 0L
      while ((!is.finite(candidate_residual) || abs(candidate_residual) > abs(fx)) &&
             backtrack < 20L) {
        step <- step / 2
        candidate <- x - step
        candidate <- min(max(candidate, left), right)
        candidate_residual <- .safe_residual(candidate)
        backtrack <- backtrack + 1L
      }

      if (!is.finite(candidate_residual)) {
        return(list(root = x, residual = fx, iter = iteration, converged = FALSE))
      }

      if (abs(candidate - x) <= tol * max(1, abs(x)) &&
          abs(candidate_residual) <= sqrt(tol)) {
        return(list(
          root = candidate,
          residual = candidate_residual,
          iter = iteration,
          converged = TRUE
        ))
      }

      x <- candidate
    }

    list(
      root = x,
      residual = .safe_residual(x),
      iter = maxiter,
      converged = FALSE
    )
  }

  roots <- numeric()
  root_methods <- character()
  root_iters <- integer()

  if (method %in% c("auto", "uniroot")) {
    roots <- scan$exact_roots
    root_methods <- rep("grid", length(roots))
    root_iters <- rep(0L, length(roots))

    if (length(scan$brackets) > 0L) {
      bracket_results <- lapply(scan$brackets, function(bracket) {
        root <- stats::uniroot(
          f = .residual,
          interval = bracket,
          tol = tol,
          maxiter = maxiter
        )
        list(root = root$root, iter = root$iter)
      })

      roots <- c(roots, vapply(bracket_results, `[[`, numeric(1), "root"))
      root_methods <- c(root_methods, rep("uniroot", length(bracket_results)))
      root_iters <- c(
        root_iters,
        vapply(bracket_results, `[[`, integer(1), "iter")
      )
    }

    if (length(roots) > 0L) {
      ordered <- order(roots)
      roots <- roots[ordered]
      root_methods <- root_methods[ordered]
      root_iters <- root_iters[ordered]

      keep_roots <- rep(TRUE, length(roots))

      if (length(roots) > 1L) {
        for (root_idx in 2:length(roots)) {
          previous_kept <- max(
            which(keep_roots[seq_len(root_idx - 1L)])
          )

          same_root <- abs(
            roots[[root_idx]] - roots[[previous_kept]]
          ) <= max(tol, 1e-12) *
            max(
              1,
              abs(roots[[root_idx]]),
              abs(roots[[previous_kept]])
            )

          if (same_root) {
            keep_roots[[root_idx]] <- FALSE
          }
        }
      }

      roots <- roots[keep_roots]
      root_methods <- root_methods[keep_roots]
      root_iters <- root_iters[keep_roots]
    }
  }

  if (length(roots) == 0L && method %in% c("auto", "newton")) {
    if (is.null(start)) {
      if (method == "newton") {
        stop("`start` must be supplied when `method = \"newton\"`.", call. = FALSE)
      }
    } else {
      newton_result <- .newton_solve(start)
      if (isTRUE(newton_result$converged)) {
        roots <- newton_result$root
        root_methods <- "newton"
        root_iters <- as.integer(newton_result$iter)
      } else if (method == "newton") {
        stop(
          "Newton--Raphson did not converge within `maxiter` iterations. ",
          "Try a different `start` or provide a bracketing `interval`.",
          call. = FALSE
        )
      }
    }
  }

  if (length(roots) == 0L) {
    stop(
      "No root was found in the searched interval. Provide a more informative ",
      "`interval`, inspect whether the target is attainable, or supply `start` ",
      "for Newton--Raphson.",
      call. = FALSE
    )
  }

  if (length(roots) > 1L) {
    if (multiple == "error") {
      stop(
        "Multiple roots were detected (", length(roots), "). Use ",
        "`multiple = \"all\"` to return them or restrict `interval`.",
        call. = FALSE
      )
    }
    if (multiple == "first") {
      roots <- roots[1L]
      root_methods <- root_methods[1L]
      root_iters <- root_iters[1L]
    }
  }

  achieved <- vapply(roots, .evaluate_value, numeric(1))
  residuals <- achieved - target

  method_for_root <- if (length(root_methods) == length(roots)) {
    root_methods
  } else {
    rep(method, length(roots))
  }

  iter_for_root <- if (length(root_iters) == length(roots)) {
    as.integer(root_iters)
  } else {
    rep(NA_integer_, length(roots))
  }

  audit <- tibble::tibble(
    root_id = seq_along(roots),
    function_name = fn_name,
    solve_for = solve_for,
    solution = roots,
    target = target,
    achieved = achieved,
    residual = residuals,
    method = method_for_root,
    converged = abs(residuals) <= max(tol, sqrt(.Machine$double.eps)),
    n_iter = iter_for_root,
    interval_left = bounds_used[1L],
    interval_right = bounds_used[2L],
    interval_inferred = interval_inferred,
    n_expansions = expansions,
    n_roots = length(roots)
  )

  if (output == "value") {
    return(audit$solution)
  }

  if (output == "summary") {
    return(
      audit |>
        dplyr::select(
          solve_for,
          solution,
          target,
          achieved,
          residual,
          method
        )
    )
  }

  audit
}
