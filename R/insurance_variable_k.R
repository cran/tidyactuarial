#' Actuarial present value of a life insurance with variable k-thly benefits
#'
#' Computes the actuarial present value of a life insurance where the death
#' benefit may vary by subperiod and is payable at the end of the subperiod of
#' death, using compact actuarial notation.
#'
#' This function is useful for level, increasing, decreasing, and credit-style
#' life insurance benefits when benefits are specified at a subannual frequency.
#'
#' @param lt A life table object or data frame containing at least columns
#'   \code{x} and \code{lx}.
#' @param x Integer actuarial age at issue.
#' @param i Annual interest-rate input.
#' @param benefit Numeric vector of benefits by subperiod, or a function of
#'   time returning the benefit at time \eqn{t}.
#' @param n Optional term in years. If \code{NULL}, the term is inferred from
#'   the length of \code{benefit} when \code{benefit} is numeric. If
#'   \code{benefit} is a function, \code{n} must be supplied.
#' @param h Nonnegative deferment period in years. Default is \code{0}.
#' @param k Positive integer. Number of subperiods per year. Default is
#'   \code{12}.
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   interest-rate inputs. Ignored for \code{"effective"} and \code{"force"}.
#' @param frac Fractional-age assumption used in survival probabilities:
#'   \code{"UDD"}, \code{"CF"}, \code{"CML"}, or \code{"Balducci"}. If not
#'   specified and \code{lt} carries a \code{frac} attribute, that value is
#'   used.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble.
#' @param check Logical. If \code{TRUE}, performs basic input checks.
#' @param tol Numeric tolerance used for integer-grid checks.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{lt} is the life table, \code{x} is the age at
#' issue, \code{i} is the interest-rate input, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' \code{n} is the term, \code{h} is the deferment period, and \code{k} is the
#' number of subperiods per year.
#'
#' Let \eqn{k} be the number of subperiods per year and \eqn{N = nk} the total
#' number of subperiods in the insurance term. With deferment \eqn{h}, the
#' actuarial present value at age \eqn{x} is:
#' \deqn{
#'   APV =
#'   \sum_{j=1}^{N}
#'   v^{h + j/k} b_j
#'   \left({}_{h + (j-1)/k}p_x - {}_{h + j/k}p_x\right).
#' }
#'
#' Here \eqn{b_j} is the benefit payable if death occurs in subperiod
#' \eqn{j}, and \eqn{v = (1+i_e)^{-1}}, where \eqn{i_e} is the annual
#' effective interest rate equivalent to the input \code{i}, \code{i_type}, and
#' \code{m}.
#'
#' If \code{benefit} is numeric of length 1, it is recycled to all subperiods.
#' If it is numeric with length greater than 1, its length must equal
#' \eqn{n k}. If \code{benefit} is a function, it is evaluated at the end of
#' each subperiod, at times \eqn{1/k, 2/k, \ldots, n}.
#'
#' Fractional survival probabilities are computed via \code{\link{t_px}} under
#' the selected fractional-age assumption.
#'
#' @return A numeric actuarial present value, or a one-row tibble if
#'   \code{tidy = TRUE}.
#'
#' @seealso \code{\link{insurance_x}} for level-benefit life insurance,
#'   \code{\link{annuity_x}} for life annuity APVs,
#'   \code{\link{t_px}} for survival probabilities.
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Monthly insurance with increasing benefits
#' insurance_variable_k(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   benefit = seq(100, 1200, length.out = 12),
#'   n = 1,
#'   k = 12
#' )
#'
#' # Credit-style insurance with declining outstanding balance
#' balance <- function(t) 2000 * exp(-0.3 * t)
#'
#' insurance_variable_k(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   benefit = balance,
#'   n = 1,
#'   k = 12
#' )
#'
#' # Level benefit with annual payments
#' insurance_variable_k(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   benefit = 1,
#'   n = 5,
#'   k = 1
#' )
#'
#' # 2-year deferred, 3-year term with monthly varying benefits
#' insurance_variable_k(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   benefit = rep(1000, 36),
#'   n = 3,
#'   h = 2,
#'   k = 12
#' )
#'
#' # Nominal annual interest convertible monthly
#' insurance_variable_k(
#'   lt = lt,
#'   x = 60,
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   benefit = rep(1000, 12),
#'   n = 1,
#'   k = 12
#' )
#'
#' # Tidy output
#' insurance_variable_k(
#'   lt = lt,
#'   x = 60,
#'   i = 0.05,
#'   benefit = rep(1000, 12),
#'   n = 1,
#'   k = 12,
#'   tidy = TRUE
#' )
#'
#' @export
insurance_variable_k <- function(
    lt,
    x,
    i,
    benefit,
    n = NULL,
    h = 0,
    k = 12,
    i_type = "effective",
    m = 1L,
    frac,
    tidy = FALSE,
    check = TRUE,
    tol = 1e-10
) {
  # --- inherit frac from lifetable attribute if not supplied ---
  if (missing(frac)) {
    lt_frac <- attr(lt, "frac")

    if (!is.null(lt_frac) && lt_frac %in% c("UDD", "CF", "Balducci")) {
      frac <- lt_frac
    } else {
      frac <- "UDD"
    }
  } else {
    frac <- match.arg(frac, c("UDD", "CF", "CML", "Balducci"))
    if (frac == "CML") frac <- "CF"
  }

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be a logical scalar.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  if (isTRUE(check)) {
    if (!is.data.frame(lt)) {
      stop("`lt` must be a data frame or tibble.", call. = FALSE)
    }

    if (!all(c("x", "lx") %in% names(lt))) {
      stop("`lt` must contain columns `x` and `lx`.", call. = FALSE)
    }

    if (missing(x) ||
        !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        abs(x - round(x)) > tol) {
      stop("`x` must be a single integer age.", call. = FALSE)
    }

    if (missing(i) ||
        !is.numeric(i) ||
        length(i) != 1L ||
        is.na(i) ||
        !is.finite(i)) {
      stop("`i` must be a single finite numeric value.", call. = FALSE)
    }

    if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
      stop("`i_type` must be a single character string.", call. = FALSE)
    }

    valid_i_type <- c(
      "effective",
      "nominal_interest",
      "nominal_discount",
      "force"
    )

    if (!i_type %in% valid_i_type) {
      stop(
        "`i_type` must be one of: ",
        paste(sprintf("'%s'", valid_i_type), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.numeric(m) ||
        length(m) != 1L ||
        is.na(m) ||
        !is.finite(m) ||
        m < 1 ||
        abs(m - round(m)) > tol) {
      stop("`m` must be a single positive integer.", call. = FALSE)
    }

    if (!is.numeric(h) ||
        length(h) != 1L ||
        is.na(h) ||
        !is.finite(h) ||
        h < 0) {
      stop("`h` must be a single finite nonnegative value.", call. = FALSE)
    }

    if (!is.numeric(k) ||
        length(k) != 1L ||
        is.na(k) ||
        !is.finite(k) ||
        k <= 0 ||
        abs(k - round(k)) > tol) {
      stop("`k` must be a single positive integer.", call. = FALSE)
    }

    if (!is.null(n)) {
      if (!is.numeric(n) ||
          length(n) != 1L ||
          is.na(n) ||
          !is.finite(n) ||
          n <= 0) {
        stop("`n` must be NULL or a single finite positive value.", call. = FALSE)
      }
    }

    if (!is.numeric(lt$x) || !is.numeric(lt$lx)) {
      stop("Columns `x` and `lx` in `lt` must be numeric.", call. = FALSE)
    }
  }

  x <- as.integer(round(x))
  k <- as.integer(round(k))
  m <- as.integer(round(m))

  h_raw <- h * k
  h_periods <- round(h_raw)

  if (abs(h_raw - h_periods) > tol) {
    stop(
      "For k-thly insurance, `h * k` must be an integer.",
      call. = FALSE
    )
  }

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  # --- determine N (total subperiods) ---
  if (is.null(n)) {
    if (is.function(benefit)) {
      stop("Provide `n` when `benefit` is a function.", call. = FALSE)
    }

    if (!is.numeric(benefit) || length(benefit) == 0L) {
      stop("`benefit` must be a numeric vector or a function.", call. = FALSE)
    }

    N <- length(benefit)
    n <- N / k
  } else {
    N_raw <- n * k
    N <- round(N_raw)

    if (abs(N_raw - N) > tol) {
      stop(
        "For k-thly insurance, `n * k` must be an integer.",
        call. = FALSE
      )
    }

    N <- as.integer(N)
  }

  if (N <= 0L) {
    stop("The total number of subperiods `n * k` must be positive.", call. = FALSE)
  }

  # --- build benefit vector ---
  if (is.function(benefit)) {
    times <- (1:N) / k
    bvec <- benefit(times)
  } else {
    if (!is.numeric(benefit) || length(benefit) == 0L) {
      stop("`benefit` must be a numeric vector or a function.", call. = FALSE)
    }

    if (length(benefit) == 1L) {
      benefit <- rep(benefit, N)
    }

    if (length(benefit) != N) {
      stop("`benefit` length must equal `n * k = ", N, "`.", call. = FALSE)
    }

    bvec <- benefit
  }

  if (!is.numeric(bvec) ||
      length(bvec) != N ||
      any(is.na(bvec)) ||
      any(!is.finite(bvec))) {
    stop("The evaluated benefit vector must contain `n * k` finite numeric values.",
         call. = FALSE)
  }

  # --- compute APV ---
  t0_vec <- h + (0:(N - 1L)) / k
  t1_vec <- h + (1:N) / k

  S_t0 <- t_px(
    lt = lt,
    x = x,
    t = t0_vec,
    frac = frac,
    check = FALSE
  )

  S_t1 <- t_px(
    lt = lt,
    x = x,
    t = t1_vec,
    frac = frac,
    check = FALSE
  )

  dq <- S_t0 - S_t1
  disc <- (1 + i_effective)^(-t1_vec)

  apv <- sum(bvec * disc * dq)

  if (!isTRUE(tidy)) {
    return(apv)
  }

  tibble::tibble(
    x = x,
    h = h,
    n = n,
    k = k,
    i = i,
    i_type = i_type,
    m = m,
    i_effective = i_effective,
    frac = frac,
    apv = apv
  )
}
