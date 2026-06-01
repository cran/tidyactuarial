#' Actuarial present value of a multi-life annuity (up to 3 independent lives)
#'
#' Computes the APV of a discrete annuity contingent on multiple independent
#' lives using compact actuarial notation.
#'
#' This implementation supports up to three lives, which covers the most common
#' practical multi-life arrangements.
#'
#' Supports status-based annuities (joint-life / last-survivor) and a
#' joint-and-survivor style annuity (`"reversionary"`) that pays 1 while all
#' lives are alive and then pays a fraction \eqn{\alpha} while at least one life
#' remains alive.
#'
#' @param lt A life table object (data frame or tibble) with column \code{x}
#'   and at least one of \code{lx}, \code{px}, or \code{qx}. For different
#'   mortality assumptions by life, pass a list of life tables of the same
#'   length as \code{ages}, for example \code{list(lt_male, lt_female)}.
#' @param ages Integer vector of actuarial ages for the lives at issue. Must
#'   have length 1, 2, or 3.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}. In
#'   \code{tidyactuarial}, \code{m} is reserved for interest conversion
#'   frequency, not deferment.
#' @param n Integer term in years after deferment. If \code{NULL}, runs to the
#'   end of the available table horizon, conservatively based on the smallest
#'   omega across life tables.
#' @param h Integer deferment period in years.
#' @param k Integer payments per year. If \code{k > 1}, Woolhouse
#'   approximations may be applied.
#' @param annuity Type of annuity logic: \code{"cohort"} uses the status
#'   defined in \code{cohort}; \code{"reversionary"} uses the \eqn{\alpha}
#'   fractional reduction.
#' @param cohort Survival status: \code{"first"} for joint-life, paying while
#'   all lives are alive, or \code{"last"} for last-survivor, paying while at
#'   least one life remains alive. Used only when \code{annuity = "cohort"}.
#' @param alpha Reversionary fraction, typically \eqn{0 \le \alpha \le 1}.
#'   Used only when \code{annuity = "reversionary"}. Note: \code{alpha = 0}
#'   matches joint-life; \code{alpha = 1} matches last-survivor.
#' @param timing \code{"immediate"} or \code{"due"}.
#' @param woolhouse \code{"none"}, \code{"first"}, or \code{"second"}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{i} is the interest rate, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' \code{n} is the term, \code{h} is the deferment period, and \code{k} is the
#' payment frequency.
#'
#' Under the assumption of independent future lifetimes, the survival
#' probability for the status is calculated as:
#' \itemize{
#'   \item \strong{Joint-life (first-death):}
#'     \eqn{{}_t p_{x_1 x_2 \dots x_n} = \prod_{j=1}^n {}_t p_{x_j}}
#'   \item \strong{Last-survivor:}
#'     \eqn{{}_t p_{\overline{x_1 x_2 \dots x_n}} = 1 - \prod_{j=1}^n (1 - {}_t p_{x_j})}
#' }
#'
#' For \code{annuity = "reversionary"}, the APV is a weighted combination of
#' the two statuses:
#' \deqn{APV = APV(\text{joint-life}) +
#' \alpha [APV(\text{last-survivor}) - APV(\text{joint-life})].}
#'
#' @return A single numeric value representing the APV.
#'
#' @seealso \code{\link{annuity_x}} for single-life annuities,
#'   \code{\link{t_px}} for survival probabilities.
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 60:90,
#'   lx = seq(100000, 0, length.out = 31)
#' )
#'
#' annuity_multi(
#'   lt = lt,
#'   ages = c(60, 62),
#'   i = 0.05,
#'   n = 5,
#'   cohort = "first",
#'   timing = "due"
#' )
#'
#' annuity_multi(
#'   lt = lt,
#'   ages = c(60, 62),
#'   i = 0.05,
#'   n = 5,
#'   cohort = "last",
#'   timing = "due"
#' )
#'
#' @export
annuity_multi <- function(
    lt,
    ages,
    i,
    i_type = "effective",
    m = 1L,
    n = NULL,
    h = 0L,
    k = 1L,
    annuity = c("cohort", "reversionary"),
    cohort = c("first", "last"),
    alpha = NULL,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second")
) {
  annuity <- match.arg(annuity)
  cohort <- match.arg(cohort)
  timing <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)

  if (missing(i) || !is.numeric(i) || length(i) != 1L || is.na(i) ||
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

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m < 1 || abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }
  m <- as.integer(round(m))

  if (!is.numeric(ages) || length(ages) < 1L) {
    stop("`ages` must be a non-empty numeric vector.", call. = FALSE)
  }
  ages <- as.integer(round(ages))

  if (length(ages) > 3L) {
    stop("`ages` must contain at most 3 lives (length <= 3).", call. = FALSE)
  }

  # --- normalize lt to a list of lifetables (one per life) ---
  if (is.data.frame(lt)) {
    lt_list <- rep(list(lt), length(ages))
  } else if (is.list(lt) && length(lt) >= 1L &&
             all(vapply(lt, is.data.frame, logical(1)))) {
    if (length(lt) == 1L) {
      lt_list <- rep(lt, length(ages))
    } else if (length(lt) == length(ages)) {
      lt_list <- lt
    } else {
      stop(
        "When `lt` is a list, it must have length 1 or length equal to `length(ages)`.",
        call. = FALSE
      )
    }
  } else {
    stop("`lt` must be a data frame or a list of data frames.", call. = FALSE)
  }

  # validate lifetable structure
  for (j in seq_along(lt_list)) {
    if (!("x" %in% names(lt_list[[j]]))) {
      stop("Each life table must contain column `x`.", call. = FALSE)
    }

    if (!("lx" %in% names(lt_list[[j]])) &&
        !("px" %in% names(lt_list[[j]])) &&
        !("qx" %in% names(lt_list[[j]]))) {
      stop("Each life table must contain `lx`, `px`, or `qx`.", call. = FALSE)
    }

    lt_list[[j]] <- lt_list[[j]][order(lt_list[[j]]$x), ]
  }

  h <- as.integer(round(h))
  k <- as.integer(round(k))

  if (h < 0) {
    stop("`h` must be a single nonnegative integer.", call. = FALSE)
  }

  if (k < 1) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }

  if (annuity == "reversionary") {
    if (is.null(alpha) || !is.numeric(alpha) || length(alpha) != 1L ||
        is.na(alpha)) {
      stop("`alpha` must be a single numeric value.", call. = FALSE)
    }
  }

  # --- Interest conversion ---
  i_effective <- standardize_interest(
    type = i_type,
    rate = i,
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

  # Conservative omega: smallest last age across lifetables
  omega_min <- min(vapply(
    lt_list,
    function(LT) max(LT$x, na.rm = TRUE),
    numeric(1)
  ))

  # term after deferment
  if (is.null(n)) {
    n <- max(0L, omega_min - max(ages) - h)
  } else {
    n <- as.integer(round(n))
    if (n < 0) {
      stop("`n` must be nonnegative or NULL.", call. = FALSE)
    }
  }

  if (n == 0L) {
    return(0)
  }

  # Ensure the lifetables support the required ages (conservative but safe)
  if (max(ages) + h + n > omega_min) {
    stop(
      "Life table horizon insufficient for requested deferment/term ",
      "(based on smallest omega).",
      call. = FALSE
    )
  }

  v <- function(t) (1 + i_effective)^(-t)

  # integer-year survival for one life with its own lifetable: {}_t p_age
  t_px1 <- function(lt1, age, t) {
    if (t == 0) return(1)

    if ("lx" %in% names(lt1)) {
      l0 <- lt1$lx[match(age, lt1$x)]
      l1 <- lt1$lx[match(age + t, lt1$x)]
      if (is.na(l0) || is.na(l1) || l0 <= 0) return(NA_real_)
      return(l1 / l0)
    }

    if ("px" %in% names(lt1)) {
      pxv <- lt1$px[match(age:(age + t - 1), lt1$x)]
      if (anyNA(pxv)) return(NA_real_)
      return(prod(pxv))
    }

    qxv <- lt1$qx[match(age:(age + t - 1), lt1$x)]
    if (anyNA(qxv)) return(NA_real_)
    prod(1 - qxv)
  }

  # Prob(all alive at time u) and Prob(any alive at time u), from issue ages
  P_all_u <- function(u) {
    p_vec <- vapply(
      seq_along(ages),
      function(j) t_px1(lt_list[[j]], ages[j], u),
      numeric(1)
    )

    if (anyNA(p_vec)) return(NA_real_)

    prod(p_vec)
  }

  P_any_u <- function(u) {
    p_vec <- vapply(
      seq_along(ages),
      function(j) t_px1(lt_list[[j]], ages[j], u),
      numeric(1)
    )

    if (anyNA(p_vec)) return(NA_real_)

    1 - prod(1 - p_vec)
  }

  # Expected payment at absolute time u (years since issue)
  E_pay_u <- function(u) {
    pa <- P_all_u(u)
    pn <- P_any_u(u)

    if (is.na(pa) || is.na(pn)) return(NA_real_)

    if (annuity == "cohort") {
      if (cohort == "first") return(pa)
      return(pn)
    }

    # reversionary: 1 while all alive; alpha while partially alive
    pa + alpha * (pn - pa)
  }

  # Annual APV, payments at u = h + t, with t depending on timing
  annual_apv <- function(n, timing) {
    tt <- if (timing == "due") 0:(n - 1L) else 1:n
    u <- h + tt

    ep <- vapply(u, E_pay_u, numeric(1))

    if (anyNA(ep)) {
      stop("Life table does not support required ages for this annuity.", call. = FALSE)
    }

    sum(v(u) * ep)
  }

  # Exact k-thly computation is not implemented here. With woolhouse = "none",
  # the annual base is returned, matching the previous implementation.
  if (k == 1L || woolhouse == "none") {
    return(annual_apv(n, timing))
  }

  # Woolhouse: adjust annual annuity-due value
  adue <- annual_apv(n, "due")

  # PV of expected payment at the first due date (time h)
  pv_first <- v(h) * E_pay_u(h)

  adj1 <- (k - 1) / (2 * k)

  if (woolhouse == "first") {
    adue_k <- adue - pv_first * adj1
  } else {
    delta <- log1p(i_effective)

    # Approximate mu of the STATUS around time h using one-year survival ratio:
    # p_status1 approx status_prob(h + 1) / status_prob(h).
    p_h  <- E_pay_u(h)
    p_h1 <- E_pay_u(h + 1L)

    if (!is.finite(p_h) || p_h <= 0 || !is.finite(p_h1) || p_h1 < 0) {
      stop(
        "Cannot compute second-order Woolhouse adjustment ",
        "(status probability at h is nonpositive).",
        call. = FALSE
      )
    }

    p_status1 <- p_h1 / p_h

    if (p_status1 <= 0) {
      stop("Cannot compute mu: derived status one-year survival <= 0.", call. = FALSE)
    }

    mu <- -log(p_status1)

    adj2 <- (k^2 - 1) / (12 * k^2) * (delta + mu)

    adue_k <- adue - pv_first * (adj1 + adj2)
  }

  if (timing == "due") return(adue_k)

  # Convert due -> immediate (approx): remove 1/k of the first due payment PV
  adue_k - pv_first * (1 / k)
}
