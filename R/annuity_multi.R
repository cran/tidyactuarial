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
#' @param frac Fractional-age assumption for exact payments: \code{"UDD"}
#'   (default), \code{"CF"}, \code{"CML"} (alias of CF), or \code{"Balducci"}.
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
#' All values are measured at issue, conditional on the stated lives being
#' alive at issue. Ages identify lives within one policy, not policy scenarios.
#' A k-thly payment is \eqn{1/k} of an annual unit: exact valuation sums
#' \eqn{F(h+j/k)/k}, with \eqn{F(t)=e^{-\delta t}w(t)} and the expected payment
#' weight \eqn{w(t)} defined above. Due uses \eqn{j=0,\ldots,nk-1}; immediate
#' uses \eqn{j=1,\ldots,nk}. The historical name \code{"reversionary"} means
#' joint-and-survivor payments; it does not mean payments only after a death.
#'
#' For \code{k > 1}, \code{woolhouse = "none"} uses the actual fractional
#' payment grid. First-order Woolhouse subtracts
#' \eqn{(1-1/k)[F(h)-F(h+n)]/2} from the annual due value. Second-order also
#' subtracts \eqn{(1-1/k^2)[(\delta+\eta_h)F(h)-
#' (\delta+\eta_{h+n})F(h+n)]/12}. Here
#' \eqn{\eta_t=-\log(w(t+1)/w(t))} is an annual forward estimate of the
#' logarithmic decay of the payment weight, not an exact instantaneous force.
#' Positive supported weights at both borders and the following year are
#' required for second-order approximation. This is a smoothness approximation,
#' not exact fractional valuation under \code{frac}; negative approximations
#' produce a warning. Immediate subtracts \eqn{[F(h)-F(h+n)]/k} from due.
#' Exact valuation requires consecutive annual ages and retains the documented
#' conservative common-age horizon; it does not extrapolate unobserved lives.
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
#' # Half-yearly payments, each equal to half an annual unit, under UDD
#' annuity_multi(lt, ages = c(60, 62), i = 0.05, n = 5, k = 2,
#'               timing = "due", frac = "UDD")
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
    lt, ages, i, i_type = "effective", m = 1L, n = NULL, h = 0L,
    k = 1L, annuity = c("cohort", "reversionary"),
    cohort = c("first", "last"), alpha = NULL,
    timing = c("immediate", "due"),
    woolhouse = c("none", "first", "second"),
    frac = c("UDD", "CF", "CML", "Balducci")
) {
  annuity <- match.arg(annuity)
  cohort <- match.arg(cohort)
  timing <- match.arg(timing)
  woolhouse <- match.arg(woolhouse)
  frac <- match.arg(frac)
  if (frac == "CML") frac <- "CF"
  integer_scalar <- function(z, minimum) {
    is.numeric(z) && length(z) == 1L && is.finite(z) &&
      z >= minimum && abs(z - round(z)) <= 1e-10
  }
  if (missing(i) || !is.numeric(i) || length(i) != 1L || !is.finite(i))
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type) ||
      !i_type %in% c("effective", "nominal_interest", "nominal_discount", "force"))
    stop("`i_type` must be a single supported rate type.", call. = FALSE)
  for (name in c("m", "h", "k")) {
    value <- get(name)
    if (!integer_scalar(value, if (name == "h") 0 else 1))
      stop("`", name, "` must be a single integer in its valid range.", call. = FALSE)
  }
  m <- round(m); k <- round(k); h <- round(h)
  if (!is.numeric(ages) || !length(ages) || length(ages) > 3L ||
      any(!is.finite(ages)) || any(abs(ages-round(ages)) > 1e-10) ||
      any(abs(ages) > .Machine$integer.max))
    stop("`ages` must contain 1 to 3 finite integer ages.", call. = FALSE)
  ages <- round(ages)
  if (annuity == "reversionary" &&
      (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha)))
    stop("`alpha` must be a single finite numeric value.", call. = FALSE)
  if (is.data.frame(lt)) {
    tables <- rep(list(lt),length(ages))
  } else if (is.list(lt) && length(lt) %in% c(1L,length(ages)) &&
             all(vapply(lt,is.data.frame,logical(1)))) {
    tables <- rep(lt,length.out=length(ages))
  } else {
    stop("`lt` must be a table or list of one table per life.", call. = FALSE)
  }
  mortality <- lapply(seq_along(ages),function(j)
    .life_payment_basis(tables[[j]],ages[j]))
  basis <- .interest_basis(i_type,i,m)
  # Preserve the explicitly documented conservative common-age horizon.
  omega <- min(vapply(mortality,function(z) z$omega,0))
  if (is.null(n)) n <- max(0,omega-max(ages)-h)
  if (!integer_scalar(n,0))
    stop("`n` must be a nonnegative integer or NULL.", call. = FALSE)
  n <- round(n)
  if (n == 0) return(0)
  if (!is.finite(h+n) || max(ages)+h+n > omega)
    stop("Life table horizon insufficient (based on smallest omega).", call. = FALSE)
  if (!is.finite(n*k) || n*k > .Machine$integer.max)
    stop("Number of payments exceeds the supported payment grid.", call. = FALSE)

  # Return sign and log absolute expected benefit, retaining rare survivors.
  weight <- function(t) {
    logs <- vapply(mortality,function(z) z$log_survival(t,frac),0)
    if (anyNA(logs)) stop("Life table does not support required ages.",call.=FALSE)
    all <- sum(logs)
    hi <- max(logs)
    any <- if (hi == -Inf) -Inf else if (hi < -36)
      hi + log(sum(exp(logs-hi))) else
      log(-expm1(sum(log1p(-exp(logs)))))
    selected <- if (cohort == "first") all else any
    if (annuity == "cohort") return(c(sign=1,log=selected))
    if (all == any || alpha == 0) return(c(sign=1,log=all))
    if (alpha == 1) return(c(sign=1,log=any))
    scaled <- alpha + (1-alpha)*exp(all-any)
    c(sign=sign(scaled),log=any+log(abs(scaled)))
  }
  payment <- function(t) {
    w <- weight(t)
    if (w[["log"]] == -Inf) return(0)
    .cash_flow_value(w[["sign"]],w[["log"]]-basis$growth*t)
  }
  annual <- function() sum(vapply(h+0:(n-1),payment,0))
  if (k == 1 || woolhouse == "none") {
    N <- as.integer(n*k)
    times <- h + (if (timing == "due") 0:(N-1L) else seq_len(N))/k
    weights <- vapply(times,weight,c(sign=0,log=0))
    logs <- vapply(seq_along(times),function(j) {
      w <- weights[,j]
      t <- times[j]
      if (w[["log"]] == -Inf) -Inf else w[["log"]]-basis$growth*t
    },0)
    signs <- weights["sign",]
    return(.cash_flow_sum(signs/k,logs))
  }
  # Euler-Maclaurin corrections at BOTH ends of the finite payment interval.
  Fh <- payment(h)
  Fn <- payment(h+n)
  loss <- .cash_flow_sum(c(Fh,-Fn),c(0,0))
  result <- annual()-(1-1/k)/2*loss
  if (woolhouse == "second") {
    force <- function(t) {
      boundary_error <- function(e)
        stop("Second-order Woolhouse requires supported positive border survival; use `none`.",call.=FALSE)
      a <- tryCatch(weight(t),error=boundary_error)
      b <- tryCatch(weight(t+1),error=boundary_error)
      if (a[["sign"]] <= 0 || b[["sign"]] <= 0 ||
          !is.finite(a[["log"]]) || !is.finite(b[["log"]]))
        stop("Second-order Woolhouse requires positive supported border survival/weights; use `none`.",call.=FALSE)
      -(b[["log"]]-a[["log"]])
    }
    correction <- .cash_flow_sum(c(basis$growth+force(h),
                                  -(basis$growth+force(h+n))),
                                c(if (Fh == 0) -Inf else log(abs(Fh)),
                                  if (Fn == 0) -Inf else log(abs(Fn))))
    result <- result-(1-1/k^2)/12*correction
  }
  if (timing == "immediate") result <- result-loss/k
  if (is.na(result)) stop("Woolhouse approximation undefined; use `none`.",call.=FALSE)
  if (result < 0) warning("Woolhouse approximation is negative; use `none`.",call.=FALSE)
  result
}
