# Log of sum(exp(rate * j)), j = 0,...,periods-1, for scalar inputs.
# Reverse increasing terms before summing to avoid overflowing the numerator.
.annuity_log_sum <- function(rate, periods) {
  if (rate == 0) return(log(periods))
  magnitude <- abs(rate)
  log(-expm1(-magnitude * periods)) - log(-expm1(-magnitude)) +
    max(rate, 0) * (periods - 1)
}

# Log of integral_0^term exp(rate * t) dt, with its exact zero-rate limit.
.annuity_log_integral <- function(rate, term) {
  magnitude <- abs(rate)
  span <- magnitude * term
  if (span == 0) return(log(term))
  if (!is.finite(span)) return(if (rate > 0) Inf else -log(magnitude))
  # Scaling by term avoids dividing a rounded subnormal span by the rate.
  log(term) + log(-expm1(-span) / span) + if (rate > 0) span else 0
}
