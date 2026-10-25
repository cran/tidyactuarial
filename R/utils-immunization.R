# Solve a full-row-rank moment system without forming normal equations.
# Scaling the constraint rows preserves the feasible set and its minimum-norm
# solution while preventing currency units from controlling the rank test.
.immunization_min_norm <- function(A, b) {
  row_scale <- apply(abs(A), 1, max)
  if (any(!is.finite(A)) || any(!is.finite(b)) ||
      any(!is.finite(row_scale)) || any(row_scale == 0)) {
    stop("System matrix is singular or non-finite; cannot compute minimum-norm solution.",
         call. = FALSE)
  }
  scaled <- A / row_scale
  rhs <- b / row_scale
  s <- svd(scaled)
  threshold <- max(dim(scaled)) * .Machine$double.eps * s$d[1L]
  if (any(s$d <= threshold)) {
    stop("System matrix is singular; cannot compute minimum-norm solution.",
         call. = FALSE)
  }
  w <- as.numeric(s$v %*% (crossprod(s$u, rhs) / s$d))
  if (any(!is.finite(w))) {
    stop("Immunization weights are non-finite. Check input magnitudes.", call. = FALSE)
  }
  w
}
