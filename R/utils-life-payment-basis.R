# Shared annual mortality basis for discrete payment grids. Keep probabilities
# in log form and obtain death mass directly, not by subtracting rounded S(t).
.life_payment_basis <- function(lt, x, tol = 1e-10, require_lx = FALSE) {
  if (!is.data.frame(lt) || !all(c("x") %in% names(lt)) ||
      !any(c("lx", "px", "qx") %in% names(lt)) ||
      (require_lx && !"lx" %in% names(lt))) {
    stop("Life table must contain numeric `x` and mortality columns",
         if (require_lx) " including `lx`." else ".", call. = FALSE)
  }
  if (!nrow(lt) || !is.numeric(lt$x) || any(!is.finite(lt$x)) ||
      any(lt$x != round(lt$x)) || any(abs(lt$x) > .Machine$integer.max) ||
      anyDuplicated(lt$x)) {
    stop("Life table ages must be unique finite integer ages.", call. = FALSE)
  }
  lt <- lt[order(lt$x), , drop = FALSE]
  if (any(diff(as.double(lt$x)) != 1))
    stop("Life table must have consecutive annual ages.", call. = FALSE)
  for (column in intersect(c("lx", "px", "qx"), names(lt))) {
    z <- lt[[column]]
    if (!is.numeric(z) || any(is.nan(z)) ||
        any(!is.finite(z) & !is.na(z)) || any(z < 0, na.rm = TRUE) ||
        (column == "lx" && anyNA(z)) ||
        (column != "lx" && any(z > 1, na.rm = TRUE)))
      stop("Invalid mortality values in `", column, "`.", call. = FALSE)
  }
  if ("lx" %in% names(lt) && any(diff(lt$lx) > 0))
    stop("`lx` must be nonincreasing.", call. = FALSE)
  if (all(c("px", "qx") %in% names(lt)) &&
      any(abs(lt$px + lt$qx - 1) > tol, na.rm = TRUE))
    stop("Inconsistent mortality probabilities `px` and `qx`.", call. = FALSE)
  origin <- match(x, lt$x)
  if (is.na(origin) || ("lx" %in% names(lt) && lt$lx[origin] <= 0))
    stop("Starting age must have positive exposure and be present in the table.", call. = FALSE)

  # Annual probabilities at each recorded age, including terminal probabilities
  # only when actually supplied. No implicit extra closing interval is added.
  lp <- lq <- rep(NA_real_, nrow(lt))
  if ("qx" %in% names(lt)) {
    lp <- log1p(-lt$qx)
    lq <- log(lt$qx)
  }
  if ("px" %in% names(lt)) {
    missing <- is.na(lp)
    lp[missing] <- log(lt$px[missing])
    lq[missing] <- log1p(-lt$px[missing])
  }
  if ("lx" %in% names(lt) && nrow(lt) > 1L) {
    j <- seq_len(nrow(lt) - 1L)
    live <- j[lt$lx[j] > 0]
    p <- lt$lx[live + 1L] / lt$lx[live]
    q <- (lt$lx[live] - lt$lx[live + 1L]) / lt$lx[live]
    if (("px" %in% names(lt) && any(abs(lt$px[live] - p) > tol, na.rm = TRUE)) ||
        ("qx" %in% names(lt) && any(abs(lt$qx[live] - q) > tol, na.rm = TRUE)))
      stop("Inconsistent `lx` and mortality probabilities.", call. = FALSE)
    lp[live] <- vapply(live, function(j)
      .single_life_log_ratio(lt$lx[j + 1L], lt$lx[j]), 0)
    # log1p of the direct decrement preserves tiny mortality near p=1.
    lp[live[p > .5]] <- log1p(-q[p > .5])
    lq[live] <- vapply(live, function(j)
      .single_life_log_ratio(lt$lx[j] - lt$lx[j + 1L], lt$lx[j]), 0)
  }
  log_integer <- function(year) {
    if (year == 0) return(0)
    target <- match(as.double(x) + year, lt$x)
    if (is.na(target)) return(NA_real_)
    if ("lx" %in% names(lt))
      return(.single_life_log_ratio(lt$lx[target], lt$lx[origin]))
    values <- lp[seq.int(origin, target - 1L)]
    # Missing mortality before extinction is unknown; after certain extinction
    # it cannot revive the cohort.
    total <- 0
    for (value in values) {
      if (is.na(value)) return(NA_real_)
      if (value == -Inf) return(-Inf)
      total <- total + value
    }
    total
  }
  log_add <- function(a, b) {
    hi <- max(a, b)
    if (hi == -Inf) return(-Inf)
    hi + log1p(exp(min(a, b) - hi))
  }
  log_survival <- function(t, frac = "UDD") {
    year <- floor(t)
    u <- t - year
    base <- log_integer(year)
    if (is.na(base) || base == -Inf || u == 0) return(base)
    j <- match(as.double(x) + year, lt$x)
    if (is.na(j) || is.na(lp[j]) || is.na(lq[j])) return(NA_real_)
    extra <- switch(frac,
      UDD = log1p(-u * exp(lq[j])),
      CF = u * lp[j],
      Balducci = lp[j] - log_add(lp[j], log(u) + lq[j]))
    base + extra
  }
  log_death <- function(t0, t1, frac = "UDD") {
    year <- floor(t0)
    u0 <- t0 - year
    u1 <- t1 - year
    if (u1 > 1 + tol || u1 <= u0)
      stop("Death interval must lie within one annual interval.", call. = FALSE)
    u1 <- min(1, u1)
    base <- log_integer(year)
    if (is.na(base) || base == -Inf) return(base)
    j <- match(as.double(x) + year, lt$x)
    if (is.na(j) || is.na(lp[j]) || is.na(lq[j])) return(NA_real_)
    if (lq[j] == -Inf) return(-Inf)
    if (frac != "UDD" && lp[j] == -Inf)
      return(if (u0 == 0) base else -Inf)
    mass <- switch(frac,
      UDD = lq[j] + log(u1 - u0),
      CF = u0 * lp[j] + log(-expm1((u1 - u0) * lp[j])),
      Balducci = lp[j] + lq[j] + log(u1 - u0) -
        log_add(lp[j], if (u0 == 0) -Inf else log(u0) + lq[j]) -
        log_add(lp[j], log(u1) + lq[j]))
    base + mass
  }
  list(log_survival = log_survival, log_death = log_death,
       omega = max(lt$x), table = lt)
}
