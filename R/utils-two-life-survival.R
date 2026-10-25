# Validate the observed annual basis without imposing a starting-life state.
# Missing or extinct initial lives retain the consumer's established policy.
.two_life_table <- function(lt) {
  if (!is.data.frame(lt) || !all(c("x", "lx") %in% names(lt)))
    stop("Life table must contain columns `x` and `lx`.", call. = FALSE)
  if (!nrow(lt) || !is.numeric(lt$x) || any(!is.finite(lt$x)) ||
      any(lt$x != round(lt$x)) || any(abs(lt$x) > .Machine$integer.max) ||
      anyDuplicated(lt$x) || any(diff(sort(as.double(lt$x))) != 1))
    stop("Life table must have unique consecutive finite integer ages.", call. = FALSE)
  ord <- order(lt$x)
  lx <- lt$lx[ord]
  if (!is.numeric(lx) || any(!is.finite(lx)) || any(lx < 0) || any(diff(lx) > 0))
    stop("`lx` must be finite, nonnegative and nonincreasing by age.", call. = FALSE)
  for (column in intersect(c("px", "qx"), names(lt))) {
    z <- lt[[column]][ord]
    if (!is.numeric(z) || any(is.nan(z)) ||
        any(!is.finite(z) & !is.na(z)) || any(z < 0 | z > 1, na.rm = TRUE))
      stop("Invalid mortality probabilities in `", column, "`.", call. = FALSE)
    if (length(lx) > 1L) {
      j <- seq_len(length(lx) - 1L)
      j <- j[lx[j] > 0]
      expected <- if (column == "px") lx[j + 1L]/lx[j] else (lx[j]-lx[j + 1L])/lx[j]
      if (any(abs(z[j] - expected) > 1e-10, na.rm = TRUE))
        stop("Inconsistent `lx` and mortality probabilities.", call. = FALSE)
    }
  }
  if (all(c("px", "qx") %in% names(lt)) &&
      any(abs(lt$px + lt$qx - 1) > 1e-10, na.rm = TRUE))
    stop("Inconsistent mortality probabilities `px` and `qx`.", call. = FALSE)
  invisible(NULL)
}

# Tables and controls are validated by callers, once per valuation.
.two_life_survival <- function(lt_x, lt_y, x, y, t, frac, status) {
  .two_life_survival_basis(lt_x, lt_y, x, y)$survival(t, frac, status)
}

# Build once per scenario. Observed lx ratios retain small survival which
# would be lost by forming p = 1-q. The payment basis supplies stable logs.
.two_life_survival_basis <- function(lt_x, lt_y, x, y) {
  marginal <- function(lt, age) {
    initial <- match(age, lt$x)
    valid <- !is.na(initial) && lt$lx[initial] > 0
    basis <- if (valid) .life_payment_basis(lt, age, require_lx=TRUE) else NULL
    horizon <- max(lt$x)-as.double(age)
    survival <- function(t, frac) {
      if (t == 0) return(1)
      if (t > horizon + 1e-10) return(0)
      if (!valid) return(t_px(lt,age,t,frac=frac,check=FALSE))
      exp(basis$log_survival(t,frac))
    }
    parameters <- function(k) {
      # At or beyond an open marginal's horizon, only the left endpoint
      # can have positive survival; its area is zero under truncation.
      if (k >= horizon) return(c(base=-Inf,p=0,q=-Inf))
      tbl <- basis$table
      j <- match(as.double(age)+k,tbl$x)
      if (tbl$lx[j] == 0) return(c(base=-Inf,p=0,q=-Inf))
      q <- (tbl$lx[j]-tbl$lx[j+1L])/tbl$lx[j]
      lp <- .single_life_log_ratio(tbl$lx[j+1L],tbl$lx[j])
      if (lp > log(.5)) lp <- log1p(-q)
      c(base=basis$log_survival(k),p=lp,
        q=.single_life_log_ratio(tbl$lx[j]-tbl$lx[j+1L],tbl$lx[j]))
    }
    list(survival=survival,parameters=parameters)
  }
  ax <- marginal(lt_x,x); by <- marginal(lt_y,y)
  survival <- function(t,frac,status) {
    px <- ax$survival(t,frac); py <- by$survival(t,frac)
    if (is.na(px) || is.na(py)) return(NA_real_)
    if (status == "joint") return(px*py)
    hi <- max(px,py); lo <- min(px,py)
    hi + (1-hi)*lo
  }
  # log(log(1+exp(z))) without intermediate overflow or underflow.
  log_log1pexp <- function(z) {
    if (z < -36) return(z)
    if (z > 36) return(log(z+log1p(exp(-z))))
    log(log1p(exp(z)))
  }
  log_add <- function(a,b) {
    hi <- max(a,b)
    if (hi == -Inf) return(-Inf)
    hi+log1p(exp(min(a,b)-hi))
  }
  cf_area <- function(lp,b) {
    if (lp == 0) return(log(b))
    if (lp == -Inf) return(-Inf)
    if (b*lp == 0) return(log(b))
    log(-expm1(b*lp))-log(-lp)
  }
  balducci_area <- function(lp,lq,b) {
    if (lp == -Inf) return(-Inf)
    if (lq == -Inf) return(log(b))
    lp-lq+log_log1pexp(log(b)+lq-lp)
  }
  balducci_joint <- function(a,c,b) {
    if (a["p"] == -Inf || c["p"] == -Inf) return(-Inf)
    if (a["p"] > c["p"]) { tmp<-a; a<-c; c<-tmp }
    den <- log_add(c["p"],log(b)+c["q"])
    leading <- log(b)+c["p"]-den
    if (a["p"] == c["p"]) return(leading)
    # Integral = b*p2/(p2+q2*b) * log1p(r)/r, p1 <= p2,
    # r = b*(p2-p1)/(p1*(p2+q2*b)). Use logs also when r overflows.
    difference <- c["p"]+log(-expm1(a["p"]-c["p"]))
    logr <- log(b)+difference-a["p"]-den
    leading+log_log1pexp(logr)-logr
  }
  area <- function(k,b,frac,status) {
    a<-ax$parameters(k); c<-by$parameters(k)
    if (frac == "UDD") {
      qx<-exp(a["q"]); qy<-exp(c["q"])
      one<-log(b-qx*b*b/2); two<-log(b-qy*b*b/2)
      joint<-log(b-(qx+qy)*b*b/2+qx*qy*b^3/3)
    } else if (frac == "CF") {
      one<-cf_area(a["p"],b); two<-cf_area(c["p"],b)
      joint<-cf_area(a["p"]+c["p"],b)
    } else {
      one<-balducci_area(a["p"],a["q"],b)
      two<-balducci_area(c["p"],c["q"],b)
      joint<-balducci_joint(a,c,b)
    }
    both<-exp(a["base"]+c["base"]+joint)
    if (status == "joint") return(unname(both))
    unname(exp(a["base"]+one)+exp(c["base"]+two)-both)
  }
  list(survival=survival,area=area)
}
