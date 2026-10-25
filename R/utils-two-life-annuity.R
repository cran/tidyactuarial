# Fractional state occupation from issue, including direct cumulative deaths.
.two_life_annuity_basis <- function(lt_x, lt_y, x, y, frac) {
  add <- function(z) {
    if (anyNA(z)) return(NA_real_)
    hi <- max(z)
    if (hi == -Inf) return(-Inf)
    hi + log(sum(exp(z-hi)))
  }
  marginal <- function(lt, age) {
    basis <- .life_payment_basis(lt, age, require_lx=TRUE)
    tbl <- basis$table
    origin <- tbl$lx[match(age,tbl$x)]
    horizon <- basis$omega-age
    function(t, strict=FALSE) {
      if (t > horizon) {
        if (strict && tail(tbl$lx,1) > 0) return(c(s=NA_real_,q=NA_real_))
        # Preserve the established truncation convention for exact annuities.
        return(c(s=-Inf,q=0))
      }
      year <- floor(t)
      count <- tbl$lx[match(as.double(age)+year,tbl$x)]
      q <- .single_life_log_ratio(origin-count,origin)
      if (t > year) q <- add(c(q,basis$log_death(year,t,frac)))
      c(s=basis$log_survival(t,frac),q=q)
    }
  }
  a <- marginal(lt_x,x)
  b <- marginal(lt_y,y)
  states <- function(t,strict=FALSE) {
    ax<-a(t,strict); by<-b(t,strict)
    c(both=unname(ax["s"]+by["s"]),
      x_only=unname(ax["s"]+by["q"]),
      y_only=unname(ax["q"]+by["s"]))
  }
  list(states=states,add=add)
}
