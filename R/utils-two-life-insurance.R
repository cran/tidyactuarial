# Disjoint annual transitions from initially living independent lives.
.two_life_insurance_basis <- function(lt_x, lt_y, x, y, status) {
  marginal <- function(lt, age) {
    .two_life_table(lt)
    lt <- lt[order(lt$x), , drop = FALSE]
    j <- match(age, lt$x)
    if (is.na(j) || lt$lx[j] <= 0)
      stop("Issue ages must have positive `lx`.", call. = FALSE)
    counts <- lt$lx[seq.int(j, nrow(lt))]
    horizon <- length(counts) - 1L
    extinct <- which(counts == 0)
    extinct <- if (length(extinct)) extinct[1L] - 1L else Inf
    survival <- function(t) {
      if (t > horizon) {
        if (extinct == Inf) stop("Cannot extrapolate an open marginal table.", call. = FALSE)
        return(-Inf)
      }
      .single_life_log_ratio(counts[t + 1L], counts[1L])
    }
    dead <- function(t) {
      if (t > horizon) {
        if (extinct == Inf) stop("Cannot extrapolate an open marginal table.", call. = FALSE)
        return(0)
      }
      .single_life_log_ratio(counts[1L] - counts[t + 1L], counts[1L])
    }
    death <- function(t) {
      if (t >= horizon) {
        if (counts[length(counts)] > 0)
          stop("Cannot infer death beyond an open marginal table horizon.", call. = FALSE)
        return(-Inf)
      }
      .single_life_log_ratio(counts[t + 1L] - counts[t + 2L], counts[1L])
    }
    list(s = survival, dead = dead, death = death, horizon = horizon,
         extinct = extinct,
         closed = counts[length(counts)] == 0)
  }
  a <- marginal(lt_x, x)
  b <- marginal(lt_y, y)
  add <- function(z) {
    hi <- max(z)
    if (hi == -Inf) return(-Inf)
    hi + log(sum(exp(z - hi)))
  }
  survival <- function(t) {
    if (status == "joint" && t >= min(a$extinct, b$extinct)) return(-Inf)
    if (status == "joint") return(a$s(t) + b$s(t))
    add(c(a$s(t) + b$s(t), a$s(t) + b$dead(t), a$dead(t) + b$s(t)))
  }
  death <- function(t) {
    if (status == "joint" && t >= min(a$extinct, b$extinct)) return(-Inf)
    dx <- a$death(t)
    dy <- b$death(t)
    if (status == "joint")
      return(add(c(dx + b$s(t), a$s(t + 1) + dy)))
    add(c(dx + b$dead(t), a$dead(t) + dy, dx + dy))
  }
  list(survival = survival, death = death,
       complete = if (status == "last") a$closed && b$closed else {
         end <- min(a$horizon, b$horizon)
         a$s(end) == -Inf || b$s(end) == -Inf
       },
       horizon = if (status == "joint") min(a$horizon, b$horizon) else max(a$horizon, b$horizon))
}

# A pure endowment is a survival payment, not the subtraction of two APVs.
.two_life_pure_endowment <- function(lt,x,y,n,h,benefit,i,i_type,m,status) {
  tables<-if(is.data.frame(lt)) list(lt,lt) else lt
  if (!is.list(tables) || length(tables)!=2L ||
      !all(vapply(tables,is.data.frame,logical(1))))
    stop("`lt` must be one life table or a list of two.",call.=FALSE)
  state<-.two_life_insurance_basis(tables[[1]],tables[[2]],x,y,status)
  start<-state$survival(h)
  if (start == -Inf) return(0)
  if (!is.finite(h+n) || h+n>state$horizon)
    stop("Pure-endowment maturity exceeds the supported status horizon.",call.=FALSE)
  survival<-state$survival(h+n)
  basis<-.interest_basis(i_type,i,m)
  if (survival == -Inf) return(0)
  .cash_flow_value(benefit,survival-basis$growth*(h+n))
}
