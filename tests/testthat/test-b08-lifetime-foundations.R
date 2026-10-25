test_that("positive annual survival retains the conditional fractional death law", {
  uniforms <- function(n,method) {
    if(method=="antithetic") {
      u<-runif(ceiling(n/2))
      as.vector(rbind(u,1-u))[seq_len(n)]
    } else runif(n)
  }
  n<-31L
  for(method in c("inverse","multinomial","antithetic")) {
    set.seed(42)
    if(method=="multinomial") invisible(sample.int(2,n,replace=TRUE,prob=c(1-1e-11,1e-11)))
    else invisible(uniforms(n,method))
    u<-uniforms(n,method)
    for(frac in c("cml","balducci")) {
      tab<-data.frame(x=60:61,qx=c(1-1e-11,1))
      out<-simulate_lifetime(tab,60,n_sim=n,frac=frac,method=method,seed=42)
      expect_true(all(out$Kx==0))
      q<-tab$qx[1]; p<-1-q
      # Evaluate the conditional CDF, independently of the inverse used to draw.
      observed_u<-if(frac=="cml") -expm1(log(p)*out$Tx)/q else
        out$Tx/(p+q*out$Tx)
      expect_equal(observed_u,u,tolerance=1e-12)
      coarse<-simulate_lifetime(tab,60,n_sim=n,frac=frac,method=method,seed=42,tol=1e-5)
      expect_equal(coarse$Tx,out$Tx,tolerance=1e-12)
    }
  }
  # Certain death retains the established zero fraction for CML/Balducci.
  for(frac in c("cml","balducci"))
    expect_equal(simulate_lifetime(data.frame(x=60,qx=1),60,n_sim=3,frac=frac,seed=1)$Tx,
                 rep(0,3))
})

test_that("simulation from survivor counts is invariant to the radix", {
  for(method in c("inverse","multinomial","antithetic")) for(frac in c("udd","cml","balducci")) {
    tab<-data.frame(x=60:63,lx=c(100,80,40,0))
    base<-simulate_lifetimes(tab,c(60,61),n_sim=25,method=method,frac=frac,seed=73)
    for(scale in c(1e-100,1e-12,1e100)) {
      scaled<-transform(tab,lx=lx*scale)
      out<-simulate_lifetimes(scaled,c(60,61),n_sim=25,method=method,frac=frac,seed=73)
      expect_equal(out,base,tolerance=1e-12)
    }
  }
  # Direct counts preserve a representable decrement which 1-l_next/l loses.
  tab<-data.frame(x=60:61,lx=c(1e16,1e16-2))
  expect_warning(out<-simulate_lifetimes(tab,60,n_sim=3,seed=1),"conditional")
  expect_equal(out$qx_at_death_year,rep(2e-16,3),tolerance=0)
})

test_that("simulation validates tolerances before using them", {
  tab<-data.frame(x=60:61,qx=c(.5,1))
  for(tol in list(NA_real_,numeric(),c(1e-10,1e-8),Inf,-1)) {
    expect_error(simulate_lifetime(tab,60,n_sim=2,tol=tol),"tol")
    expect_error(simulate_lifetimes(tab,60,n_sim=2,tol=tol),"tol")
  }
})

test_that("simulation validates integer ages before coercion", {
  tab<-data.frame(x=60:61,qx=c(.5,1))
  for(age in c(3e9,.Machine$integer.max+1)) {
    expect_error(simulate_lifetime(tab,age,n_sim=2),"integer|range")
    expect_error(simulate_lifetimes(tab,age,n_sim=2),"integer|range")
    bad<-data.frame(x=age,qx=1)
    expect_error(simulate_lifetime(bad,60,n_sim=2),"integer|range")
    expect_error(simulate_lifetimes(bad,60,n_sim=2),"integer|range")
  }
})

test_that("death ages retain values beyond integer storage", {
  # death_age may exceed integer storage even though the input age is valid.
  edge<-data.frame(x=.Machine$integer.max,qx=1)
  out<-simulate_lifetime(edge,.Machine$integer.max,n_sim=2,seed=1)
  expect_equal(out$death_age,rep(as.double(.Machine$integer.max)+1,2))
  expect_false(anyNA(out$death_age))
})

test_that("rounded simulation counts cannot become zero", {
  tab<-data.frame(x=60:61,qx=c(.5,1))
  expect_error(simulate_lifetime(tab,60,n_sim=1e-12),"positive|integer")
  expect_error(simulate_lifetimes(tab,60,n_sim=1e-12),"positive|integer")
})

test_that("seeded generation restores RNG even when a later life fails", {
  tab<-data.frame(x=60:61,qx=c(.5,0))
  set.seed(827)
  state<-.Random.seed
  expect_error(simulate_lifetimes(tab,c(60,61),n_sim=5,seed=3),"zero|probability")
  expect_identical(.Random.seed,state)
  expect_error(simulate_lifetime(tab,60,n_sim=5,tol=NA_real_,seed=3),"tol")
  expect_identical(.Random.seed,state)
})

test_that("distribution diagnostics agree with independent annual death masses", {
  q<-c(.2,.5,1)
  tab<-data.frame(x=60:62,qx=q)
  for(method in c("inverse","multinomial","antithetic")) {
    out<-simulate_lifetime(tab,60,n_sim=7,method=method,frac="udd",seed=12,
                           include_distribution=TRUE)
    dist<-out$distribution[[1]]
    expect_equal(dist$survival_to_k,c(1,.8,.4))
    expect_equal(dist$prob_unconditional,c(.2,.4,.4))
    expect_equal(dist$prob,c(.2,.4,.4))
    expect_equal(dist$cdf,c(.2,.6,1))
    expect_true(all(out$Kx==floor(out$Tx)))
  }
})

test_that("seeded streams restore absence of RNG state and unseeded calls advance it", {
  check<-function() {
    had<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
    old<-if(had) get(".Random.seed",envir=.GlobalEnv)
    on.exit({
      if(had) assign(".Random.seed",old,envir=.GlobalEnv)
      else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
        rm(list=".Random.seed",envir=.GlobalEnv)
    })
    tab<-data.frame(x=60:61,qx=c(.5,1))
    if(had) rm(list=".Random.seed",envir=.GlobalEnv)
    invisible(simulate_lifetime(tab,60,n_sim=5,seed=1))
    expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
    invisible(simulate_lifetimes(tab,c(60,60),n_sim=5,seed=1))
    expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
    set.seed(321)
    invisible(runif(10)) # Five inverse draws and five UDD fractions.
    expected<-get(".Random.seed",envir=.GlobalEnv)
    set.seed(321)
    invisible(simulate_lifetime(tab,60,n_sim=5))
    expect_identical(get(".Random.seed",envir=.GlobalEnv),expected)
    set.seed(321)
    invisible(runif(20))
    expected<-get(".Random.seed",envir=.GlobalEnv)
    set.seed(321)
    invisible(simulate_lifetimes(tab,c(60,60),n_sim=5))
    expect_identical(get(".Random.seed",envir=.GlobalEnv),expected)
  }
  check()
})
