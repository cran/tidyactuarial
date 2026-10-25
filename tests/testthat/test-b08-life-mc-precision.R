test_that("lx simulations retain finite CML forces after annual ratio underflow", {
  lt<-data.frame(x=60:61,lx=c(1e300,1e-300))
  n_sim<-11L
  set.seed(42)
  invisible(sample.int(2,n_sim,replace=TRUE,prob=c(1,0)))
  u<-runif(n_sim)
  expected<--log1p(-u)/(600*log(10))
  out<-simulate_annuity_x(lt,60,i=0,n=1,n_sim=n_sim,k=1000,frac="cml",seed=42)
  expect_equal(out$Tx,expected,tolerance=1e-12)
  expect_true(all(out$Tx>0))
  dates<-(0:999)/1000
  count<-vapply(expected,function(t) sum(dates<t),0L)
  expect_equal(out$n_payments,count)
  expect_equal(out$present_value,as.double(count))
})

test_that("fractional death CDFs and RNG streams agree with independent laws", {
  n_sim<-17L
  for(p in c(.8,1e-11,1e-200)) for(frac in c("udd","cml","balducci")) {
    lt<-data.frame(x=60:61,lx=c(1,p))
    # Retain exactly which paths die, without calling a package engine.
    set.seed(17)
    index<-sample.int(2,n_sim,replace=TRUE,prob=c(1-p,p))
    u<-runif(sum(index==1L))
    out<-simulate_annuity_x(lt,60,i=0,n=1,n_sim=n_sim,k=2,frac=frac,seed=17)
    tt<-out$Tx[index==1L]
    cdf<-switch(frac,udd=tt,cml=-expm1(log(p)*tt)/(1-p),
                balducci=tt/(p+(1-p)*tt))
    expect_equal(cdf,u,tolerance=1e-11)
    expect_identical(out$died_within_horizon,index==1L)
    expect_true(all(is.na(out$Tx[index==2L])))
  }
})

test_that("simulated cash flows retain original rates and monetary range", {
  dies<-data.frame(x=60:61,lx=c(1,0))
  for(force in c(-1000,720)) {
    benefit<-if(force<0) 1e-300 else 1e300
    out<-simulate_insurance_x(dies,60,i=force,i_type="force",benefit=benefit,
                              n_sim=3,seed=1)
    expected<-exp(log(benefit)-force)
    expect_equal(out$present_value/expected,rep(1,3),tolerance=1e-12)
    expect_equal(out$benefit_time,rep(1,3))
  }
  alive<-data.frame(x=60:62,lx=c(1,1,0))
  for(timing in c("due","immediate")) {
    out<-simulate_annuity_x(alive,60,i=720,i_type="force",n=1,k=2,
                            timing=timing,payment=1e300,n_sim=3,seed=1)
    dates<-if(timing=="due") c(0,.5) else c(.5,1)
    expected<-sum(exp(log(1e300)-720*dates))
    expect_equal(out$present_value/expected,rep(1,3),tolerance=1e-12)
  }
  zero<-simulate_annuity_x(alive,60,i=-1000,i_type="force",payment=0,n_sim=3,seed=1)
  expect_identical(zero$present_value,rep(0,3))
  expect_identical(simulate_insurance_x(dies,60,i=-1000,i_type="force",benefit=0,
                                       n_sim=3,seed=1)$present_value,rep(0,3))
  # Even an unrepresentable annual log growth cannot alter a sole payment at 0.
  edge<-simulate_annuity_x(alive,60,i=-1e308*(1-1e-15),i_type="nominal_interest",
                           m=1e308,n=1,k=1,payment=7,n_sim=2,seed=1)
  expect_equal(edge$present_value,c(7,7))
})

test_that("small per-period interest is not replaced by zero", {
  lt<-data.frame(x=60:62,lx=c(1,1,0))
  out<-simulate_annuity_x(lt,60,i=1e-12,n=1,k=4,payment=1e12,n_sim=2,seed=1)
  expected<-sum(1e12*exp(-log1p(1e-12)*(0:3)/4))
  expect_lt(max(abs(out$present_value-expected)),.02)
})

test_that("ordinary simulated paths match explicit payment and benefit ledgers", {
  lt<-data.frame(x=60:64,lx=c(100,80,40,10,0))
  for(rate in c(-.02,0,.05)) for(k in c(1,2,4)) for(timing in c("due","immediate")) {
    out<-simulate_annuity_x(lt,60,i=rate,n=3,k=k,timing=timing,payment=7,
                            n_sim=13,seed=42)
    grid<-(if(timing=="due") 0:(3*k-1) else 1:(3*k))/k
    expected<-vapply(seq_len(nrow(out)),function(j) {
      eligible<-if(!out$died_within_horizon[j]) rep(TRUE,length(grid)) else if(k==1)
        if(timing=="due") grid<=out$Kx[j] else grid<=out$Kx[j] else
          if(timing=="due") grid<out$Tx[j] else grid<=out$Tx[j]
      sum(7*(1+rate)^(-grid[eligible]))
    },0)
    expect_equal(out$present_value,expected,tolerance=1e-12)
  }
  for(type in c("whole","term","endowment")) {
    out<-simulate_insurance_x(lt,60,i=.05,type=type,n=3,benefit=10,n_sim=13,seed=42)
    time<-if(type=="endowment") ifelse(out$died_within_horizon,out$Kx+1,3) else
      ifelse(out$died_within_horizon,out$Kx+1,NA_real_)
    expect_equal(out$present_value,ifelse(is.na(time),0,10*1.05^(-time)))
  }
})

for(name in c("simulate_annuity_x","simulate_insurance_x")) {
  cases<-list(list(x=3e9),list(n=3e9),list(n_sim=3e9),list(n_sim=1e-12),list(n=1e-12))
  patterns<-c("integer|range","horizon|table|support","n_sim|range","positive|integer","positive|integer")
  for(j in seq_along(cases)) test_that(paste(name,"validates range case",j), {
    args<-list(lt=data.frame(x=60:62,lx=c(1,1,0)),x=60,i=.05,n_sim=2)
    args[names(cases[[j]])]<-cases[[j]]
    expect_error(do.call(get(name),args),patterns[j])
  })
  test_that(paste(name,"preserves large conversion frequency"), {
    large<-do.call(get(name),list(lt=data.frame(x=60:62,lx=c(1,1,0)),x=60,
                                 i=.05,m=3e9,n_sim=2,seed=1))
    expect_equal(large$m,rep(3e9,2))
  })
}

test_that("lx simulation validates payment frequency and observed annual tables", {
  lt<-data.frame(x=60:62,lx=c(1,1,0))
  expect_error(simulate_annuity_x(lt,60,i=.05,k=3e9,n_sim=2),"k|range")
  invalid<-data.frame(x=60:62,lx=c(1e-100,1e-100*(1+1e-12),0))
  expect_error(simulate_annuity_x(invalid,60,i=0,n_sim=2),"nonincreasing")
  duplicated<-data.frame(x=c(60,60+5e-11,61),lx=c(100,80,0))
  expect_error(simulate_annuity_x(duplicated,60,i=0,n=1,n_sim=2),"duplicated")
  bad_age<-data.frame(x=c(60,3e9),lx=c(1,0))
  expect_error(simulate_insurance_x(bad_age,60,i=0,n=1,n_sim=2),"integer")
})

test_that("seeded lx simulations preserve presence and absence of RNG state", {
  lt<-data.frame(x=60:62,lx=c(1,1,0))
  check<-function() {
    had<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
    old<-if(had) get(".Random.seed",envir=.GlobalEnv)
    on.exit({if(had) assign(".Random.seed",old,envir=.GlobalEnv) else
      if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(list=".Random.seed",envir=.GlobalEnv)})
    for(fun in list(simulate_annuity_x,simulate_insurance_x)) {
      set.seed(73); before<-get(".Random.seed",envir=.GlobalEnv)
      invisible(fun(lt,60,i=.05,n_sim=3,seed=1))
      expect_identical(get(".Random.seed",envir=.GlobalEnv),before)
      rm(list=".Random.seed",envir=.GlobalEnv)
      invisible(fun(lt,60,i=.05,n_sim=3,seed=1))
      expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
    }
  }
  check()
})

test_that("lx scaling and unseeded random streams retain their contracts", {
  lt<-data.frame(x=60:62,lx=c(100,80,0))
  for(frac in c("udd","cml","balducci")) {
    base<-simulate_annuity_x(lt,60,i=.05,n=1,k=4,frac=frac,n_sim=17,seed=1)
    for(scale in c(1e-100,1e100)) {
      out<-simulate_annuity_x(transform(lt,lx=lx*scale),60,i=.05,n=1,k=4,
                              frac=frac,n_sim=17,seed=1)
      expect_equal(out,base,tolerance=1e-12)
    }
  }
  set.seed(23)
  index<-sample.int(2,17,replace=TRUE,prob=c(.2,.8))
  invisible(runif(sum(index==1)))
  expected<-.Random.seed
  set.seed(23)
  invisible(simulate_annuity_x(lt,60,i=.05,n=1,k=4,n_sim=17))
  expect_identical(.Random.seed,expected)
  set.seed(23)
  invisible(sample.int(2,17,replace=TRUE,prob=c(.2,.8)))
  expected<-.Random.seed
  set.seed(23)
  invisible(simulate_insurance_x(lt,60,i=.05,n=1,n_sim=17))
  expect_identical(.Random.seed,expected)
})
