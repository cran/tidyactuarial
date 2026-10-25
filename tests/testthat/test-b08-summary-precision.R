test_that("summary moments retain representable deviation and standard error", {
  for (a in c(1e-300,1e-160,1e308,1.7e308)) {
    out<-summary_mc(data.frame(present_value=c(-a,a)),probs=.5,var_probs=.5)
    if(a<1.3e308) expect_equal(out$sd/(sqrt(2)*a),1,tolerance=1e-12) else
      expect_true(is.infinite(out$sd))
    expect_equal(out$se_mean/a,1,tolerance=1e-12)
    expect_equal(out$mean,0)
  }
})

test_that("empirical tails avoid overflow and cancellation in excess losses", {
  for(a in c(1e308,1.7e308)) {
    out<-summary_mc(data.frame(present_value=c(-a,a)),probs=.5,var_probs=c(0,.25,.5,1))
    expect_equal(out$TVaR_000,0)
    expect_equal(out$TVaR_250/(a/3),1,tolerance=1e-12)
    expect_equal(out$TVaR_500,a)
    expect_equal(out$TVaR_1000,a)
  }
  out<-summary_mc(data.frame(present_value=c(-1e16,1e16-2)),var_probs=0)
  expect_equal(out$TVaR_000,-1)
  tiny<-summary_mc(data.frame(present_value=c(-1e16,1e16-2)),var_probs=1e-20)
  expect_equal(tiny$TVaR_000,(-1+1e-4)/(1-1e-20),tolerance=1e-12)
})

test_that("subnormal tail masses preserve their representable mean", {
  a<-.Machine$double.xmin*.Machine$double.eps
  out<-summary_mc(data.frame(present_value=c(0,a,a,a,a)),var_probs=c(0,.5,1))
  expect_equal(out$TVaR_000,a)
  expect_equal(out$TVaR_500,a)
  expect_equal(out$TVaR_1000,a)
})

test_that("TVaR matches independent probability-mass integration", {
  samples<-list(c(-8,-2,0,5,13),c(0,100,100,200),c(-10,-10,-10,-2),rep(40,5))
  ps<-c(0,.1,.25,.5,.8,1)
  for(x in samples) {
    out<-summary_mc(data.frame(present_value=x),probs=c(0,.5,1),var_probs=ps)
    ordered<-sort(x)
    expected<-vapply(ps,function(p) {
      if(p==1) return(max(x))
      right<-seq_along(x)/length(x);left<-(seq_along(x)-1)/length(x)
      mass<-pmax(0,right-pmax(left,p))
      sum(ordered*mass)/(1-p)
    },0)
    expect_equal(as.numeric(unlist(out[1,paste0("TVaR_",sprintf("%03.0f",ps*1000))],
                                  use.names=FALSE)),expected,
                 tolerance=1e-12)
    expect_equal(out$q500,stats::median(x))
  }
})

test_that("extended tails do not subtract infinity from itself", {
  for(x in list(c(0,Inf),c(Inf,Inf),c(-Inf,0))) {
    out<-summary_mc(data.frame(present_value=x),var_probs=c(.5,1))
    expect_equal(out$TVaR_500,max(x))
    expect_equal(out$TVaR_1000,max(x))
    expect_true(is.na(out$variance))
    expect_true(is.na(out$se_mean))
  }
})

test_that("grouping, NA policies and empty samples retain their schemas", {
  d<-data.frame(g=c("a","a",NA,NA),present_value=c(1,NA,3,5))
  out<-summary_mc(d,by="g",probs=numeric(),var_probs=numeric())
  expect_equal(out$n_sim,c(1L,2L))
  expect_equal(out$mean,c(1,4))
  expect_equal(out,summary_mc(dplyr::group_by(d,g),probs=numeric(),var_probs=numeric()))
  kept<-summary_mc(d,by="g",na_rm=FALSE)
  expect_true(is.na(kept$mean[1]))
  expect_equal(kept$n_sim,c(2L,2L))
  empty<-summary_mc(data.frame(present_value=numeric()))
  expect_equal(empty$n_sim,0L)
  expect_true(is.na(empty$mean))
  one<-summary_mc(data.frame(present_value=7))
  expect_true(is.na(one$sd))
  expect_equal(one$mean,7)
  empty_groups<-summary_mc(data.frame(g=character(),present_value=numeric()),by="g")
  expect_equal(nrow(empty_groups),0L)
  expect_true(all(c("g","n_sim","mean","TVaR_950") %in% names(empty_groups)))
})

test_that("iid error is descriptive and does not silently infer antithetic pairs", {
  d<-data.frame(present_value=c(1,9,2,8,4,6),method="antithetic")
  raw<-summary_mc(d)
  expect_equal(raw$se_mean,sd(d$present_value)/sqrt(6))
  pairs<-data.frame(present_value=c(5,5,5))
  paired<-summary_mc(pairs)
  expect_equal(paired$mean,5)
  expect_equal(paired$se_mean,0)
  expect_gt(raw$se_mean,0)
  d$method<-"inverse"
  expect_equal(summary_mc(d),raw)
  set.seed(37);before<-.Random.seed
  invisible(summary_mc(d))
  expect_identical(.Random.seed,before)
})

test_that("B08 status products premiums losses reserves and summaries reconcile", {
  lives<-data.frame(sim_id=rep(1:3,each=2),life_id=rep(c("x","y"),3),
                     Tx=c(.3,2.7,1.2,.8,3.2,4.1))
  lives$Kx<-floor(lives$Tx)
  for(status in c("joint","last_survivor")) {
    d<-mc_multilife_status(lives,status=status)
    T<-if(status=="joint") c(.3,.8,3.2) else c(2.7,1.2,4.1)
    expect_equal(d$T_status,T)
    z<-100*1.04^(-(floor(T)+1))
    factors<-vapply(T,function(t) sum(1.04^(-(0:20/2)[(0:20/2)<t])),0)
    quota<-sum(z)/sum(factors)
    d<-mc_insurance(d,i=.04,benefit=100,col_K="K_status",col_T="T_status")
    d<-mc_annuity(d,i=.04,payment=.5,k=2,timing="due",col_K="K_status",col_T="T_status")
    expect_equal(d$pv_benefit,z,tolerance=1e-12)
    expect_equal(d$pv_annuity,.5*factors,tolerance=1e-12)
    d<-mc_premium(d)
    expect_equal(d$premium_per_payment,rep(quota,3),tolerance=1e-12)
    loss<-mc_loss(d)
    expected<-z-quota*factors
    expect_equal(loss$L,expected,tolerance=1e-12)
    reserve<-mc_reserve(d,t=0,i=.04,benefit=100,col_K="K_status",col_T="T_status")
    expect_equal(reserve$L_t,expected,tolerance=1e-12)
    s<-summary_mc(loss,col_value="L")
    expect_equal(s$mean,mean(expected),tolerance=1e-12)
    expect_equal(s$variance,sum((expected-mean(expected))^2)/2,tolerance=1e-12)
    expect_equal(s$n_sim,3L)
  }
})

test_that("B08 wrappers retain a common finite-horizon observed cash-flow book", {
  lt<-data.frame(x=60:63,lx=c(100,70,20,0))
  a<-simulate_annuity_x(lt,60,i=.04,n=2,k=1,payment=.5,n_sim=9,seed=13)
  z<-simulate_insurance_x(lt,60,i=.04,n=2,type="endowment",benefit=100,n_sim=9,seed=13)
  expect_identical(a$Kx,z$Kx)
  expect_identical(a$died_within_horizon,z$died_within_horizon)
  expected<-vapply(seq_len(nrow(a)),function(j) {
    dates<-c(0,1)
    dates<-if(a$died_within_horizon[j]) dates[dates<=a$Kx[j]] else dates
    sum(.5*1.04^(-dates))
  },0)
  expect_equal(a$present_value,expected,tolerance=1e-12)
  date<-ifelse(z$died_within_horizon,z$Kx+1,2)
  expect_equal(z$present_value,100*1.04^(-date),tolerance=1e-12)
})

test_that("B08 consumers accept only coherent overflowing annual basis metadata", {
  d<-mc_annuity(data.frame(Kx=1,Tx=1.5),i=1000,i_type="force",payment=1e308,k=2,
                 timing="immediate")
  expect_true(is.finite(d$pv_annuity))
  expect_true(is.infinite(d$payment_annualized))
  d$pv_benefit<-1e90
  priced<-mc_premium(d)
  expected_quota<-1e90/(exp(-500)+exp(-1000))
  expect_equal(priced$premium_per_payment/expected_quota,1,tolerance=1e-12)
  loss<-mc_loss(priced)
  expect_equal(loss$pv_premiums/1e90,1,tolerance=1e-12)
  expect_lt(abs(loss$L)/1e90,1e-12)
  reserve<-mc_reserve(priced,i=1000,i_type="force",benefit=0,premium_timing="immediate")
  expect_equal(reserve$future_pv_premiums/1e90,1,tolerance=1e-12)
  bad<-data.frame(pv_benefit=1,pv_annuity=1,Kx=1,Tx=1.5,
                   payment_per_payment=1,payment_annualized=Inf,payments_per_year=1)
  expect_error(mc_premium(bad),"payment_annualized|inconsistent")
  expect_error(mc_loss(bad,P=1),"payment_annualized|inconsistent")
  expect_error(mc_reserve(bad,i=0,P=1),"payment_annualized|inconsistent")
  only_annual<-data.frame(pv_benefit=1,pv_annuity=1,payment_annualized=Inf)
  expect_error(mc_premium(only_annual),"finite")
  expect_error(mc_loss(only_annual,P=1),"finite")
  expect_error(mc_reserve(transform(only_annual,Kx=1),i=0,P=1),"finite")
})
