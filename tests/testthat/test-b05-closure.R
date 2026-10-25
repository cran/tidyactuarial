test_that("pure endowments preserve discounted survival instead of a probability cap", {
  lt <- lifetable(0:2,lx=c(100,100,100))
  expect_equal(t_Ex(lt,0,1,i=-.5),2)
  expect_equal(t_Ex(lt,0,1e20,i=1e-20),0) # outside the documented horizon
  # Use a valid annual table for original-force products, with binary oracles.
  rare <- lifetable(0:2,lx=c(2^1000,1,0))
  expect_equal(t_Ex(rare,0,1,i=-1100*log(2),i_type="force") / 2^100,c(1),tolerance=1e-12)
  expect_equal(t_Ex(lt,0,0,i=1000,i_type="force"),1)
  expect_no_warning(t_Ex(lt,0,1,i=.05,m=2^31))
})

test_that("expectancy integrals retain small horizons and extinction limits", {
  tiny <- lifetable(0:2,qx=c(1e-12,0,1),radix=1)
  for(f in c("UDD","CF","Balducci")) {
    expect_equal(e_x(tiny,0,t=1e-8,type="complete",frac=f)/1e-8,1,tolerance=1e-10)
    dead <- lifetable(0:3,lx=c(100,0,0,0))
    expect_equal(e_x(dead,0,type="complete",frac=f),if(f=="UDD") .5 else 0)
    expect_equal(t_qx(tiny,0,t=1e-8,frac=f)/1e-20,1,tolerance=1e-10)
  }
  for(f in c("CF","Balducci"))
    expect_equal(e_x(tiny,0,t=1e-310,type="complete",frac=f)/1e-310,1,tolerance=1e-10)
  expect_equal(t_qx(tiny,1e-12,1e-8,frac="CF")/1e-20,1,tolerance=1e-10)
  irregular <- data.frame(x=c(0,2,4),lx=c(100,80,60))
  expect_error(t_px(irregular,0,1),"annual")
  expect_error(e_x(irregular,0),"annual")
})

test_that("commutation retains original log discounts and annual assumptions", {
  out <- commutation_table(data.frame(x=0:1,lx=c(2^1000,0)),
    i=1100*log(2),i_type="force")
  expect_equal(out$Cx[1]/2^-100,1,tolerance=1e-12)
  expect_equal(out$Dx,c(2^1000,0))
  expect_equal(out$Mx, c(2^-100,0), tolerance=1e-12)
  expect_error(commutation_table(data.frame(x=c(0,2),lx=c(100,50)),i=.05),"consecutive")
  expect_error(commutation_table(data.frame(x=0:1,lx=c(100,110)),i=.05),"nonincreasing")
})

test_that("KM exposure and lived-time ax respect radix and interval width", {
  for(r in c(1e-250,1000,1e250)) for(a in c("UDD","CF","Balducci")) {
    tab <- km_lifetable(c(.5,3),c(1,0),breaks=c(0,2),radix=r,assumption=a)$lifetable
    # Integrate survival from 1 to .5 over width 2, independent of monetary scaling.
    sf <- switch(a, UDD=function(u)1-u/4, CF=function(u)exp(-log(2)*u/2),
      Balducci=function(u)1/(1+u/2))
    exposure <- integrate(sf,0,2,rel.tol=1e-12)$value
    expect_equal(tab$Lx/r,exposure,tolerance=1e-12)
    expect_equal(tab$ax,(exposure-1)/.5,tolerance=1e-12)
    expect_equal(tab$mx,.5/exposure,tolerance=1e-12)
  }
  for(a in c("CF","Balducci")) {
    tab <- km_lifetable(c(.5,.5),c(1,1),breaks=0:1,assumption=a)$lifetable
    expect_equal(tab$Lx,0)
    expect_identical(tab$mx,Inf)
  }
  expect_error(km_lifetable(c(1,2),c(1,0),radix=Inf),"radix")
  expect_error(km_lifetable(c(1,2),c(1,0),conf_level=1),"conf_level")
  expect_error(km_lifetable(numeric(),numeric()),"nonempty")
})

test_that("decrement tables reject destructive names and retain rare causes", {
  expect_error(md_table(data.frame(x=0:1,q_a=c(.1,1),d_a=c(.2,0))),"collid|reserved")
  expect_error(md_table(data.frame(x=0:1,q_total=c(.2,1))),"reserved")
  expect_error(md_table(data.frame(lx=0:1,q_a=c(.2,1)),age_col="lx"),"reserved")
  expect_error(md_table(data.frame(x=0:1,q_a=c(.1,1),a=c(.1,0),d_a.1=c(.1,0))),"collid")
  expect_error(md_table(data.frame(x=0:1,q_a=c(.2,1)),radix=Inf),"radix")
  md <- md_table(data.frame(x=0:2,q_a=c(1e-12,.1,1)),radix=100)
  expect_equal(t_qxj(md,0,1e-8,cause="q_a",frac="CF")/1e-20,1,tolerance=1e-10)
  expect_error(t_qxj(md,0,Inf,cause="q_a"),"finite")
})

test_that("KM graphics sort step data and validate estimates", {
  p <- plot_km(data.frame(time=c(3,1,2),S=c(.4,.8,.6),
    ci_low=c(.2,.6,.4),ci_high=c(.6,.9,.8)))
  expect_equal(p$data$time,c(1,2,3))
  expect_s3_class(p+ggplot2::theme_bw(),"ggplot")
  expect_error(plot_km(data.frame(time=1:2,S=c(.9,1.1))),"\\[0, ?1\\]")
})

test_that("unknown survivor counts propagate without indexing errors", {
  lt <- data.frame(x=0:2,lx=c(100,NA_real_,0))
  expect_true(is.na(t_px(lt,0,1)))
  expect_true(is.na(t_qx(lt,0,1)))
  expect_true(is.na(e_x(lt,1)))
  expect_true(is.na(e_x(lt,0,type="complete")))
})

test_that("delayed-entry KM confidence limits retain small survival estimates", {
  # One continuing censored life plus one new life at each event: n_j=2,d_j=1.
  n <- 600
  fit <- km_lifetable(c(seq_len(n),n+1),c(rep(1,n),0),
    entry=c(seq_len(n)-.5,0),breaks=c(0,n+1),assumption="CF")
  out <- fit$km
  survival <- 2^-n
  greenwood <- n/2
  z <- stats::qnorm(.975)
  loglog <- log(n*log(2)); se <- sqrt(greenwood)/(n*log(2))
  expect_equal(out$seS[n]/(survival*sqrt(greenwood)),1,tolerance=1e-12)
  expect_equal(out$ci_low[n]/exp(-exp(loglog+z*se)),1,tolerance=1e-10)
  expect_equal(out$ci_high[n]/exp(-exp(loglog-z*se)),1,tolerance=1e-10)
  expect_equal(fit$lifetable$px/survival,1,tolerance=1e-12)
  expect_equal(fit$lifetable$Lx/1e5,(n+1)/(n*log(2)),tolerance=1e-12)
})
