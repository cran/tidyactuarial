test_that("approved C07 leaves unknown open terminal intervals missing", {
  lt <- lifetable(60:62, lx = c(1000,800,600), close = FALSE)
  expect_equal(lt$dx[1:2], c(200,200))
  expect_equal(lt$qx[1:2], c(.2,.25))
  for (name in c("dx","qx","px","mx")) expect_true(is.na(lt[[name]][3]))
  for (probability in list(list(qx=c(.2,.25,.5)),list(px=c(.8,.75,.5)))) {
    out <- do.call(lifetable,c(list(x=60:62,lx=c(1000,800,600),close=FALSE),probability))
    expect_equal(out$dx[3],300)
    expect_equal(out$qx[3],.5)
    expect_equal(out$mx[3],2/3)
  }
  closed <- lifetable(60:62,lx=c(1000,800,600))
  expect_equal(closed$dx[3],600)
  expect_equal(closed$qx[3],1)
})

test_that("contradictory probabilities and forbidden NaN inputs are rejected", {
  expect_error(lifetable(0:2,lx=c(100,90,80),qx=c(.1,1/9,1),px=c(.2,.2,0)),
    "Inconsistent.*px|Inconsistent.*probabil")
  for (name in c("lx","qx","px","mx")) {
    values <- c(NaN,.1,1)
    expect_error(do.call(lifetable,c(list(x=0:2,radix=100),setNames(list(values),name))),
      "Inf/NaN")
  }
  expect_error(lifetable(0:2,lx=c(100,90,80),tol=Inf),"tol")
  expect_error(lifetable(0:2,lx=c(100,90,80),check=NA),"check")
  expect_error(lifetable(0:2,lx=c(100,90,80),close=NA),"close")
  expect_error(lifetable(0:2,lx=c(100,90,80),ax=c(.2,.3)),"ax")
  expect_error(lifetable(0:2,lx=c(100,90,80),radix=c(100,200)),"radix")
  expect_error(lifetable((2^31)+(0:2),lx=c(100,90,80)),"representable|integer range")
  expect_error(mortality_law_table('Exponential',2^31,2^31+1,lambda=.1),"representable|integer range")
  expect_error(lifetable(0:3,lx=c(100,90,80),omega=2),"same length")
  expect_error(lifetable(0:2,lx=c(100,90,80),radix=NaN),"radix")
  expect_error(lifetable(0:2,lx=c(100,90,80),ax=-1e-12),"ax")
  expect_error(mortality_law_table('Exponential',0,2,lambda=.1,tol=NA),"tol")
})

test_that("rounded annual ages and explicit omega share a terminal index", {
  lt <- lifetable((0:3)+1e-12,lx=c(100,90,80,70),omega=2)
  expect_identical(lt$x,0:2)
  expect_equal(lt$qx[3],1)
  expect_equal(lt$dx[3],80)
})

test_that("continuous-law evaluation avoids intermediate power overflow", {
  # At integer ages 1024 and 1025, B*2^x is finite although 2^x overflows.
  for (law in c("Gompertz","Makeham")) {
    pars <- list(B=2^-1000,c=2)
    if(law=="Makeham") pars$A <- .01
    out <- mortality_law_table(law,1024,1025,params=pars,close=FALSE)
    force <- 2^c(24,25) + if(law=="Makeham") .01 else 0
    expect_equal(out$mu_x/force,c(1,1),tolerance=1e-12)
    expect_equal(out$qx,c(1,1))
  }
  # Logistic tends to B/C; evaluate an algebraically independent reciprocal form.
  out <- mortality_law_table("Logistic",1100,1102,A=.01,B=.1,c=2,C=.5,close=FALSE)
  expect_equal(out$mu_x,rep(.2,3),tolerance=1e-12)
  expect_equal(out$qx,rep(-expm1(-.2),3),tolerance=1e-12)
  edge <- mortality_law_table("Logistic",2^31-2,2^31-1,
    A=0,B=.1,c=1e308,C=.5,close=FALSE)
  expect_equal(edge$mu_x/.2,c(1,1),tolerance=1e-12)
  expect_equal(edge$qx/(-expm1(-.2)),c(1,1),tolerance=1e-12)
  for (law in c("Makeham","Logistic")) {
    pars <- list(A=.1,B=0,c=2)
    if(law=="Logistic") pars$C <- 0
    tab <- mortality_law_table(law,1100,1101,params=pars,close=FALSE)
    expect_equal(tab$mu_x,c(.1,.1))
    expect_equal(tab$qx,rep(-expm1(-.1),2))
  }
})

test_that("Weibull preserves small integrated hazards independently of frac", {
  # Integral of mu is bounded by endpoint forces for this decreasing force.
  for(frac in c("CF","UDD","Balducci","CML")) {
    out <- mortality_law_table("Weibull",1e9,1e9+2,shape=1e-8,scale=1e9,
      frac=frac,close=FALSE)
    left <- 1e-17 * (out$x/1e9)^(1e-8-1)
    right <- 1e-17 * ((out$x+1)/1e9)^(1e-8-1)
    hazard <- -log1p(-out$qx)
    expect_true(all(hazard >= right*(1-1e-12)))
    expect_true(all(hazard <= left*(1+1e-12)))
    expect_true(all(out$qx>0))
    expect_equal(out$qx/1e-17,rep(1,3),tolerance=1e-8)
  }
})

test_that("Heligman-Pollard unit childhood base is not zero times infinity", {
  tab <- mortality_law_table("HeligmanPollard",2,3,A=1,B=1e308,C=2,
    D=0,E=1,F_hp=20,G=0,H=1,close=FALSE)
  expect_equal(tab$qx,c(.5,.5))
})
