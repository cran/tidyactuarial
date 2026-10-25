test_that("C12 exact multi-life annuities match independent payment sums", {
  mu <- c(.03,.08,.12)
  tables <- lapply(mu, function(z) data.frame(x=40:50,lx=1000*exp(-z*(0:10))))
  for (lives in 1:3) for (timing in c("due","immediate")) {
    k <- 4; h <- 1; n <- 3
    tt <- h + (if (timing == "due") 0:(n*k-1) else 1:(n*k))/k
    for (frac in c("UDD","CF","Balducci")) {
      probs <- sapply(mu[seq_len(lives)],function(z) {
        p <- exp(-z); year <- floor(tt); u <- tt-year
        exp(-z*year) * switch(frac, UDD=1-u*(1-p),CF=p^u,
                              Balducci=p/(p+u*(1-p)))
      })
      all <- apply(probs,1,prod)
      any <- 1-apply(1-probs,1,prod)
      for (status in c("first","last","reversionary")) {
        w <- switch(status,first=all,last=any,reversionary=.6*all+.4*any)
        actual <- annuity_multi(tables[seq_len(lives)],rep(40,lives),.04,
                                 i_type="force",n=n,h=h,k=k,timing=timing,
                                 annuity=if(status=="reversionary") status else "cohort",
                                 cohort=if(status=="last") "last" else "first",
                                 alpha=.4,frac=frac)
        expect_equal(actual,sum(w*exp(-.04*tt))/k)
      }
    }
  }
})

test_that("C12 Woolhouse includes both finite and deferred endpoints", {
  mu <- .03; delta <- .04; k <- 12; n <- 3; h <- 2
  lt <- data.frame(x=40:50,lx=1000*exp(-mu*(0:10)))
  F <- function(t) exp(-(mu+delta)*t)
  D <- F(h)-F(h+n)
  annual <- sum(F(h+0:(n-1)))
  first <- annual-(1-1/k)/2*D
  second <- first-(1-1/k^2)/12*(mu+delta)*D
  for(order in c("first","second")) for(timing in c("due","immediate")) {
    ref <- if(order=="first") first else second
    if(timing=="immediate") ref <- ref-D/k
    expect_equal(annuity_multi(lt,40,delta,i_type="force",n=n,h=h,k=k,
                               timing=timing,woolhouse=order),ref)
  }
  # Immortal lives at zero interest must receive precisely one annual unit.
  lt$lx <- 1000
  for(order in c("none","first","second")) for(timing in c("due","immediate"))
    expect_equal(annuity_multi(lt,40,0,n=1,k=2,timing=timing,woolhouse=order),1)
})

test_that("multi annuities validate ages, grids and positive Woolhouse borders", {
  lt <- data.frame(x=40:44,lx=c(100,80,60,40,0))
  for(arg in c("ages","n","h","k")) {
    args <- list(lt=lt,ages=40,i=0,n=1)
    args[[arg]] <- 1.2
    expect_error(do.call(annuity_multi,args),"integer")
  }
  expect_error(annuity_multi(lt,NA_real_,0,n=1),"ages")
  expect_error(annuity_multi(lt,40,0,n=1,k=3e9),"payments")
  expect_error(annuity_multi(lt,40,0,n=4,k=2,woolhouse="second"),"Woolhouse|border|survival")
  expect_error(annuity_multi(lt,40,0,n=1,annuity="reversionary",alpha=Inf),"alpha")
  expect_equal(annuity_multi(lt,40,0,n=0),0)
  expect_equal(annuity_multi(lt,40,.03,n=1,m=3e9),.8/1.03)
})

test_that("multi annuities retain rare last survivors and original rates", {
  lt <- data.frame(x=40:42,lx=c(1e300,1e-100,1e-100))
  ref <- exp(log(2)-400*log(10)+1000)
  expect_equal(annuity_multi(lt,c(40,40),-1000,i_type="force",n=1,
                             cohort="last"),ref)
  lt <- data.frame(x=40:42,lx=c(100,0,0))
  expect_equal(annuity_multi(lt,c(40,40),-1e308,i_type="force",n=1,
                             timing="due"),1)
  lt <- data.frame(x=40:44,lx=c(100,80,60,40,20))
  px <- data.frame(x=40:44,px=c(.8,.75,2/3,.5,NA))
  qx <- data.frame(x=40:44,qx=1-px$px)
  for(frac in c("UDD","CF","Balducci")) for(status in c("first","last")) {
    ref <- annuity_multi(lt,c(40,41),.04,n=2,k=12,cohort=status,frac=frac)
    expect_equal(annuity_multi(px,c(40,41),.04,n=2,k=12,cohort=status,frac=frac),ref)
    expect_equal(annuity_multi(qx,c(40,41),.04,n=2,k=12,cohort=status,frac=frac),ref)
  }
  expect_error(annuity_multi(lt,40,0,n=1,h=c(0,1)),"integer")
  expect_error(annuity_multi(lt,40.2,0,n=1),"integer")
  expect_error(annuity_multi(lt[c(1,3,4),],40,0,n=1),"consecutive")
})
