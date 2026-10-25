test_that("two-life fractional survival preserves small observed annual survival", {
  for (p in c(1e-20, 1e-100, 1e-200, 1e-300)) {
    lt <- data.frame(x=60:61,lx=c(1,p))
    # Independent closed forms, compared relatively even near underflow.
    for (u in c(.25,.5,.75)) {
      for (frac in c("CF","Balducci")) {
        s <- if(frac=="CF") p^u else p/(p+(1-p)*u)
        expect_equal(t_pxy(lt,60,60,u,frac=frac,status="last")/(2*s-s*s),1)
        if (s*s > 0)
          expect_equal(t_pxy(lt,60,60,u,frac=frac)/(s*s),1)
      }
    }
  }
  # A marginal ratio underflows, but its fractional survival need not.
  lt <- data.frame(x=60:61,lx=c(1e300,1e-300))
  expect_equal(t_pxy(lt,60,60,.25,frac="CF")/1e-300,1,tolerance=1e-12)
})

test_that("complete two-life expectations resolve steep CF and Balducci survival", {
  for (p in c(.8,.01,1e-20,1e-100,1e-200,1e-300)) {
    lt <- data.frame(x=60:61,lx=c(1,p))
    for (b in c(.25,1)) {
      mu <- -log(p)
      cf_one <- -expm1(-mu*b)/mu
      cf_both <- -expm1(-2*mu*b)/(2*mu)
      # Integral of [p/(p+(1-p)u)]^2, not numerical quadrature.
      ba_both <- b*p/(p+(1-p)*b)
      ba_one <- p/(1-p)*log1p((1-p)*b/p)
      for (cohort in c("first","last")) {
        cf <- if(cohort=="first") cf_both else 2*cf_one-cf_both
        ba <- if(cohort=="first") ba_both else 2*ba_one-ba_both
        expect_equal(e_xy(lt,60,60,t=b,type="complete",frac="CF",cohort=cohort)/cf,
                     1,tolerance=1e-11)
        expect_equal(e_xy(lt,60,60,t=b,type="complete",frac="Balducci",cohort=cohort)/ba,
                     1,tolerance=1e-11)
      }
    }
  }
})

test_that("annual area formulas handle unequal, identical and certain lives", {
  for (p in c(0,.2,.2+1e-15,.9,1)) {
    lt <- data.frame(x=60:62,lx=c(1,.2,.2*p))
    for(frac in c("UDD","CF","Balducci")) {
      marginal <- function(u,p) switch(frac,
        UDD=1-(1-p)*u, CF=if(p==0) rep(0,length(u)) else p^u,
        Balducci=if(p==0) rep(0,length(u)) else p/(p+(1-p)*u))
      both <- integrate(function(u) marginal(u,.2)*marginal(u,p),0,.7,
                        rel.tol=1e-12)$value
      last <- integrate(function(u) {
        a<-marginal(u,.2); b<-marginal(u,p); a+b-a*b
      },0,.7,rel.tol=1e-12)$value
      expect_equal(e_xy(lt,60,61,t=.7,type="complete",frac=frac),both,tolerance=1e-11)
      expect_equal(e_xy(lt,60,61,t=.7,type="complete",frac=frac,cohort="last"),last,
                   tolerance=1e-11)
    }
  }
  certain<-data.frame(x=60:62,lx=c(1,1,1))
  for(frac in c("UDD","CF","Balducci"))
    expect_equal(e_xy(certain,60,60,t=1.5,type="complete",frac=frac),1.5)
})

test_that("last-survivor probability stays in the unit interval", {
  p <- c(.999999999999925,.9999999999999998)
  lt <- lapply(p,function(z) data.frame(x=60:61,lx=c(1,z)))
  actual <- t_pxy(lt,60,60,1,status="last")
  expect_lte(actual,1)
  expect_gte(actual,max(p))
})

test_that("the six two-life functions reconcile with independent state ledgers", {
  counts<-c(100,80,40,10,0)
  lt<-data.frame(x=60:64,lx=counts)
  for(status in c("joint","last")) for(rate in c(-.02,0,.05)) {
    sx<-counts/100; sy<-c(counts[-1]/80,0)
    s<-if(status=="joint") sx*sy else sx+sy-sx*sy
    cohort<-if(status=="joint") "first" else "last"
    death<-sum((s[1:3]-s[2:4])*(1+rate)^(-(1:3)))
    a<-sum(s[1:3]*(1+rate)^(-(0:2)))
    for(j in 0:3)
      expect_equal(t_pxy(lt,60,61,j,status=status),s[j+1])
    expect_equal(e_xy(lt,60,61,t=3,cohort=cohort),sum(s[2:4]))
    expect_equal(annuity_xy(lt,60,61,rate,n=3,status=status,timing="due"),a)
    expect_equal(insurance_xy(lt,60,61,rate,type="term",n=3,status=status),death)
    expect_equal(premium_xy(lt,60,61,rate,type="term",n=3,status=status),death/a)
    expect_equal(unname(reserve_xy(lt,60,61,rate,type="term",n=3,status=status,
                                  P=death/a,t=0,output="value")),0,tolerance=1e-12)
    # Annual end-of-year death versus due annuity, including the maturity term.
    expect_equal(death,1-rate/(1+rate)*a-(1+rate)^(-3)*s[4])
  }
})

test_that("truncation and nonrepresentable premium factors remain explicit limits", {
  open<-data.frame(x=60:62,lx=c(100,80,40))
  expect_equal(e_xy(open,60,60,t=10),.8^2+.4^2)
  expect_equal(t_pxy(open,60,60,3),0)
  expect_error(insurance_xy(open,60,60,0,type="whole"),"extinction|horizon|whole")
  closed<-data.frame(x=60:62,lx=c(1,.5,0))
  expect_error(premium_xy(closed,60,60,1000,i_type="force",type="term",n=1,
                          timing="immediate"),"annuity|positive|finite")
})
