test_that("pure endowments retain small and underflowed joint survival", {
  lt<-data.frame(x=60:62,lx=c(1e16,2,0))
  out<-premium_xy(lt,60,60,0,type="pure_endowment",n=1,benefit=1e16,output="summary")
  expect_equal(out$apv_benefits/4e-16,1,tolerance=1e-12)
  expect_equal(out$premium_annualized/4e-16,1,tolerance=1e-12)
  for(method in c("prospective","recursive"))
    expect_equal(unname(reserve_xy(lt,60,60,0,type="pure_endowment",n=1,benefit=1e16,
                   P=0,t=0,method=method,output="value"))/4e-16,1,tolerance=1e-12)
  rare<-data.frame(x=60:62,lx=c(1e300,1e-100,0))
  # Each survival ratio is 1e-400; combine with benefit and negative force.
  expected<-exp(log(1e300)+2*(log(1e-100)-log(1e300))+700)
  for(method in c("prospective","recursive"))
    expect_equal(unname(reserve_xy(rare,60,60,-700,i_type="force",type="pure_endowment",
                   n=1,benefit=1e300,P=0,t=0,method=method,output="value"))/expected,1,
                   tolerance=1e-12)
})

test_that("zero premiums do not require an overflowing unit annuity", {
  lt<-data.frame(x=60:81,lx=c(rep(100,21),0))
  expected<-exp(log(1e-300)+800)
  for(method in c("prospective","recursive"))
    expect_equal(unname(reserve_xy(lt,60,60,-40,i_type="force",type="pure_endowment",
                   n=20,benefit=1e-300,P=0,t=0,method=method,output="value"))/expected,1,
                   tolerance=1e-12)
})

test_that("recursive reserves retain original rates and tiny first-second deaths", {
  lt<-data.frame(x=60:63,lx=c(100,80,40,0))
  for(status in c("joint","last")) {
    q<-if(status=="joint") .36 else .04
    expected<-exp(log(1e300)+log(q)-720)
    expect_equal(unname(reserve_xy(lt,60,60,720,i_type="force",type="term",status=status,
                  n=2,benefit=1e300,P=0,t=0,method="recursive",output="value"))/expected,1,
                  tolerance=1e-12)
    # At t=1 the reserve is conditioned on BOTH lives alive, not on an
    # unspecified survivor at the first death.
    q1<-if(status=="joint") .75 else .25
    expected1<-exp(log(1e300)+log(q1)-720)
    expect_equal(unname(reserve_xy(lt,60,60,720,i_type="force",type="term",status=status,
                  n=2,benefit=1e300,P=0,t=1,method="recursive",output="value"))/expected1,1,
                  tolerance=1e-12)
  }
  tiny<-data.frame(x=60:62,lx=c(1e16,1e16-2,0))
  survivor<-data.frame(x=60:62,lx=c(100,100,0))
  expected<-1e16*2e-16
  expect_equal(unname(reserve_xy(list(tiny,survivor),60,60,0,type="term",n=1,
                  benefit=1e16,P=0,t=0,method="recursive",output="value")),expected)
  expect_equal(unname(reserve_xy(tiny,60,60,0,type="term",n=1,status="last",
                  benefit=1e300,P=0,t=0,method="recursive",output="value"))/4e268,1,
                  tolerance=1e-12)
  rare<-data.frame(x=60:62,lx=c(1e16,2,0))
  expect_equal(unname(reserve_xy(rare,60,60,0,type="term",n=2,benefit=1e16,
                  P=0,t=1,method="recursive",output="value")),1e16,tolerance=1e-12)
})

test_that("approved immediate-premium reserves include the original date at t", {
  lx<-c(100,80,40,10,0); ly<-c(200,140,100,30,0)
  lt<-list(data.frame(x=60:64,lx=lx),data.frame(x=65:69,lx=ly))
  survival<-function(counts,time,t) {
    j<-floor(time); u<-time-j
    (counts[j+1]*(1-u)+counts[j+2]*u)/counts[t+1]
  }
  for(status in c("joint","last")) for(k in c(1,2)) for(h in 0:1)
    for(type in c("term","endowment","pure_endowment")) {
      end<-h+2
      dates<-h+(1:(2*k))/k
      expected<-vapply(0:end,function(t) {
        s<-function(u) {
          a<-survival(lx,u,t); b<-survival(ly,u,t)
          if(status=="joint") a*b else a+b-a*b
        }
        r<-if(t<end) seq.int(max(t,h),end-1) else numeric()
        deaths<-if(type=="pure_endowment") 0 else
          sum(vapply(r,function(j) 250*1.05^(-(j+1-t))*(s(j)-s(j+1)),0))
        maturity<-if(type=="term") 0 else 250*1.05^(-(end-t))*s(end)
        premium<-10/k*sum(vapply(dates[dates>=t],function(j) 1.05^(-(j-t))*s(j),0))
        deaths+maturity-premium
      },0)
      actual<-reserve_xy(lt,60,65,.05,type=type,status=status,n=2,h=h,benefit=250,
                         P=10,n_prem=2,k=k,timing="immediate",premium_start="deferred",
                         t=0:end,output="value")
      expect_equal(unname(actual),expected,tolerance=1e-12)
    }
  # A short fractional premium period: there is a payment at t=1, but no
  # payment at the insurance maturity t=2.
  expect_equal(unname(reserve_xy(lt,60,65,0,type="pure_endowment",n=2,benefit=250,
              P=10,k=2,n_prem=1.5,timing="immediate",t=2,output="value")),250)
})

test_that("pure-endowment premiums and conditional reserves match explicit UDD books", {
  lx<-c(100,80,40,10,0); ly<-c(200,140,100,30,0)
  lt<-list(data.frame(x=60:64,lx=lx),data.frame(x=65:69,lx=ly))
  marginal<-function(counts,t,origin=0) {
    j<-floor(t); u<-t-j
    (counts[j+1]*(1-u)+counts[j+2]*u)/counts[origin+1]
  }
  for(status in c("joint","last")) for(k in c(1,2,4)) for(h in 0:1) {
    for(start in c("issue","deferred")) {
      pstart<-if(start=="issue") 0 else h
      n<-2; duration<-h+n
      sx<-lx[duration+1]/100; sy<-ly[duration+1]/200
      s<-if(status=="joint") sx*sy else sx+sy-sx*sy
      B<-250*1.05^(-duration)*s
      times<-pstart+(0:(2*k-1))/k
      ax<-vapply(times,function(t) marginal(lx,t),0)
      ay<-vapply(times,function(t) marginal(ly,t),0)
      ss<-if(status=="joint") ax*ay else ax+ay-ax*ay
      a<-sum(1.05^(-times)*ss)/k
      out<-premium_xy(lt,60,65,.05,type="pure_endowment",status=status,n=n,h=h,
                    benefit=250,k=k,n_prem=2,premium_start=start,output="summary")
      expect_equal(out$apv_benefits,B,tolerance=1e-12)
      expect_equal(out$premium_annualized,B/a,tolerance=1e-12)
      expect_equal(out$premium_per_payment,B/a/k,tolerance=1e-12)
      # P=10 fixed: conditional book before the due premium at t=1.
      t<-1
      cx<-lx[duration+1]/lx[t+1]; cy<-ly[duration+1]/ly[t+1]
      sc<-if(status=="joint") cx*cy else cx+cy-cx*cy
      future<-times[times>=t]
      ax<-vapply(future,function(u) marginal(lx,u,t),0)
      ay<-vapply(future,function(u) marginal(ly,u,t),0)
      ss<-if(status=="joint") ax*ay else ax+ay-ax*ay
      expected<-250*1.05^(-(duration-t))*sc-10*sum(1.05^(-(future-t))*ss)/k
      got<-reserve_xy(lt,60,65,.05,type="pure_endowment",status=status,n=n,h=h,
                     benefit=250,P=10,k=k,n_prem=2,premium_start=start,t=t,output="value")
      expect_equal(unname(got),expected,tolerance=1e-12)
    }
  }
})

test_that("premium and reserve storage reject bad grids without integer coercion", {
  lt<-data.frame(x=60:62,lx=c(100,80,0))
  expected<-.64*exp(-.03)
  expect_equal(premium_xy(lt,60,60,.03,i_type="force",m=3e9,type="pure_endowment",
                         n=1,benefit=1),expected,tolerance=1e-12)
  expect_equal(unname(reserve_xy(lt,60,60,.03,i_type="force",m=3e9,type="pure_endowment",
                         n=1,benefit=1,P=0,t=0,output="value")),expected,tolerance=1e-12)
  expect_error(premium_xy(lt,60,60,0,type="pure_endowment",n=1,n_prem=1e308,k=2),
               "n_prem|integer")
  expect_error(reserve_xy(lt,60,60,0,type="pure_endowment",n=1,n_prem=1e308,k=2,P=0),
               "n_prem|integer")
  expect_error(reserve_xy(lt[c(1,3),],60,60,0,type="pure_endowment",n=2,P=0,t=2),
               "annual|consecutive")
  expect_error(reserve_xy(transform(lt,lx=c(100,110,0)),60,60,0,type="pure_endowment",
                         n=2,P=0,t=2),"lx")
})
