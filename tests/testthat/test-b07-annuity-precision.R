test_that("designated reversions retain tiny annual and fractional death masses", {
  a<-data.frame(x=60:62,lx=c(1e16,1e16-2,0))
  b<-data.frame(x=60:62,lx=c(100,100,0))
  for(frac in c("UDD","CF","Balducci")) {
    expect_equal(annuity_xy(list(a,b),60,60,0,n=1,
                 benefit=list(both=0,x_only=0,y_only=1e16),frac=frac),2)
    # At a half-year the limiting death probability is q/2 for each assumption.
    expect_equal(annuity_xy(list(a,b),60,60,0,n=.5,k=2,
                 benefit=list(both=0,x_only=0,y_only=1e16),frac=frac),.5,
                 tolerance=1e-12)
  }
})

test_that("state payments preserve monetary range and original interest", {
  lt<-data.frame(x=60:63,lx=c(100,80,40,0))
  expect_equal(annuity_xy(lt,60,60,720,i_type="force",n=1,
                 benefit=list(both=1e300,x_only=0,y_only=0))/
                 exp(log(.64)+log(1e300)-720),1,tolerance=1e-12)
  extreme<-data.frame(x=60:62,lx=c(1e300,1e-100,0))
  immortal<-data.frame(x=60:62,lx=c(100,100,100))
  expect_equal(annuity_xy(list(extreme,immortal),60,60,0,n=1,
                 benefit=list(both=1e300,x_only=0,y_only=0))/1e-100,1,
                 tolerance=1e-12)
  extinct<-data.frame(x=60:62,lx=c(100,0,0))
  expect_equal(annuity_xy(extinct,60,60,-700,i_type="force",n=2),0)
  expect_equal(annuity_xy(lt,60,60,0,n=1,
                 benefit=list(both=-2,x_only=3,y_only=4)), -.64*2+.16*7)
})

test_that("Woolhouse endpoints agree with analytic constant-force corrections", {
  mu<-.2; g<-.03
  lt<-data.frame(x=60:66,lx=100*exp(-mu*(0:6)))
  for(status in c("joint","last")) for(h in 0:1) for(timing in c("due","immediate")) {
    s<-function(t) if(status=="joint") exp(-2*mu*t) else 2*exp(-mu*t)-exp(-2*mu*t)
    f<-function(t) exp(-g*t)*s(t)
    eta<-function(t) -log(s(t+1)/s(t))
    loss<-f(h)-f(h+2)
    annual<-f(h)+f(h+1)
    for(order in c("first","second")) {
      expected<-annual-(1-1/4)/2*loss
      if(order=="second") expected<-expected-(1-1/16)/12*
        ((g+eta(h))*f(h)-(g+eta(h+2))*f(h+2))
      if(timing=="immediate") expected<-expected-loss/4
      expect_equal(annuity_xy(lt,60,60,g,i_type="force",status=status,n=2,h=h,k=4,
                  timing=timing,woolhouse=order,frac="CF"),expected,tolerance=1e-12)
    }
  }
  short<-lt[1:3,]
  expect_error(annuity_xy(short,60,60,g,i_type="force",n=2,k=2,
                  woolhouse="second"),"border|Woolhouse|second")
})

test_that("annuity tables and payment grids fail explicitly when incompatible", {
  lt<-data.frame(x=60:62,lx=c(100,80,0))
  expect_error(annuity_xy(transform(lt,lx=c(100,110,0)),60,60,0,n=1),"lx")
  expect_error(annuity_xy(lt[c(1,3),],60,60,0,n=1),"annual|consecutive")
  expect_error(annuity_xy(lt,60,60,0,n=1e308,k=2),"grid|payments|finite")
  expect_error(annuity_xy(lt,60,60,0,n=1,k=3e9),"grid|payments")
  expect_equal(annuity_xy(lt,60,60,0,n=0,tidy=TRUE)$apv,0)
})

test_that("all fractional state annuities match independent payment sums", {
  lx<-c(100,85,60,30,10,0)
  ly<-c(200,150,110,70,20,0)
  lt<-list(data.frame(x=60:65,lx=lx),data.frame(x=65:70,lx=ly))
  patterns<-list(c(1,0,0),c(1,1,1),c(0,1,0),c(0,0,1),c(.7,.25,-.3))
  survival<-function(counts,t,frac) {
    j<-floor(t); u<-t-j
    p<-counts[j+2]/counts[j+1]
    f<-switch(frac,UDD=1-u*(1-p),CF=p^u,Balducci=p/(p+u*(1-p)))
    counts[j+1]/counts[1]*f
  }
  growth<-2*log1p(.06/2)
  for(frac in c("UDD","CF","Balducci")) for(k in c(2,4))
    for(timing in c("due","immediate")) for(b in patterns) {
      j<-if(timing=="due") 0:(1.5*k-1) else 1:(1.5*k)
      times<-1+j/k
      sx<-vapply(times,function(t) survival(lx,t,frac),0)
      sy<-vapply(times,function(t) survival(ly,t,frac),0)
      payment<-b[1]*sx*sy+b[2]*sx*(1-sy)+b[3]*sy*(1-sx)
      expected<-sum(exp(-growth*times)*payment)/k
      actual<-annuity_xy(lt,60,65,.06,i_type="nominal_interest",m=2,n=1.5,h=1,k=k,
                         timing=timing,frac=frac,
                         benefit=list(both=b[1],x_only=b[2],y_only=b[3]))
      expect_equal(actual,expected,tolerance=1e-12)
    }
})

test_that("Woolhouse reports a negative approximation without clipping", {
  immortal<-data.frame(x=60:64,lx=rep(100,5))
  expect_warning(value<-annuity_xy(immortal,60,60,100,i_type="force",n=1,k=2,
                      timing="due",woolhouse="second"),"negative")
  expect_equal(value,1-.25*(1-exp(-100))-.0625*100*(1-exp(-100)))
  expect_gt(annuity_xy(immortal,60,60,100,i_type="force",n=1,k=2,
                      timing="due",woolhouse="none"),0)
})
