test_that("state annuities and designated reversions match explicit UDD payments", {
  lx <- c(100,80,40,10,0)
  ly <- c(200,140,100,30,0)
  lt <- list(data.frame(x=60:64,lx=lx),data.frame(x=65:69,lx=ly))
  survival <- function(counts,t) {
    j <- floor(t)
    u <- t-j
    (counts[j+1]*(1-u) + counts[j+2]*u)/counts[1]
  }
  for (k in c(1,2,4)) for (timing in c("due","immediate")) {
    times <- 1 + (if (timing == "due") 0:(2*k-1) else 1:(2*k))/k
    sx <- vapply(times,function(t) survival(lx,t),0)
    sy <- vapply(times,function(t) survival(ly,t),0)
    weights <- list(both=sx*sy,x_only=sx*(1-sy),y_only=sy*(1-sx))
    values <- vapply(names(weights),function(name) {
      b <- list(both=0,x_only=0,y_only=0)
      b[[name]] <- 1
      actual <- annuity_xy(lt,60,65,.05,benefit=b,n=2,h=1,k=k,timing=timing)
      expect_equal(actual,sum(1.05^(-times)*weights[[name]])/k,tolerance=1e-12)
      actual
    },0)
    joint <- annuity_xy(lt,60,65,.05,status="joint",n=2,h=1,k=k,timing=timing)
    last <- annuity_xy(lt,60,65,.05,status="last",n=2,h=1,k=k,timing=timing)
    expect_equal(joint,values[1],ignore_attr=TRUE)
    expect_equal(last,sum(values))
    # Reversion to y is the y-life annuity less the both-alive annuity.
    expect_equal(values[3],sum(1.05^(-times)*sy)/k-joint,ignore_attr=TRUE)
  }
})

test_that("premiums and reserves link to an independent annual death ledger", {
  lx <- c(100,80,40,10,0)
  ly <- c(200,140,100,30,0)
  lt <- list(data.frame(x=60:64,lx=lx),data.frame(x=65:69,lx=ly))
  for (status in c("joint","last")) {
    sx <- lx/100
    sy <- ly/200
    s <- if (status == "joint") sx*sy else sx+sy-sx*sy
    benefit_apv <- sum(250*1.05^(-(1:3))*(s[1:3]-s[2:4]))
    for (k in c(1,2)) {
      times <- (0:(3*k-1))/k
      marginal <- function(counts,t) {
        j<-floor(t); u<-t-j
        (counts[j+1]*(1-u)+counts[j+2]*u)/counts[1]
      }
      sxk<-vapply(times,function(t) marginal(lx,t),0)
      syk<-vapply(times,function(t) marginal(ly,t),0)
      sk<-if(status=="joint") sxk*syk else sxk+syk-sxk*syk
      annualized <- benefit_apv/(sum(1.05^(-times)*sk)/k)
      out <- premium_xy(lt,60,65,.05,type="term",status=status,n=3,benefit=250,k=k)
      # The scalar premium output is annualized by the public contract.
      expect_equal(out,annualized,tolerance=1e-12)
      expected_reserves <- vapply(0:2,function(t) {
        j<-0:(3-t)
        sx<-lx[t+j+1]/lx[t+1]
        sy<-ly[t+j+1]/ly[t+1]
        ss<-if(status=="joint") sx*sy else sx+sy-sx*sy
        death<-sum(250*1.05^(-j[-1])*(-diff(ss)))
        grid<-(0:((3-t)*k-1))/k
        conditional <- function(counts,u) {
          jj<-floor(t+u); f<-t+u-jj
          (counts[jj+1]*(1-f)+counts[jj+2]*f)/counts[t+1]
        }
        ax<-vapply(grid,function(u) conditional(lx,u),0)
        ay<-vapply(grid,function(u) conditional(ly,u),0)
        sp<-if(status=="joint") ax*ay else ax+ay-ax*ay
        death-10*sum(1.05^(-grid)*sp)/k
      },0)
      methods<-if(k==1) c("prospective","recursive") else "prospective"
      for(method in methods) {
        reserves<-reserve_xy(lt,60,65,.05,type="term",status=status,n=3,
                              benefit=250,P=10,n_prem=3,k=k,t=0:2,
                              method=method,output="value")
        expect_equal(unname(reserves),expected_reserves,tolerance=1e-12)
      }
    }
  }
})
