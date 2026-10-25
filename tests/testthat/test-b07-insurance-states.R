# Enumerate all pairs of annual death times; no valuation engine is the oracle.
two_death_ledger <- function(lx, ly, g, h, n, status, type, benefit = 1) {
  dx <- -diff(lx)/lx[1]
  dy <- -diff(ly)/ly[1]
  total <- 0
  for (tx in seq_along(dx)) for (ty in seq_along(dy)) {
    trigger <- if (status == "joint") min(tx, ty) else max(tx, ty)
    payment <- if (trigger > h && trigger <= h+n) exp(-g*trigger) else 0
    if (type == "endowment" && trigger > h+n)
      payment <- payment + exp(-g*(h+n))
    total <- total + benefit * dx[tx]*dy[ty]*payment
  }
  total
}

test_that("two-life insurance matches independent death-pair enumeration", {
  lx <- c(100, 80, 40, 10, 0)
  ly <- c(200, 140, 100, 30, 0)
  lt <- list(data.frame(x = 60:64, lx = lx), data.frame(x = 65:69, lx = ly))
  for (status in c("joint", "last")) for (type in c("term", "endowment", "whole")) {
    for (h in 0:2) for (g in c(0, .03, -.1)) {
      n <- if (type == "whole") 4-h else 2
      expected <- two_death_ledger(lx, ly, g, h, n, status, type, 250)
      out <- insurance_xy(lt, 60, 65, g, i_type = "force", status = status,
                          type = type, n = if (type == "whole") Inf else n,
                          h = h, benefit = 250, tidy = TRUE)
      expect_equal(out$apv, expected, tolerance = 1e-12)
      expect_equal(out$unit_apv, expected/250, tolerance = 1e-12)
    }
  }
})

test_that("insurance preserves tiny death masses and discount range", {
  a <- data.frame(x = 60:62, lx = c(1e16, 1e16-2, 0))
  b <- data.frame(x = 60:62, lx = c(100, 100, 0))
  expect_equal(insurance_xy(list(a,b),60,60,0,type="term",n=1,benefit=1e16), 2)
  expect_equal(insurance_xy(list(a,a),60,60,0,type="term",status="last",
                            n=1,benefit=1e32), 4)
  simple <- data.frame(x=60:62,lx=c(100,80,0))
  expect_equal(insurance_xy(simple,60,60,1e20,type="term",n=1)/(.36/(1+1e20)), 1,
               tolerance = 1e-12)
  expect_equal(insurance_xy(simple,60,60,700,i_type="force",type="term",n=2,
                            benefit=1e300), .36*exp(log(1e300)-700),
               tolerance = 1e-12)
  # The status survives h with probability .64, although its discounted factor
  # underflows; the monetary APV remains representable.
  expect_equal(insurance_xy(simple,60,60,700,i_type="force",type="term",n=1,
                            h=1,benefit=1e300)/exp(log(.64)+log(1e300)-1400), 1,
               tolerance = 1e-12)
  extreme<-data.frame(x=60:62,lx=c(1e300,1e-100,0))
  expect_equal(insurance_xy(list(extreme,b),60,60,0,type="term",n=1,h=1,
                            benefit=1e300)/1e-100,1,tolerance=1e-12)
})

test_that("annual state transitions distinguish first and second death", {
  # x dies in year 1, y in year 3: 11 -> 01 -> 01 -> 00.
  a <- data.frame(x=60:63,lx=c(100,0,0,0))
  b <- data.frame(x=60:63,lx=c(100,100,100,0))
  for (status in c("joint","last")) {
    trigger <- if (status == "joint") 1 else 3
    expect_equal(insurance_xy(list(a,b),60,60,.05,status=status), 1.05^(-trigger))
    expect_equal(insurance_xy(list(b,a),60,60,.05,status=status), 1.05^(-trigger))
    expect_equal(insurance_xy(list(a,b),60,60,.05,status=status,type="term",n=2),
                 if (trigger <= 2) 1.05^(-trigger) else 0)
  }
  # Both die in the same annual interval: 11 -> 00 pays once for either status.
  for (status in c("joint","last"))
    expect_equal(insurance_xy(list(a,a),60,60,.05,status=status), 1/1.05)
})

test_that("approved zero-term endowments pay at the deferred maturity", {
  lx <- c(100,80,40,0)
  ly <- c(100,50,20,0)
  lt <- list(data.frame(x=60:63,lx=lx),data.frame(x=60:63,lx=ly))
  for (status in c("joint","last")) for (h in 0:2) {
    sx <- lx[h+1]/100
    sy <- ly[h+1]/100
    s <- if (status == "joint") sx*sy else sx+sy-sx*sy
    expect_equal(insurance_xy(lt,60,60,.05,type="endowment",status=status,
                              n=0,h=h,benefit=250), 250*1.05^(-h)*s)
    expect_equal(insurance_xy(lt,60,60,.05,type="term",status=status,n=0,h=h),0)
  }
  expect_equal(insurance_xy(lt,60,60,700,i_type="force",type="endowment",
                            n=0,h=2,benefit=1e300)/exp(log(.08)+log(1e300)-1400),1,
               tolerance=1e-12)
})

test_that("approved whole-life support never implies a terminal payment", {
  open <- data.frame(x=60:62,lx=c(100,80,60))
  closed <- data.frame(x=60:62,lx=c(100,50,0))
  expect_error(insurance_xy(open,60,60,0,status="last"),"extinction|open")
  expect_error(insurance_xy(list(closed,open),60,60,0,status="last"),"extinction|open")
  expect_equal(insurance_xy(list(closed,open),60,60,0,status="joint"),1)
  expect_equal(insurance_xy(open,60,60,0,status="last",type="term",n=2),.16)
  expect_equal(insurance_xy(open,60,60,0,status="last",type="endowment",n=2),1)
  short <- data.frame(x=60:61,lx=c(100,80))
  expect_error(insurance_xy(list(short,closed),60,60,0,status="last"),"extinction|open")
  expect_error(insurance_xy(list(short,closed),60,60,0,status="last",type="term",n=2),
               "open marginal")
})
