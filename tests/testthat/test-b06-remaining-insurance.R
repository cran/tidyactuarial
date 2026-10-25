test_that("variable insurance retains small decrements and original interest", {
  lt <- data.frame(x = 40:42, lx = c(1e16, 1e16 - 2, 1e16 - 4))
  for (frac in c("UDD", "CF", "CML", "Balducci")) {
    expect_equal(insurance_variable_k(lt, 40, 0, 1e16, n = 1, k = 1,
                                     frac = frac), 2)
    expect_equal(insurance_variable_k(lt, 40, 0, 1e16, n = 1, k = 2,
                                     frac = frac), 2)
  }
  lt <- data.frame(x = 40:42, lx = c(100, 50, 25))
  expect_equal(insurance_variable_k(lt, 40, -1000, 1e-300, n = 1, k = 1,
                                   i_type = "force"), exp(1000 + log(.5) + log(1e-300)))
  expect_equal(insurance_variable_k(lt, 40, -1000, 0, n = 1, k = 1,
                                   i_type = "force"), 0)
  expect_equal(insurance_variable_k(lt, 40, .03, 1, n = 1, k = 1, m = 3e9),
               .5 / 1.03)
  tiny <- data.frame(x=40:42,lx=c(1e300,1e-100,5e-101))
  for(frac in c("UDD","CF","Balducci"))
    expect_equal(insurance_variable_k(tiny,40,0,1e300,n=1,h=1,k=2,frac=frac),
                 5e-101,tolerance=1e-12)
  extinct <- data.frame(x=40:42,lx=c(100,0,0))
  expect_equal(insurance_variable_k(extinct,40,-1e308,100,n=1,h=1,k=2,
                                   i_type="force"),0)
})

test_that("variable insurance matches independent fractional death ledgers", {
  lt <- data.frame(x = 40:44, lx = c(100, 80, 56, 28, 0))
  h <- .5
  n <- 1.5
  k <- 2
  times <- h + (1:3) / k
  benefits <- c(120, 80, -20)
  for (frac in c("UDD", "CF", "Balducci")) {
    survival <- function(t) {
      year <- floor(t)
      u <- t - year
      l <- lt$lx[year + 1]
      if (u == 0) return(l / 100)
      p <- lt$lx[year + 2] / l
      f <- switch(frac, UDD = 1 - u * (1 - p), CF = p^u,
                  Balducci = p / (p + u * (1 - p)))
      l / 100 * f
    }
    mass <- vapply(times - 1/k, survival, 0) - vapply(times, survival, 0)
    for (rate in c(-.03, 0, .06)) {
      ref <- sum(benefits * exp(-log1p(rate) * times) * mass)
      expect_equal(insurance_variable_k(lt, 40, rate, benefits, n, h, k,
                                       frac = frac), ref)
    }
    captured <- NULL
    fun <- function(t) { captured <<- t; benefits }
    actual <- insurance_variable_k(lt, 40, .06, fun, n, h, k, frac = frac)
    expect_identical(captured, (1:3) / k)
    expect_equal(actual, sum(benefits * exp(-log1p(.06) * times) * mass))
  }
})

test_that("variable insurance validates the actual table and payment grid", {
  lt <- data.frame(x = 40:43, lx = c(100, 80, 40, 0))
  expect_error(insurance_variable_k(lt, 40, 0, 1, n = 1e308, k = 12), "payments|subperiod|finite")
  expect_error(insurance_variable_k(lt, 40, 0, 1, n = 1, h = 1e308, k = 12), "finite|grid|horizon|deferment")
  expect_error(insurance_variable_k(lt, 40, 0, 1, n = 1, tol = Inf), "tol")
  expect_error(insurance_variable_k(lt[c(1,3,4), ], 40, 0, 1, n = 1), "annual|consecutive")
  expect_error(insurance_variable_k(transform(lt, lx = c(100, 110, 40, 0)),
                                   40, 0, 1, n = 1), "nonincreasing|non-increasing")
  expect_error(insurance_variable_k(lt, 40, 0, 1, n = 5), "horizon|support")
  inherited <- lt
  attr(inherited, "frac") <- "CML"
  expect_identical(insurance_variable_k(inherited, 40, 0, 1, n = 1,
                                       tidy = TRUE)$frac, "CF")
})

test_that("cause insurance is a probability of exit from the active state", {
  md <- md_table(data.frame(x = 40:43, q_death = c(.1,.2,.3,1),
                            q_lapse = c(.2,.1,.2,0)))
  for (rate in c(-.05, 0, .04)) {
    p <- c(1, .7, .49, .245)
    times <- 1:4
    for (cause in c("q_death", "q_lapse")) {
      ref <- sum(300 * p[2:4] * md[[cause]][2:4] * (1 + rate)^(-times[2:4]))
      expect_equal(insurance_xj(md, 40, rate, cause, "term", 300, 3, 1), ref)
    }
    total <- sum(vapply(c("q_death", "q_lapse"), function(cause)
      insurance_xj(md, 40, rate, cause), 0))
    expect_equal(total, sum(p * (1 - md$p_total) * (1 + rate)^(-times)))
  }
  expect_identical(A_xj, insurance_xj)
  out <- insurance_xj(md, x = c(40,41), i = c(.03,.04), "q_death",
                       type = "term", n = 2, tidy = TRUE)
  expect_identical(out$x, c(40L,41L))
  expect_equal(out$insurance_xj, c(insurance_xj(md,40,.03,"q_death","term",n=2),
                                 insurance_xj(md,41,.04,"q_death","term",n=2)))
})

test_that("cause insurance validates intervals, probabilities and controls", {
  md <- data.frame(x = 40:42, p_total = c(.8,.7,0), q_death = c(.1,.2,1))
  expect_error(insurance_xj(md[c(1,3), ], 40, 0, "q_death"), "annual|consecutive")
  expect_error(insurance_xj(transform(md, q_death = c(.5,.2,1)),40,0,"q_death"), "probabil|total")
  expect_error(insurance_xj(md,40,0,"p_total"), "cause|reserved")
  expect_error(insurance_xj(md,40,0,"x"), "cause|reserved")
  expect_error(insurance_xj(md,40,0,NA_character_), "cause")
  expect_error(insurance_xj(md,40,0,"q_death",tol=Inf), "tol")
  expect_error(insurance_xj(md,numeric(),0,"q_death"), "empty|length|numeric")
  expect_error(insurance_xj(md,40,0,"q_death",n=1e20,type="term"), "horizon|available")
  expect_equal(insurance_xj(md,40,.03,"q_death",m=3e9),
               insurance_xj(md,40,.03,"q_death"))
})

test_that("cause insurance retains underflowing survival and monetary products", {
  md <- data.frame(x = 40:42, p_total = c(1e-200,1e-200,0),
                   q_death = c(0,0,1))
  expect_equal(insurance_xj(md,40,0,"q_death","term",benefit=1e300,n=1,h=2),1e-100)
  md <- data.frame(x=40, p_total=.5, q_death=.5)
  expect_equal(insurance_xj(md,40,-1000,"q_death",benefit=1e-300,i_type="force"),
               exp(1000 + log(.5) + log(1e-300)))
  expect_equal(insurance_xj(md,40,-1000,"q_death",benefit=0,i_type="force"),0)
  expect_error(insurance_xj(md,x=c(40,40),i=c(.01,.02,.03),"q_death"),"length")
})
