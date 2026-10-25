test_that("net premiums match independent annual and fractional payment ledgers", {
  lt <- data.frame(x = 60:66, lx = c(1000, 800, 500, 200, 80, 20, 0))
  for (i in c(-.03, 0, .05)) for (h in 0:1) for (start in c("issue", "deferred")) {
    n <- 3
    benefit <- 500
    j <- h + 0:2
    b <- benefit * sum((lt$lx[j + 1] - lt$lx[j + 2]) / 1000 * (1 + i)^-(j + 1))
    for (k in c(1, 4)) {
      t <- (if (start == "issue") 0 else h) + (0:(2 * k - 1)) / k
      f <- floor(t)
      u <- t - f
      survival <- ((1 - u) * lt$lx[f + 1] + u * lt$lx[f + 2]) / 1000
      a <- sum(survival * (1 + i)^-t) / k
      out <- premium_x(lt, 60, i, type = "term", benefit = benefit, n = n,
                       h = h, k = k, premium_start = start, n_prem = 2,
                       output = "summary")
      expect_equal(out$premium_annualized, b / a)
      expect_equal(out$premium_per_payment, b / (a * k))
      expect_equal(out$apv_benefits, b)
      expect_equal(out$apv_premium_annuity, a)
      expect_lt(abs(out$equivalence_residual) / max(1, b), 1e-14)
    }
  }
})

test_that("variable-benefit premium terms use their actual inferred coverage", {
  lt <- data.frame(x = 60:65, lx = c(100, 80, 50, 20, 10, 0))
  expect_error(premium_x(lt, 60, .05, type = "variable_k", benefit = c(100, 100),
                         k = 2, n_prem = 2), "beyond.*coverage")
  out <- premium_x(lt, 60, .05, type = "variable_k", benefit = c(100, 100),
                   k = 2, output = "audit")
  expect_equal(out$value[out$component == "insurance_term_years"], 1)
  # Two quarters of UDD deaths, discounted at the ends of those quarters.
  b <- 100 * .05 * sum(1.05^-(1:2 / 4))
  a <- sum(c(1, .95) * 1.05^-(0:1 / 4)) / 4
  out <- premium_x(lt, 60, .05, type = "variable_k", benefit = c(100, 100),
                   n = .5, k = 4, n_prem = .5, output = "summary")
  expect_equal(out$premium_annualized, b / a)
})

test_that("premium inputs reject invalid controls and preserve large valid m", {
  lt <- data.frame(x = 60:63, lx = c(100, 80, 50, 0))
  expect_error(premium_x(lt, 60, .05, check = NA), "check")
  expect_error(premium_x(lt, 1e20, .05), "x")
  expect_error(premium_x(lt, 60, .05, h = 1e20), "h")
  expect_error(premium_x(lt, 60, .05, n_prem = 1e308, k = 12), "n_prem")
  expect_no_warning(p <- premium_x(lt, 60, .05, type = "term", n = 2, m = 3e9))
  expect_true(is.finite(p))
})

test_that("gross premiums enforce monetary consistency without a unit floor", {
  for (scale in c(1e-100, 1, 1e100)) {
    prem <- data.frame(premium_annualized = 1200 * scale,
                       premium_per_payment = 90 * scale,
                       payments_per_year = 12, apv_premium_annuity = 10)
    expect_error(premium_gross(prem), "inconsistent")
    prem$premium_per_payment <- 100 * scale
    out <- premium_gross(prem, alpha = .5, beta = .05, gamma = 20 * scale,
                         output = "summary")
    expected <- (1200 + 240) / (.95 - .5 / 120)
    expect_equal(out$gross_premium_annualized / scale, expected)
  }
  prem <- data.frame(premium_annualized = 1200, premium_per_payment = NA_real_,
                     payments_per_year = 12, apv_premium_annuity = 10)
  expect_error(premium_gross(prem), "payment|premium")
  prem$premium_per_payment <- 100
  prem$payments_per_year <- NA_real_
  expect_error(premium_gross(prem, k = 12), "frequency|payments_per_year")
})

test_that("gross expense shares and loadings retain representable values", {
  prem <- data.frame(premium_annualized = 1, payments_per_year = 12,
                     apv_premium_annuity = 1e308)
  out <- premium_gross(prem, alpha = 1e308, output = "summary")
  expect_equal(out$gross_premium_annualized, 12 / 11)
  prem <- data.frame(premium_annualized = 1e20, payments_per_year = 1,
                     apv_premium_annuity = 1)
  out <- premium_gross(prem, beta = 1e-20, output = "summary")
  expect_equal(out$loading_annualized, 1)
  prem <- data.frame(premium_annualized = 2, payments_per_year = 12,
                     apv_premium_annuity = 1e308)
  out <- premium_gross(prem, output = "summary")
  expect_equal(out$gross_premium_annualized, 2)
  expect_equal(out$equivalence_residual, 0)
  prem <- data.frame(premium_annualized = 0, payments_per_year = 3e9,
                     apv_premium_annuity = 1)
  expect_no_warning(g <- premium_gross(prem, gamma = 1e-9))
  expect_equal(g, 3)
})

test_that("gross expense APVs match an independent normalized ledger", {
  for (k in c(1, 4, 12)) for (scale in c(1e-100, 1, 1e100)) {
    a <- 2.5
    net <- 150 * scale
    gamma <- 7 * scale
    alpha <- .4
    beta <- .06
    g <- (net * a + gamma * k * a) / ((1 - beta) * a - alpha / k)
    out <- premium_gross(data.frame(premium_annualized = net,
      payments_per_year = k, apv_premium_annuity = a), alpha = alpha,
      beta = beta, gamma = gamma, output = "audit")
    val <- function(nm) out$value[out$component == nm] / scale
    expect_equal(val("gross_premium_annualized"), g / scale)
    expect_equal(val("apv_initial_expense"), alpha * g / k / scale)
    expect_equal(val("apv_proportional_expense"), beta * g * a / scale)
    expect_equal(val("apv_fixed_per_payment_expense"), gamma * k * a / scale)
    expect_lt(abs(val("equivalence_residual")), 1e-10)
  }
})

test_that("prospective reserves match conditional future cash flows at every duration", {
  lt <- data.frame(x = 60:66, lx = c(1000, 800, 500, 200, 80, 20, 0))
  for (type in c("term", "endowment")) for (h in 0:1) for (k in c(1, 4))
    for (start in c("issue", "deferred")) {
    i <- .04
    n <- 3
    np <- 2
    benefit <- 500
    P <- 35
    durations <- 0:(h + n)
    ref <- vapply(durations, function(t) {
      l0 <- lt$lx[t + 1]
      deaths <- (h + 0:(n - 1))[h + 1:n > t]
      b <- benefit * sum((lt$lx[deaths + 1] - lt$lx[deaths + 2]) /
                          l0 * (1 + i)^-(deaths + 1 - t))
      if (type == "endowment") b <- b + benefit * lt$lx[h + n + 1] /
        l0 * (1 + i)^-(h + n - t)
      times <- (if (start == "issue") 0 else h) + (0:(np * k - 1)) / k
      times <- times[times >= t]
      f <- floor(times)
      u <- times - f
      p <- ((1 - u) * lt$lx[f + 1] + u * lt$lx[f + 2]) / l0
      b - P / k * sum(p * (1 + i)^-(times - t))
    }, numeric(1))
    out <- reserve_x(lt, 60, i, type = type, n = n, h = h, benefit = benefit,
                     P = P, k = k, n_prem = np, premium_start = start,
                     t = durations, output = "value")
    expect_equal(unname(out), ref)
    if (k == 1) {
      rec <- reserve_x(lt, 60, i, type = type, n = n, h = h, benefit = benefit,
                       P = P, n_prem = np, premium_start = start,
                       t = durations, method = "recursive",
                       output = "value")
      expect_equal(unname(rec), ref)
    }
  }
})

test_that("recursive reserves are stable for rapidly decreasing survival", {
  lt <- data.frame(x = 60:71, lx = 1e6 * .01^(0:11))
  for (P in c(70, 90, 110)) {
    durations <- 0:10
    # q=.99, v=1/1.1, so each unit of remaining due-annuity has net cost 90-P.
    ref <- (90 - P) * (1 - (.01 / 1.1)^(10 - durations)) / (1 - .01 / 1.1)
    rec <- reserve_x(lt, 60, .1, type = "term", n = 10, benefit = 100,
                     P = P, t = durations, method = "recursive", output = "value")
    expect_equal(unname(rec), ref, tolerance = 1e-10)
  }
})

test_that("recursive reserves preserve rare decrements and original growth", {
  lt <- data.frame(x = 60:62, lx = c(1e16, 1e16 - 2, 1e16 - 4))
  out <- reserve_x(lt, 60, .05, type = "term", n = 1, benefit = 1e16, P = 0,
                   t = 0:1, method = "recursive", output = "value")
  expect_equal(unname(out), c(2 / 1.05, 0), tolerance = 1e-13)
  lt <- data.frame(x = 60:62, lx = c(100, 80, 50))
  for (force in c(-40, 1000)) {
    benefit <- if (force == 1000) 1e300 else 1e-20
    out <- reserve_x(lt, 60, force, i_type = "force", type = "term", n = 1,
                     benefit = benefit, P = 0, t = 0:1, method = "recursive", output = "value")
    ref <- exp(log(benefit) + log(.2) - force)
    expect_equal(unname(out) / ref, c(1, 0))
  }
  lt <- data.frame(x = 60:62, lx = c(100, 1e-18, 1e-20))
  expect_no_error(reserve_x(lt, 60, .05, type = "term", n = 2, P = 0,
                            method = "recursive"))
})

test_that("recursive computation does not require unrequested extinct states", {
  lt <- data.frame(x = 60:64, lx = c(100, 80, 0, 0, 0))
  rec <- reserve_x(lt, 60, .05, type = "term", n = 3, P = 0, t = 0,
                   method = "recursive", output = "value")
  expect_equal(unname(rec), .2 / 1.05 + .8 / 1.05^2)
})

test_that("reserve controls and extreme inputs fail without coercion warnings", {
  lt <- data.frame(x = 60:63, lx = c(100, 80, 50, 0))
  expect_error(reserve_x(lt, 1e20, .05), "x")
  expect_error(reserve_x(lt, 60, .05, t = 1e20), "t")
  expect_error(reserve_x(lt, 60, .05, h = 1e20), "h")
  expect_error(reserve_x(lt, 60, .05, n_prem = 1e308, k = 12), "n_prem")
  expect_error(reserve_x(transform(lt, lx = c(100, 110, 50, 0)),
                         60, .05, P = 0, t = 0, n_prem = 1), "nonincreasing")
  expect_no_warning(r <- reserve_x(lt, 60, .05, P = 0, m = 3e9, t = 0))
  expect_true(is.finite(r$reserve))
})

test_that("whole-life recursion keeps the nonzero final supported reserve", {
  lt <- data.frame(x = 60:63, lx = c(100, 80, 50, 0))
  ref <- c(0, 0, 100 / 1.05 - 10)
  ref[2] <- (.625 * ref[3] + .375 * 100) / 1.05 - 10
  ref[1] <- (.8 * ref[2] + .2 * 100) / 1.05 - 10
  for (method in c("prospective", "recursive")) {
    out <- reserve_x(lt, 60, .05, P = 10, benefit = 100,
                     method = method, output = "summary")
    expect_equal(out$reserve, ref)
    expect_equal(out$t, 0:2) # Preserve the documented supported positive cohorts.
  }
})

test_that("fractional deferred premium terms retain the actual future installments", {
  lt <- data.frame(x = 60:65, lx = c(100, 80, 50, 20, 10, 0))
  times <- 1 + (0:5) / 4
  for (t in 0:2) {
    future <- times[times >= t]
    f <- floor(future)
    u <- future - f
    l <- (1 - u) * lt$lx[f + 1] + u * lt$lx[f + 2]
    a <- sum(l / lt$lx[t + 1] * 1.05^-(future - t)) / 4
    # Use a unit benefit at final endowment maturity; enumerate annual deaths.
    deaths <- (1:3)[2:4 > t]
    b <- sum((lt$lx[deaths + 1] - lt$lx[deaths + 2]) / lt$lx[t + 1] *
               1.05^-(deaths + 1 - t)) + lt$lx[5] / lt$lx[t + 1] * 1.05^-(4 - t)
    out <- reserve_x(lt, 60, .05, type = "endowment", n = 3, h = 1,
                     P = 2, k = 4, n_prem = 1.5, premium_start = "deferred",
                     t = t, output = "value")
    expect_equal(unname(out), b - 2 * a)
  }
})

test_that("gross residuals report expenses below the rounded premium increment", {
  prem <- data.frame(premium_annualized = 1e20, payments_per_year = 1,
                     apv_premium_annuity = 1)
  out <- premium_gross(prem, beta = 1e-20, output = "summary")
  expect_equal(out$equivalence_residual, -1)
  prem <- data.frame(premium_annualized = 1e300, payments_per_year = 12,
                     apv_premium_annuity = 1e300)
  out <- premium_gross(prem, alpha = 1e-300, output = "summary")
  expect_equal(out$loading_annualized / (1e-300 / 12), 1)
  expect_equal(out$equivalence_residual, -1 / 12)
})

test_that("a zero premium stays zero when its diagnostic annuity factor overflows", {
  lt <- data.frame(x = 60:62, lx = c(100, 80, 80))
  out <- reserve_x(lt, 60, -1000, i_type = "force", type = "term", n = 2,
                   benefit = 1e-300, P = 0, n_prem = 2, t = 0, output = "audit")
  ref <- exp(log(1e-300) + log(.2) + 1000)
  expect_equal(out$value[out$component == "apv_future_premiums"], 0)
  expect_equal(out$value[out$component == "reserve"] / ref, 1)
})
