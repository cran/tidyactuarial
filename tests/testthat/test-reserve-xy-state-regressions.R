test_that("last-survivor recursion preserves the both-alive conditioning state", {
  lx <- c(100, 80, 40, 10, 0)
  ly <- c(100, 50, 30, 5, 0)
  tables <- list(data.frame(x = 60:64, lx = lx),
                 data.frame(x = 60:64, lx = ly))
  v <- 1 / 1.05
  for (type in c("term", "endowment", "pure_endowment", "whole")) {
    horizon <- if (type == "whole") 4 else 3
    for (n_prem in c(2, horizon)) {
      expected <- vapply(0:2, function(t) {
        j <- 0:(horizon - t)
        # Independent survival from the valuation state: both lives alive.
        px <- lx[t + j + 1] / lx[t + 1]
        py <- ly[t + j + 1] / ly[t + 1]
        s <- px + py - px * py
        death <- if (type == "pure_endowment") 0 else
          100 * sum(v^j[-1] * (-diff(s)))
        maturity <- if (type %in% c("endowment", "pure_endowment"))
          100 * v^max(j) * tail(s, 1) else 0
        payable <- head(j, -1) + t < n_prem
        premiums <- 10 * sum(v^head(j, -1) * head(s, -1) * payable)
        death + maturity - premiums
      }, numeric(1))
      for (method in c("prospective", "recursive")) {
        out <- reserve_xy(tables, 60, 60, .05, status = "last", type = type,
          n = if (type == "whole") Inf else 3, benefit = 100, P = 10,
          n_prem = n_prem, t = 0:2, method = method, output = "value")
        expect_equal(unname(out), expected, tolerance = 1e-10)
      }
    }
  }
})
