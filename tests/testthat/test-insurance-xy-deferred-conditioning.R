test_that("deferred insurance retains issue-age survival conditioning", {
  lx <- c(100, 80, 40, 10, 0)
  ly <- c(100, 50, 30, 5, 0)
  lt <- list(data.frame(x = 60:64, lx = lx),
             data.frame(x = 60:64, lx = ly))
  for (status in c("joint", "last")) {
    px <- lx / lx[1]
    py <- ly / ly[1]
    s <- if (status == "joint") px * py else px + py - px * py
    for (i in c(0, .05, -.02)) {
      for (type in c("term", "endowment", "whole")) {
        h <- 1
        n <- if (type == "whole") 3 else 2
        years <- h + seq_len(n)
        death <- sum((1 + i)^(-years) * (s[years] - s[years + 1]))
        maturity <- if (type == "endowment")
          (1 + i)^(-max(years)) * s[max(years) + 1] else 0
        out <- insurance_xy(lt, 60, 60, i, type = type, status = status,
          h = h, n = if (type == "whole") Inf else n, benefit = 100,
          tidy = TRUE)
        expect_equal(out$apv, 100 * (death + maturity), tolerance = 1e-10)
        # Audit columns describe continuation conditional on status alive at h,
        # including single-survivor paths for the last-survivor status.
        times <- h + seq_len(n) - 1
        ann <- sum((1 + i)^(-(times - h)) * s[times + 1] / s[h + 1])
        expect_equal(out$annuity_due_value, ann, tolerance = 1e-12)
        if (type == "term") {
          expect_equal(out$survival_to_maturity,
            s[h + n + 1] / s[h + 1], tolerance = 1e-12)
        }
      }
    }
  }
})
