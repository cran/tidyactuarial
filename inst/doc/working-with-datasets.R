## ----setup--------------------------------------------------------------------
library(tidyactuarial)

## ----wpp----------------------------------------------------------------------
col <- subset(mortality_regions, location_code == 170L & sex == "both",
              select = c(x, lx))
t_px(col, x = 65, t = 5)

## ----melanoma, eval=requireNamespace("MASS", quietly = TRUE)------------------
d <- MASS::Melanoma
times <- sort(unique(d$time[d$status != 2L]))
risk <- vapply(times, function(t) sum(d$time >= t), 0L)
death <- sapply(c(melanoma = 1L, other = 3L), function(cause) {
  vapply(times, function(t) sum(d$time == t & d$status == cause), 0L)
})
increment <- death / risk
survival <- cumprod(1 - rowSums(increment))
before <- c(1, head(survival, -1))
incidence <- apply(increment * before, 2, cumsum)
tail(data.frame(days = times, survival, incidence), 1)

