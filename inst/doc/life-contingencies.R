## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(tidyactuarial)
library(dplyr)
library(tibble)

## ----life-table---------------------------------------------------------------
lt <- lifetable(
  x = 60:66,
  lx = c(
    100000,
    99000,
    97500,
    95500,
    93000,
    90000,
    86000
  ),
  close = FALSE,
  frac = "UDD"
)

lt

## ----survival-probabilities---------------------------------------------------
survival_audit <- tibble(
  duration = 0:5
) |>
  mutate(
    survival = t_px(
      lt = lt,
      x = 60,
      t = duration,
      frac = "UDD"
    ),
    death = t_qx(
      lt = lt,
      x = 60,
      t = duration,
      frac = "UDD"
    )
  )

survival_audit

## ----pure-endowment-----------------------------------------------------------
i <- 0.05
v <- 1 / (1 + i)

pure_endowment_manual <- v^5 * 90000 / 100000

pure_endowment_package <- t_Ex(
  lt = lt,
  x = 60,
  t = 5,
  i = i,
  frac = "UDD"
)

tibble(
  method = c(
    "Manual",
    "tidyactuarial"
  ),
  value = c(
    pure_endowment_manual,
    pure_endowment_package
  )
)

## ----life-annuity-------------------------------------------------------------
survival_0_4 <- lt$lx[1:5] / lt$lx[1]

annuity_manual <- sum(
  v^(0:4) * survival_0_4
)

annuity_package <- annuity_x(
  lt = lt,
  x = 60,
  i = i,
  n = 5,
  timing = "due",
  frac = "UDD"
)

tibble(
  method = c(
    "Manual",
    "tidyactuarial"
  ),
  apv = c(
    annuity_manual,
    annuity_package
  )
)

## ----term-insurance-----------------------------------------------------------
death_0_4 <- (
  lt$lx[1:5] - lt$lx[2:6]
) / lt$lx[1]

insurance_manual <- sum(
  v^(1:5) * death_0_4
)

insurance_package <- insurance_x(
  lt = lt,
  x = 60,
  i = i,
  n = 5,
  type = "term",
  benefit = 1
)

tibble(
  method = c(
    "Manual",
    "tidyactuarial"
  ),
  apv = c(
    insurance_manual,
    insurance_package
  )
)

## ----insurance-benefit--------------------------------------------------------
benefit <- 100000

insurance_x(
  lt = lt,
  x = 60,
  i = i,
  n = 5,
  type = "term",
  benefit = benefit
)

## ----identity-audit-----------------------------------------------------------
d <- i / (1 + i)

identity_value <- 1 -
  d * annuity_package -
  pure_endowment_package

tibble(
  direct_insurance = insurance_package,
  identity_value = identity_value,
  difference = insurance_package - identity_value
)

## ----net-premium--------------------------------------------------------------
premium_manual <- benefit *
  insurance_manual /
  annuity_manual

premium_package <- premium_x(
  lt = lt,
  x = 60,
  i = i,
  type = "term",
  benefit = benefit,
  n = 5,
  k = 1,
  n_prem = 5,
  timing = "due",
  output = "value"
)

tibble(
  method = c(
    "Manual",
    "tidyactuarial"
  ),
  premium = c(
    premium_manual,
    premium_package
  )
)

## ----premium-summary----------------------------------------------------------
premium_x(
  lt = lt,
  x = 60,
  i = i,
  type = "term",
  benefit = benefit,
  n = 5,
  k = 1,
  n_prem = 5,
  timing = "due",
  output = "summary"
)

## ----reserve-manual-----------------------------------------------------------
survival_62 <- lt$lx[3:5] / lt$lx[3]

annuity_62 <- sum(
  v^(0:2) * survival_62
)

death_62 <- (
  lt$lx[3:5] - lt$lx[4:6]
) / lt$lx[3]

insurance_62 <- sum(
  v^(1:3) * death_62
)

reserve_2_manual <- benefit *
  insurance_62 -
  premium_package *
  annuity_62

reserve_2_manual

## ----reserve-schedule---------------------------------------------------------
reserve_schedule <- reserve_x(
  lt = lt,
  x = 60,
  i = i,
  type = "term",
  n = 5,
  benefit = benefit,
  P = premium_package,
  k = 1,
  n_prem = 5,
  timing = "due",
  t = 0:5,
  method = "prospective",
  output = "summary"
)

reserve_schedule

## ----reserve-recursive--------------------------------------------------------
reserve_x(
  lt = lt,
  x = 60,
  i = i,
  type = "term",
  n = 5,
  benefit = benefit,
  P = premium_package,
  k = 1,
  n_prem = 5,
  timing = "due",
  t = 0:5,
  method = "recursive",
  output = "summary"
)

## ----contract-pipeline--------------------------------------------------------
contract <- life_contract(
  lt = lt,
  lives = "single",
  x = 60,
  i = i
) |>
  add_insurance(
    type = "term",
    benefit = benefit,
    n = 5
  ) |>
  add_premium_schedule(
    k = 1,
    n_prem = 5,
    timing = "due"
  )

contract |>
  premium_x(
    output = "summary"
  )

## ----soa-benchmark------------------------------------------------------------
data("soa08lt")

soa08lt |>
  filter(
    x >= 60,
    x <= 65
  ) |>
  select(
    x,
    lx,
    dx,
    qx,
    px
  )

## ----soa-values---------------------------------------------------------------
tibble(
  quantity = c(
    "Whole-life annuity-due",
    "Whole-life insurance, benefit 100000"
  ),
  value = c(
    annuity_x(
      lt = soa08lt,
      x = 65,
      i = 0.05,
      timing = "due"
    ),
    insurance_x(
      lt = soa08lt,
      x = 65,
      i = 0.05,
      type = "whole",
      benefit = 100000
    )
  )
)

