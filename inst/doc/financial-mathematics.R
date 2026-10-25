## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

library(tidyactuarial)

## ----general-accumulation-law-------------------------------------------------
accumulation_law <- function(t) {
  exp(
    0.03 * t +
      0.002 * t^2
  )
}

accumulation <-
  accumulation_factor(
    s = 2,
    t = 5,
    a = accumulation_law
  )

discount <-
  discount_factor(
    s = 2,
    t = 5,
    a = accumulation_law
  )

tibble::tibble(
  quantity = c(
    "Accumulation factor",
    "Discount factor"
  ),
  value = c(
    accumulation,
    discount
  )
)

## ----interest-equivalents-----------------------------------------------------
annual_rate <- 0.08
conversion_frequency <- 12

rate_equivalents <-
  interest_equivalents(
    i_type = "effective",
    i = annual_rate,
    m = conversion_frequency
  )

rate_equivalents

## ----future-value-------------------------------------------------------------
initial_capital <- 1000
interest_rate <- 0.06
term <- 4

future_amount <-
  future_value(
    C = initial_capital,
    i = interest_rate,
    t = term
  )

future_amount

## ----present-value------------------------------------------------------------
future_payment <- 1500

present_amount <-
  present_value(
    C = future_payment,
    i = interest_rate,
    t = term
  )

present_amount

## ----replacement-payment------------------------------------------------------
replacement_rate <- 0.06

replacement_payment <-
  future_value(
    C = 800,
    i = replacement_rate,
    t = 1
  ) +
  present_value(
    C = 1200,
    i = replacement_rate,
    t = 1
  )

replacement_payment

## ----replacement-audit--------------------------------------------------------
original_present_value <-
  pv_flow(
    cf = c(800, 1200),
    i = replacement_rate,
    t = c(1, 3)
  )

replacement_present_value <-
  present_value(
    C = replacement_payment,
    i = replacement_rate,
    t = 2
  )

tibble::tibble(
  valuation = c(
    "Original obligations",
    "Replacement payment"
  ),
  present_value = c(
    original_present_value,
    replacement_present_value
  )
)

## ----general-cash-flow--------------------------------------------------------
cash_flow <-
  tibble::tibble(
    time = c(1, 3, 5, 6),
    amount = c(800, 1200, -500, 1500)
  )

cash_flow_rate <- 0.07

present_value_flow <-
  cash_flow |>
  dplyr::summarise(
    value = pv_flow(
      cf = amount,
      i = cash_flow_rate,
      t = time
    )
  ) |>
  dplyr::pull(value)

future_value_flow <-
  cash_flow |>
  dplyr::summarise(
    value = fv_flow(
      cf = amount,
      i = cash_flow_rate,
      t = time
    )
  ) |>
  dplyr::pull(value)

present_value_flow
future_value_flow

## ----unknown-date-------------------------------------------------------------
unknown_date_rate <- 0.08
replacement_amount <- 12000

obligations_present_value <-
  pv_flow(
    cf = c(4000, 7000),
    i = unknown_date_rate,
    t = c(1, 4)
  )

unknown_date <-
  solve_value(
    fn = present_value,
    solve_for = "t",
    target = obligations_present_value,
    args = list(
      C = replacement_amount,
      i = unknown_date_rate
    ),
    interval = c(3, 4),
    method = "uniroot"
  )

unknown_date

## ----annuity-immediate--------------------------------------------------------
level_payment <- 1000
annuity_term <- 10
annuity_rate <- 0.05

immediate_factor <-
  a_angle(
    n = annuity_term,
    i = annuity_rate,
    i_type = "effective",
    timing = "immediate"
  )

immediate_present_value <-
  level_payment *
  immediate_factor

immediate_present_value

## ----annuity-due--------------------------------------------------------------
due_factor <-
  a_angle(
    n = annuity_term,
    i = annuity_rate,
    i_type = "effective",
    timing = "due"
  )

due_present_value <-
  level_payment *
  due_factor

tibble::tibble(
  timing = c(
    "Immediate",
    "Due"
  ),
  factor = c(
    immediate_factor,
    due_factor
  ),
  present_value = c(
    immediate_present_value,
    due_present_value
  )
)

## ----accumulated-annuity------------------------------------------------------
immediate_accumulation_factor <-
  s_angle(
    n = annuity_term,
    i = annuity_rate,
    i_type = "effective",
    timing = "immediate"
  )

immediate_accumulated_value <-
  level_payment *
  immediate_accumulation_factor

present_value_accumulated <-
  future_value(
    C = immediate_present_value,
    i = annuity_rate,
    t = annuity_term
  )

tibble::tibble(
  route = c(
    "s_angle()",
    "Present value accumulated"
  ),
  accumulated_value = c(
    immediate_accumulated_value,
    present_value_accumulated
  )
)

