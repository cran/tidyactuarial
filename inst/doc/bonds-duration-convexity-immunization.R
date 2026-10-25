## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

library(tidyactuarial)

## ----bond-cash-flows----------------------------------------------------------
face_value <- 1000000
redemption_value <- 1000000
annual_coupon_rate <- 0.09
years_to_maturity <- 6
coupons_per_year <- 2

bond_flows <-
  bond_cash_flows(
    face = face_value,
    c = annual_coupon_rate,
    n = years_to_maturity,
    k = coupons_per_year,
    R = redemption_value
  )

bond_flows

## ----bond-price---------------------------------------------------------------
semiannual_yield <- 0.05

bond_value <-
  bond_price(
    face = face_value,
    c = annual_coupon_rate,
    n = years_to_maturity,
    k = coupons_per_year,
    y_effective_per_period = semiannual_yield,
    R = redemption_value
  )

bond_value

## ----bond-ytm-----------------------------------------------------------------
market_price <- 935000

ytm_result <-
  bond_ytm(
    P = market_price,
    face = face_value,
    c = annual_coupon_rate,
    n = years_to_maturity,
    k = coupons_per_year,
    R = redemption_value,
    interval = c(0.05, 0.06)
  )

semiannual_ytm <-
  ytm_result |>
  dplyr::pull(i_period) |>
  dplyr::first()

ytm_result

## ----bond-book-value----------------------------------------------------------
valuation_time <- 3

book_value <-
  bond_book_value(
    face = face_value,
    c = annual_coupon_rate,
    n = years_to_maturity,
    t = valuation_time,
    k = coupons_per_year,
    y_effective_per_period = semiannual_ytm,
    R = redemption_value
  )

book_value

## ----callable-bond------------------------------------------------------------
call_dates <- c(4, 6)
call_values <- c(1040000, 1020000)

target_yield <- 0.04

callable_price <-
  bond_callable_price(
    face = 1000000,
    c = 0.10,
    n = 8,
    k = 2,
    call_t = call_dates,
    call_R = call_values,
    y_effective_per_period = target_yield,
    R = 1000000
  )

callable_price

## ----cash-flow-duration-convexity---------------------------------------------
cash_flows <- c(120, 160, 220, 300)
payment_times <- c(1, 2, 4, 6)
valuation_rate <- 0.07

macaulay_duration <-
  duration_cash_flow(
    cf = cash_flows,
    t = payment_times,
    rate = valuation_rate,
    duration = "macaulay",
    output = "value"
  )

modified_duration <-
  duration_cash_flow(
    cf = cash_flows,
    t = payment_times,
    rate = valuation_rate,
    duration = "modified",
    output = "value"
  )

cash_flow_convexity <-
  convexity_cash_flow(
    cf = cash_flows,
    t = payment_times,
    rate = valuation_rate,
    output = "value"
  )

tibble::tibble(
  measure = c("Macaulay duration", "Modified duration", "Convexity"),
  value = c(macaulay_duration, modified_duration, cash_flow_convexity)
)

## ----bond-duration-convexity--------------------------------------------------
bond_duration_result <-
  bond_duration(
    face = 1000,
    c = 0.08,
    n = 10,
    k = 2,
    y = 0.09,
    y_type = "effective"
  )

bond_convexity_result <-
  bond_convexity(
    face = 1000,
    c = 0.08,
    n = 10,
    k = 2,
    y = 0.09,
    y_type = "effective"
  )

bond_duration_result
bond_convexity_result

## ----portfolio-risk-----------------------------------------------------------
portfolio <-
  tibble::tibble(
    asset = c("A", "B", "C"),
    units = c(450, 300, 180),
    price = c(96.20, 102.75, 88.50),
    duration = c(2.80, 5.60, 9.40),
    convexity = c(11.40, 43.20, 104.60)
  ) |>
  dplyr::mutate(
    market_value = units * price
  )

portfolio_duration_result <-
  portfolio_duration(
    portfolio,
    col_P = "market_value",
    col_D = "duration",
    col_portfolio = NULL
  )

portfolio_convexity_result <-
  portfolio_convexity(
    portfolio,
    col_P = "market_value",
    col_C = "convexity",
    col_portfolio = NULL
  )

portfolio_duration_result
portfolio_convexity_result

## ----immunization-inputs------------------------------------------------------
liabilities <- c(30000, 40000, 50000)
liability_times <- c(3, 5, 8)
reference_rate <- 0.05

asset_maturities <- c(2, 5, 10)
asset_redemption <- 1000

asset_prices <-
  asset_redemption *
  (1 + reference_rate)^(-asset_maturities)

asset_durations <-
  asset_maturities

asset_convexities <-
  asset_maturities *
  (asset_maturities + 1) /
  (1 + reference_rate)^2

## ----duration-match-----------------------------------------------------------
duration_match <-
  immunize_duration(
    L = liabilities,
    t = liability_times,
    P = asset_prices[c(1, 3)],
    D = asset_durations[c(1, 3)],
    i = reference_rate,
    i_type = "effective"
  )

duration_match

## ----duration-convexity-match-------------------------------------------------
duration_convexity_match <-
  immunize_duration_convexity(
    L = liabilities,
    t = liability_times,
    P = asset_prices,
    D = asset_durations,
    C = asset_convexities,
    i = reference_rate,
    i_type = "effective"
  )

duration_convexity_match

## ----immunization-gap, fig.cap="Immunization gap around the reference rate."----
asset_cash_flows <- list(
  list(
    cf = 444998.22,
    t = 3
  ),
  list(
    cf = 561800.00,
    t = 7
  )
)

plot_immunization_gap(
  L = 1000000,
  t = 5,
  asset_cashflows = asset_cash_flows,
  w = c(1, 1),
  i = 0.06,
  i_type = "effective",
  delta = 0.02
)

