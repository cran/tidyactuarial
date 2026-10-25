## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

library(tidyactuarial)

## ----level-payment------------------------------------------------------------
loan_principal <- 50000
loan_term_years <- 5
payments_per_year <- 12
nominal_rate <- 0.12
number_payments <- 60

annuity_factor <-
  a_angle(
    n = loan_term_years,
    k = payments_per_year,
    i = nominal_rate,
    i_type = "nominal_interest",
    m = 12,
    timing = "immediate"
  )

level_payment <-
  loan_principal /
  annuity_factor

level_payment

## ----level-schedule-----------------------------------------------------------
level_schedule <-
  amort_schedule(
    principal = loan_principal,
    n = number_payments,
    i = nominal_rate,
    i_type = "nominal_interest",
    m = 12,
    k = payments_per_year,
    timing = "immediate"
  )

level_schedule |>
  dplyr::filter(
    period %in% c(1, 2, 12, 24, 60)
  ) |>
  dplyr::select(
    period,
    ob_start,
    interest,
    payment,
    principal,
    ob_end
  )

## ----extra-principal----------------------------------------------------------
extra_principal <- c(
  "12" = 5000,
  "24" = 3000
)

schedule_term <-
  amort_schedule(
    principal = loan_principal,
    n = number_payments,
    i = nominal_rate,
    i_type = "nominal_interest",
    m = 12,
    k = payments_per_year,
    timing = "immediate",
    extra_principal = extra_principal,
    adjust = "term"
  )

schedule_payment <-
  amort_schedule(
    principal = loan_principal,
    n = number_payments,
    i = nominal_rate,
    i_type = "nominal_interest",
    m = 12,
    k = payments_per_year,
    timing = "immediate",
    extra_principal = extra_principal,
    adjust = "payment"
  )

## ----extra-principal-summary--------------------------------------------------
summary_original <-
  level_schedule |>
  dplyr::summarise(
    scenario = "No extra principal",
    number_payments = dplyr::n(),
    total_interest = sum(interest),
    total_paid = sum(cashflow),
    ending_balance = dplyr::last(ob_end)
  )

summary_term <-
  schedule_term |>
  dplyr::summarise(
    scenario = "Shorter term",
    number_payments = dplyr::n(),
    total_interest = sum(interest),
    total_paid = sum(cashflow),
    ending_balance = dplyr::last(ob_end)
  )

summary_payment <-
  schedule_payment |>
  dplyr::summarise(
    scenario = "Lower payment",
    number_payments = dplyr::n(),
    total_interest = sum(interest),
    total_paid = sum(cashflow),
    ending_balance = dplyr::last(ob_end)
  )

dplyr::bind_rows(
  summary_original,
  summary_term,
  summary_payment
)

## ----variable-rates-----------------------------------------------------------
variable_rates <- c(
  0.04,
  0.05,
  0.06,
  0.07,
  0.08,
  0.09
)

variable_rate_schedule <-
  amort_schedule_general(
    principal = 120000,
    n = 6,
    i = variable_rates,
    i_type = "effective",
    k = 1,
    payment = NULL,
    timing = "immediate"
  )

variable_rate_schedule |>
  dplyr::select(
    period,
    ob_start,
    interest,
    payment,
    principal,
    ob_end
  )

## ----constant-principal-rule--------------------------------------------------
fixed_principal <-
  loan_principal /
  number_payments

constant_principal_rule <- function(
  period,
  ob_start,
  i_effective_period
) {
  fixed_principal +
    i_effective_period * ob_start
}

constant_principal_schedule <-
  amort_schedule_general(
    principal = loan_principal,
    n = number_payments,
    i = nominal_rate,
    i_type = "nominal_interest",
    m = 12,
    k = payments_per_year,
    payment = constant_principal_rule,
    timing = "immediate"
  )

constant_principal_schedule |>
  dplyr::filter(
    period %in% c(1, 2, 12, 24, 60)
  ) |>
  dplyr::select(
    period,
    ob_start,
    interest,
    payment,
    principal,
    ob_end
  )

