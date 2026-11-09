# tidyactuarial

**Actuarial mathematics in R with transparent formulas, reproducible workflows, and tidy outputs where they add value.**

`tidyactuarial` provides tools for financial mathematics, life contingencies, and collective risk in R. The package is designed so that the actuarial model remains visible: formulas, timing conventions, mortality assumptions, and valuation dates come first; computation is then used to reproduce, audit, and scale the calculation.

The current stable release is **tidyactuarial 0.2.0**, published on CRAN
on October 25, 2026. This repository and its reference website document
the released features, including collective risk, insurance coverage,
aggregate loss distributions, and tail risk measures.

## Installation

Install the stable version from CRAN:

```r
install.packages("tidyactuarial")
```

Load the package with:

```r
library(tidyactuarial)
```

## What tidyactuarial covers

The package currently includes tools for:

- **Financial mathematics and interest rates:** accumulation, discounting, equivalent rates, present and future values, equations of value, cash-flow valuation, and IRR.
- **Annuities and loans:** level and varying annuities, amortization schedules, extra principal, variable rates, and general payment rules.
- **Bonds and interest-rate risk:** bond cash flows, price, yield, book value, callable bonds, duration, convexity, and portfolio measures.
- **Term structure and immunization:** spot discounting, forward rates, yield curves, duration matching, and duration-convexity immunization.
- **Life tables and survival:** life-table construction, fractional-age assumptions, survival and death probabilities, life expectancy, Kaplan-Meier estimation, and commutation functions.
- **Life contingencies:** life annuities, life insurance, premiums, reserves, and contract-based workflows.
- **Multiple-life and multiple-decrement models:** joint-life and last-survivor quantities, premiums, reserves, and cause-specific decrement calculations.
- **Simulation and actuarial risk:** lifetime simulation and Monte Carlo tools for annuities, insurance, premiums, reserves, losses, and multiple-life status models.
- **Collective risk and coverage:** frequency and severity models, deductibles and limits, aggregate loss distributions and simulation, VaR, TVaR, and stop-loss premiums.

The full API is organized by actuarial area in the **Reference** section of the package website.

## Collective risk: from a model to an annual loss

In the collective risk model, the annual aggregate loss is the sum of
individual claim amounts:

$$
S = \sum_{j=1}^{N} X_j.
$$

For this example, suppose the number of annual claims follows a Poisson
distribution with mean 3, while individual claim amounts follow a Gamma
distribution with shape 2 and rate 0.005. We can calculate aggregate
moments from the model and then use simulation to examine its tail.

```r
frequency <- frequency_spec("poisson", lambda = 3)

risk_model <- frequency |>
  collective_risk_spec(
    severity = severity_spec(
      "gamma", shape = 2, rate = 0.005
    )
  )

# Expected aggregate loss and second moment
risk_model |>
  aggregate_moments(order = 1:2)

# Simulate annual aggregate losses
losses <- risk_model |>
  simulate_aggregate(n_sim = 5000, seed = 2026)

# Empirical tail risk at 95% and 99%
losses |>
  value_at_risk(level = c(0.95, 0.99))

losses |>
  tail_value_at_risk(level = c(0.95, 0.99))
```

The theoretical expected annual loss is $E[S] = 3(2/0.005) = 1,200$.
The simulated VaR and TVaR are estimates, and will not necessarily match
their theoretical values exactly. The distinction between analytical
calculations and simulation is important when interpreting actuarial results.

## A simple design principle

The package does not try to replace actuarial reasoning with software.

For a payment of amount $C$ at time $t$, financial valuation begins with

$$
PV = C(1+i)^{-t}.
$$

For a payment contingent on survival of a life aged $x$, the corresponding actuarial structure becomes

$$
{}_tE_x
=
v^t\,{}_tp_x.
$$

`tidyactuarial` keeps these quantities explicit and then provides functions that make the same calculations easier to reproduce across scenarios, portfolios, contracts, and datasets.

## Quick start

### Interest-rate conversion and valuation

Different rate conventions can first be placed on a common annual-effective basis.

```r
library(tidyactuarial)
library(dplyr)
library(tibble)

rates <- tibble(
  product = c(
    "effective",
    "nominal interest",
    "force",
    "nominal discount"
  ),
  i = c(
    0.0620,
    0.0605,
    0.0605,
    0.0590
  ),
  i_type = c(
    "effective",
    "nominal_interest",
    "force",
    "nominal_discount"
  ),
  m = c(
    1,
    12,
    1,
    4
  )
)

rates |>
  mutate(
    i_effective = standardize_interest(
      i_type = i_type,
      i = i,
      m = m
    ),
    value_6y = future_value(
      C = 25000,
      i = i_effective,
      t = 6
    )
  )
```

This pattern is useful when actuarial assumptions arrive in heterogeneous forms but the valuation model requires a common basis.

### General loan schedules

For a loan balance $B_{t-1}$, interest $I_t$, scheduled payment $R_t$, and extra principal $E_t$,

$$
B_t
=
B_{t-1}
+
I_t
-
R_t
-
E_t.
$$

A standard amortization schedule can be generated directly:

```r
amort_schedule(
  principal = 50000,
  n = 5,
  i = 0.12,
  i_type = "nominal_interest",
  m = 12,
  k = 12
)
```

When rates or payment rules vary over time, `amort_schedule_general()` provides the more flexible engine.

```r
principal <- 120000
rates <- c(
  0.04,
  0.05,
  0.06,
  0.07,
  0.08,
  0.09
)

amort_schedule_general(
  principal = principal,
  n = 6,
  i = rates,
  payment = NULL
)
```

### Bond valuation and interest-rate risk

```r
price <- bond_price(
  face = 1000,
  c = 0.05,
  n = 10,
  k = 2,
  y = 0.06,
  y_type = "effective"
)

duration <- bond_duration(
  face = 1000,
  c = 0.05,
  n = 10,
  k = 2,
  y = 0.06,
  y_type = "effective"
)

convexity <- bond_convexity(
  face = 1000,
  c = 0.05,
  n = 10,
  k = 2,
  y = 0.06,
  y_type = "effective"
)

tibble(
  price = price,
  duration = duration,
  convexity = convexity
)
```

The package uses the same actuarial specification to move from bond price to duration and convexity.

### From a life table to a premium

```r
data("soa08lt")

premium_x(
  lt = soa08lt,
  x = 40,
  i = 0.05,
  type = "term",
  benefit = 100000,
  n = 20,
  k = 12,
  n_prem = 10,
  output = "summary"
)
```

The same contract can be assembled progressively:

```r
soa08lt |>
  life_contract(
    lives = "single",
    x = 40,
    i = 0.05
  ) |>
  add_insurance(
    type = "term",
    benefit = 100000,
    n = 20
  ) |>
  add_premium_schedule(
    k = 12,
    n_prem = 10,
    timing = "due"
  ) |>
  premium_x(
    output = "summary"
  )
```

The pipe records the actuarial construction of the contract; it does not replace the underlying equation of equivalence.

## Datasets

The package includes selected 2023 period life tables from UN World
Population Prospects 2024 in `mortality_regions`: 67 countries or areas
and 28 geographic aggregates. To compute five-year survival at age 65
using the published survivor counts for Colombia:

```r
col <- subset(mortality_regions, location_code == 170L & sex == "both",
              select = c(x, lx))
t_px(col, x = 65, t = 5)
```

The 100+ interval is open; it does not supply a one-year probability or
an age of extinction. See `?mortality_regions` for units, source and
licence. The existing mortality datasets are unchanged. The
`working-with-datasets` vignette also illustrates competing risks using
the external `MASS::Melanoma` dataset.

`pet_lifetables` includes study-specific life tables for US dogs and cats
and UK dogs. Each table retains its own denominator and open last age;
see `?pet_lifetables` for the published quantities and limitations.

## Articles

Four introductory articles develop complete workflows rather than isolated function calls:

- **Financial Mathematics with tidyactuarial**
- **Loans, Amortization, and General Payment Schedules**
- **Bonds, Duration, Convexity, and Immunization**
- **Life Contingencies with tidyactuarial**

They are available from the **Articles** menu of the package website.

## Output philosophy

Not every actuarial result needs to be a tibble.

`tidyactuarial` uses several output patterns depending on the task:

- direct numerical values for simple actuarial quantities;
- tidy one-row summaries when intermediate quantities are useful;
- schedules for loans, bonds, reserves, and related time-dependent calculations;
- portfolio-level summaries for grouped financial positions;
- pipe-friendly contract objects for life-contingency workflows.

This keeps simple calculations simple while allowing larger analyses to remain inspectable and reproducible.

## Reliability

The package is developed with automated unit and regression tests. Continuous
integration is configured for Linux, Windows, and macOS, including the declared
minimum R version on Ubuntu. A successful run applies to the commit it checks;
it does not certify later development changes.

For released work where stability is the priority, use the CRAN version.

## Documentation

Function-level documentation is available from R:

```r
help(package = "tidyactuarial")
?standardize_interest
?amort_schedule_general
?bond_price
?annuity_x
?premium_x
?reserve_x
```

The package website provides:

- **Reference:** functions organized by actuarial area;
- **Articles:** worked actuarial workflows;
- **Changelog:** package development history.

## Reporting problems

Bug reports and reproducible examples are welcome through the project issue tracker.

For numerical issues, include whenever possible:

- the function call;
- all relevant inputs;
- the expected actuarial result;
- the observed result;
- `sessionInfo()`.

## References

The package follows standard actuarial and financial mathematics. Useful background references include:

- Dickson, D. C. M., Hardy, M. R., & Waters, H. R. *Actuarial Mathematics for Life Contingent Risks*.
- Kellison, S. G. *The Theory of Interest*.
- Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., & Nesbitt, C. J. *Actuarial Mathematics*.
- Klugman, S. A., Panjer, H. H., & Willmot, G. E. (2019). *Loss Models: From Data to Decisions* (5th ed.). Wiley. A primary reference for the collective risk module.

## License

MIT License. See `LICENSE` for details.
