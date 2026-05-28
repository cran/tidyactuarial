# tidyactuarial 0.1.3

## Major changes

* Added `life_contract()` as a lightweight contract specification object for pipe-friendly life-contingency workflows.
* Added support for using `life_contract()` with single-life functions such as `annuity_x()`, `insurance_x()`, `premium_x()`, and `reserve_x()`.
* Updated two-life functions to use a unified API based on `mortality_table`, `age_x`, `age_y`, `rate`, `term_years`, `deferment_years`, and `output`.
* Updated `annuity_xy()`, `insurance_xy()`, `premium_xy()`, and `reserve_xy()` to support pipe-friendly workflows and two independent life tables through `mortality_table = list(table_x, table_y)`.
* Added `summary_mc()` for tidy summaries of Monte Carlo simulation output, including means, variances, standard errors, quantiles, VaR, and TVaR.
* Added `soa08lt`, a tidy version of the SOA Illustrative Life Table, for examples, validation, and benchmark tests.

## Monte Carlo simulation

* Added or improved support for Monte Carlo simulation workflows for single-life annuities and life insurance.
* Added compatibility between simulation functions and `life_contract()` objects.
* Improved Monte Carlo summaries through `summary_mc()`, with support for grouped summaries using `by` or `group_cols`.

## Cash-flow visualization

* Redesigned `plot_cash_flow()` as a more professional deterministic cash-flow diagram tool.
* Added support for both classical vector input and tidyverse-style pipe workflows.
* Added support for numeric times and calendar dates.
* Added real-scale cash-flow diagrams by default, with optional normalized schematic diagrams.
* Added automatic aggregation of cash flows occurring at the same time or date.
* Added improved labels, subtitles, present-value display, and currency formatting.

## Data

* Added `soa08lt`, a tidy actuarial life table with columns `x`, `lx`, `dx`, `qx`, and `px`.
* Added tests validating the structure of `soa08lt` and its compatibility with basic single-life and two-life actuarial calculations.

## Internal improvements

* Updated tests to reflect the new unified API.
* Removed obsolete documentation and examples for retired geometric annuity helper functions.
* Replaced old references to `summarise_mc()` with `summary_mc()`.
* Improved consistency of argument names across life-contingency functions.
* Improved documentation for life-contingency workflows, Monte Carlo summaries, and cash-flow diagrams.

## Breaking changes

* The two-life functions now use the unified argument names `mortality_table`, `age_x`, `age_y`, `rate`, `term_years`, `deferment_years`, and `output`.
* Older argument names such as `lt`, `x`, `y`, `i`, `n`, `m`, and `tidy` are no longer used in the updated two-life API.
* Retired old geometric annuity helper documentation in favor of the unified annuity interface.
