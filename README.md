
## Overview

**tidyactuarial** is a comprehensive R package for financial mathematics
and life contingencies. It bridges the gap between classical actuarial
notation (SOA/CAS standards) and the modern **tidyverse** workflow.

The package is designed for actuaries, students, and researchers who
need reproducible, vectorized, and "tidy" actuarial calculations.

## Key Features

- **Life Contingencies:** Valuation of life insurance, annuities, and
  premiums for single and multiple lives.
- **Flexible Assumptions:** Supports UDD (Uniform Distribution of
  Deaths), Constant Force, and Balducci assumptions for fractional ages.
- **Financial Math:** Bond pricing (including callable bonds), yield
  curve interpolation (Spot/Forward rates), and immunization
  (Redington/Full).
- **Tidy Design:** All main functions return `tibble` objects and are
  compatible with the native R pipe `|>` or `%>%`.

## Installation

Once the package is on CRAN, you can install it with:

``` r
install.packages("tidyactuarial")
#> Installing package into 'C:/Users/LENOVO/AppData/Local/Temp/RtmpG8Qiyo/temp_libpath3c505756476c'
#> (as 'lib' is unspecified)
#> Warning: package 'tidyactuarial' is not available for this version of R
#> 
#> A version of this package for your version of R might be available elsewhere,
#> see the ideas at
#> https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
```

Alternatively, you can install the development version from
[GitHub](https://www.google.com/search?q=https://github.com/jfajardo1908/tidyactuarial)
with:

``` r
# install.packages("devtools")
devtools::install_github("julian.fajardo1908/tidyactuarial")
```

## Quick Start

### 1. Life Tables and Premiums

Building a life table and calculating a net level premium for a 20-year
endowment insurance:

``` r
library(tidyactuarial)

# 1. Create a life table (ensure qx has the same length as x)
ages <- 0:110
probs <- rep(0.002, length(ages)) 

lt <- lifetable(x = ages, qx = probs, radix = 100000)

# 2. Calculate net premium for a 35-year-old (20-year endowment)
# Added benefit = 100000
premium_x(
  lt = lt, 
  x = 35, 
  n = 20, 
  i = 0.05, 
  product = "endowment", 
  benefit = 100000,
  tidy = TRUE
)
#> # A tibble: 1 x 15
#>       x     m     n product benefit     k frac  premium_timing prem_start n_prem
#>   <int> <int> <int> <chr>     <dbl> <int> <chr> <chr>          <chr>       <int>
#> 1    35     0    20 endowm...  100000     1 UDD   due            issue          20
#> # i 5 more variables: woolhouse <chr>, premium <dbl>, premium_annual <dbl>,
#> #   apv_benefits <dbl>, apv_premiums <dbl>
```

### 2. Bond Cash Flows

Generate a tidy schedule for a coupon bond:

``` r
bond_cash_flows(
  face = 1000, 
  coupon_rate = 0.06, 
  years_to_maturity = 5, 
  coupons_per_year = 2
)
#> # A tibble: 11 x 5
#>    cashflow_id period  time cash_flow type      
#>          <int>  <int> <dbl>     <dbl> <chr>     
#>  1           1      1   0.5        30 coupon    
#>  2           2      2   1          30 coupon    
#>  3           3      3   1.5        30 coupon    
#>  4           4      4   2          30 coupon    
#>  5           5      5   2.5        30 coupon    
#>  6           6      6   3          30 coupon    
#>  7           7      7   3.5        30 coupon    
#>  8           8      8   4          30 coupon    
#>  9           9      9   4.5        30 coupon    
#> 10          10     10   5          30 coupon    
#> 11          11     10   5        1000 redemption
```

## References

Mathematical formulas and actuarial notation follow:

- Finan, M. B. *A Reading of the Theory of Life Contingencies*.
- Kellison, S. G. *The Theory of Interest*.
- Dickson, D. C., Hardy, M. R., & Waters, H. R. *Actuarial Mathematics
  for Life Contingent Risks*.
