# tidyactuarial 0.1.6

## Bug fixes

* Fixed the construction of the cumulative distribution function used by
  `simulate_lifetime()` to prevent platform-dependent numerical failures
  when floating-point accumulation produces values slightly greater than one.

* Added a regression test for the multiple-life simulation workflow using
  `simulate_lifetimes()`, `mc_multilife_status()`, and `mc_annuity()`.

* This release addresses the errors reported by CRAN on macOS ARM64.

# tidyactuarial 0.1.5

* Fixed the multiple-life Monte Carlo workflow used in the `mc_annuity()` examples.
* `simulate_lifetimes()` now constructs lifetime simulations defensively when mortality inputs contain invalid terminal values.
* `mc_multilife_status()` now uses defaults compatible with `simulate_lifetimes()`: `sim_id`, `life_id`, `Kx`, and `Tx`.
* Added a regression test for the full multiple-life simulation workflow.
