# tidyactuarial 0.1.4

* Fixed the multiple-life Monte Carlo workflow used in the `mc_annuity()` examples.
* `simulate_lifetimes()` now constructs lifetime simulations defensively when mortality inputs contain invalid terminal values.
* `mc_multilife_status()` now uses defaults compatible with `simulate_lifetimes()`: `sim_id`, `life_id`, `Kx`, and `Tx`.
* Added a regression test for the full multiple-life simulation workflow.
