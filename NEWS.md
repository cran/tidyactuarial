# tidyactuarial 0.2.0

## Life contingencies and Monte Carlo

* Single-life fractional survival and complete life expectancy retain small
  annual survival derived from survivor counts. CF and Balducci avoid
  cancellation through mortality complements and use logarithmic ratios
  when annual survival underflows; integer horizons and table truncation
  conventions are preserved.

* Monte Carlo summaries retain representable standard deviations and standard
  errors when variance loses range. Empirical TVaR integrates ordered tail
  masses without overflowing excess-loss differences, preserving fractional
  atoms, signed values and endpoint conventions. Empty probability selectors
  omit their columns and empty grouped summaries retain their schema.
  The independence assumption behind `se_mean` and antithetic limits are
  documented explicitly; no pairing is inferred from a method label.

* Monte Carlo premium, loss and reserve consumers accept overflowing annualized
  annuity-basis metadata only when they match a finite per-payment basis and
  frequency. Inconsistent metadata and inference from an infinite annualized
  value alone remain errors, preserving valid chained monetary calculations.

* `mc_reserve` moment-of-death pure endowments require
  survival strictly beyond maturity. Internally estimated premiums retain
  their scale interpretation regardless of the input premium-unit selector;
  supplied premiums retain their declared units.

* Monte Carlo reserves value future monetary flows with the original log
  interest rate. Known uncovered benefits remain zero, missing issue
  annuities are excluded as complete estimation pairs, and actual premium
  income retains range independently of its displayed annuity basis.
  Coincident flows can cancel despite separate present-value overflow.
  Existing in-force and before/after-payment conventions are clarified.

* Monte Carlo premiums and losses retain representable monetary amounts when
  an intermediate annuity-scale coefficient underflows or overflows. Scaled
  means preserve positive subnormal premium exposure, and large valid payment
  frequencies retain numeric storage. Existing paired-NA, grouping and premium
  unit conventions are retained; displayed unrepresentable coefficients and
  amounts retain their double-precision limits.

* Moment-of-death pure endowments in `mc_insurance` require
  survival strictly beyond maturity. A death exactly at maturity belongs
  to the term component, making term plus pure endowment equal endowment
  on each trajectory. The annual curtate convention is retained.

* Pathwise Monte Carlo annuities and insurance retain the original log
  interest rate and representable monetary cash flows at extreme rates.
  Zero amounts, payments at time zero, unknown lifetimes and existing payment
  calendars retain their meanings. Large conversion frequencies preserve
  numeric storage; invalid payment frequencies fail before integer coercion.

* Single-life simulation wrappers retain finite constant-force lifetimes when
  annual survivor ratios underflow. Valuation uses the original interest
  convention and stable monetary annuity sums, including small period rates.
  Invalid ranges and incompatible survivor counts fail before coercion or
  sampling. Finite-horizon censoring and existing random-stream conventions
  are retained and clarified.

* Multiple-life simulated statuses require the exact lifetime
  relation K <= T < K+1 before aggregation; tolerance still checks integer K.
  Complete lifetimes are rejected rather than repaired. Large curtate status
  values retain numeric storage, and identifier/output-name collisions give
  an explicit error. Existing valid min/max results and schemas are retained.

* Lifetime simulation retains positive annual survival under CML and Balducci
  even below the checking tolerance. Simulation from survivor counts uses
  direct decrements and positive exposure independently of the radix.
  Invalid tolerances and out-of-range ages fail before integer coercion;
  valid terminal death ages retain values beyond integer storage. Seeded
  calls preserve their existing random-stream and restoration conventions.

* Two-life fractional survival retains tiny observed annual survival in log
  form. Complete two-life expectancy uses analytic annual UDD, constant-force
  and Balducci integrals, avoiding quadrature loss near steep mortality.
  Existing conditioning, zero-duration and open-table truncation are retained.

* Two-life pure-endowment premiums and reserves use direct discounted status
  survival, retaining small benefits lost by subtracting insurance APVs.
  Joint reserves now use backward conditional recursion; last-survivor reserves
  retain their three-state model. Both use direct annual decrements and original
  log interest. Immediate premiums include those due exactly at the
  reserve date, including a scheduled last installment at maturity.

* Two-life annuities retain tiny designated-reversion probabilities, original
  log interest and representable signed monetary payments on fractional grids.
  Annual-table and payment-grid consistency are checked explicitly.
  Second-order Woolhouse corrections require supported positive borders and following
  years; negative approximations warn without clipping. State meanings and
  open-table truncation limits are documented without changing the signature.

* Two-life insurance sums disjoint annual first/second-death transition masses
  from issue, retaining tiny mortality and representable monetary values with
  extreme interest. State conditioning and existing audit columns are retained.
  Whole insurance without observed status extinction is rejected;
  a zero-term endowment pays at its deferred maturity.

* Two-life survival and expectancy reject incompatible annual ages, increasing
  or invalid survivor counts and contradictory mortality probabilities.
  Invalid controls and out-of-range age coercions fail explicitly; existing
  initial-state, zero-duration and table-horizon conventions are preserved.

* Two-life survival and expectancy normalize inherited `CML` attributes to
  constant force, including equivalent `CF`/`CML` attributes on two tables.

* `annuity_multi()` adds trailing `frac`, values actual payments
  of 1/k under UDD, CF/CML or Balducci with `woolhouse = "none"`, and includes
  both deferred and terminal endpoints in temporary Woolhouse corrections.
  The historical joint-and-survivor meaning of `reversionary` is clarified.
* Variable subperiod death benefits use direct mortality masses and original
  log interest, retaining small decrements and representable extreme APVs.
  Coverage-relative functional benefit times and full payment amounts are
  documented. Invalid tables, horizons and payment grids fail explicitly.
* Cause-specific annual insurance validates competing-exit probabilities and
  consecutive ages, retains cumulative log survival and original interest,
  and preserves its vectorized output and `A_xj` alias. Available-table whole
  insurance is explicitly described as truncated when the table is open.

* Life-contract specifications preserve large integer-valued inputs without
  coercing them to missing values and validate the financial domain of stored
  rates. Premium schedules reject overflowing or empty rounded payment grids;
  standard positive insurance terms cannot round to zero. Component schemas,
  ordinary integer types and original rate specifications are retained.

* Net premiums use the actual term inferred from a variable benefit schedule,
  reject premiums beyond that coverage and accept aligned fractional variable
  terms. Invalid controls and integer coercion failures are checked explicitly.
* Gross premiums validate monetary consistency independently of currency units,
  retain acquisition shares and small expense loadings, and factor the
  equivalence residual to avoid overflowing intermediate APV subtraction.
* Single-life recursive reserves use the equivalent backward identity with a
  supported prospective boundary condition, preserving original rate growth,
  rare death/survival probabilities and nonzero whole-life terminal values.
  Unrequested extinct states are not traversed. Public signatures and output
  schemas remain unchanged.

* `annuity_x()` and `insurance_x()` now value parallel single-life scenarios
  on a shared mortality table, recycling only length-one arguments. Numeric
  vectors work in `mutate()`; tidy output has one row per scenario. Missing or
  incompatible scenario inputs fail explicitly. Existing formals are retained.
* Exact single-life valuations retain original log interest and survival ratios,
  value deferred payments from issue, and handle extinct tails. Annual death
  probabilities use survivor decrements before division. Explicit interest
  overrides on insurance contracts are respected. Woolhouse reports negative
  approximations and rejects undefined positive-survival boundary forces.
* Zero-term endowments pay on survival at deferred maturity.
  Empty annuities retain zero APV, with a pure endowment factor
  of one and the ordinary `annuity_value_at_start` diagnostic column.

* Life cash-flow APVs retain original log interest rates and combine monetary,
  survival and discount factors before rounding the product. Rare last-survivor
  probabilities avoid subtractive cancellation. Missing numeric times remain
  missing rows; mortality queries are bounded by table size and preserve known
  extinction. Invalid mortality values, duplicate ages and infinite dates are
  rejected explicitly. The existing nonmissing APV subtotal is now documented;
  signatures, output columns and row order remain unchanged.

## New features

* `plot_contingent_cash_flow()` plots payment-level contractual amounts,
  unconditional event probabilities, expected cash flows and APV contributions.
  It reuses `apv_life_flow()` outputs and accepts explicit discrete death-event
  tables. Reversionary benefit weights are labeled separately from probabilities;
  `apv_life_flow()` now attaches additive contingency metadata for this purpose.

## Bug fixes

* `t_Ex()` returns discounted survival without an upper cap of
  one, allowing actuarial values above one under negative interest rates.
  Pure endowments and commutation values retain original log discount rates.
* Survival/death probabilities preserve small fractional decrement losses;
  complete expectancy integrals retain small horizons and treat extinct
  cohorts and unknown survivor counts without spurious NaN/index errors.
* Kaplan-Meier `ax` measures time lived by interval decedents
  in years. Empirical exposures are stable across radix scales and use the
  CF/Balducci zero-survival limits; positive deaths at zero exposure yield Inf.
  Standard errors and confidence limits retain very small survival estimates.
* Multiple-decrement tables reject cause names that overwrite totals or other
  causes; fractional cause probabilities retain rare losses. Annual consumers
  validate annual age grids, and KM plots sort and validate step/band data.


* Open life tables built only from `lx` now leave
  terminal deaths/probabilities/central rates missing instead of implying certain
  death. Supplied terminal probabilities and closed-table conventions remain.
* Annual life tables reject contradictory `qx`/`px`, forbidden NaN values,
  invalid controls/radices and unsupported integer ages. Input lengths are
  checked before truncation; accepted near-integer ages are rounded consistently.
* Mortality laws avoid intermediate power overflow in Gompertz/Makeham and
  Logistic forces. Weibull annual increments retain tiny positive hazards;
  Heligman-Pollard handles a unit childhood base with very large exponents.


* Cash-flow diagrams retain finite graphical coordinates at extreme monetary
  magnitudes, use scientific labels for tiny nonzero/very large amounts, reject
  infinite dates explicitly and avoid unnecessary integer coercion of `m`.
  Overlapping labels may be omitted while all arrows and payment times remain;
  existing argument names, raw cash-flow columns and financial valuation remain.


* Immunization retains original log liability rates and forms relative moments
  and scaled constraints before monetary products overflow. The existing
  minimum norm over asset units, SVD rank checks and sensitivity conventions
  are preserved; resulting portfolio measures must be finite.
* Immunization gap plots combine signed weighted cash flows with log discounts,
  preserving finite surplus values and zero-weight assets when individual
  valuations lose range. The grid respects representable annual effective
  rates above -1; invalid grids and nonfinite gaps are rejected explicitly.

* Cash-flow duration and convexity form relative contributions directly,
  avoiding overflowing monetary numerators and recovering representable
  discounted products. Signed near-cancellation checks remain unchanged.
* Bond duration and convexity retain original log yields and normalized
  present-value moments. Calendar tolerances no longer reject small positive
  prices solely because of currency units; unsupported indices and truly
  nonrepresentable risk measures fail explicitly.
* Portfolio duration and convexity avoid monetary products that overflow or
  underflow before taking a finite weighted average. Yield-coordinate and
  time-unit requirements for aggregation are now explicit in documentation.
* With `.na = "propagate"`, portfolio positions
  with missing identifiers now produce a final row with an NA identifier and
  missing sensitivity, value and position count, instead of silently vanishing.
  Missing-value `drop` and `error` policies are preserved.

* Spot discount and forward calculations retain original log growth rates,
  avoiding overflowed power ratios and preserving small rates and flat curves
  across intervals with large common origins. Linear interpolation remains
  on annual effective spot rates.
* Par coupon rates separate monetary scale from discount sums and preserve
  small positive or negative coupons instead of rounding them to zero.
* Missing spot curves keep their plot list positions; forward intervals with
  missing start times now follow the requested missing-value policy.

* Bond book values retain original log yields when valuing remaining payments.
  Callable scenarios pass the original yield specification to pricing rather
  than an effective equivalent that may have rounded to -1 or overflowed.
* Callable dates must identify distinct coupon periods strictly inside
  maturity after alignment; nonrepresentable period counts are rejected before
  integer conversion in book-value and callable schedules.
* Callable `is_worst_case` now marks exactly the
  computed minimum price and all exact computed ties. Calendar tolerance no
  longer determines monetary near ties. Arguments and columns are unchanged.

* Bond prices retain original log yields and recover representable discounted
  amounts when intermediate powers overflow or underflow. Coupon schedules and
  accrued interest also recover finite coupons from overflowing products.
* Bond yields use a stable log price equation, avoiding bracket failures on
  ordinary long bonds and dependence on currency scale. Automatic brackets
  cover representable yields above -1; exhausted iterations now cause an error.
  Nonfinite tolerances, intervals and unsupported integer indices are rejected.
  Public arguments and result columns are preserved.

* General amortization now identifies negative amortization under advance
  payments when accrued interest exceeds regular and extra payments combined.
* Named extra-payment periods are checked before integer conversion, preventing
  fractional names from silently assigning prepayments to the wrong period.
* Loan schedules require finite balance tolerances and reject nonrepresentable
  calculated payments or monetary ledgers explicitly. General level payments
  use a log discount sum to preserve finite payments when discounts overflow.

* Arithmetic annuities retain original log rates and recover representable
  monetary products when standalone factors lose range, preserving signed
  and known zero payments.
* Geometric annuities use stable log sums instead of subtracting nearly equal
  rates or accumulation levels. Small interest-growth spreads are retained
  before period conversion, and convergent perpetuities no longer face an
  arbitrary rate-gap cutoff. Missing finite terms propagate missing values.

* Level annuities retain original log interest rates, small-rate effects over
  long terms and representable monetary values when standalone factors lose
  range. Strictly positive perpetuity rates no longer face an arbitrary cutoff.
* Loan payments and sinking-fund deposits no longer use monetary balance
  tolerances to classify rates as zero. Period rates and variable-rate discount
  factors use original log growth. Public signatures and columns are preserved.
* s_angle documentation now identifies the existing horizon as the end of
  the term, one payment interval after the final due payment.

* `solve_value()` now has a trailing argument `residual_tol = tol`,
  separating result-unit residual accuracy from unknown-coordinate accuracy.
  Audit convergence uses this explicit threshold without a hidden floor.
  Exhausting uniroot iterations is now an error. Existing positional arguments
  and output columns are preserved; an audit flag can become more stringent.

* Present/future cash-flow totals and inferred plot present values recover
  cancellation of overflowing terms with exactly equal finite log factors.
  Representable ordinary sums retain their existing evaluation; unequal-factor
  extreme cancellation remains outside this recovery rule.

* Basic cash-flow valuations recover representable signed capital-factor
  products even when a separate factor overflows, underflows or is subnormal.
  Known zero payments stay zero without converting missing inputs to zero.
* Newton solving now requires a small correction as well as a small residual;
  tiny residual units or a boundary-clipped step no longer imply convergence.
  Nonzero derivatives are no longer discarded by an absolute scale threshold.
* Cash-flow plot present values now use the corrected financial basis.
  Unrepresentable aggregate amounts and inferred present values are rejected;
  tabular pipes, negative numeric timelines and editable ggplots are preserved.

* IRR searches no longer classify small monetary NPVs as grid roots.
  Positive rescaling preserves roots across monetary units, sign checks
  avoid residual products, and reported NPVs retain the original currency.
  Scaling that would erase a nonzero cash flow is rejected explicitly.
* solve_value interval scans now refine sign-changing brackets rather than
  treating all residuals below an absolute tolerance as distinct grid roots.

* Interval factors, equivalent-rate tables and basic cash-flow valuations
  retain the original log accumulation rate when an annual effective
  equivalent rounds to -1 or overflows. Interval discounts no longer require
  a representable reciprocal accumulation. Future cash-flow ratios avoid
  overflow of their separate accumulation levels. Public schemas are retained;
  extreme annual-equivalent audit fields may still display -1 or Inf.

* Rate conversions and accumulation/spot discount factors retain precision
  near zero using log1p and expm1. Small positive or negative rates no longer
  collapse to zero before conversion or lose their effect over long terms.
  Missing rate conventions now propagate missing equivalents as documented.

* Deferred two-life insurance now retains survival conditioning from issue.
  Last-survivor continuation includes paths with only one survivor at the
  end of deferment, instead of restarting both lives at advanced ages.
  Its annuity and maturity-survival audit values use the same conditioning.

* Fractional deferred Monte Carlo annuities now anchor their payment grid
  at the deferment rather than issue time. Non-grid-aligned deferments no
  longer omit eligible payments; the shared reserve calendar is corrected too.
* Negative-binomial and geometric frequency thinning now avoids intermediate
  overflow for tiny original `prob`. Conversions that would round a positive
  payment frequency to zero now fail explicitly and retain the per-loss
  model as the available representation.

* Internal development instructions and audit records are excluded from
  source-package archives.

* Life-contingent Monte Carlo insurance, annuity and reserve calculations
  now reject noninteger or infinite curtate lifetimes, infinite complete
  lifetimes, and observed pairs inconsistent with K <= T < K + 1.
  Missing values retain their existing meaning. Certain annuities remain
  independent of lifetime values. A shared internal validator prevents
  chained valuations from accepting contradictory versions of the same life.

* Monte Carlo life annuities now preserve missing lifetime information as
  missing present values and payment counts, instead of zero payments.
  Certain annuities remain independent of lifetime data.
* Fractional Monte Carlo annuity schedules now compare payment grid dates
  directly with the complete lifetime. A fixed epsilon no longer removes
  payments just before death or introduces payments just after death.
  Existing exact-death endpoint conventions are preserved.
* Guaranteed fractional annuities retain their guaranteed payments even
  when death occurs at time zero. The shared payment schedule also serves
  Monte Carlo reserve and fractional annuity simulation calculations.

* Last-survivor recursive reserves now use backward recursion across the
  both-alive and each single-survivor state. The previous scalar recursion
  mixed these conditioning states and disagreed with prospective reserves.
  Documentation explicitly conditions reserve_xy() on both lives being alive
  at valuation. Joint-life recursion and public output columns are preserved.
* Constant-force fractional lifetime simulation now evaluates log(1-q) using
  log1p(-q). Tiny positive mortality probabilities no longer round the
  denominator to zero and collapse conditional death times to the year start.

* Discrete VaR now compares the requested level directly with cumulative
  probabilities, without shifting it by an absolute epsilon. Levels just
  above an atom or the represented aggregate mass no longer select the
  previous atom or slip into an uncomputed tail.
* VaR from FFT aggregates now checks whether omitted claim-count probability
  could change the selected quantile, using the minimum possible loss of
  omitted counts. Ambiguous quantiles are rejected with advice to reduce
  `tol`. Documentation distinguishes unrepresented mass from an upper loss
  tail: omitted large counts can still produce small aggregate losses.

* Aggregate distributions now infer decimal lattices relative to the severity
  magnitude and merge support points using relative monetary tolerances.
  Small monetary units no longer collapse distinct losses into zero, and
  convolution/FFT `max_loss` cutoffs respect the same scale. Strictly positive
  masses below 1e-18 are retained rather than silently discarded.
* Stop-loss premiums for covered exponential per-loss severities now use
  closed forms for ordinary and franchise deductibles, inflation, coinsurance,
  and finite or infinite limits. Very large limits no longer produce spurious
  zero premiums through quadrature. Their reported method changes from
  `numerical_survival` to `closed_form`; TVaR uses the corrected premiums too.

* Immunization functions now use relative checks for liability present values
  and two-asset duration degeneracy, avoiding rejection of valid portfolios
  solely because of monetary units or small asset prices. Non-finite
  liability valuations and sensitivities are rejected explicitly.
* Immunization systems with more assets than moment constraints now use row-scaled singular
  value decomposition for the minimum-norm asset units, instead of normal
  equations that could falsely classify distinct sensitivities as singular.
* Clarified that `immunize_duration_convexity()` requires convexity with
  respect to the annual effective yield. The coupon-frequency-scaled
  `bond_convexity()` output generally uses a different yield variable;
  documentation now gives the chain-rule conversion. Bond output meanings
  are preserved. Examples now use convexities on the liability yield basis.

* Corrected complete `e_xy(..., cohort = "last")` calculations to preserve
  conditioning on both lives being alive at the initial ages. Annual
  integration now uses survival from those original ages, including paths
  with only one survivor, instead of restarting both lives each year.
  Annual integrals use analytic UDD, constant-force, or Balducci formulas
  for both cohorts, replacing numerical quadrature.
  Missing initial ages and zero initial survivors return `NA` consistently.
* `t_px()` now obtains integer survival directly from `lx`, including at
  the last recorded age, without requiring a next-year probability.
  An extinct attained cohort has zero survival even when its fractional
  probability is undefined. This also prevents spurious `NA` in `e_xy()`.

* Corrected `mortality_law_table("HeligmanPollard", ...)` to the classical
  eight-parameter odds model. `H` is now required; `C` is the childhood
  exponent, `D`/`E`/`F_hp` describe the accident hump, and `G`/`H` describe
  old-age mortality. Previous seven-parameter calls must be respecified;
  adding `H` to the old values is insufficient. The existing `F` alias
  denotes the same parameter as `F_hp`. Log-space evaluation avoids
  undefined probabilities when odds overflow.
* Beta-law conditional probabilities now use log upper-tail survival,
  avoiding premature `qx = 1` when survival is small.
* Mortality-law tables now compute `mx = dx / Lx` and `ex = Tx / lx`
  without artificial denominator floors. Zero-exposure positive deaths
  give `Inf`; undefined `0/0` ratios for extinct cohorts give `NaN`.

* `mortality_law_table()` now derives annual death probabilities from the
  integrated continuous force of mortality, independently of `frac`.
  This changes Gompertz, Makeham, Weibull and Logistic results and the former
  UDD/Balducci conversions. `mu_x` remains the force at integer age;
  `a_x` still controls the exposure approximation. Direct `qx` laws and
  forced terminal closure retain their existing behavior.

* Corrected `lifetable()` conversions between `mx` and `qx` for `ax != 0.5`,
  using the mean fraction of the year lived by those who die. This correction
  can change previous actuarial results; results for `ax = 0.5` are preserved.

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
