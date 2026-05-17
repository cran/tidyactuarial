#' Sample cash flows for interest theory examples
#'
#' @description
#' A small pedagogical dataset containing cash-flow scenarios for present value,
#' future value, net present value, internal rate of return, equations of value,
#' and cash-flow diagrams.
#'
#' @format A tibble with 19 rows and 5 variables:
#' \describe{
#'   \item{scenario_id}{Scenario identifier.}
#'   \item{time}{Payment time.}
#'   \item{amount}{Cash-flow amount. Negative values represent outflows and
#'   positive values represent inflows.}
#'   \item{cashflow_type}{Type of cash flow.}
#'   \item{description}{Short description of the cash flow.}
#' }
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(cash_flows_sample)
#'
#' cash_flows_sample |>
#'   dplyr::filter(scenario_id == "investment_project")
#'
#' @docType data
#' @keywords datasets
"cash_flows_sample"


#' Sample bond contracts for fixed-income examples
#'
#' @description
#' A small pedagogical dataset containing bond contracts for pricing,
#' yield-to-maturity, duration, and convexity examples.
#'
#' @format A tibble with 5 rows and 8 variables:
#' \describe{
#'   \item{bond_id}{Bond identifier.}
#'   \item{face_value}{Face value of the bond.}
#'   \item{coupon_rate}{Annual coupon rate.}
#'   \item{coupon_frequency}{Number of coupon payments per year.}
#'   \item{maturity_years}{Maturity in years.}
#'   \item{yield_rate}{Annual nominal yield rate consistent with coupon frequency.}
#'   \item{bond_type}{Short bond type label.}
#'   \item{price}{Theoretical bond price computed from the listed yield rate.}
#' }
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(bonds_sample)
#'
#' bonds_sample |>
#'   dplyr::select(bond_id, face_value, coupon_rate, maturity_years, price)
#'
#' @docType data
#' @keywords datasets
"bonds_sample"


#' Sample loan contracts for amortization examples
#'
#' @description
#' A small pedagogical dataset containing level-payment loan contracts for
#' amortization schedule and outstanding balance examples.
#'
#' @format A tibble with 4 rows and 7 variables:
#' \describe{
#'   \item{loan_id}{Loan identifier.}
#'   \item{principal}{Initial loan principal.}
#'   \item{annual_effective_rate}{Annual effective interest rate.}
#'   \item{term_months}{Loan term in months.}
#'   \item{payments_per_year}{Number of payments per year.}
#'   \item{loan_type}{Short loan type label.}
#'   \item{payment}{Level payment amount.}
#' }
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(loans_sample)
#'
#' loans_sample |>
#'   dplyr::select(loan_id, principal, annual_effective_rate, term_months, payment)
#'
#' @docType data
#' @keywords datasets
"loans_sample"


#' Sample multiple decrement probabilities
#'
#' @description
#' A small pedagogical annual multiple decrement dataset with three causes:
#' death, disability, and withdrawal. It is intended for examples involving
#' multiple decrement tables, total-decrement life tables, cause-specific
#' decrement probabilities, and cause-specific insurance benefits.
#'
#' @format A tibble with 7 rows and 6 variables:
#' \describe{
#'   \item{x}{Integer age.}
#'   \item{q_death}{One-year death decrement probability.}
#'   \item{q_disability}{One-year disability decrement probability.}
#'   \item{q_withdrawal}{One-year withdrawal decrement probability.}
#'   \item{q_total}{Total one-year decrement probability.}
#'   \item{p_total}{Total one-year survival probability.}
#' }
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(multiple_decrement_sample)
#'
#' md <- md_table(
#'   qx_df = multiple_decrement_sample |>
#'     dplyr::select(x, q_death, q_disability, q_withdrawal),
#'   radix = 100000,
#'   close = FALSE
#' )
#'
#' md
#'
#' @docType data
#' @keywords datasets
"multiple_decrement_sample"
