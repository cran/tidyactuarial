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
#'   \item{t}{Payment time.}
#'   \item{C}{Cash-flow amount. Negative values represent outflows and
#'   positive values represent inflows.}
#'   \item{cashflow_type}{Type of cash flow.}
#'   \item{description}{Short description of the cash flow.}
#' }
#'
#' @details
#' This dataset uses the compact financial-actuarial notation used throughout
#' \code{tidyactuarial}: \code{t} denotes time and \code{C} denotes the cash-flow
#' amount at that time.
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(cash_flows_sample)
#'
#' cash_flows_sample |>
#'   dplyr::filter(scenario_id == "investment_project") |>
#'   dplyr::select(scenario_id, t, C, cashflow_type)
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
#'   \item{face}{Face value of the bond.}
#'   \item{c}{Annual coupon rate.}
#'   \item{k}{Number of coupon payments per year.}
#'   \item{n}{Maturity in years.}
#'   \item{y}{Annual nominal yield rate consistent with coupon frequency.}
#'   \item{bond_type}{Short bond type label.}
#'   \item{P}{Theoretical bond price computed from the listed yield rate.}
#' }
#'
#' @details
#' This dataset uses the compact bond notation adopted in \code{tidyactuarial}:
#' \code{face} is the face value, \code{c} is the annual coupon rate,
#' \code{k} is the coupon frequency, \code{n} is the maturity, \code{y} is the
#' yield input, and \code{P} is the bond price.
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(bonds_sample)
#'
#' bonds_sample |>
#'   dplyr::select(bond_id, face, c, k, n, y, P)
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
#'   \item{L}{Initial loan principal.}
#'   \item{i}{Annual effective interest rate.}
#'   \item{n_months}{Loan term in months.}
#'   \item{k}{Number of payments per year.}
#'   \item{loan_type}{Short loan type label.}
#'   \item{R}{Level payment amount per payment period.}
#' }
#'
#' @details
#' This dataset uses compact loan notation: \code{L} is the initial loan
#' principal, \code{i} is the annual effective interest rate, \code{n_months}
#' is the loan term in months, \code{k} is the payment frequency, and \code{R}
#' is the level payment.
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(loans_sample)
#'
#' loans_sample |>
#'   dplyr::select(loan_id, L, i, n_months, k, R)
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
#'   \item{x}{Integer actuarial age.}
#'   \item{q_death}{One-year death decrement probability.}
#'   \item{q_disability}{One-year disability decrement probability.}
#'   \item{q_withdrawal}{One-year withdrawal decrement probability.}
#'   \item{q_total}{Total one-year decrement probability.}
#'   \item{p_total}{Total one-year survival probability.}
#' }
#'
#' @details
#' This dataset already follows the compact actuarial convention \code{x} for
#' age and \code{q} for decrement probabilities.
#'
#' @source Synthetic pedagogical data created for tidyactuarial examples.
#'
#' @examples
#' data(multiple_decrement_sample)
#'
#' multiple_decrement_sample |>
#'   dplyr::select(x, q_death, q_disability, q_withdrawal, q_total, p_total)
#'
#' @docType data
#' @keywords datasets
"multiple_decrement_sample"
