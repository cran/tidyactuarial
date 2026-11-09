#' UN World Population Prospects 2024 period life tables for 2023
#'
#' A selection of official 2023 period estimates, with single-year ages
#' 0 through 99 and an open 100+ interval, for 67 countries or areas and
#' 28 geographic aggregates. These are model-based estimates for a
#' reference period, not observed individual lifetimes or cohort forecasts.
#'
#' @format A tibble with 28,785 rows and 21 columns:
#' \describe{
#'   \item{location_code, location}{UN WPP location identifier and name.
#'     For countries and areas the identifier is the numeric M49 code;
#'     aggregates use WPP identifiers, which need not be M49 codes.}
#'   \item{iso3}{ISO alpha-3 code for countries and areas; missing for aggregates.}
#'   \item{level, region, subregion}{Selection level (country, world, major
#'     region or subregion) and geographic grouping. Missing grouping labels
#'     denote levels where the grouping does not apply.}
#'   \item{year, sex}{Calendar year (2023) and sex category
#'     (\code{"female"}, \code{"male"} or \code{"both"}).}
#'   \item{x, age_interval, age_span, open_interval}{Starting age in years,
#'     original interval label, interval width in years (\code{-1} is the
#'     original open-interval marker), and an explicit open-interval flag.}
#'   \item{mx}{Central mortality rate, deaths per person-year.}
#'   \item{qx, px}{Death and survival probabilities over the interval, in
#'     base 1. At 100+, 1 and 0 refer to the open interval.}
#'   \item{lx, dx}{Survivors at the start and published deaths in the
#'     interval, per initial 100,000 lives.}
#'   \item{Lx, Tx}{Person-years in the interval and total remaining
#'     person-years, per initial 100,000 lives.}
#'   \item{ex, ax}{Remaining life expectancy and mean years lived in
#'     the interval by those who die, in years.}
#' }
#'
#' @details
#' Keys are \code{location_code}, \code{year}, \code{sex} and \code{x}.
#' No numeric source values were rescaled, interpolated or replaced.
#' Rounding can cause small differences between the published probabilities
#' and ratios of published survivor counts. For valuation from survivor
#' counts, supply only \code{x} and \code{lx} from a single table and
#' choose a supported finite horizon. Do not interpret 100+ as a one-year
#' interval or infer an age of extinction from it. In particular, the
#' published open-interval \code{dx} in the 84 aggregate rows differs
#' from \code{lx * qx}; it is retained as published, without treating it
#' as a one-year decrement or forcing the identity at that boundary.
#'
#' @source United Nations, Department of Economic and Social Affairs,
#'   Population Division (2024). World Population Prospects 2024, complete
#'   period life tables, 2023 estimates. \url{https://population.un.org/wpp/}.
#'   Data attribution: United Nations, CC BY 3.0 IGO
#'   (\url{https://creativecommons.org/licenses/by/3.0/igo/}).
#'   This package subsets and renames columns, adds geographic selection
#'   labels and an open-interval flag, and preserves the numeric values.
#'
#' @examples
#' col <- subset(mortality_regions, location_code == 170L & sex == "both",
#'               select = c(x, lx))
#' t_px(col, x = 65, t = 5)
"mortality_regions"
