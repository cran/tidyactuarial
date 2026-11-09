#' Published life tables for companion dogs and cats
#'
#' Three tables from veterinary studies: US dogs and cats (Montoya et
#' al., 2023, tables 3 and 4) and UK dogs (Teng et al., 2022, table 2).
#' These are study-specific life tables, not mortality bases for insured
#' pets or data from a common population.
#'
#' @format A tibble with 57 rows and 17 columns:
#' \describe{
#'   \item{table_id, species, country}{Source table identifier, species,
#'     and study country. Use one table at a time.}
#'   \item{x, age_interval, open_interval}{Age at interval start in years,
#'     interval label and open-interval flag.}
#'   \item{n_deaths, n_midyear}{Deaths in the source study and mid-interval
#'     population. Mid-interval population is unavailable for UK dogs.}
#'   \item{mx, qx, ax}{Published central death rate, probability of death
#'     in the stated interval (base 1) and mean years lived in the interval
#'     by those who die. \code{mx} is unavailable for UK dogs and
#'     \code{ax} for the US tables.}
#'   \item{lx, Lx, Tx}{Published survivors at interval start, animal-years
#'     lived in the interval and remaining animal-years. The US tables
#'     use a hypothetical radix of 100,000; UK dog counts start at 30,563.}
#'   \item{ex, ex_ci_low, ex_ci_high}{Remaining life expectancy and
#'     published 95 percent confidence limits, in years.}
#' }
#'
#' @details
#' Values are reproduced at their published precision. US \code{n_deaths}
#' is a study count and is not the decrement of the hypothetical
#' \code{lx} column. Do not compare raw survivor counts across studies
#' or mix cohorts and denominators. The US tables end at 17+ and the
#' UK dog table at 20+. The terminal \code{qx} is reproduced as published
#' and is not a one-year probability; no age of extinction is inferred.
#' For a supported finite survival calculation, use only \code{x} and
#' \code{lx} from the selected table.
#'
#' @source Montoya et al. (2023), \emph{Life expectancy tables for dogs
#'   and cats derived from clinical data}, tables 3 and 4,
#'   \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC9989186/}.
#'   Teng et al. (2022), \emph{Life tables of annual life expectancy and
#'   mortality for companion dogs in the United Kingdom}, table 2,
#'   \doi{10.1038/s41598-022-10341-6}.
#'   Both articles are licensed under CC BY 4.0
#'   (\url{https://creativecommons.org/licenses/by/4.0/}).
#'   This package selects and renames columns and normalizes interval
#'   punctuation; it does not change the numeric values.
#'
#' @examples
#' dogs <- subset(pet_lifetables, table_id == "montoya2023_us_dog",
#'                select = c(x, lx))
#' t_px(dogs, x = 0, t = 5)
"pet_lifetables"
