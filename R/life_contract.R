#' Create a life-contingency contract specification
#'
#' Creates a lightweight contract object that stores common life-contingency
#' inputs for use in pipe workflows.
#'
#' This function does not compute actuarial values. It validates and stores
#' common inputs such as the mortality table, life status, ages, and interest
#' rate specification. Calculation functions such as [annuity_x()],
#' [insurance_x()], [premium_x()], [reserve_x()], [annuity_xy()],
#' [insurance_xy()], [premium_xy()], and simulation functions can then consume
#' this object.
#'
#' @param mortality_table A life table or a list of two life tables. For
#'   single-life contracts, provide one data.frame or tibble. For two-life
#'   contracts, provide either one data.frame used for both lives, or
#'   `list(table_x, table_y)` with one table for each life. Each life table must
#'   contain columns `x` and `lx`.
#' @param lives Character string. Use `"single"` for a single-life contract,
#'   `"joint"` for a joint-life two-life contract, or `"last_survivor"` for a
#'   last-survivor two-life contract.
#' @param age Numeric scalar. Age for a single-life contract.
#' @param age_x Numeric scalar. First age for two-life contracts.
#' @param age_y Numeric scalar. Second age for two-life contracts.
#' @param rate Numeric scalar. Annual interest-rate input.
#' @param rate_type Character string indicating the rate type. Allowed values
#'   are `"effective"`, `"nominal_interest"`, `"nominal_discount"`, and `"force"`.
#' @param m Positive integer. Compounding frequency for nominal rates. Ignored
#'   for `rate_type = "effective"` and `rate_type = "force"`.
#' @param ... Reserved for future extensions.
#'
#' @return An object of class `"tidyact_life_contract"`.
#'
#' @family life-contingencies
#'
#' @export
life_contract <- function(
    mortality_table,
    lives = c("single", "joint", "last_survivor"),
    age = NULL,
    age_x = NULL,
    age_y = NULL,
    rate,
    rate_type = "effective",
    m = 1L,
    ...
) {
  lives <- match.arg(lives)

  validate_one_life_table <- function(tab, label = "mortality_table") {
    if (!is.data.frame(tab)) {
      stop(
        "`", label, "` must be a data.frame or tibble.",
        call. = FALSE
      )
    }

    if (!all(c("x", "lx") %in% names(tab))) {
      stop(
        "`", label, "` must contain columns `x` and `lx`.",
        call. = FALSE
      )
    }

    if (!is.numeric(tab$x)) {
      stop(
        "Column `x` in `", label, "` must be numeric.",
        call. = FALSE
      )
    }

    if (!is.numeric(tab$lx)) {
      stop(
        "Column `lx` in `", label, "` must be numeric.",
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  if (lives == "single") {
    validate_one_life_table(mortality_table, "mortality_table")
  } else {
    if (is.data.frame(mortality_table)) {
      validate_one_life_table(mortality_table, "mortality_table")
    } else if (
      is.list(mortality_table) &&
      length(mortality_table) == 2L &&
      all(vapply(mortality_table, is.data.frame, logical(1L)))
    ) {
      validate_one_life_table(mortality_table[[1L]], "mortality_table[[1]]")
      validate_one_life_table(mortality_table[[2L]], "mortality_table[[2]]")
    } else {
      stop(
        "For two-life contracts, `mortality_table` must be either one life ",
        "table or a list of two life tables `list(table_x, table_y)`.",
        call. = FALSE
      )
    }
  }

  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  if (!is.numeric(rate) ||
      length(rate) != 1L ||
      is.na(rate) ||
      !is.finite(rate)) {
    stop("`rate` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(rate_type) ||
      length(rate_type) != 1L ||
      is.na(rate_type)) {
    stop("`rate_type` must be a single character string.", call. = FALSE)
  }

  valid_rate_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  if (!rate_type %in% valid_rate_type) {
    stop(
      "`rate_type` must be one of: ",
      paste(sprintf("'%s'", valid_rate_type), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.numeric(m) ||
      length(m) != 1L ||
      is.na(m) ||
      !is.finite(m) ||
      m < 1 ||
      abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  m <- as.integer(round(m))

  if (lives == "single") {
    if (is.null(age)) {
      stop("`age` must be provided when `lives = 'single'`.", call. = FALSE)
    }

    if (!is.numeric(age) ||
        length(age) != 1L ||
        is.na(age) ||
        !is.finite(age) ||
        abs(age - round(age)) > 1e-10) {
      stop("`age` must be a single integer age.", call. = FALSE)
    }

    age <- as.integer(round(age))
    age_x <- NULL
    age_y <- NULL
  } else {
    if (is.null(age_x) || is.null(age_y)) {
      stop(
        "`age_x` and `age_y` must be provided for two-life contracts.",
        call. = FALSE
      )
    }

    if (!is.numeric(age_x) ||
        length(age_x) != 1L ||
        is.na(age_x) ||
        !is.finite(age_x) ||
        abs(age_x - round(age_x)) > 1e-10) {
      stop("`age_x` must be a single integer age.", call. = FALSE)
    }

    if (!is.numeric(age_y) ||
        length(age_y) != 1L ||
        is.na(age_y) ||
        !is.finite(age_y) ||
        abs(age_y - round(age_y)) > 1e-10) {
      stop("`age_y` must be a single integer age.", call. = FALSE)
    }

    age <- NULL
    age_x <- as.integer(round(age_x))
    age_y <- as.integer(round(age_y))
  }

  out <- list(
    mortality_table = mortality_table,
    lives = lives,
    age = age,
    age_x = age_x,
    age_y = age_y,
    rate = rate,
    rate_type = rate_type,
    m = m
  )

  class(out) <- c("tidyact_life_contract", "list")

  out
}


#' @export
print.tidyact_life_contract <- function(x, ...) {
  cat("<tidyact_life_contract>\n")
  cat("  lives:     ", x$lives, "\n", sep = "")

  if (identical(x$lives, "single")) {
    cat("  age:       ", x$age, "\n", sep = "")
  } else {
    cat("  age_x:     ", x$age_x, "\n", sep = "")
    cat("  age_y:     ", x$age_y, "\n", sep = "")
  }

  cat("  rate:      ", x$rate, "\n", sep = "")
  cat("  rate_type: ", x$rate_type, "\n", sep = "")
  cat("  m:         ", x$m, "\n", sep = "")

  if (is.list(x$mortality_table) && !is.data.frame(x$mortality_table)) {
    cat("  tables:    ", length(x$mortality_table), "\n", sep = "")
  } else {
    cat("  tables:    1\n", sep = "")
  }

  invisible(x)
}


.as_life_contract <- function(x) {
  inherits(x, "tidyact_life_contract")
}
