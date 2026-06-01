#' Create a life-contingency contract specification
#'
#' Creates a lightweight actuarial contract object that stores common
#' life-contingency inputs for use in pipe workflows.
#'
#' This function does not compute actuarial values. It validates and stores
#' common actuarial inputs such as the life table, life status, ages, and
#' interest-rate specification. Calculation functions such as [annuity_x()],
#' [insurance_x()], [premium_x()], [reserve_x()], [annuity_xy()],
#' [insurance_xy()], [premium_xy()], [reserve_xy()], and simulation functions
#' can then consume this object.
#'
#' @param lt A life table or a list of two life tables. For single-life
#'   contracts, provide one data frame or tibble. For two-life contracts,
#'   provide either one data frame used for both lives, or `list(lt_x, lt_y)`
#'   with one table for each life. Each life table must contain columns `x`
#'   and `lx`.
#' @param lives Character string. Use `"single"` for a single-life contract,
#'   `"joint"` for a joint-life two-life contract, or `"last_survivor"` for a
#'   last-survivor two-life contract.
#' @param x Numeric scalar. Age of the single life, or age of the first life in
#'   a two-life contract.
#' @param y Numeric scalar. Age of the second life in a two-life contract.
#'   Required when `lives` is `"joint"` or `"last_survivor"`.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are `"effective"`, `"nominal_interest"`, `"nominal_discount"`,
#'   and `"force"`.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for `i_type = "effective"` and `i_type = "force"`.
#' @param ... Reserved for future extensions. Deprecated argument names such as
#'   `mortality_table`, `age`, `age_x`, `age_y`, `rate`, and `rate_type` are
#'   not accepted.
#'
#' @details
#' `life_contract()` follows the compact actuarial notation used throughout
#' `tidyactuarial`:
#'
#' * `lt`: life table;
#' * `x`: age of the first or single life;
#' * `y`: age of the second life;
#' * `i`: interest rate;
#' * `i_type`: type of interest rate;
#' * `m`: conversion frequency for nominal rates.
#'
#' The object stores actuarial fields using the compact names above. During the
#' 0.1.4 API transition, it also stores internal compatibility fields so that
#' functions not yet migrated can continue to read the contract. These internal
#' fields are not part of the preferred user-facing notation.
#'
#' @return An object of class `"tidyact_life_contract"`.
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x = 40:90,
#'   lx = round(100000 * exp(-0.018 * (0:50)^1.35))
#' )
#' lt$lx[nrow(lt)] <- 0
#'
#' life_contract(
#'   lt = lt,
#'   lives = "single",
#'   x = 40,
#'   i = 0.05
#' )
#'
#' life_contract(
#'   lt = lt,
#'   lives = "joint",
#'   x = 60,
#'   y = 58,
#'   i = 0.05
#' )
#'
#' @export
life_contract <- function(
    lt,
    lives = c("single", "joint", "last_survivor"),
    x = NULL,
    y = NULL,
    i,
    i_type = "effective",
    m = 1L,
    ...
) {
  lives <- match.arg(lives)

  extra <- list(...)

  if (length(extra) > 0L) {
    deprecated <- intersect(
      names(extra),
      c("mortality_table", "age", "age_x", "age_y", "rate", "rate_type")
    )

    if (length(deprecated) > 0L) {
      stop(
        "Deprecated argument name(s): ",
        paste(sprintf("`%s`", deprecated), collapse = ", "),
        ". Use compact actuarial notation: `lt`, `x`, `y`, `i`, and `i_type`.",
        call. = FALSE
      )
    }

    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", names(extra)), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (missing(lt)) {
    stop("`lt` must be provided.", call. = FALSE)
  }

  validate_one_life_table <- function(tab, label = "lt") {
    if (!is.data.frame(tab)) {
      stop(
        "`", label, "` must be a data frame or tibble.",
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

  if (identical(lives, "single")) {
    validate_one_life_table(lt, "lt")
  } else {
    if (is.data.frame(lt)) {
      validate_one_life_table(lt, "lt")
    } else if (
      is.list(lt) &&
      length(lt) == 2L &&
      all(vapply(lt, is.data.frame, logical(1L)))
    ) {
      validate_one_life_table(lt[[1L]], "lt[[1]]")
      validate_one_life_table(lt[[2L]], "lt[[2]]")
    } else {
      stop(
        "For two-life contracts, `lt` must be either one life table or ",
        "a list of two life tables `list(lt_x, lt_y)`.",
        call. = FALSE
      )
    }
  }

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  if (!is.numeric(i) ||
      length(i) != 1L ||
      is.na(i) ||
      !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) ||
      length(i_type) != 1L ||
      is.na(i_type)) {
    stop("`i_type` must be a single character string.", call. = FALSE)
  }

  valid_i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  if (!i_type %in% valid_i_type) {
    stop(
      "`i_type` must be one of: ",
      paste(sprintf("'%s'", valid_i_type), collapse = ", "),
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

  if (identical(lives, "single")) {
    if (is.null(x)) {
      stop("`x` must be provided when `lives = 'single'`.", call. = FALSE)
    }

    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        abs(x - round(x)) > 1e-10) {
      stop("`x` must be a single integer age.", call. = FALSE)
    }

    x <- as.integer(round(x))
    y <- NULL
  } else {
    if (is.null(x) || is.null(y)) {
      stop(
        "`x` and `y` must be provided for two-life contracts.",
        call. = FALSE
      )
    }

    if (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        abs(x - round(x)) > 1e-10) {
      stop("`x` must be a single integer age.", call. = FALSE)
    }

    if (!is.numeric(y) ||
        length(y) != 1L ||
        is.na(y) ||
        !is.finite(y) ||
        abs(y - round(y)) > 1e-10) {
      stop("`y` must be a single integer age.", call. = FALSE)
    }

    x <- as.integer(round(x))
    y <- as.integer(round(y))
  }

  out <- list(
    lt = lt,
    lives = lives,
    x = x,
    y = y,
    i = i,
    i_type = i_type,
    m = m,

    # Internal compatibility fields for the 0.1.4 transition.
    # These keep older contract-consuming functions working while their
    # signatures are migrated to compact actuarial notation.
    mortality_table = lt,
    age = if (identical(lives, "single")) x else NULL,
    age_x = if (!identical(lives, "single")) x else NULL,
    age_y = if (!identical(lives, "single")) y else NULL,
    rate = i,
    rate_type = i_type
  )

  class(out) <- c("tidyact_life_contract", "list")

  out
}


#' @export
print.tidyact_life_contract <- function(x, ...) {
  cat("<tidyact_life_contract>\n")
  cat("  lives:  ", x$lives, "\n", sep = "")

  if (identical(x$lives, "single")) {
    cat("  x:      ", x$x, "\n", sep = "")
  } else {
    cat("  x:      ", x$x, "\n", sep = "")
    cat("  y:      ", x$y, "\n", sep = "")
  }

  cat("  i:      ", x$i, "\n", sep = "")
  cat("  i_type: ", x$i_type, "\n", sep = "")
  cat("  m:      ", x$m, "\n", sep = "")

  if (is.list(x$lt) && !is.data.frame(x$lt)) {
    cat("  tables: ", length(x$lt), "\n", sep = "")
  } else {
    cat("  tables: 1\n", sep = "")
  }

  invisible(x)
}


.as_life_contract <- function(x) {
  inherits(x, "tidyact_life_contract")
}
