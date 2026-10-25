# Shared display helpers only: no financial valuation or mutation of amounts.
.cash_flow_display_scale <- function(x) {
  magnitude <- max(abs(x))
  exponent <- if (magnitude > 0 && (magnitude > 1e150 || magnitude < 1e-150)) {
    max(-300, min(308, floor(log10(magnitude))))
  } else 0
  list(divisor = 10^exponent, exponent = exponent)
}

.cash_flow_label <- function(x, digits, format, currency = "") {
  magnitude <- abs(x)
  valid <- is.finite(x)
  scientific <- valid & format != "full" & magnitude != 0 &
    (magnitude >= 1e15 | magnitude < 10^(-digits))
  out <- character(length(x))
  out[scientific] <- formatC(magnitude[scientific], format = "e", digits = min(digits, 6L))
  ordinary <- valid & !scientific
  compact <- ordinary & (format == "compact" |
    (format == "auto" & magnitude >= 1e6))
  fixed <- ordinary & !compact
  out[fixed] <- formatC(magnitude[fixed], format = "f", digits = digits)
  if (any(compact)) {
    exponent <- pmin(4, pmax(0, floor(log10(pmax(magnitude[compact], 1)) / 3)))
    text <- formatC(magnitude[compact] / 1000^exponent, format = "f", digits = digits)
    text <- sub("(\\.[0-9]*?)0+$", "\\1", text)
    text <- sub("\\.$", "", text)
    out[compact] <- paste0(text, c("", "K", "M", "B", "T")[exponent + 1L])
  }
  result <- paste0(ifelse(x < 0, "-", ""), currency, out)
  result[!valid] <- NA_character_
  result
}
