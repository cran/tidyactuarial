# Reproducible financial and actuarial figures. No random sampling required.
# Run after installing/loading tidyactuarial; optional first argument is output folder.
library(tidyactuarial)
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1L]] else file.path(tempdir(), "cash-flow-visualizations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

financial <- tibble::tibble(t = c(0, 1, 2, 2, 3, 4),
  C = c(-10000, 2500, 4000, -500, 4500, 3000))
p_financial <- financial |> plot_cash_flow(C = C, t = t, i = .06, currency = "$",
  title = "Investment cash flows", subtitle = NULL)

dated <- tibble::tibble(date = as.Date("2026-01-01") + 30 * (0:59),
  C = c(-50000, rep(c(1500, 2000, 2500), length.out = 59)))
p_dense <- dated |> plot_cash_flow(C = C, date = date, i = .05, show_labels = FALSE,
  currency = "$", title = "A dense dated payment stream", date_labels = "%b %Y")

lt <- tibble::tibble(x = 60:65, lx = c(1000, 800, 600, 300, 100, 0))
annuity <- apv_life_flow(lt, ages = 60, t = 1:4, cf = rep(1200, 4), i = .04)
p_annuity <- annuity |> plot_contingent_cash_flow(currency = "$",
  title = "Temporary life annuity: four annual payments")

# Death during (j-1,j], paid at j, relative to a life alive at age 60.
death <- tibble::tibble(t = 1:4, cf = rep(10000, 4),
  probability = c(200, 200, 300, 200) / 1000,
  discount = 1.04^-(1:4), event = "death_interval")
p_death <- death |> plot_contingent_cash_flow(currency = "$",
  title = "Term insurance: benefits paid at year-end of death")

reversionary <- apv_life_flow(lt, ages = c(60, 61), t = 1:3,
  cf = rep(1200, 3), i = .04, status = "reversionary", alpha = .5)
p_reversionary <- reversionary |> plot_contingent_cash_flow(currency = "$",
  title = "Joint-and-survivor annuity: 50% after first death")

mixed <- tibble::tibble(t = c(1, 1, 2, 2), cf = c(1200, 10000, 1200, 10000),
  probability = c(.8, .2, .6, .2), discount = 1.04^-c(1, 1, 2, 2),
  event = rep(c("survival", "death_interval"), 2))
p_mixed <- mixed |> plot_contingent_cash_flow(currency = "$",
  title = "Survival and death benefits on the same payment dates")

p_extreme <- plot_cash_flow(C = c(-1.7e308, 1e308, 8e307), t = 0:2,
  title = "Extreme monetary amounts: finite graphical coordinates")

plots <- list(financial = p_financial, dense_dates = p_dense,
  annuity = p_annuity, death = p_death, reversionary = p_reversionary,
  simultaneous = p_mixed, extreme = p_extreme)
for (name in names(plots)) {
  for (extension in c("png", "pdf")) {
    ggplot2::ggsave(file.path(output_dir, paste0(name, ".", extension)),
      plots[[name]], width = 11, height = if (name %in% c("financial", "dense_dates", "extreme"))
        5.5 else 8, units = "in", dpi = 150, bg = "white")
  }
}
message("Figures saved to: ", normalizePath(output_dir, winslash = "/"))
