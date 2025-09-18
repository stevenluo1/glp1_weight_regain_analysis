library(tidyverse)
library(ggplot2)
library(metafor)
source("R/calculate_weight_regain.R")

# Import data and derive weight regain
data_raw <- read_csv("data/data_sheet.csv")
data_raw <- data_raw %>% mutate(study = as.factor(study))
data_regain <- calculate_weight_regain(data_raw)

get_closest_within_window <- function(data, target_week, window) {
  data %>%
    filter(abs(wks_after_cessation - target_week) <= window) %>%
    group_by(study) %>%
    slice_min(abs(wks_after_cessation - target_week), with_ties = FALSE) %>%
    ungroup() %>%
    mutate(target = target_week)
}

dat12 <- get_closest_within_window(data_regain, 12, 2)
dat26 <- get_closest_within_window(data_regain, 26, 2)
dat52 <- get_closest_within_window(data_regain, 52, 4)

res12 <- rma(yi = weight_regain_pct, sei = weight_regain_sem, data = dat12, method = "REML")
res26 <- rma(yi = weight_regain_pct, sei = weight_regain_sem, data = dat26, method = "REML")
res52 <- rma(yi = weight_regain_pct, sei = weight_regain_sem, data = dat52, method = "REML")

png("figures/forest_12weeks.png", width = 2000, height = 1200, res = 300)
forest(res12,
       slab = dat12$study,
       alim = c(0, 100),     # x-axis range for data
       xlim = c(-100, 200),      # pull plot region leftward
       xlab = "% weight regain",
       main = "12 weeks",
       refline = res12$b)
dev.off()

png("figures/forest_26weeks.png", width = 2000, height = 1200, res = 300)
forest(res26,
       slab = dat26$study,
       alim = c(0, 100),     # x-axis range for data
       xlim = c(-100, 200),      # pull plot region leftward
       xlab = "% weight regain",
       main = "26 weeks",
       refline = res26$b)
dev.off()

png("figures/forest_52weeks.png", width = 2000, height = 1200, res = 300)
forest(res52,
       slab = dat52$study,
       alim = c(0, 100),     # x-axis range for data
       xlim = c(-100, 200),      # pull plot region leftward
       xlab = "% weight regain",
       main = "52 weeks",
       refline = res52$b)
dev.off()
