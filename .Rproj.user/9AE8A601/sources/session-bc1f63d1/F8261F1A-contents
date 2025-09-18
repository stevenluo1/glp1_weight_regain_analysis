library(tidyverse)
library(ggplot2)
source("R/calculate_weight_regain.R")
source("R/fit_nlme_exp_recovery.R")
source("R/plot_weight_regain.R")

# Import data and derive weight regain
data_raw <- read_csv("data/data_sheet.csv")
data_raw <- data_raw %>% mutate(study = as.factor(study))
data_regain <- calculate_weight_regain(data_raw)

# Fit model
model <- fit_nlme_exp_recovery(data_regain, method="REML")
summary(model)

# Draw graph

p <- create_weight_regain_plot(data_regain, model, "study")

ggsave("figures/graph_main.png", plot = p, width = 8, height = 6, dpi = 300)

p
