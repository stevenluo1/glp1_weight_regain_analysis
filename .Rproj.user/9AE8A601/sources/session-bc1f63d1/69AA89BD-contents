library(tidyverse)
library(nlme)
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

tau <- as.numeric(VarCorr(model)[1, "StdDev"])
A_hat <- unname(fixef(model)["A"])
k_hat <- unname(fixef(model)["k"])

max_weeks <- 60
x_grid <- seq(0, max_weeks, by = 0.1)


pred_df <- data.frame(
  wks_after_cessation = x_grid,
  mean = A_hat * (1 - exp((-k_hat) * x_grid)),
  lo   = A_hat * (1 - exp( -(k_hat + qnorm(0.025) * tau) * x_grid) ),
  hi   = A_hat * (1 - exp( -(k_hat + qnorm(0.975) * tau) * x_grid) )
)

p <- create_weight_regain_plot(data_regain, model, "study")

p <- p + geom_ribbon(data = pred_df,
              aes(x = wks_after_cessation, ymin = lo, ymax = hi),
              inherit.aes = FALSE,          # <-- key line
              fill = "lightgray", alpha = 0.3) 

p
  
ggsave("figures/graph_95pc_population_level_prediction_band.png", plot = p, width = 8, height= 6, dpi = 300)
