library(tidyverse)
library(ggplot2)
library(nlme)
source("R/calculate_weight_regain.R")
source("R/fit_nlme_exp_recovery.R")
source("R/plot_weight_regain.R")

# Import data and derive weight regain
data_raw <- read_csv("data/data_sheet.csv")
data_raw <- data_raw %>% mutate(study = as.factor(study))
data_regain <- calculate_weight_regain(data_raw)

# Fit model
model <- fit_nlme_exp_recovery(data_regain, method="ML")
summary(model)

# Convert from proportion regained (use 0–1 scale, not %) to logit(proportion)
data_regain$proportion <- data_regain$weight_regain_pct / 100
data_regain$logit_proportion <- qlogis(data_regain$proportion)

model_logit <- nlme(
  model = logit_proportion ~ qlogis(A * (1 - exp(-k * wks_after_cessation))),
  data = data_regain,
  fixed = A + k ~ 1,
  random = k ~ 1 | study,
  start = c(A = 0.75, k = 0.03),  # use estimates from your identity model
  #  correlation = corAR1(form = ~ wks_after_cessation | study),
  weights = varFixed(~ weight_regain_sem^2)
)

fixed_curve_logit  <- make_fixed_curve(model_logit)
random_curves_logit <- make_random_curves(model_logit, "study", levels(data_regain$study))

fixed_curve_backtransformed <- fixed_curve_logit
random_curves_backtransformed <- random_curves_logit

fixed_curve_backtransformed$pred_fixed <- plogis(fixed_curve_logit$pred_fixed) * 100
random_curves_backtransformed$pred_random <- plogis(random_curves_logit$pred_random) * 100

p <- ggplot_weight_regain(data_regain, fixed_curve_backtransformed, random_curves_backtransformed, "study")

ggsave("figures/graph_logit_prop_backtransformed.png", plot = p, width = 8, height = 6, dpi = 300)

p

