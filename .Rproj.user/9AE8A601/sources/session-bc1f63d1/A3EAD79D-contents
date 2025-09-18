library(tidyverse)
library(ggplot2)
source("R/calculate_weight_regain.R")
source("R/plot_weight_regain.R")

# Import data and derive weight regain
data_raw <- read_csv("data/data_sheet.csv")
data_raw <- data_raw %>% mutate(study = as.factor(study))
data_regain <- calculate_weight_regain(data_raw)

model_rand_k <- nlme(
  model = weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
  data = data_regain,
  fixed = A + k ~ 1,
  random = k ~ 1 | study,
  weights = varFixed(~ weight_regain_sem^2),
  start = c(A = 100, k = 0.1)
)

model_rand_A <- nlme(
  model = weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
  data = data_regain,
  fixed = A + k ~ 1,
  random = A ~ 1 | study,
  weights = varFixed(~ weight_regain_sem^2),
  start = c(A = 100, k = 0.1)
)

model_rand_Ak <- nlme(
  model = weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
  data = data_regain,
  fixed = A + k ~ 1,
  random = A + k ~ 1 | study,
  weights = varFixed(~ weight_regain_sem^2),
  start = c(A = 100, k = 0.1)
)

summary(model_rand_k)
summary(model_rand_A)
summary(model_rand_Ak) 

ranef(model_rand_Ak) # random effects on A collapse to near zero


# Draw graph

p_Ak <- create_weight_regain_plot(data_regain, model_rand_Ak, "study")
p_A <- create_weight_regain_plot(data_regain, model_rand_A, "study")
p_k <- create_weight_regain_plot(data_regain, model_rand_k, "study")


ggsave("figures/graph_random_effects_A.png", plot = pA, width = 8, height = 6, dpi = 300)

