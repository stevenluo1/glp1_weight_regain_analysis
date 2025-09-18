library(tidyverse)
library(nlme)
library(ggplot2)
library(ggsci)

study_data <- read_csv("data_set_export.csv")
s <- "SURMOUNT-4"

idx_leave <- which(study_data$study == s)
study_data <- study_data[-idx_leave, ]

# Get weight lost at end of treatment (this should be the nadir)
eot_weight_loss <- study_data %>%
  filter(wks_after_cessation == 0) %>%
  select(
    study,
    drug,
    dose,
    eot_weight_delta = delta_weight_mean,
    eot_weight_delta_sem = delta_weight_sem,
  )

# Append EOT weight loss to main table 
study_data <- study_data %>%
  left_join(eot_weight_loss, by = c("study", "drug", "dose"))

# Error Propagation for calculating percentage weight regain (X/Y)
study_data <- study_data %>%
  mutate(
    # Weight regain and SEM
    weight_regain_pct = (1 - delta_weight_mean / eot_weight_delta) * 100,
    weight_regain_sem = if_else(
      delta_weight_mean == 0,
      sqrt((delta_weight_sem^2) / (eot_weight_delta^2)) * 100,
      (delta_weight_mean / eot_weight_delta) *
        sqrt((delta_weight_sem / delta_weight_mean)^2 + 
               (eot_weight_delta_sem / eot_weight_delta)^2) * 100
    )
  )

# Remove redundant (0,0) points
study_data <- study_data %>%
  filter(weight_regain_pct != 0 | wks_after_cessation != 0)

# Define the model function
exp_recovery <- function(week, A, k) {
  A * (1 - exp(-k * week))
}


# Fit the non-linear mixed-effects model
model <- nlme(
  model = weight_regain_pct ~ exp_recovery(wks_after_cessation, A, k),
  data = study_data,
  fixed = A + k ~ 1,
  random = k ~ 1 | study,
  weights = varFixed(~ weight_regain_sem^2),
  start = c(A = 100, k = 0.1)
)

summary(model)