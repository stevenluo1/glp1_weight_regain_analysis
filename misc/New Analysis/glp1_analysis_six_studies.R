library(tidyverse)
library(nlme)
library(ggplot2)
library(ggsci)

study_data <- read_csv("data_set_export.csv")

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


# Add predictions to the dataset
study_data <- study_data %>%
  mutate(predicted = predict(model))


# Generate a fine sequence of 'week' values
week_grid <- data.frame(
  wks_after_cessation = seq(0, 60, length.out = 600)
)


# Fixed effects predictions (smooth curve)
fixed_effects_grid <- week_grid %>%
  mutate(predicted_fixed = predict(model, newdata = week_grid, level = 0))



# Random effects predictions (study-specific)
random_effects_grid <- expand.grid(
  wks_after_cessation = seq(0, 60, by = 0.1),
  study = unique(study_data$study)
) %>%
  mutate(predicted_random = predict(model, newdata = ., level = 1))





# Draw graph

okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#E6B800", # yellow; darker than default
  "#0072B2", # blue
  "#D55E00", # vermilion
  "#CC79A7", # reddish purple
  "#000000"  # black
)


ggplot(study_data, aes(x = wks_after_cessation, y = weight_regain_pct, color = study)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), breaks = seq(0, 100, by = 10)) +
  scale_x_continuous(limits = c(0, 60), expand = c(0, 0), breaks = seq(0, 100, by = 4)) +
  
  # Fixed effect curve
  geom_line(data = fixed_effects_grid,
            aes(x = wks_after_cessation, y = predicted_fixed),
            color = "black", size = 1.6, alpha = 1) +
  
  
  # Observed data points
  geom_point(size = 2.5, alpha = 0.5, stroke = 1, aes(shape = study, fill = study)) +
  
  # Error bars
  geom_errorbar(aes(
    ymin = pmax(weight_regain_pct - weight_regain_sem * 1.96, 0),
    ymax = weight_regain_pct + weight_regain_sem * 1.96
  ), width = 0.4, alpha = 0.6, show.legend = FALSE) +
  
  
  # Study-specific curves (random effects)
  geom_line(data = random_effects_grid,
            aes(x = wks_after_cessation, y = predicted_random, group = study, linetype = study), size = 0.8, alpha = 0.8) +
  

  labs(
    x = "Weeks after cessation",
    y = "Weight regained (% of weight lost)",
  ) +
  guides(
    color = guide_legend("Trial"),
    shape = guide_legend("Trial"),
    fill = guide_legend("Trial"),
    linetype = guide_legend("Trial"),
  ) +
  scale_color_manual(values = okabe_ito, breaks = c("STEP 1", "SCALE Diabetes", "SCALE Obesity", "STEP 10", "SURMOUNT-4", "STEP 4", "SURMOUNT-1")) +
  scale_fill_manual(values = okabe_ito, breaks = c("STEP 1", "SCALE Diabetes", "SCALE Obesity", "STEP 10", "SURMOUNT-4", "STEP 4", "SURMOUNT-1")) +
  scale_shape_manual(values = c(0, 1, 21, 22, 24, 23, 25), breaks = c("STEP 1", "SCALE Diabetes", "SCALE Obesity", "STEP 10", "SURMOUNT-4", "STEP 4", "SURMOUNT-1")) + 
  scale_linetype_manual(values = c("longdash", "dotted", "dotdash", "31", "twodash", "dashed"), breaks = c("STEP 1", "SCALE Obesity", "STEP 10", "SURMOUNT-4", "STEP 4", "SURMOUNT-1")) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.width = unit(2.5, "cm"),  # Adjust width to taste
    legend.position = c(0.16, 0.82),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 10)  # top, right, bottom, left
  )


# RMSE for fixed-effect predictions
study_data <- study_data %>%
  mutate(predicted_fixed = predict(model, level = 0))

rmse <- sqrt(mean((study_data$weight_regain_pct - study_data$predicted_fixed)^2, na.rm = TRUE))
print(rmse)


# Extract tau (standard deviation of random effect)
tau <- as.numeric(VarCorr(model)[1, "StdDev"])

# Get tau-squared
tau_sq <- tau^2
print(tau_sq)
