library(tidyverse)
library(nlme)
library(ggplot2)
library(ggsci)

study_data <- read_csv("data_sheet.csv")

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
  wks_after_cessation = seq(0, 60, length.out = 600),
  study = unique(study_data$study)
) %>%
  mutate(predicted_random = predict(model, newdata = ., level = 1))




loso_nlme <- function(fit,
                      data,
                      study_col = "study",
                      sem_col   = "weight_regain_sem",
                      response  = all.vars(formula(fit))[1],
                      pred_fun  = function(fit, newdata) predict(fit, newdata, level = 0),
                      parallel  = FALSE) {
  # helper to refit after dropping a set of rows
  refit_after_drop <- function(idx_drop) {
    try(update(fit, data = data[-idx_drop, ], start = fixef(fit)), silent = TRUE)
  }
  
  # studies and book-keeping
  studies <- unique(data[[study_col]])
  full_fx <- fixef(fit)
  
  
  results <- lapply(studies, function(s) {
    idx <- which(data[[study_col]] == s)
    fit_s <- refit_after_drop(idx)
    
    if (inherits(fit_s, "try-error")) {
      return(tibble(
        study = s, ok = FALSE,
        A = NA_real_, k = NA_real_, dA = NA_real_, dk = NA_real_,
        logLik = NA_real_, AIC = NA_real_,
        n_left_out = length(idx),
        cv_err = NA_real_
      ))
    }
    
    fx <- fixef(fit_s)
    # prediction for left-out rows (population level, no BLUPs)
    preds <- pred_fun(fit_s, data[idx, , drop = FALSE])
    
    # weighted squared error using SEM (matches your varFixed(~ SEM^2))
    sem  <- data[[sem_col]][idx]
    obs  <- data[[response]][idx]
    werr <- sum(((obs - preds)^2) / (sem^2))
    
    tibble(
      study      = s,
      ok         = TRUE,
      A          = unname(fx["A"]),
      k          = unname(fx["k"]),
      dA         = unname(fx["A"] - full_fx["A"]),
      dk         = unname(fx["k"] - full_fx["k"]),
      logLik     = as.numeric(logLik(fit_s)),
      AIC        = AIC(fit_s),
      n_left_out = length(idx),
      cv_err     = werr
    )
  }) %>% list_rbind() %>% arrange(study)
  
  # overall PRESS-like metric across left-out studies
  press <- sum(results$cv_err[results$ok], na.rm = TRUE)
  
  list(
    full_fixef = unclass(full_fx),
    table = results,
    press = press
  )
}


loso <- loso_nlme(model, study_data, study_col = "study", sem_col = "weight_regain_sem")
loso$table         # per-study estimates, deltas, and weighted CV error
loso$press         # overall PRESS-like metric
loso$full_fixef    # baseline A, k


library(ggplot2)

ggplot(loso$table, aes(reorder(study, dA), dA)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() + coord_flip() +
  labs(x = "Study left out", y = "ΔA (omit study)", title = "LOSO influence on A")

ggplot(loso$table, aes(reorder(study, dk), dk)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() + coord_flip() +
  labs(x = "Study left out", y = "Δk (omit study)", title = "LOSO influence on k")
