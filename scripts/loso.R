library(tidyverse)
library(nlme)
library(ggplot2)
source("R/calculate_weight_regain.R")
source("R/fit_nlme_exp_recovery.R")

loso_nlme <- function(fit,
                      data,
                      study_col = "study",
                      sem_col   = "weight_regain_sem",
                      response  = all.vars(formula(fit))[1],
                      pred_fun  = function(fit, newdata) predict(fit, newdata, level = 0),
                      parallel  = FALSE) {

  # studies and book-keeping
  studies <- unique(data[[study_col]])
  full_fx <- fixef(fit)


  results <- lapply(studies, function(s) {
    idx <- which(data[[study_col]] == s)
    fit_s <- fit_nlme_exp_recovery(data[-idx, ])

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
    confidence_intervals <- intervals(fit_s)
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
      A_ci_lo    = confidence_intervals$fixed["A", "lower"],
      A_ci_hi    = confidence_intervals$fixed["A", "upper"],
      k          = unname(fx["k"]),
      k_ci_lo    = confidence_intervals$fixed["k", "lower"],
      k_ci_hi    = confidence_intervals$fixed["k", "upper"],
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

# Import data and derive weight regain
data_raw <- read_csv("data/data_sheet.csv")
data_raw <- data_raw %>% mutate(study = as.factor(study))
data_regain <- calculate_weight_regain(data_raw)

# Fit model
model <- fit_nlme_exp_recovery(data_regain)


loso <- loso_nlme(model, data_regain, study_col = "study", sem_col = "weight_regain_sem")
loso$table         # per-study estimates, deltas, and weighted CV error
loso$press         # overall PRESS-like metric
loso$full_fixef    # baseline A, k



A_hat <- unname(loso$full_fixef[1])
k_hat <- unname(loso$full_fixef[2])


loso_forest_plot_A <- ggplot(loso$table, aes(A, reorder(study, A))) +
  theme_classic(base_size = 12) +
  geom_vline(xintercept = A_hat, color="darkred") +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbarh(aes(xmin = A_ci_lo, xmax = A_ci_hi),
                 height = 0, color = "darkgreen") +
  labs(x = "A with study omitted", y = "Omitted study")

ggsave("figures/loso_forest_A.png", plot = loso_forest_plot_A, width = 8, height = 6, dpi = 300)



loso_forest_plot_k <- ggplot(loso$table, aes(k, reorder(study, k))) +
  theme_classic(base_size = 12) +
  geom_vline(xintercept = k_hat, color="darkred") +
  geom_point(size = 3, color = "darkgreen") +
  geom_errorbarh(aes(xmin = k_ci_lo, xmax = k_ci_hi),
                 height = 0, color = "darkgreen") +
  labs(x = "k with study omitted", y = "Omitted study")

ggsave("figures/loso_forest_k.png", plot = loso_forest_plot_k, width = 8, height = 6, dpi = 300)
