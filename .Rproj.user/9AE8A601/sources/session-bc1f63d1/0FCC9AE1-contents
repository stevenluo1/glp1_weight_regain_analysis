library(tidyverse)
library(ggplot2)
source("R/calculate_weight_regain.R")
source("R/fit_nlme_exp_recovery.R")
source("R/plot_weight_regain.R")

monte_carlo <- function(fit, dat, B = 100, digitize_error_sd) {
  fe <- fixef(fit); A_hat <- unname(fe["A"]); k_hat <- unname(fe["k"])
  
  dat <- transform(dat, study = as.factor(study))
  
  draws <- matrix(NA_real_, nrow=B, ncol=2, dimnames=list(NULL, c("A","k")))
  start_vals <- c(A=A_hat, k=k_hat)
  
  skipped <- integer(0)   # store failed iteration numbers
  
  for (b in seq_len(B)) {
    
    
    # Add jitter to the raw extracted points
    delta_weight_mean_observed <- dat$delta_weight_mean
    delta_weight_mean_jitter <- delta_weight_mean_observed + rnorm(length(delta_weight_mean_observed), mean = 0, sd = digitize_error_sd)
    
    
    # Put jittered points in data set
    dat_jitter <- dat
    dat_jitter$delta_weight_mean <- delta_weight_mean_jitter
    
    dat_calc <- calculate_weight_regain(dat_jitter)
    
    #    print(head(dat_calc))
    
    # attempt refit
    fit_b <- tryCatch(
      nlme(
        weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
        data    = dat_calc,
        fixed   = A + k ~ 1,
        random  = k ~ 1 | study,
        weights = varFixed(~ weight_regain_sem^2),
        start   = start_vals
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_b)) {
      skipped <- c(skipped, b)
      next
    }
    
    fe_b <- fixef(fit_b)
    if (anyNA(fe_b)) {
      skipped <- c(skipped, b)
      next
    }
    
    cat("Iter", b, ":  A =", fe_b["A"], " k =", fe_b["k"], "\n")
    
    draws[b, ] <- unname(fe_b)
    start_vals <- fe_b
  }
  
  
  if (length(skipped) > 0) {
    message("Skipped iterations: ", paste(skipped, collapse = ", "))
  } else {
    message("No iterations skipped")
  }  
  
  
  as.data.frame(draws)
}

# Import data and derive weight regain
data_raw <- read_csv("data/data_sheet.csv")
data_raw <- data_raw %>% mutate(study = as.factor(study))
data_regain <- calculate_weight_regain(data_raw)

# Fit model
model <- fit_nlme_exp_recovery(data_regain, method="ML")
summary(model)

digitize_error_sd <- 0.05 # 0.05 is a conservative overestimate

set.seed(42)
monte_carlo_data <- monte_carlo(model, data_raw, B = 1000, digitize_error_sd)
summary(monte_carlo_data)

apply(monte_carlo_data, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)


# Graphing prediction bands (simulation band not visible at SD = 0.05)

x_grid <- seq(0, 60, by = 0.1)

# Create matrix where rows = bootstrap draws, columns = x_grid
pred_matrix <- sapply(seq_len(nrow(monte_carlo_data)), function(i) {
  A_i <- monte_carlo_data$A[i]; k_i <- monte_carlo_data$k[i]
  A_i * (1 - exp(-k_i * x_grid))
})

# rows = x_grid, cols = bootstrap draws
pred_matrix <- t(pred_matrix)

# 95% confidence interval
pred_mean <- apply(pred_matrix, 2, mean, na.rm = TRUE)
pred_lo   <- apply(pred_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
pred_hi   <- apply(pred_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)

pred_df <- data.frame(
  x   = x_grid,
  mean = pred_mean,
  lo   = pred_lo,
  hi   = pred_hi
)

p <- create_weight_regain_plot(data_regain, model, "study")

p <- p + geom_ribbon(data = pred_df,
                  aes(x = x_grid, ymin = lo, ymax = hi),
                  inherit.aes = FALSE,          # <-- key line
                  fill = "lightgray", alpha = 0.3)

ggsave(paste0("figures/graph_monte_carlo_error_sd_", digitize_error_sd, ".png"), plot = p, width = 8, height = 6, dpi = 300)


