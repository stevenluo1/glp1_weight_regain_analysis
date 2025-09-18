library(tidyverse)
library(nlme)
library(ggplot2)
source("R/calculate_weight_regain.R")
source("R/fit_nlme_exp_recovery.R")
source("R/plot_weight_regain.R")


# Import data and derive weight regain
data_raw <- read_csv("data/data_sheet.csv")

# AR model requires only one data point at each time point per study; SURMOUNT-1 violates this so must regroup by arm
data_raw <- data_raw %>%
  mutate(
    study = factor(study),
    arm   = interaction(study, drug, dose, drop = TRUE)  # one cohort/arm per group
  )

levels(data_raw$arm) <- c("SURMOUNT-1 10 mg", "SURMOUNT-1 15 mg", "SURMOUNT-4", "STEP 1", "STEP 10", "STEP 4", "SCALE Obesity", "SURMOUNT-1 5 mg")

data_regain <- calculate_weight_regain(data_raw)

m_iid <- nlme(
  model  = weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
  data   = data_regain,
  fixed  = A + k ~ 1,
  random = k ~ 1 | arm,                       
  weights = varFixed(~ weight_regain_sem^2),
  start  = c(A = 100, k = 0.1),
  control = nlmeControl(msMaxIter = 200, pnlsTol = 1e-6)
)

m_car1 <- nlme(
  model  = weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
  data   = data_regain,
  fixed  = A + k ~ 1,
  random = k ~ 1 | arm,                       
  weights = varFixed(~ weight_regain_sem^2),
  correlation = corCAR1(form = ~ wks_after_cessation | arm),  # <- key fix
  start  = c(A = 100, k = 0.1),
  control = nlmeControl(msMaxIter = 200, pnlsTol = 1e-6)
)

m_cs <- nlme(
  model  = weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
  data   = data_regain,
  fixed  = A + k ~ 1,
  random = k ~ 1 | arm,
  weights = varFixed(~ weight_regain_sem^2),
  correlation = corCompSymm(form = ~ 1 | arm),
  start  = c(A = 100, k = 0.1)
)


# Helper to pull fixed effects, CIs, sigma, AIC
summ <- function(fit) {
  fe  <- fixef(fit)
  ci  <- intervals(fit)$fixed[, c("lower","est.","upper")]
  data.frame(
    A_est = fe["A"], A_lo = ci["A","lower"], A_hi = ci["A","upper"],
    k_est = fe["k"], k_lo = ci["k","lower"], k_hi = ci["k","upper"],
    sigma = sigma(fit), AIC = AIC(fit)
  )
}

# Correlation parameters: CAR(1) phi → rho(1wk), half-life; CS rho
car1_phi_ci <- tryCatch(as.data.frame(intervals(m_car1)$corStruct)[1,], error = function(e) NULL)
phi  <- if (!is.null(car1_phi_ci)) car1_phi_ci$`est.` else NA_real_
rho1 <- if (!is.na(phi)) exp(-phi*1) else NA_real_
hl   <- if (!is.na(phi)) log(2)/phi else NA_real_

cs_ci <- tryCatch(as.data.frame(intervals(m_cs)$corStruct)[1,], error = function(e) NULL)
rho_cs <- if (!is.null(cs_ci)) cs_ci$`est.` else NA_real_

# Sensitivity table
sens_tbl <- bind_rows(
  cbind(model = "iid",          summ(m_iid)),
  cbind(model = "CAR1 (time)",  summ(m_car1)),
  cbind(model = "CS",           summ(m_cs))
) %>% mutate(
  across(c(A_est:AIC), ~ signif(., 3))
)

print(sens_tbl)

# Brief correlation summary to report (optional print)
car1_summary <- data.frame(
  model = "CAR1 (time)",
  phi_est = signif(phi, 3),
  rho_1wk = signif(rho1, 3),
  half_life_weeks = signif(hl, 3)
)
cs_summary <- data.frame(
  model = "CS",
  rho_cs = signif(rho_cs, 3)
)

print(car1_summary)
print(cs_summary)

group_levels = c("STEP 1","SCALE Obesity","STEP 10","SURMOUNT-4","STEP 4","SURMOUNT-1 5 mg","SURMOUNT-1 10 mg","SURMOUNT-1 15 mg")

data_regain$arm

p <- create_weight_regain_plot(data=data_regain, m_car1, "arm", group_levels)
ggsave("figures/graph_car1.png", plot = p, width = 8, height = 6, dpi = 300)
p

p <- create_weight_regain_plot(data_regain, m_cs, "arm", group_levels)
ggsave("figures/graph_compound_symmetry.png", plot = p, width = 8, height = 6, dpi = 300)
