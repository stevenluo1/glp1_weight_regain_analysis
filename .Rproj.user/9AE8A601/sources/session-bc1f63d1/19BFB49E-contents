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
model <- fit_nlme_exp_recovery(data_regain)


idx <- which(data_regain$study == "SURMOUNT-4")
data_drop_study <- data_regain[-idx, ]
fit_s <- fit_nlme_exp_recovery(data_regain[-idx, ])


summary(fit_s)
intervals(fit_s)