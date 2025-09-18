fit_nlme_exp_recovery <- function(dat, start = c(A = 100, k = 0.1), method = "REML") {
  nlme::nlme(
    model   = weight_regain_pct ~ A * (1 - exp(-k * wks_after_cessation)),
    data    = dat,
    fixed   = A + k ~ 1,
    random  = k ~ 1 | study,
    weights = nlme::varFixed(~ weight_regain_sem^2),
    start   = start,
    method  = method
  )
}