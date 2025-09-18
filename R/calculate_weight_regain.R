calculate_weight_regain <- function(dat) {
  
  stopifnot(all(c("study","drug","dose","wks_after_cessation",
                  "delta_weight_mean","delta_weight_sem") %in% names(dat)))
  
  eot <- dat |>
    dplyr::filter(wks_after_cessation == 0) |>
    dplyr::select(study, drug, dose,
                  eot_weight_delta = delta_weight_mean,
                  eot_weight_delta_sem = delta_weight_sem) |>
    dplyr::distinct()
  
  dat <- dat |>
    dplyr::left_join(eot, by = c("study","drug","dose")) |>
    dplyr::mutate(
      weight_regain_pct = (1 - delta_weight_mean / eot_weight_delta) * 100,
      weight_regain_sem = dplyr::if_else(
        delta_weight_mean == 0,
        sqrt((delta_weight_sem^2) / (eot_weight_delta^2)) * 100,
        (delta_weight_mean / eot_weight_delta) *
          sqrt((delta_weight_sem / delta_weight_mean)^2 +
                 (eot_weight_delta_sem / eot_weight_delta)^2) * 100
      )
    )

  out <- dplyr::filter(dat, weight_regain_pct != 0 | wks_after_cessation != 0)
  
  out
}