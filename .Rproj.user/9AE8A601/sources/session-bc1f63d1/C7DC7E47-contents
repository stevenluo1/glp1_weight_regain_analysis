# Make week grid for use in fixed and random effect curves
make_week_grid <- function(max_week=60, step=0.1) {
  tibble::tibble(wks_after_cessation = seq(0, max_week, by = step))
}

# Make curve for fixed effect
make_fixed_curve <- function(model, max_week=60, step=0.1) {
  grid <- make_week_grid(max_week, step)
  grid$pred_fixed <- predict(model, newdata=grid, level=0)
  grid
}

# Make individual group (study or arm) curves
make_random_curves <- function(model, group_var, groups, max_week=60, step=0.1) {
  grid <- expand.grid(
    wks_after_cessation = seq(0, max_week, by=step),
    group = groups
  )
  
  grid$group <- factor(grid$group, levels = groups)
  
  # rename to the actual grouping variable name
  names(grid)[names(grid) == "group"] <- group_var
  
  grid$pred_random <- predict(model, newdata=grid, level=1)
  
  tibble::as_tibble(grid)
}


default_group_levels <- c("STEP 1","SCALE Obesity","STEP 10","SURMOUNT-4", "STEP 4", "SURMOUNT-1")
default_pal <- c("#E69F00","#009E73","#E6B800","#0072B2","#D55E00","#CC79A7","#56B4E9","#000000") # okabe
default_shapes <- c(0,1,21,22,24,23,20,19)
default_linetypes <- c("longdash","dotted","dotdash","31","twodash","dashed","21","13","73")

ggplot_weight_regain <- function(data_regain, fixed_curve, random_curves, 
                                 group_var,   
                                 group_levels = default_group_levels,
                                 pal = default_pal,
                                 shapes = default_shapes,
                                 ltypes = default_linetypes,
                                 x_lim = c(0, 60), y_lim = c(0, 100), show_legend=TRUE) {

  # ensure factor order
  data_regain$group   <- factor(data_regain[[group_var]],   levels = group_levels)
  random_curves$group <- factor(random_curves[[group_var]], levels = group_levels)
  
  p <- ggplot(data_regain,
              aes(wks_after_cessation, weight_regain_pct,
                  colour = group, shape = group, fill = group, linetype = group)) +
    # fixed-effect curve
    geom_line(data = fixed_curve,
              aes(wks_after_cessation, pred_fixed),
              inherit.aes = FALSE, linewidth = 1.4, colour = "black") +
    # study-specific curves
    geom_line(data = random_curves,
              aes(wks_after_cessation, pred_random),
              linewidth = 0.8, alpha = 0.8) +
    # points
    geom_point(size = 2.5, alpha = 0.7, stroke = 1) +
    # error bars
    geom_errorbar(aes(ymin = pmax(weight_regain_pct + qnorm(0.025) * weight_regain_sem, 0),
                      ymax = weight_regain_pct + qnorm(0.975) * weight_regain_sem),
                  linetype = "solid", width = 0.6, alpha = 0.6, show.legend = FALSE) +

    scale_x_continuous(limits = x_lim, expand = c(0, 0), breaks = seq(x_lim[1], x_lim[2], by = 4)) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0), breaks = seq(y_lim[1], y_lim[2], by = 10)) +
    scale_color_manual(values = pal, limits = group_levels) +
    scale_fill_manual(values  = pal, limits = group_levels) +
    scale_shape_manual(values = shapes, na.value = 1, limits = group_levels) +
    scale_linetype_manual(values = ltypes, limits = group_levels) +
    labs(x = "Weeks after cessation", y = "Weight regained (% of weight lost)") +
    guides(colour = guide_legend("Trial"),
           shape  = guide_legend("Trial"),
           fill   = guide_legend("Trial"),
           linetype = guide_legend("Trial")) +
    theme_classic(base_size = 12)

    if (show_legend) {
      p <- p + theme(legend.position = c(0.16, 0.82),
                     legend.key.width = grid::unit(2.5, "cm"),
                     plot.margin = margin(10, 20, 10, 10))   
    }
    else {
      p <- p + theme(legend.position = "none")   
    }      
 

  p
}

create_weight_regain_plot <- function(data, model, group_var,
                                      group_levels=default_group_levels,
                                      pal=default_pal,
                                      shapes=default_shapes,
                                      ltypes=default_linetypes,
                                      max_week=60,
                                      show_legend=TRUE) {
  
  
  fixed_curve  <- make_fixed_curve(model, max_week=max_week)
  random_curves <- make_random_curves(model, group_var, levels(data[[group_var]]), max_week=max_week)
  
  p <- ggplot_weight_regain(data, fixed_curve, random_curves, group_var, group_levels, pal, shapes, ltypes, x_lim = c(0, max_week), show_legend=show_legend)
  
  p
  
}

