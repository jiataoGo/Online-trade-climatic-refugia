# Figure 3: saved-model predictions.

suppressPackageStartupMessages({
  library(lme4)
  library(ggeffects)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
})

# Set the working directory to the repository root before running.
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)
input_dir <- "data"
model_path <- file.path("models", "Step4_best_model_Results.rds")
model_checksum_before <- tools::md5sum(model_path)
best_model <- readRDS(model_path)
df_agg <- model.frame(best_model)
stopifnot(nrow(df_agg) == 3333L)

# Generate the original prediction grid and point estimates, then explicitly
# calculate fixed-effect Wald confidence intervals on the model's logit scale.
# This avoids dependence on package versions that return response-scale SEs
# but subsequently use them as logit-scale SEs. It reproduces the original CIs.
predict_for_figure <- function(model, terms, condition = NULL) {
  pred <- ggeffects::ggpredict(model, terms = terms, condition = condition,
                              type = "fixed", ci_level = 0.95,
                              bias_correction = FALSE, verbose = FALSE)
  frame <- model.frame(model)
  grid <- frame[rep(1L, nrow(pred)), , drop = FALSE]
  focal_names <- sub("[[:space:]\\[].*$", "", terms)
  set_predictor <- function(name, values) {
    if (!name %in% names(grid)) stop("Unknown predictor: ", name)
    reference <- frame[[name]]
    if (is.factor(reference)) {
      grid[[name]] <<- factor(rep(as.character(values), length.out = nrow(pred)),
                              levels = levels(reference))
    } else {
      grid[[name]] <<- rep(as.numeric(as.character(values)), length.out = nrow(pred))
    }
  }
  constants <- attr(pred, "constant.values")
  fixed_names <- all.vars(delete.response(terms(model)))
  for (name in intersect(names(constants), fixed_names)) {
    set_predictor(name, constants[[name]])
  }
  if (length(condition)) {
    for (name in names(condition)) set_predictor(name, condition[[name]])
  }
  focal_columns <- c("x", "group", "facet")
  for (i in seq_along(focal_names)) set_predictor(focal_names[i], pred[[focal_columns[i]]])
  X <- model.matrix(delete.response(terms(model)), data = grid,
                    contrasts.arg = attr(lme4::getME(model, "X"), "contrasts"))
  X <- X[, names(lme4::fixef(model)), drop = FALSE]
  eta <- as.vector(X %*% lme4::fixef(model))
  stopifnot(max(abs(plogis(eta) - pred$predicted)) < 1e-7)
  V <- as.matrix(vcov(model))
  se_link <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  pred$std.error <- se_link
  pred$conf.low <- plogis(eta - qnorm(0.975) * se_link)
  pred$conf.high <- plogis(eta + qnorm(0.975) * se_link)
  attr(pred, "std.error") <- se_link
  attr(pred, "source_grid") <- grid[, fixed_names, drop = FALSE]
  stopifnot(all(pred$conf.low <= pred$predicted),
            all(pred$conf.high >= pred$predicted))
  pred
}

my_theme <- theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "white", linetype = "solid"),
    text = element_text(family = "serif"),
    plot.title = element_text(family = "serif", colour = "black", size = 14),
    plot.subtitle = element_text(family = "serif", colour = "black", size = 14),
    axis.title = element_text(family = "serif", colour = "black", size = 14),
    axis.text = element_text(family = "serif", colour = "black", size = 14),
    legend.position = c(1.0, 0.01),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(family = "serif", size = 14),
    legend.text = element_text(family = "serif", size = 14)
  )
dat_c1_lat <- predict_for_figure(
  best_model,
  terms = c("wcg_lat_scaled", "C1_NatStatus [-0.6666667, 0.3333333]"),
  condition = c(C2_Invasiveness = 0)
)
dat_c1_lat$group <- factor(dat_c1_lat$group,
                           labels = c("Non-naturalized", "Naturalized"))

p1 <- plot(dat_c1_lat,
           line_size = 1) +
  labs(title = "Latitude \u00D7 Naturalization: ***",
       y = "Proportion of\nrefugium records (%)", x = "Latitude (\u00B0; scaled)",
       color = "Naturalization status", fill = "Naturalization status") +
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_y_continuous(
    labels = function(x) x * 100,
    limits = c(0, 1)
  ) +
  my_theme

dat_c1_lng <- predict_for_figure(
  best_model,
  terms = c("wcg_lng_scaled", "C1_NatStatus [-0.6666667, 0.3333333]"),
  condition = c(C2_Invasiveness = 0)
)
dat_c1_lng$group <- factor(dat_c1_lng$group,
                           labels = c("Non-naturalized", "Naturalized"))

p2 <- plot(dat_c1_lng,
           line_size = 1) +
  labs(title = "Longitude \u00D7 Naturalization: ***",
       y = "Proportion of\nrefugium records (%)", x = "Longitude (\u00B0; scaled)",
       color = "Naturalization status", fill = "Naturalization status") +
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_y_continuous(
    labels = function(x) x * 100,
    limits = c(0, 1)
  ) +
  my_theme

dat_c2_lat <- predict_for_figure(
  best_model,
  terms = c("wcg_lat_scaled", "C2_Invasiveness [-0.5, 0.5]"),
  condition = c(C1_NatStatus = 1/3)
)
dat_c2_lat$group <- as.factor(dat_c2_lat$group)

p3 <- plot(dat_c2_lat, line_size = 1) +
  labs(title = "Latitude \u00D7 Invasiveness: ***",
       y = "Proportion of\nrefugium records (%)", x = "Latitude (\u00B0; scaled)",
       color = "Invasiveness status", fill = "Invasiveness status") +
  scale_color_aaas(labels = c("-0.5" = "Non-invasive", "0.5" = "Invasive")) +
  scale_fill_aaas(labels = c("-0.5" = "Non-invasive", "0.5" = "Invasive")) +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  my_theme

dat_prop_lat <- predict_for_figure(best_model, terms = c("wcg_lat_scaled", "propagule_type"))
dat_prop_lat$group <- factor(dat_prop_lat$group,
                             levels = c("seed", "adult"),
                             labels = c("Seeds", "Live plants"))

p4 <- plot(dat_prop_lat, line_size = 1) +
  labs(title = "Latitude \u00D7 Sale form: ***",
       y = "Proportion of\nrefugium records (%)", x = "Latitude (\u00B0; scaled)",
       color = "Sale form", fill = "Sale form") +
  scale_color_aaas() +
  scale_fill_aaas() +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  my_theme

plot_interaction <- wrap_plots(p1, p2, p3, p4, nrow = 2) +
  plot_annotation(tag_levels = 'A')

save_figure <- function(plot, basename, width, height, expected_axes) {
  prefix <- file.path(output_dir, basename)
  ggsave(paste0(prefix, ".png"), plot, device = ragg::agg_png,
         width = width, height = height, units = "in", dpi = 600, bg = "white")
  ggsave(paste0(prefix, ".pdf"), plot, device = grDevices::cairo_pdf,
         width = width, height = height, units = "in", family = "serif", bg = "white")

}

save_figure(plot_interaction, "Figure_3", 8, 8, 4L)

