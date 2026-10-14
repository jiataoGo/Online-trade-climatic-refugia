# Figure 2: saved-model predictions.

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
my_main_theme <- my_theme + theme(legend.position = "none")
dat_m1 <- predict_for_figure(
  best_model,
  terms = "C1_NatStatus [-0.6666667, 0.3333333]",
  condition = c(C2_Invasiveness = 0)
)
dat_m1$x <- factor(dat_m1$x,
                   labels = c("Non-naturalized", "Naturalized"))
p6 <- ggplot(dat_m1, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete() +
  scale_fill_aaas() +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  labs(title = "Naturalization: ***", x = NULL, y = "Proportion of\nrefugium records (%)") +
  my_main_theme

dat_m2 <- predict_for_figure(
  best_model,
  terms = "C2_Invasiveness [-0.5, 0.5]",
  condition = c(C1_NatStatus = 1/3)
)
p7 <- ggplot(dat_m2, aes(x = factor(x), y = predicted, fill = factor(x))) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete(labels = c("-0.5" = "Non-invasive", "0.5" = "Invasive")) +
  scale_fill_aaas() +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  labs(title = "Invasiveness: ***", x = NULL, y = "Proportion of\nrefugium records (%)") +
  my_main_theme

dat_m3 <- predict_for_figure(best_model, terms = "propagule_type")
dat_m3$x <- factor(dat_m3$x, levels = c("seed", "adult"))

p8 <- ggplot(dat_m3, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete(labels = c("seed" = "Seeds", "adult" = "Live plants")) +
  scale_fill_aaas() +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  labs(title = "Sale form: ns", x = NULL, y = "Proportion of\nrefugium records (%)") +
  my_main_theme

dat_m4 <- predict_for_figure(best_model, terms = "year_label")
p9 <- ggplot(dat_m4, aes(x = x, y = predicted,  fill = x)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_fill_aaas() +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  labs(title = "Year: ns", x = NULL, y = "Proportion of\nrefugium records (%)") +
  my_main_theme

aaas_blue <- pal_aaas("default")(10)[1]
dat_m5 <- predict_for_figure(best_model, terms = "wcg_lat_scaled")
p10 <- plot(dat_m5) +
  geom_line(color = aaas_blue, linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = aaas_blue, alpha = 0.15) +
  labs(title = "Latitude: ***",
       x = "Standardized latitude",
       y = "Proportion of\nrefugium records (%)") +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  my_main_theme

dat_m6 <- predict_for_figure(best_model, terms = "wcg_lng_scaled")
p11 <- plot(dat_m6) +
  geom_line(color = aaas_blue, linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = aaas_blue, alpha = 0.15) +
  labs(title = "Longitude: ***",
       x = "Standardized longitude",
       y = "Proportion of\nrefugium records (%)") +
  scale_y_continuous(labels = function(x) x * 100, limits = c(0, 1)) +
  my_main_theme

plot_mian <- wrap_plots(p6, p7, p8, p9, p10, p11, nrow = 3) +
  plot_annotation(tag_levels = 'A')




save_figure <- function(plot, basename, width, height, expected_axes) {
  prefix <- file.path(output_dir, basename)
  ggsave(paste0(prefix, ".png"), plot, device = ragg::agg_png,
         width = width, height = height, units = "in", dpi = 600, bg = "white")
  ggsave(paste0(prefix, ".pdf"), plot, device = grDevices::cairo_pdf,
         width = width, height = height, units = "in", family = "serif", bg = "white")

}

save_figure(plot_mian, "Figure_2", 8, 10, 6L)

