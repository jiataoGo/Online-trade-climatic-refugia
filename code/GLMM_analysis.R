# GLMM fitting, AICc selection, diagnostics and likelihood-ratio tests.
# Run from the repository root. This script fits models from the trade CSV.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lme4)
  library(MuMIn)
  library(performance)
})
dir.create("models", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
options(na.action = "na.fail")

# 1. Aggregate records and construct the original orthogonal contrasts.
df <- read_csv("data/Step4_Refugia_Full_Filtered_NotVeg.csv", show_col_types = FALSE)
glimpse(df)
df_agg <- df %>%
  group_by(shop_id, wcg_lat, wcg_lng, year_label, invasion_label, propagule_type) %>%
  summarise(
    refugia_count = sum(is_refugia),
    total_count = n(),
    ratio = refugia_count / total_count,
    not.refugia_count = total_count - refugia_count,
    .groups = 'drop'
  ) %>%
  mutate(
    shop_id = as.factor(shop_id),
    year_label = as.factor(year_label),
    invasion_label = as.factor(invasion_label),
    propagule_type = as.factor(propagule_type)
  )
cat("Aggregated Data Summary:\n")
print(head(df_agg))
write_csv(df_agg, "results/Step4_df_agg.csv")
mat_comp <- cbind(
  c(1/3, 1/3, 1/3),
  c(-1, 0.5, 0.5),
  c(0, -1, 1)
)
mymat_comp <- solve(t(mat_comp))
contrast_comp <- as.data.frame(mymat_comp[, 2:3])
contrast_comp$invasion_label <- c("nonnat", "nat", "inv")
names(contrast_comp)[1:2] <- c("C1_NatStatus", "C2_Invasiveness")
df_agg <- df_agg %>% left_join(contrast_comp, by = "invasion_label")
glimpse(df_agg)
df_agg$wcg_lat_scaled <- scale(df_agg$wcg_lat)
df_agg$wcg_lng_scaled <- scale(df_agg$wcg_lng)
strict_control <- glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 2e5))
model_interaction_full <- glmer(cbind(refugia_count, not.refugia_count) ~
                                  wcg_lat_scaled +
                                  wcg_lng_scaled +
                                  C1_NatStatus +
                                  C2_Invasiveness +
                                  propagule_type +
                                  year_label +
                                  wcg_lat_scaled:C1_NatStatus +
                                  wcg_lat_scaled:C2_Invasiveness +
                                  wcg_lat_scaled:propagule_type +
                                  wcg_lat_scaled:year_label +
                                  wcg_lng_scaled:C1_NatStatus +
                                  wcg_lng_scaled:C2_Invasiveness +
                                  wcg_lng_scaled:propagule_type +
                                  wcg_lng_scaled:year_label +
                                  (1 | shop_id),
                                data = df_agg,
                                family = binomial,
                                control = strict_control)

stopifnot(nrow(df) == 21351L, nrow(df_agg) == 3333L)
write_csv(mutate(df_agg, across(ends_with("_scaled"), as.numeric)),
          "results/Step4_model_data.csv")

# 2. Fit every hierarchical candidate model and rank by AICc.
# No saved model or selection table is read during this workflow.
available_cores <- parallel::detectCores()
if (is.na(available_cores)) available_cores <- 1L
n_cores <- max(1L, min(6L, available_cores - 2L))
worker_cluster <- parallel::makeCluster(n_cores)
parallel::clusterExport(worker_cluster, c("df_agg", "strict_control"), envir = .GlobalEnv)
parallel::clusterEvalQ(worker_cluster, {library(lme4); options(na.action = "na.fail")})
message("Fitting candidate models with ", n_cores, " workers.")
dd <- tryCatch(
  MuMIn::dredge(model_interaction_full, rank = "AICc", cluster = worker_cluster),
  finally = parallel::stopCluster(worker_cluster)
)
best_model <- MuMIn::get.models(dd, subset = 1)[[1]]
saveRDS(model_interaction_full, "models/Step4_Full_Model_Results.rds")
saveRDS(dd, "models/Step4_Dredge_Results.rds")
saveRDS(best_model, "models/Step4_best_model_Results.rds")
write_csv(as.data.frame(dd), "results/Step4_Model_Selection.csv")

# 3. Report estimates and diagnostics for the full and selected models.
coefficient_table <- as.data.frame(coef(summary(best_model)))
coefficient_table$term <- rownames(coefficient_table)
write_csv(coefficient_table, "results/Step4_Model_Coefficients.csv")
sink("results/Step4_Model_Diagnostics.txt")
print(summary(model_interaction_full))
print(model_interaction_full@optinfo$conv)
print(performance::check_overdispersion(model_interaction_full))
print(performance::check_collinearity(model_interaction_full))
print(car::vif(model_interaction_full))
print(summary(best_model))
print(best_model@optinfo$conv)
cat("Singular full model:", lme4::isSingular(model_interaction_full), "\n")
cat("Singular selected model:", lme4::isSingular(best_model), "\n")
print(performance::check_overdispersion(best_model))
print(performance::check_collinearity(best_model))
sink()
write_csv(as.data.frame(performance::model_performance(best_model)),
          "results/Step4_Model_Performance_Metrics.csv")

# 4. Test each selected interaction against the selected model.
# Test main effects against an additive model, as in the original analysis.
# These are two sets of drop-one tests, not sequential backward elimination.
selected_terms <- attr(terms(best_model), "term.labels")
interaction_terms <- selected_terms[grepl(":", selected_terms)]
main_terms <- selected_terms[!grepl(":", selected_terms)]
drop_term <- function(model, term) {
  update(model, as.formula(paste(". ~ . -", term)))
}
lrt_row <- function(full, reduced, term, reference) {
  test <- anova(reduced, full, test = "Chisq")
  data.frame(term = term, reference_model = reference,
             AIC_full = AIC(full), AIC_reduced = AIC(reduced),
             Chisq = test$Chisq[2], Df = test$Df[2],
             P_value = test$`Pr(>Chisq)`[2])
}
interaction_tests <- lapply(interaction_terms, function(term) {
  lrt_row(best_model, drop_term(best_model, term), term, "selected")
})
additive_formula <- as.formula(paste(". ~ .", paste(paste("-", interaction_terms), collapse = " ")))
m_additive <- update(best_model, additive_formula)
main_tests <- lapply(main_terms, function(term) {
  lrt_row(m_additive, drop_term(m_additive, term), term, "additive")
})
write_csv(bind_rows(c(interaction_tests, main_tests)), "results/Step4_Model_LRT.csv")
capture.output(sessionInfo(), file = "results/GLMM_sessionInfo.txt")
message("GLMM workflow complete. Run code/Figure_2.R and code/Figure_3.R next.")

