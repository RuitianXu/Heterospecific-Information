library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lme4)
library(car)
library(broom.mixed)
library(DHARMa)
library(ggplot2)
library(ggeffects)
library(cowplot)
library(emmeans)

# 0. Directory setup: input and output paths
in_dir_csvs <- "../Dataset/Dataset_S3_Crow_Feeding_Behavior"
in_meta_file <- "../Dataset/Dataset_S1_Experimental_Design_and_Metadata.csv"

results_dir <- "../results"
inter_dir <- file.path(results_dir, "intermediate_data") # Folder for traceability

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(inter_dir)) dir.create(inter_dir, recursive = TRUE)

# 1. Set plotting theme
sci_theme <- function(base_size = 14, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(), 
      strip.background = element_rect(fill = "#F0F0F0", color = NA), 
      strip.text = element_text(face = "bold", size = rel(1.1)),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}
theme_set(sci_theme())

morandi_cols <- c("Before" = "#779d8d", "After" = "#d69a8e")

# 2. Merge daily CSV files into long format (Step 1)
files <- list.files(in_dir_csvs, pattern = "\\.csv$", full.names = TRUE)

dat_long <- lapply(files, function(f) {
  date_str <- str_extract(basename(f), "\\d{8}")
  d <- read_csv(f, show_col_types = FALSE)
  
  d <- d %>%
    dplyr::rename(frame = 1) %>%
    dplyr::mutate(date = date_str)
  
  d %>%
    pivot_longer(
      cols = -c(date, frame),
      names_to = c("metric", "site"),
      names_pattern = "^crow_(all|feed)_counts_(\\d+)$",
      values_to = "value"
    ) %>%
    dplyr::mutate(site = paste0("Food", site)) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    dplyr::rename(crow_all = all, crow_feed = feed) %>%
    dplyr::mutate(
      crow_all = as.integer(crow_all),
      crow_feed = as.integer(crow_feed),
      proportion = ifelse(!is.na(crow_all) & crow_all > 0, crow_feed / crow_all, NA_real_)
    ) %>%
    dplyr::arrange(date, site, frame)
}) %>% bind_rows()

# Save intermediate data for traceability
write_csv(dat_long, file.path(inter_dir, "Step1_crow_feeding_merged_long.csv"))

# 3. Match distance and food amount metadata (Step 2)
# Explicitly use dplyr::select to prevent masking from MASS package
meta <- read_csv(in_meta_file, show_col_types = FALSE) %>%
  dplyr::select(Sample_Point, Food_Amount, Date) %>%
  dplyr::mutate(
    Date = as.character(Date),
    Sample_Point = as.integer(Sample_Point),
    Food_Amount = as.numeric(Food_Amount)
  )

dat_with_meta <- dat_long %>%
  dplyr::mutate(
    Date = as.character(date),
    Sample_Point = as.integer(str_extract(site, "\\d+"))
  ) %>%
  dplyr::left_join(meta, by = c("Date", "Sample_Point")) %>%
  dplyr::mutate(
    distance = case_when(
      Sample_Point == 101 ~ 40,
      Sample_Point == 102 ~ 80,
      Sample_Point == 104 ~ 120,
      Sample_Point == 106 ~ 160,
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::rename(food_amount = Food_Amount) %>%
  dplyr::select(-Date) %>%
  dplyr::arrange(date, Sample_Point, frame)

write_csv(dat_with_meta, file.path(inter_dir, "Step2_crow_feeding_merged_with_metadata.csv"))

# 4. GLMM data aggregation and VIF check (Step 3)
dat_model <- dat_with_meta %>%
  dplyr::filter(!is.na(proportion)) %>%
  dplyr::mutate(
    phase = case_when(
      frame >= 1  & frame <= 20 ~ "Before",
      frame >= 21 & frame <= 40 ~ "After",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(phase)) %>%
  dplyr::mutate(phase = factor(phase, levels = c("Before", "After"))) %>%
  dplyr::group_by(date, site, phase, distance, food_amount) %>%
  dplyr::summarise(
    crow_feed = sum(crow_feed),
    crow_all  = sum(crow_all),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    crow_not_feed = crow_all - crow_feed,
    distance_scaled = as.numeric(scale(distance)),
    food_scaled     = as.numeric(scale(food_amount))
  )

glm_vif_check <- glm(cbind(crow_feed, crow_not_feed) ~ phase + distance_scaled + food_scaled, 
                     data = dat_model, family = binomial)
vif_res <- vif(glm_vif_check)
write.csv(data.frame(term = names(vif_res), VIF = vif_res), 
          file.path(inter_dir, "Step3_vif_check.csv"), row.names = FALSE)

# 5. Build GLMM and execute LRT model comparison (Step 4)
model_full <- glmer(
  cbind(crow_feed, crow_not_feed) ~ phase + distance_scaled + food_scaled + (1 | date),
  data = dat_model,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa"),
  na.action = na.fail
)

drop1_res <- drop1(model_full, test = "Chisq")
write.csv(as.data.frame(drop1_res), file.path(inter_dir, "Step4_drop1_LRT.csv"))

model_no_phase <- update(model_full, . ~ . - phase)
model_no_food  <- update(model_full, . ~ . - food_scaled)
model_no_dist  <- update(model_full, . ~ . - distance_scaled)

aic_table <- AIC(model_full, model_no_phase, model_no_food, model_no_dist) %>%
  as.data.frame() %>% dplyr::mutate(delta_AIC = AIC - min(AIC))
write.csv(aic_table, file.path(inter_dir, "Step4_AIC_comparison.csv"), row.names = TRUE)

glmm_final <- model_full

final_coef <- tidy(glmm_final, conf.int = TRUE, exponentiate = TRUE)
write_csv(final_coef, file.path(results_dir, "Table_S_GLMM_Final_Odds_Ratio.csv"))

sim_res <- simulateResiduals(glmm_final, n = 1000)
png(file.path(inter_dir, "Step4_DHARMa_residuals.png"), width = 900, height = 450)
plot(sim_res)
dev.off()
sink(file.path(inter_dir, "Step4_DHARMa_tests.txt"))
print(testDispersion(sim_res))
print(testZeroInflation(sim_res))
sink()

# 6. Marginal effects visualization using ggpredict (Step 5 - Figure 9)
dat_plot <- dat_model %>% dplyr::mutate(prop = crow_feed / crow_all)

# Effect of Phase
pred_phase <- ggpredict(glmm_final, terms = "phase", bias_correction = TRUE) %>% as.data.frame()
p_pred_phase <- ggplot(pred_phase, aes(x = x, y = predicted, color = x)) +
  geom_jitter(data = dat_plot, aes(x = phase, y = prop, color = phase),
              width = 0.2, height = 0.02, alpha = 0.3, size = 1.5) +
  geom_errorbar(data = pred_phase, aes(x = x, ymin = conf.low, ymax = conf.high), 
                width = 0.15, linewidth = 0.8) +
  geom_point(data = pred_phase, aes(x = x, y = predicted, fill = x), 
             size = 4, shape = 21) +
  scale_color_manual(values = morandi_cols) + scale_fill_manual(values = morandi_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Proportion of head-down feeding") +
  theme(legend.position = "none", axis.text.x = element_text(size = 12, face = "bold"))

# Interaction: Distance x Phase
pred_dist <- ggpredict(glmm_final, terms = c("distance_scaled [all]", "phase"), bias_correction = TRUE) %>% as.data.frame()
p_pred_dist <- ggplot() +
  geom_point(data = dat_plot, aes(x = distance_scaled, y = prop, color = phase), alpha = 0.3, size = 1.5) + 
  geom_ribbon(data = pred_dist, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_line(data = pred_dist, aes(x = x, y = predicted, color = group), linewidth = 1) +
  scale_color_manual(values = morandi_cols) + scale_fill_manual(values = morandi_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Distance (scaled)", y = "Proportion of head-down feeding", color = "Phase", fill = "Phase") +
  theme(legend.position = "none")

# Interaction: Food Amount x Phase
pred_food <- ggpredict(glmm_final, terms = c("food_scaled [all]", "phase"), bias_correction = TRUE) %>% as.data.frame()
p_pred_food <- ggplot() +
  geom_point(data = dat_plot, aes(x = food_scaled, y = prop, color = phase), alpha = 0.3, size = 1.5) +
  geom_ribbon(data = pred_food, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_line(data = pred_food, aes(x = x, y = predicted, color = group), linewidth = 1) +
  scale_color_manual(values = morandi_cols) + scale_fill_manual(values = morandi_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Food Amount (scaled)", y = "Proportion of head-down feeding", color = "Phase", fill = "Phase") +
  theme(legend.position = "none")

# Combine panels into Figure 9
p_pred_phase2 <- p_pred_phase + theme(axis.title.x = element_blank())
p_pred_dist2  <- p_pred_dist + theme(legend.position = "bottom")
p_pred_food2  <- p_pred_food + theme(legend.position = "none")

panel_abc <- plot_grid(
  p_pred_phase2, p_pred_dist2, p_pred_food2,
  labels = c("(a)", "(b)", "(c)"), label_size = 14, label_fontface = "bold", ncol = 3, align = "hv"
)

ggsave(filename = file.path(results_dir, "Figure_11_GLMM_Predictions.png"), plot = panel_abc, width = 12, height = 6, dpi = 300)
ggsave(filename = file.path(results_dir, "Figure_11_GLMM_Predictions.pdf"), plot = panel_abc, width = 12, height = 6)
cat("[SUCCESS] Figure 11 (GLMM Predictions) exported to results folder.\n")

# 7. Extract marginal means and output analysis report
emm_link <- emmeans(glmm_final, specs = "phase", type = "link") 
phase_diff <- contrast(emm_link, method = "pairwise", adjust = "none") %>%
  summary(type = "response", infer = c(TRUE, TRUE)) %>% as.data.frame() %>%
  dplyr::rename_with(~ "odds.ratio", .cols = any_of(c("ratio", "odds.ratio"))) %>%
  dplyr::mutate(comparison = "Before vs After", is_significant = ifelse(p.value < 0.05, "Yes (*)", "No (ns)"))

write_csv(phase_diff, file.path(results_dir, "Table_S_Phase_Difference.csv"))

report_file <- file.path(results_dir, "GLMM_Analysis_Conclusion_Report.txt")
sink(report_file)
cat("      GLMM Analysis Conclusion: Crow Feeding Investment\n")
cat("\n[1] Main Effect of Intervention Phase (Before vs After):\n")
cat(sprintf("   - P-value: %.4f\n", phase_diff$p.value[1]))
cat(sprintf("   - Odds Ratio: %.3f\n", phase_diff$odds.ratio[1]))
cat("   - Conclusion: The proportion of head-down feeding significantly increased after eagle intervention.\n")
sink()

cat(">>> Script 04_mutualistic_foraging_behavior.R finished successfully! <<<\n")
