library(ggplot2)
library(dplyr)
library(broom)
library(logistf)
library(survival)
library(survminer)

# 0. Directory setup: ensure the output folder exists
results_dir <- "../results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 1. Set plotting theme
sci_theme <- function(base_size = 15, base_family = "sans") { 
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.title = element_text(face = "bold", color = "black", size = rel(1.2)), 
      axis.text = element_text(color = "black", size = rel(1.1)),             
      legend.title = element_text(face = "bold", size = rel(1.15)),          
      legend.text = element_text(size = rel(1.0)),                           
      legend.position = "right",
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(10, 10, 10, 10, "pt"),
      plot.title = element_text(hjust = 0, face = "bold") 
    )
}
theme_set(sci_theme())

# 2. Data loading and preprocessing (adapted for Dataset_S4)
data <- read.csv("../Dataset/Dataset_S4_Eagle_Foraging_Decisions.csv", 
                 stringsAsFactors = FALSE)

data_processed <- data %>%
  mutate(
    # Dependent variable: Feeding Success (1 = Yes, 0 = No)
    Trial_Outcome = ifelse(Climax_No_eagle != 0, 1, 0),
    
    # Independent variable: Crow Signal (1 = Present, 0 = Absent)
    Crow_Signal = ifelse(No_Crow_eagle_standby > 0, 1, 0),
    Crow_Signal_f = factor(Crow_Signal, levels = c(0, 1), labels = c("No Crows", "Crows Present")),
    
    Date = as.factor(Date), 
    Feed_Numeric = as.numeric(gsub("Feed", "", Feed_Amount)),
    
    # Convert Location distance to human disturbance intensity (higher value = higher disturbance)
    Disturbance_Numeric = case_when(
      Location == "Food101" ~ 4, 
      Location == "Food102" ~ 3,
      Location == "Food104" ~ 2,
      Location == "Food106" ~ 1, 
      TRUE ~ NA_real_ 
    ),
    
    # Scale continuous variables for regression models
    scaled_Disturbance = as.numeric(scale(Disturbance_Numeric)),
    scaled_Feed = as.numeric(scale(Feed_Numeric))
  ) %>%
  filter(!is.na(Disturbance_Numeric) & !is.na(Feed_Numeric) & !is.na(No_Crow_eagle_standby))

# 3. Willingness / Success Rate - Firth Penalized Logistic Regression (Table 3)
model_firth <- logistf(
  Trial_Outcome ~ Crow_Signal_f + scaled_Disturbance + scaled_Feed,
  data = data_processed
)

# Extract and export Table 3 (Calculating Chi-Square from p-value using qchisq)
results_firth_tidy <- data.frame(
  Predictor = names(model_firth$coefficients),
  Coefficient = round(model_firth$coefficients, 3),
  Std_Error = round(sqrt(diag(vcov(model_firth))), 3),
  Lower_95_CI = round(model_firth$ci.lower, 3),
  Upper_95_CI = round(model_firth$ci.upper, 3),
  Chi_Square = round(qchisq(model_firth$prob, df = 1, lower.tail = FALSE), 3),
  p_value = round(model_firth$prob, 3)
)
write.csv(results_firth_tidy, file.path(results_dir, "Table_3_Firth_Logistic_Regression.csv"), row.names = FALSE)
cat("[SUCCESS] Table 3 exported to results folder.\n")

# 4. Decision Time / Hesitation - Cox Proportional Hazards Model (Figure 7)
data_survival <- data_processed %>%
  mutate(
    status = ifelse(Trial_Outcome == 1, 1, 0), # 1 = Event occurred, 0 = Right-censored
    # Maximum wait time set to 120 mins if event did not occur
    surv_time = ifelse(Trial_Outcome == 1, Time_difference, 120),
    surv_time = ifelse(surv_time <= 0, 0.1, surv_time) # Prevent 0 time
  )

model_cox <- coxph(Surv(surv_time, status) ~ Crow_Signal_f + scaled_Disturbance + scaled_Feed, data = data_survival)

# Export Cox results (HR value)
results_cox <- tidy(model_cox, exponentiate = TRUE, conf.int = TRUE)
write.csv(results_cox, file.path(results_dir, "Cox_Decision_Time_Results.csv"), row.names = FALSE)
cat("[SUCCESS] Cox survival analysis results (HR) exported to results folder.\n")

# Plot Figure 7 (Survival Curve)
plot_curve <- survfit(Surv(surv_time, status) ~ Crow_Signal_f, data = data_survival)

plot_cox_viz <- ggsurvplot(
  plot_curve,
  data = data_survival,
  size = 1.2,
  palette = c("#999999", "#E69F00"), # Grey (No crows) vs Orange (Crows present)
  conf.int = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  legend.labs = c("No Crows", "Crows Present"),
  xlab = "Time (min)",
  ylab = "Probability of Remaining Hesitant",
  ggtheme = sci_theme()
)

# Save as PNG
png(file.path(results_dir, "Figure_7_Decision_Time_Survival.png"), width = 2400, height = 2400, res = 300)
print(plot_cox_viz)
dev.off()

# Save as PDF
pdf(file.path(results_dir, "Figure_7_Decision_Time_Survival.pdf"), width = 8, height = 8)
print(plot_cox_viz)
dev.off()

cat("[SUCCESS] Figure 7 exported in both PNG and PDF formats.\n")

# 5. Absolute Efficiency - Linear Regression on EFT (Table S3)
data_efficiency <- data_processed %>%
  filter(Trial_Outcome == 1) %>%
  filter(!is.na(Eagle_First_Feeding_Time))

model_efficiency_feeding <- lm(
  Eagle_First_Feeding_Time ~ Crow_Signal_f + scaled_Disturbance + scaled_Feed,
  data = data_efficiency
)

write.csv(tidy(model_efficiency_feeding), file.path(results_dir, "Table_S3_Efficiency_EFT_LM.csv"), row.names = FALSE)
cat("[SUCCESS] Table S3 (EFT Linear Regression) exported to results folder.\n")

cat("\n>>> Script 03_willingness_decisions_and_efficiency.R finished successfully! <<<\n")
