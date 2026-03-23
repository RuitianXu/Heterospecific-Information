library(lme4)
library(ggplot2)
library(dplyr)
library(MASS)
library(Hmisc)
library(reshape2)

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
      plot.title = element_blank()
    )
}
theme_set(sci_theme())

# 2. Data loading and preprocessing (adapted for Table_S1.csv)
# Load dataset using relative path
data <- read.csv("../Dataset/Table_S1.csv", stringsAsFactors = FALSE)

# Extract numeric values from strings to use as continuous variables in linear models
# Map Food101, 102, 104, 106 to distances 40, 80, 120, 160
data <- data %>%
  mutate(
    Location_Num = gsub("Food", "", Site),
    Distance = case_when(
      Location_Num == "101" ~ 40,
      Location_Num == "102" ~ 80,
      Location_Num == "104" ~ 120,
      Location_Num == "106" ~ 160,
      TRUE ~ NA_real_
    ),
    # Remove 'Feed' prefix and convert to numeric
    Food_Amount_Num = as.numeric(gsub("Feed", "", Food_Amount))
  )

# 3. Log transformation of behavioral times
data$log_CAT <- log(data$CAT)
data$log_CFT <- log(data$CFT)
data$log_EAT <- log(data$EAT)
data$log_EFT <- log(data$EFT)

# 4. Extract residuals via linear models (controlling for environmental factors)
# Formula: log_time ~ Distance + Food_Amount_Num
fit_CAT <- lm(log_CAT ~ Distance + Food_Amount_Num, data = data)
fit_CFT <- lm(log_CFT ~ Distance + Food_Amount_Num, data = data)
fit_EAT <- lm(log_EAT ~ Distance + Food_Amount_Num, data = data)
fit_EFT <- lm(log_EFT ~ Distance + Food_Amount_Num, data = data)

resid_CAT <- residuals(fit_CAT)
resid_CFT <- residuals(fit_CFT)
resid_EAT <- residuals(fit_EAT)
resid_EFT <- residuals(fit_EFT)

# 5. Normality test of residuals & Correlation method selection
sw_test <- list(
  sw_CAT = shapiro.test(resid_CAT),
  sw_CFT = shapiro.test(resid_CFT),
  sw_EAT = shapiro.test(resid_EAT),
  sw_EFT = shapiro.test(resid_EFT)
)

# Print P-values
p_values <- sapply(sw_test, function(x) x$p.value)
is_normal <- all(p_values > 0.05)

cat("\n>>> Shapiro-Wilk Test P-values (P > 0.05 indicates normality) <<<\n")
print(p_values)

# Automatically determine the correlation method
corr_method <- if(is_normal) "pearson" else "spearman"

cat(paste0("\n>>> Statistical Recommendation: <<<\nSince ", 
           if(is_normal) "all residuals are normally distributed" else "non-normal residuals exist", 
           ", ", toupper(corr_method), " correlation analysis is recommended.\n"))

# 6. Calculate correlation
residual_data <- data.frame(
  CAT = resid_CAT,
  CFT = resid_CFT,
  EAT = resid_EAT,
  EFT = resid_EFT
)

res_corr <- rcorr(as.matrix(residual_data), type = corr_method)

cat(paste0("\n>>> Correlation Matrix (Method: ", corr_method, ") <<<\n"))
print(res_corr$r)
cat("\n>>> P-value Matrix <<<\n")
print(res_corr$P)

# 7. Visualization and result export
corr_matrix <- res_corr$r
p_matrix <- res_corr$P
corr_melt <- melt(corr_matrix)
p_melt <- melt(p_matrix)

plot_data <- corr_melt
plot_data$p_value <- p_melt$value

# Add significance asterisks
plot_data$significance <- ifelse(plot_data$p_value < 0.001, "***",
                                 ifelse(plot_data$p_value < 0.01, "**",
                                        ifelse(plot_data$p_value <= 0.05, "*", " ")))

# Keep only the lower triangle
plot_data <- plot_data %>% filter(as.numeric(Var1) > as.numeric(Var2))

# Format variable labels
eng_labels <- c(
  "CAT" = "Crow Arr (Resid)",
  "CFT" = "Crow Feed (Resid)",
  "EAT" = "WTSE Arr (Resid)",
  "EFT" = "WTSE Feed (Resid)"
)
plot_data <- plot_data %>%
  mutate(Var1 = recode(Var1, !!!eng_labels), 
         Var2 = recode(Var2, !!!eng_labels))

my_colors <- c("#c7e0ed", "#bfafd2", "#327db7", "#134687", "#053061")

p1 <- ggplot(plot_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = my_colors, limit = c(0, 1), 
                       name = paste0(toupper(substr(corr_method,1,1)), substr(corr_method,2,nchar(corr_method)), "\nCorrelation")) +
  geom_text(aes(label = paste0(round(value, 2), significance)), color = "white", size = 5, family = "sans") +
  scale_x_discrete(position = "bottom") +
  labs(x = NULL, y = NULL, title = "Correlation of Residuals (Controlled for Dist & Food)") +
  theme(axis.text.x = element_text(size = 12, hjust = 0.5), 
        axis.text.y = element_text(size = 12), 
        panel.grid = element_blank())

print(p1)

# Output to results folder
output_path <- file.path(results_dir, "Figure_4_temporal_correlation.png")
ggsave(filename = output_path, 
       plot = p1, 
       width = 8, 
       height = 6, 
       bg = "white",
       dpi = 300)

cat(paste0("\n[SUCCESS] Correlation heatmap successfully generated and saved to: ", output_path, "\n"))
