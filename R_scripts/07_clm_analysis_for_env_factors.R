library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(Hmisc)    
library(patchwork)
library(lubridate)
library(broom)
library(ordinal)  
library(brms)
library(tidybayes)
library(posterior)

# 0. Directory setup: input and output paths
data_file <- "../Dataset/Dataset_S5_Temporal_Sequences.csv"
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

# 2. Data loading and preprocessing
data <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(
    Location = as.factor(Location),
    Distance = case_when(
      as.character(Location) == "101" ~ 40,
      as.character(Location) == "102" ~ 80,
      as.character(Location) == "104" ~ 120,
      as.character(Location) == "106" ~ 160,
      TRUE ~ NA_real_
    ),
    Food_Amount_Num = as.numeric(Food_Amount),
    
    # Numeric variables for Y-axis plotting
    EOS_num = as.numeric(Eagle_Peak_Sequence_EOS),
    AOS_E_num = as.numeric(Eagle_Arrival_Sequence_AOS_E),
    AOS_C_num = as.numeric(Crow_Arrival_Sequence_AOS_C),
    
    # Ordered factors for ordinal regression models (CLM / Bayesian)
    EOS_fac = factor(EOS_num, levels = c("1", "2", "3"), ordered = TRUE),
    AOS_E_fac = factor(AOS_E_num, levels = c("1", "2", "3"), ordered = TRUE),
    AOS_C_fac = factor(AOS_C_num, levels = c("1", "2", "3"), ordered = TRUE),
    
    # Standardization for modeling
    Dist_Std = as.numeric(scale(Distance)),
    Food_Std = as.numeric(scale(Food_Amount_Num)),
    Dist40 = Distance / 40
  )

# 3. Non-parametric Spearman Correlation
cat(">>> Calculating Spearman Correlation...\n")
correlation_vars <- c("Distance", "EOS_num", "AOS_E_num", "AOS_C_num", "Food_Amount_Num")
correlation_results <- rcorr(as.matrix(data[, correlation_vars]), type = "spearman")

write.csv(correlation_results$r, file.path(results_dir, "Spearman_Correlation_r.csv"))
write.csv(correlation_results$P, file.path(results_dir, "Spearman_Correlation_P_values.csv"))

# 4. Ordinal Regression Analysis (CLM & Bayesian)
cat(">>> Running Bayesian Ordinal Model for EOS (This may take 1-2 minutes to compile)...\n")
# Bayesian model for EOS
fit_eos40 <- brm(
  EOS_fac ~ Dist40 + scale(Food_Amount_Num),
  data = data,
  family = cumulative("logit"),
  prior = c(prior(normal(0, 1), class = "b")),
  chains = 4, cores = 4, iter = 4000,
  refresh = 0, silent = 2, seed = 123
)

res_eos <- list(
  model = fit_eos40,
  type = "Bayes Ordinal",
  label = "Bayes Ordinal\nβ(Dist) = -2.48 [-3.62, -1.44]\nβ(Food) = -0.04 [-0.95, 0.82]")

# Wrapper function for Cumulative Link Models (CLM)
analyze_sequence <- function(y_fac_col, y_num_col, data_in) {
  model_fit <- tryCatch({
    clm(data_in[[y_fac_col]] ~ Dist_Std + Food_Std, data = data_in)
  }, error = function(e) return(NULL), warning = function(w) return(NULL))
  
  if (!is.null(model_fit)) {
    model_type <- "CLM (Ordinal)"
    final_model <- model_fit
    coefs <- summary(model_fit)$coefficients
    p_dist <- coefs[grep("Dist_Std", rownames(coefs)), 4]
    p_food <- coefs[grep("Food_Std", rownames(coefs)), 4]
  } else {
    model_type <- "LM (Linear)"
    final_model <- lm(data_in[[y_num_col]] ~ Dist_Std + Food_Std, data = data_in)
    coefs <- summary(final_model)$coefficients
    p_dist <- coefs["Dist_Std", "Pr(>|t|)"]
    p_food <- coefs["Food_Std", "Pr(>|t|)"]
  }
  
  fmt_p <- function(p) { if(is.na(p)) return("NA"); if(p < 0.001) return("< 0.001"); return(paste0("= ", round(p, 3))) }
  label_text <- paste0(model_type, "\n", "P(Dist) ", fmt_p(p_dist), "\n", "P(Food) ", fmt_p(p_food))
  return(list(model = final_model, label = label_text, type = model_type))
}

cat(">>> Running CLM Models for AOS_E and AOS_C...\n")
res_aos_e <- analyze_sequence("AOS_E_fac", "AOS_E_num", data)
res_aos_c <- analyze_sequence("AOS_C_fac", "AOS_C_num", data)

# 5. Visualization 1: Six-Panel Bubble Plot (Figure 11)
plot_configs <- list(
  list(x_var="Distance", y_var="EOS_num", res=res_eos, xlab="Distance (m)", ylab="Eagle Peak Seq (EOS)"),
  list(x_var="Distance", y_var="AOS_E_num", res=res_aos_e, xlab="Distance (m)", ylab="Eagle Arrival Seq (AOS_E)"),
  list(x_var="Distance", y_var="AOS_C_num", res=res_aos_c, xlab="Distance (m)", ylab="Crow Arrival Seq (AOS_C)"),
  list(x_var="Food_Amount_Num", y_var="EOS_num", res=res_eos, xlab="Food Amount (items)", ylab="Eagle Peak Seq (EOS)"),
  list(x_var="Food_Amount_Num", y_var="AOS_E_num", res=res_aos_e, xlab="Food Amount (items)", ylab="Eagle Arrival Seq (AOS_E)"),
  list(x_var="Food_Amount_Num", y_var="AOS_C_num", res=res_aos_c, xlab="Food Amount (items)", ylab="Crow Arrival Seq (AOS_C)")
)

individual_plots <- list()
food_discrete_colors <- c("6" = "#F0CA8C", "18" = "#DD8800", "36" = "#783124")
location_colors <- c("101" = "#355483", "102" = "#5494BE", "104" = "#97D2F0", "106" = "#CFE3EF")

for (i in seq_along(plot_configs)) {
  cfg <- plot_configs[[i]]
  
  if (cfg$xlab == "Distance (m)") {
    p <- ggplot(data, aes(x = .data[[cfg$x_var]], y = .data[[cfg$y_var]], color = as.factor(Food_Amount_Num), fill = as.factor(Food_Amount_Num))) + 
      geom_count(shape = 1, stroke = 1.5, position = "identity", alpha = 0.8) +
      scale_color_manual(values = food_discrete_colors) + scale_fill_manual(values = food_discrete_colors) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
      labs(color = "Food (items)", fill = "Food (items)", size = "Sample Count") 
  } else {
    p <- ggplot(data, aes(x = .data[[cfg$x_var]], y = .data[[cfg$y_var]], color = Location, fill = Location)) +
      geom_count(shape = 1, stroke = 1.5, position = "identity", alpha = 0.8) +
      scale_color_manual(values = location_colors) + scale_fill_manual(values = location_colors) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
      labs(color = "Location", fill = "Location", size = "Sample Count")
  }
  
  p <- p + scale_size_area(max_size = 8, breaks = c(1, 2, 4, 6, 8)) + 
    labs(x = cfg$xlab, y = cfg$ylab) +
    annotate("text", x = -Inf, y = Inf, label = cfg$res$label, hjust = -0.1, vjust = 1.2, size = 3.5, fontface = "italic")
  
  if (i == 3 || i == 6) { p <- p + theme(legend.position = "right") } else { p <- p + theme(legend.position = "none") }
  individual_plots[[i]] <- p
}

final_plot <- (individual_plots[[1]] | individual_plots[[2]] | individual_plots[[3]]) /
  (individual_plots[[4]] | individual_plots[[5]] | individual_plots[[6]]) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(results_dir, "Figure_9_Sequence_Analysis_Bubble.png"), final_plot, width = 18, height = 10, dpi = 300)
ggsave(file.path(results_dir, "Figure_9_Sequence_Analysis_Bubble.pdf"), final_plot, width = 18, height = 10)
cat("[SUCCESS] Figure 9 exported.\n")

cat("\n>>> Script 07_clm_analysis_for_env_factors.R finished successfully! <<<\n")
