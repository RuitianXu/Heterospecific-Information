library(tidyverse)
library(lubridate)
library(patchwork)
library(broom)

# 0. Directory setup: input and output paths
counts_dir <- "../Dataset/Dataset_S2_Time_Series_Individual_Counts"
metadata_path <- "../Dataset/Dataset_S1_Experimental_Design_and_Metadata.csv"

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

# 2. Read and process count data to extract maximum aggregation scales
count_files <- list.files(counts_dir, pattern = "\\.csv$", full.names = TRUE)
all_counts_data <- list()

for (file_path in count_files) {
  file_name <- basename(file_path)
  date_str <- str_extract(file_name, "\\d{8}") 
  current_date <- ymd(date_str)
  
  # Use skip = 1 to bypass the messy original headers
  df_raw <- read_csv(file_path, skip = 1, show_col_types = FALSE) 
  
  file_max_counts <- tibble(
    Date = as.Date(character()),
    Sample_Point = numeric(),
    Max_Eagle_Count = numeric(),
    Max_Crow_Count = numeric()
  )
  
  # Dynamically extract sample points (e.g., 101, 102, 104, 106)
  sample_points_in_file <- unique(str_extract(names(df_raw), "\\d+$"))
  sample_points_in_file <- as.numeric(na.omit(sample_points_in_file)) %>% sort()
  
  for (sp in sample_points_in_file) {
    eagle_col_name <- paste0("Feed_Eagle_", sp)
    crow_col_name <- paste0("Crow_", sp)
    
    if (eagle_col_name %in% names(df_raw) && crow_col_name %in% names(df_raw)) {
      max_eagle <- max(df_raw[[eagle_col_name]], na.rm = TRUE)
      max_crow <- max(df_raw[[crow_col_name]], na.rm = TRUE)
      file_max_counts <- rbind(file_max_counts, 
                               data.frame(Date = current_date, 
                                          Sample_Point = sp,
                                          Max_Eagle_Count = max_eagle, 
                                          Max_Crow_Count = max_crow))
    }
  }
  all_counts_data[[file_name]] <- file_max_counts
}

final_counts_df <- bind_rows(all_counts_data)

# 3. Read and process metadata
metadata_df <- read_csv(metadata_path, show_col_types = FALSE) %>%
  rename_with(~ str_trim(.) %>% str_replace_all("\\s+", "_")) %>%
  mutate(Date = ymd(Date)) %>%
  dplyr::select(Sample_Point, Food_Amount, Date) %>%
  distinct()

# 4. Merge dataframes and convert distance to numeric
combined_df <- left_join(final_counts_df, metadata_df, by = c("Date", "Sample_Point")) %>%
  mutate(
    Sample_Point_Str = as.character(Sample_Point),
    # Map sample points to actual distance from road
    Distance = case_when(
      Sample_Point_Str == "101" ~ 40,
      Sample_Point_Str == "102" ~ 80,
      Sample_Point_Str == "104" ~ 120,
      Sample_Point_Str == "106" ~ 160,
      TRUE ~ NA_real_
    ),
    Food_Amount = as.numeric(Food_Amount)
  ) %>%
  drop_na(Food_Amount, Distance)

# 5. Statistical Analysis: Generalized Linear Models (GLM - Poisson distribution)
model_eagle <- glm(Max_Eagle_Count ~ Food_Amount * Distance, 
                   data = combined_df, 
                   family = poisson(link = "log"))

model_crow <- glm(Max_Crow_Count ~ Food_Amount * Distance, 
                  data = combined_df, 
                  family = poisson(link = "log"))

# Export tidy results for Table S4 and Table S5
eagle_results <- tidy(model_eagle) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(eagle_results, file.path(results_dir, "Table_S4_Eagle_GLM.csv"), row.names = FALSE)

crow_results <- tidy(model_crow) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(crow_results, file.path(results_dir, "Table_S5_Crow_GLM.csv"), row.names = FALSE)
cat("[SUCCESS] Table S4 (WTE GLM) and Table S5 (Crow GLM) exported.\n")

# 6. Visualization
# Helper function to extract overall model significance (Likelihood Ratio Test)
get_model_stats <- function(model_obj) {
  null_model <- update(model_obj, . ~ 1)
  lr_test <- anova(null_model, model_obj, test = "Chisq")
  p_value <- lr_test$`Pr(>Chi)`[2]
  
  if(is.na(p_value)) { p_text <- "P = NA" } 
  else if (p_value < 0.001) { p_text <- "P < 0.001" } 
  else { p_text <- paste0("P = ", round(p_value, 3)) }
  
  return(p_text)
}

dist_colors <- c(
  "40"  = "#355483",
  "80"  = "#5494BE",
  "120" = "#97D2F0",
  "160" = "#CFE3EF"
)

# Plotting function (X = Food_Amount, Color = Distance)
create_food_plot <- function(data, y_col, y_lab, model_obj) {
  ggplot(data, aes(x = Food_Amount, y = .data[[y_col]], 
                   color = as.factor(Distance), 
                   fill = as.factor(Distance))) + 
    geom_point(shape = 21, size = 3.5, stroke = 1.2, alpha = 0.6) +
    # Crucial: Use GLM Poisson for the smooth lines to correctly model count data
    geom_smooth(method = "glm", 
                method.args = list(family = poisson(link = "log")), 
                se = TRUE, alpha = 0.2) +
    scale_color_manual(values = dist_colors, name = "Distance (m)") +
    scale_fill_manual(values = dist_colors, name = "Distance (m)") +
    scale_x_continuous(breaks = c(6, 18, 36)) +
    labs(x = "Food Amount (items)", y = y_lab) +
    annotate("text", x = Inf, y = Inf, label = get_model_stats(model_obj), 
             hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold", color = "black")
}

p_eagle <- create_food_plot(combined_df, "Max_Eagle_Count", "Maximum WTE Count", model_eagle)
p_crow  <- create_food_plot(combined_df, "Max_Crow_Count", "Maximum Crow Count", model_crow)

# 7. Combine plots using patchwork (Figure 10)
combined_counts_plot <- (p_eagle | p_crow) +
  plot_annotation(
    tag_levels = 'a',
    tag_prefix = '(',   
    tag_suffix = ')',   
    theme = theme(
      plot.tag = element_text(face = 'bold', size = 18),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

# 8. Export plots in both PNG and PDF formats
ggsave(filename = file.path(results_dir, "Figure_8_Max_Counts_GLM.png"),
       plot = combined_counts_plot, width = 16, height = 7, dpi = 300)

ggsave(filename = file.path(results_dir, "Figure_8_Max_Counts_GLM.pdf"),
       plot = combined_counts_plot, width = 16, height = 7)

cat("[SUCCESS] Figure 8 (GLM Max Counts) exported in both PNG and PDF formats.\n")
cat(">>> Script 06_max_count_analysis_for_env_factors.R finished successfully! <<<\n")
