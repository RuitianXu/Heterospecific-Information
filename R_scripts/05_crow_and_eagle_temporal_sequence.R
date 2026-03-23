library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)

# 0. Directory setup: input and output paths
counts_folder <- "../Dataset/Dataset_S2_Time_Series_Individual_Counts"
metadata_file <- "../Dataset/Dataset_S1_Experimental_Design_and_Metadata.csv"

results_dir <- "../results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 1. Set plotting theme
sci_theme <- function(base_size = 16, base_family = "sans") { 
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
      plot.title = element_blank(),
      strip.background = element_rect(fill = "#F0F0F0", color = NA),
      strip.text = element_text(face = "bold", size = rel(1.1))
    )
}
theme_set(sci_theme())

color_palette <- c(
  "Crow" = "#22577A",            # Dark blue
  "Watching Eagle" = "#E86A33",  # Bright orange
  "Feeding Eagle" = "#FFB200"    # Warm yellow
)

# 2. Data loading and preprocessing
# Read all daily CSVs and merge into one dataframe
count_files <- list.files(path = counts_folder, pattern = "\\.csv$", full.names = TRUE)

all_counts_raw <- map_dfr(count_files, function(file_path) {
  # Add skip = 1 to skip the messy first row ("Unnamed:...") and use the actual headers
  df <- read_csv(file_path, skip = 1, show_col_types = FALSE)
  date_extracted <- str_extract(basename(file_path), "\\d{8}")
  df <- df %>% mutate(Date = date_extracted)
  return(df)
})

# Read metadata
metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    Date = as.character(Date),
    Sample_Point = as.character(Sample_Point)
  ) 

# Convert raw wide data to cleaned long format
data_long <- all_counts_raw %>%
  group_by(Date) %>% 
  mutate(Time_minutes = (row_number() - 1) * 0.5) %>% # 1 frame = 30 seconds = 0.5 mins
  ungroup() %>%
  pivot_longer(
    cols = starts_with(c("Crow_", "Watch_Eagle_", "Feed_Eagle_")),
    names_to = "Variable",
    values_to = "Count"
  ) %>%
  separate(Variable, into = c("Animal_Type", "Sample_Point"), sep = "_(?=[0-9]+$)") %>%
  mutate(
    Count = replace_na(as.numeric(Count), 0),
    Animal_Type = case_when(
      Animal_Type == "Watch_Eagle" ~ "Watching Eagle",
      Animal_Type == "Feed_Eagle"  ~ "Feeding Eagle",
      TRUE ~ "Crow"
    )
  )

# Merge with metadata
all_data_final <- left_join(data_long, metadata, by = c("Date", "Sample_Point"))

# Factorize variables for correct facet ordering
food_levels <- sort(unique(metadata$Food_Amount))
distance_levels <- c("40m", "80m", "120m", "160m")

all_data_final <- all_data_final %>%
  filter(!is.na(Food_Amount)) %>% # Remove unmapped records
  mutate(
    Food_Amount_Label = factor(paste(Food_Amount, "items"), levels = paste(food_levels, "items")),
    Distance_Label = case_when(
      Sample_Point == "101" ~ "40m",
      Sample_Point == "102" ~ "80m",
      Sample_Point == "104" ~ "120m",
      Sample_Point == "106" ~ "160m",
      TRUE ~ Sample_Point
    ),
    Distance_Label = factor(Distance_Label, levels = distance_levels)
  )

# 3. Split data into First 4 days (Figure 8) and Last 4 days (Figure S1)
unique_dates <- sort(unique(all_data_final$Date))
dates_first_4 <- unique_dates[1:4]
dates_last_4 <- unique_dates[5:8]

data_fig8 <- all_data_final %>% filter(Date %in% dates_first_4)
data_figS1 <- all_data_final %>% filter(Date %in% dates_last_4)

# 4. Create a reusable plotting function for smoothed lines
plot_smoothed_sequence <- function(df_subset) {
  ggplot(df_subset, aes(x = Time_minutes, y = Count, color = Animal_Type, fill = Animal_Type)) +
    stat_smooth(geom = "area", method = "loess", span = 0.1, alpha = 0.2) +
    stat_smooth(geom = "line", method = "loess", span = 0.1, linewidth = 1) +
    geom_point(size = 1, alpha = 0.4) + 
    facet_grid(rows = vars(Distance_Label), cols = vars(Food_Amount_Label), scales = "free_y") +
    coord_cartesian(ylim = c(0, NA)) +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    labs(
      x = "Time (Minutes from 7:00 AM)",
      y = "Smoothed Count Trend",
      color = "Species Activity",
      fill = "Species Activity"
    )
}

# 5. Generate plots
p_fig8 <- plot_smoothed_sequence(data_fig8)
p_figS1 <- plot_smoothed_sequence(data_figS1)

# 6. Export plots in both PNG and PDF formats
# Figure 8 (First 4 days / Replicate 1)
ggsave(file.path(results_dir, "Figure_10_Temporal_Sequence_Smooth_Rep1.png"), plot = p_fig8, width = 16, height = 12, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Figure_10_Temporal_Sequence_Smooth_Rep1.pdf"), plot = p_fig8, width = 16, height = 12, bg = "white")

# Figure S1 (Last 4 days / Replicate 2)
ggsave(file.path(results_dir, "Figure_S1_Temporal_Sequence_Smooth_Rep2.png"), plot = p_figS1, width = 16, height = 12, dpi = 300, bg = "white")
ggsave(file.path(results_dir, "Figure_S1_Temporal_Sequence_Smooth_Rep2.pdf"), plot = p_figS1, width = 16, height = 12, bg = "white")

cat("[SUCCESS] Script 05_crow_and_eagle_temporal_sequence.R finished successfully!\n")
cat("          Figure 10 (First 4 days) and Figure S1 (Last 4 days) exported to results folder.\n")
