library(mediation)
library(writexl)
library(cowplot)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(dplyr)
library(magick)

# 0. Directory setup: ensure the output folder exists
results_dir <- "../results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 1. Data loading and preprocessing (adapted for Table_S1.csv)
data <- read.csv("../Dataset/Table_S1.csv", stringsAsFactors = FALSE)

# Extract and clean covariates as factors
data <- data %>%
  mutate(
    Distance_factor = as.factor(gsub("Food", "", Site)),
    Food_Amount_factor = as.factor(gsub("Feed", "", Food_Amount))
  )

# 2. Build mediation models (controlling for Distance and Food Amount)
# Model 1: CAT -> CFT -> EAT
med_model1 <- lm(CFT ~ CAT + Distance_factor + Food_Amount_factor, data = data)
med_model1_1 <- lm(EAT ~ CAT + CFT + Distance_factor + Food_Amount_factor, data = data)

# Model 2: CFT -> EAT -> EFT
med_model2 <- lm(EAT ~ CFT + Distance_factor + Food_Amount_factor, data = data)
outcome_model <- lm(EFT ~ EAT + CFT + Distance_factor + Food_Amount_factor, data = data)

# 3. Perform mediation analysis and export Table 2
set.seed(123) # Set seed for reproducibility of simulations
med_analysis1 <- mediate(med_model1, med_model1_1, treat = "CAT", mediator = "CFT", robustSE = TRUE, sims = 1000)
med_analysis2 <- mediate(med_model2, outcome_model, treat = "CFT", mediator = "EAT", robustSE = TRUE, sims = 1000)

# Extract specific metrics for Table 2
# d0: ACME, z0: ADE, tau.coef: Total Effect, n0: Proportion Mediated, d0.p: p-value of ACME
table_2 <- data.frame(
  Model = c("Model 1: CAT-CFT-EAT", "Model 2: CFT-EAT-EFT"),
  ACME = round(c(med_analysis1$d0, med_analysis2$d0), 3),
  ADE = round(c(med_analysis1$z0, med_analysis2$z0), 3),
  Total_Effect = round(c(med_analysis1$tau.coef, med_analysis2$tau.coef), 3),
  Proportion_mediated = round(c(med_analysis1$n0, med_analysis2$n0), 3),
  p_value = round(c(med_analysis1$d0.p, med_analysis2$d0.p), 3)
)

# Export Table 2
write_xlsx(table_2, path = file.path(results_dir, "Table_2_Mediation_Results.xlsx"))
cat("[SUCCESS] Table 2 exported to results folder.\n")

# 4. Sensitivity analysis and export Figure 6
sens.out1 <- medsens(med_analysis1, rho.by = 0.1, effect.type = "indirect", sims = 1000)
sens.out2 <- medsens(med_analysis2, rho.by = 0.1, effect.type = "indirect", sims = 1000)

# Function to save high-resolution base R plots
save_sens_png <- function(sens_obj, file, width_in = 4.25, height_in = 3.5, dpi = 600) {
  px_w <- width_in * dpi
  px_h <- height_in * dpi
  ragg::agg_png(filename = file, width = px_w, height = px_h, units = "px", res = dpi)
  op <- par(family = "serif", mar = c(5.1, 5.1, 3.5, 2.1), cex = 0.95)
  on.exit({par(op); dev.off()}, add = TRUE)
  plot(sens_obj, main = "")
}

save_sens_png(sens.out1, file.path(results_dir, "sens_A.png"))
save_sens_png(sens.out2, file.path(results_dir, "sens_B.png"))

# Combine the two sensitivity plots into Figure 6
gA <- cowplot::ggdraw() + cowplot::draw_image(magick::image_read(file.path(results_dir, "sens_A.png")))
gB <- cowplot::ggdraw() + cowplot::draw_image(magick::image_read(file.path(results_dir, "sens_B.png")))

fig_6 <- cowplot::plot_grid(
  gA, gB, ncol = 2, labels = c("a", "b"), label_size = 18, 
  label_fontface = "bold", label_fontfamily = "serif", 
  label_x = 0.02, label_y = 0.98, hjust = 0, vjust = 1
)

ggsave(file.path(results_dir, "Figure_6_Sensitivity_Analysis.png"), plot = fig_6, width = 8.8, height = 3.8, dpi = 600)
cat("[SUCCESS] Figure 6 (Sensitivity Analysis) exported to results folder.\n")

# 5. Path diagram plotting (DiagrammeR) - Figure 5
# Note: Path coefficients here are hardcoded for visual layout purposes.
graph_code <- "
digraph mediation {
  graph [layout = dot, rankdir = LR]
  X [label = 'Crows First Arrival Time\\n(CAT)', shape = ellipse, fontname = 'serif']
  M1 [label = 'Crow First Feeding Time\\n(CFT)', shape = ellipse, fontname = 'serif']
  M2 [label = 'WTE First Arrival Time\\n(EAT)', shape = ellipse, fontname = 'serif']
  Y [label = 'WTE First Feeding Time\\n(EFT)', shape = ellipse, fontname = 'serif']
  
  X -> M1 [label = 'a1 = 0.9105', fontsize = 10, fontname = 'serif']
  M1 -> M2 [label = 'b1 = 0.7309', fontsize = 10, fontname = 'serif']
  X -> M2 [label = 'a2 = 0.2017', fontsize = 10, style = 'dashed', fontname = 'serif']
  M2 -> Y [label = 'b2 = 0.9269', fontsize = 10, fontname = 'serif']
  M1 -> Y [label = 'c = 0.0363', fontsize = 10, style = 'dashed', fontname = 'serif']
  X -> Y [label = 'direct effect\\n(c\\' = 0.0363)', fontsize = 10, fontname = 'serif']
  
  {rank=same; X; Y}
  {rank=same; M1; M2}
}
"
g <- DiagrammeR::grViz(graph_code)
svg_txt <- export_svg(g)
rsvg_png(charToRaw(svg_txt), file = file.path(results_dir, "Figure_5_Mediation_Pathways.png"), width = 10*600, height = 5*600)
cat("[SUCCESS] Figure 5 (Mediation Pathways) exported to results folder.\n")

cat("\n>>> Script 02_parallel_analysis.R finished successfully! <<<\n")
