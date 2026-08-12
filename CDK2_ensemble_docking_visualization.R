#!/usr/bin/env Rscript

# ============================================================
# CDK2 ENSEMBLE DOCKING — PUBLICATION-READY VISUALIZATION
# ============================================================
#
# Input:
#   ensemble_docking_ranked.csv
#
# Receptors:
#   1VYW, 7RWF, 5A14, 5IF1, 1W98
#
# The script generates:
#
#   Figure 1 — Top-20 consensus ligands
#   Figure 2 — Top-20 ligand × receptor affinity heatmap
#   Figure 3 — Receptor-wise affinity distributions
#   Figure 4 — Ensemble mean affinity vs variability
#   Figure 5 — Consensus rank vs mean affinity
#   Figure 6 — Top-10 ligand affinity profiles across states
#   Figure 7 — Distribution of ensemble mean affinities
#
# All figures are saved as high-resolution PNG and PDF files.
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})


# ============================================================
# 2. SETTINGS
# ============================================================

INPUT_FILE <- "ensemble_docking_ranked.csv"

OUTPUT_DIR <- "ensemble_figures"

TOP_N <- 20

TOP_PROFILE_N <- 10

RECEPTORS <- c(
  "1VYW",
  "7RWF",
  "5A14",
  "5IF1",
  "1W98"
)

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 3. READ DATA
# ============================================================

if (!file.exists(INPUT_FILE)) {
  stop(
    paste0(
      "Input file not found: ", INPUT_FILE,
      "\nPlace ensemble_docking_ranked.csv in the working directory ",
      "or change INPUT_FILE in this script."
    )
  )
}

df <- read_csv(
  INPUT_FILE,
  show_col_types = FALSE
)


# ============================================================
# 4. VALIDATE DATA
# ============================================================

required_columns <- c(
  "Rank",
  "Ligand",
  "Initial_1VYW_Affinity",
  "Mean_Affinity",
  "Median_Affinity",
  "Best_Affinity",
  "Worst_Affinity",
  "SD_Affinity",
  "Successful_Receptors",
  RECEPTORS
)

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}


# ============================================================
# 5. CLEAN DATA
# ============================================================

df <- df %>%
  mutate(
    Rank = as.numeric(Rank),
    Mean_Affinity = as.numeric(Mean_Affinity),
    Median_Affinity = as.numeric(Median_Affinity),
    Best_Affinity = as.numeric(Best_Affinity),
    Worst_Affinity = as.numeric(Worst_Affinity),
    SD_Affinity = as.numeric(SD_Affinity),
    Successful_Receptors = as.numeric(Successful_Receptors)
  )


# Keep only ligands with complete 5-state ensemble docking
complete_df <- df %>%
  filter(
    Successful_Receptors == length(RECEPTORS)
  )

message(
  "Total ligands: ",
  nrow(df)
)

message(
  "Complete 5/5 ensemble ligands: ",
  nrow(complete_df)
)


# ============================================================
# 6. TOP-N DATA
# ============================================================

top20 <- complete_df %>%
  arrange(Rank) %>%
  slice_head(n = TOP_N) %>%
  mutate(
    Ligand = factor(
      Ligand,
      levels = rev(Ligand)
    )
  )


top10 <- complete_df %>%
  arrange(Rank) %>%
  slice_head(n = TOP_PROFILE_N)


# ============================================================
# 7. PUBLICATION THEME
# ============================================================

theme_publication <- theme_classic(
  base_size = 12
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      color = "grey30"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


# ============================================================
# 8. HELPER FUNCTION
# ============================================================

save_figure <- function(
    plot_object,
    filename,
    width = 8,
    height = 6
) {

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(filename, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    filename = file.path(
      OUTPUT_DIR,
      paste0(filename, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
}


# ============================================================
# FIGURE 1
# TOP-20 CONSENSUS LIGANDS
# ============================================================

fig1 <- ggplot(
  top20,
  aes(
    x = Ligand,
    y = Mean_Affinity
  )
) +

  geom_col(
    width = 0.72
  ) +

  geom_errorbar(
    aes(
      ymin = Mean_Affinity - SD_Affinity,
      ymax = Mean_Affinity + SD_Affinity
    ),
    width = 0.20,
    linewidth = 0.35
  ) +

  coord_flip() +

  scale_y_continuous(
    name = "Mean docking affinity (kcal/mol)"
  ) +

  labs(
    title = "Top 20 CDK2 ensemble-docking ligands",
    subtitle = "Mean affinity across five CDK2 structural states; error bars represent SD",
    x = "Ligand"
  ) +

  theme_publication


save_figure(
  fig1,
  "Figure_1_Top20_Consensus_Ligands",
  width = 8,
  height = 7
)


# ============================================================
# FIGURE 2
# TOP-20 LIGAND × RECEPTOR HEATMAP
# ============================================================

heatmap_data <- top20 %>%
  select(
    Ligand,
    all_of(RECEPTORS)
  ) %>%
  pivot_longer(
    cols = all_of(RECEPTORS),
    names_to = "Receptor",
    values_to = "Affinity"
  ) %>%
  mutate(
    Ligand = factor(
      Ligand,
      levels = rev(top20$Ligand)
    ),
    Receptor = factor(
      Receptor,
      levels = RECEPTORS
    )
  )


fig2 <- ggplot(
  heatmap_data,
  aes(
    x = Receptor,
    y = Ligand,
    fill = Affinity
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Affinity
      )
    ),
    size = 3.0
  ) +

  scale_fill_gradient2(
    name = "Affinity\n(kcal/mol)",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = -10
  ) +

  labs(
    title = "State-dependent docking profile of the top 20 ligands",
    subtitle = "More negative values indicate more favorable predicted binding",
    x = "CDK2 structural state",
    y = "Ligand"
  ) +

  theme_publication +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )


save_figure(
  fig2,
  "Figure_2_Top20_State_Heatmap",
  width = 8,
  height = 8
)


# ============================================================
# FIGURE 3
# RECEPTOR-WISE AFFINITY DISTRIBUTIONS
# ============================================================

distribution_data <- complete_df %>%
  select(
    Ligand,
    all_of(RECEPTORS)
  ) %>%
  pivot_longer(
    cols = all_of(RECEPTORS),
    names_to = "Receptor",
    values_to = "Affinity"
  )


# ------------------------------------------------------------
# Define colors for CDK2 structural states
# ------------------------------------------------------------

receptor_colors <- c(
  "1VYW" = "#4C78A8",
  "1W98" = "#F58518",
  "5A14" = "#54A24B",
  "5IF1" = "#E45756",
  "7RWF" = "#B279A2"
)


# ============================================================
# FIGURE
# ============================================================

fig3 <- ggplot(
  distribution_data,
  aes(
    x = Receptor,
    y = Affinity
  )
) +
  
  # ----------------------------------------------------------
# Violin distribution
# ----------------------------------------------------------

geom_violin(
  fill = "grey92",
  color = "grey35",
  trim = FALSE,
  linewidth = 0.45
) +
  
  # ----------------------------------------------------------
# Colored boxplots
# ----------------------------------------------------------

geom_boxplot(
  aes(fill = Receptor),
  width = 0.20,
  outlier.shape = NA,
  color = "grey20",
  linewidth = 0.45
) +
  
  # ----------------------------------------------------------
# Screening threshold
# ----------------------------------------------------------

geom_hline(
  yintercept = -10,
  linetype = "dashed",
  linewidth = 0.5,
  color = "black"
) +
  
  # ----------------------------------------------------------
# Manual colors
# ----------------------------------------------------------

scale_fill_manual(
  values = receptor_colors
) +
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  title = "Distribution of docking affinities across CDK2 states",
  subtitle = "Dashed line indicates the initial screening threshold of −10 kcal/mol",
  x = "CDK2 structural state",
  y = "Docking affinity (kcal/mol)"
) +
  
  # ----------------------------------------------------------
# Publication theme
# ----------------------------------------------------------

theme_publication +
  
  theme(
    legend.position = "none"
  )


# ============================================================
# SAVE FIGURE
# ============================================================

save_figure(
  fig3,
  "Figure_3_Receptor_Affinity_Distributions",
  width = 8,
  height = 6
)


# ============================================================
# FIGURE 4
# MEAN AFFINITY VS ENSEMBLE VARIABILITY
# ============================================================

library(ggplot2)
library(ggrepel)

fig4 <- ggplot(
  complete_df,
  aes(
    x = Mean_Affinity,
    y = SD_Affinity
  )
) +
  
  # ----------------------------------------------------------
# All ligands
# ----------------------------------------------------------

geom_point(
  size = 2.2,
  alpha = 0.70
) +
  
  # ----------------------------------------------------------
# Highlight top 20 ligands
# ----------------------------------------------------------

geom_point(
  data = top20,
  size = 3.0
) +
  
  # ----------------------------------------------------------
# Labels for top 10 ligands
# ----------------------------------------------------------

geom_text_repel(
  data = top10,
  aes(
    label = Ligand
  ),
  
  # Label positioning
  size = 3.2,
  
  # Avoid excessive overlap
  box.padding = 0.6,
  point.padding = 0.35,
  
  # Strength of repulsion
  force = 2,
  
  # Allow labels to move in all directions
  direction = "both",
  
  # Connecting lines
  min.segment.length = 0,
  segment.size = 0.3,
  segment.alpha = 0.6,
  
  # White background around text improves readability
  bg.color = "white",
  bg.r = 0.15,
  
  # Do not automatically remove labels
  max.overlaps = Inf,
  
  # Keep labels inside plotting area as much as possible
  xlim = c(
    min(complete_df$Mean_Affinity) - 0.5,
    max(complete_df$Mean_Affinity) + 1.5
  ),
  ylim = c(
    min(complete_df$SD_Affinity) - 0.2,
    max(complete_df$SD_Affinity) + 0.8
  )
) +
  
  # ----------------------------------------------------------
# Axis labels and title
# ----------------------------------------------------------

labs(
  title = "Ensemble affinity versus state-to-state variability",
  subtitle = "Top-ranked ligands are highlighted and labeled",
  x = "Mean affinity across five CDK2 states (kcal/mol)",
  y = "Standard deviation of affinity (kcal/mol)"
) +
  
  # ----------------------------------------------------------
# Publication theme
# ----------------------------------------------------------

theme_publication +
  
  # Give labels more room around the plot
  coord_cartesian(
    clip = "off"
  ) +
  
  theme(
    plot.margin = margin(
      t = 10,
      r = 80,
      b = 10,
      l = 10
    )
  )


# ============================================================
# SAVE FIGURE
# ============================================================

save_figure(
  fig4,
  "Figure_4_Mean_Affinity_vs_Variability",
  width = 8,
  height = 6
)

# ============================================================
# FIGURE 5
# CONSENSUS RANK VS MEAN AFFINITY
# ============================================================

rank_df <- complete_df %>%
  arrange(Rank)


fig5 <- ggplot(
  rank_df,
  aes(
    x = Rank,
    y = Mean_Affinity
  )
) +

  geom_point(
    size = 1.8,
    alpha = 0.65
  ) +

  geom_point(
    data = top20,
    size = 2.7
  ) +

  geom_hline(
    yintercept = -10,
    linetype = "dashed",
    linewidth = 0.5
  ) +

  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +

  labs(
    title = "Consensus ranking of CDK2 ensemble-docking ligands",
    subtitle = "Ranking is based on mean affinity across five structural states",
    x = "Consensus rank",
    y = "Mean affinity (kcal/mol)"
  ) +

  theme_publication


save_figure(
  fig5,
  "Figure_5_Consensus_Rank",
  width = 8,
  height = 6
)


# ============================================================
# FIGURE 6
# TOP-10 LIGAND AFFINITY PROFILES ACROSS STATES
# ============================================================

profile_data <- top10 %>%
  select(
    Rank,
    Ligand,
    all_of(RECEPTORS)
  ) %>%
  pivot_longer(
    cols = all_of(RECEPTORS),
    names_to = "Receptor",
    values_to = "Affinity"
  ) %>%
  mutate(
    Receptor = factor(
      Receptor,
      levels = RECEPTORS
    )
  )


fig6 <- ggplot(
  profile_data,
  aes(
    x = Receptor,
    y = Affinity,
    group = Ligand
  )
) +

  geom_line(
    linewidth = 0.65,
    alpha = 0.60
  ) +

  geom_point(
    size = 1.8
  ) +

  facet_wrap(
    ~ Ligand,
    ncol = 2
  ) +

  geom_hline(
    yintercept = -10,
    linetype = "dashed",
    linewidth = 0.4
  ) +

  labs(
    title = "State-dependent affinity profiles of the top 10 ligands",
    subtitle = "Each panel tracks one ligand across the five CDK2 structural states",
    x = "CDK2 structural state",
    y = "Docking affinity (kcal/mol)"
  ) +

  theme_publication +
  theme(
    strip.background = element_rect(
      fill = "grey90",
      color = NA
    ),
    strip.text = element_text(
      face = "bold"
    )
  )


save_figure(
  fig6,
  "Figure_6_Top10_State_Profiles",
  width = 9,
  height = 12
)


# ============================================================
# FIGURE 7
# DISTRIBUTION OF ENSEMBLE MEAN AFFINITIES
# ============================================================

fig7 <- ggplot(
  complete_df,
  aes(
    x = Mean_Affinity
  )
) +

  geom_histogram(
    bins = 30,
    fill = "grey65",
    color = "white"
  ) +

  geom_vline(
    xintercept = -10,
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_vline(
    xintercept = mean(
      complete_df$Mean_Affinity,
      na.rm = TRUE
    ),
    linetype = "dotted",
    linewidth = 0.7
  ) +

  labs(
    title = "Distribution of ensemble mean docking affinities",
    subtitle = "Dashed line: −10 kcal/mol; dotted line: dataset mean",
    x = "Mean affinity across five CDK2 states (kcal/mol)",
    y = "Number of ligands"
  ) +

  theme_publication


save_figure(
  fig7,
  "Figure_7_Ensemble_Affinity_Distribution",
  width = 8,
  height = 6
)


# ============================================================
# 9. EXPORT TOP-LIGAND TABLE
# ============================================================

top20_export <- complete_df %>%
  arrange(Rank) %>%
  slice_head(n = TOP_N)

write_csv(
  top20_export,
  file.path(
    OUTPUT_DIR,
    "top20_consensus_ligands.csv"
  )
)


# ============================================================
# 10. EXPORT LONG-FORMAT DATA FOR FUTURE ANALYSIS
# ============================================================

long_export <- complete_df %>%
  select(
    Rank,
    Ligand,
    Initial_1VYW_Affinity,
    Mean_Affinity,
    Median_Affinity,
    Best_Affinity,
    Worst_Affinity,
    SD_Affinity,
    Successful_Receptors,
    all_of(RECEPTORS)
  ) %>%
  pivot_longer(
    cols = all_of(RECEPTORS),
    names_to = "Receptor",
    values_to = "Affinity"
  )

write_csv(
  long_export,
  file.path(
    OUTPUT_DIR,
    "ensemble_docking_long_format.csv"
  )
)


# ============================================================
# 11. SUMMARY STATISTICS
# ============================================================

summary_table <- complete_df %>%
  summarise(
    Ligands = n(),
    Mean_of_Mean_Affinity = mean(
      Mean_Affinity,
      na.rm = TRUE
    ),
    Median_of_Mean_Affinity = median(
      Mean_Affinity,
      na.rm = TRUE
    ),
    Best_Ensemble_Affinity = min(
      Mean_Affinity,
      na.rm = TRUE
    ),
    Mean_SD = mean(
      SD_Affinity,
      na.rm = TRUE
    )
  )

write_csv(
  summary_table,
  file.path(
    OUTPUT_DIR,
    "ensemble_summary_statistics.csv"
  )
)


# ============================================================
# 12. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CDK2 ENSEMBLE DOCKING VISUALIZATION COMPLETED\n")
cat("============================================================\n")
cat("\n")

cat(
  "Input file: ",
  INPUT_FILE,
  "\n",
  sep = ""
)

cat(
  "Total ligands: ",
  nrow(df),
  "\n",
  sep = ""
)

cat(
  "Complete 5/5 ensemble ligands: ",
  nrow(complete_df),
  "\n",
  sep = ""
)

cat(
  "Best consensus ligand: ",
  complete_df$Ligand[
    which.min(complete_df$Mean_Affinity)
  ],
  "\n",
  sep = ""
)

cat(
  "Best mean affinity: ",
  sprintf(
    "%.3f kcal/mol",
    min(
      complete_df$Mean_Affinity,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat("\n")
cat("Figures saved to:\n")
cat(
  normalizePath(
    OUTPUT_DIR,
    mustWork = FALSE
  ),
  "\n"
)

cat("\n")
cat("Generated figures:\n")
cat("  Figure_1_Top20_Consensus_Ligands\n")
cat("  Figure_2_Top20_State_Heatmap\n")
cat("  Figure_3_Receptor_Affinity_Distributions\n")
cat("  Figure_4_Mean_Affinity_vs_Variability\n")
cat("  Figure_5_Consensus_Rank\n")
cat("  Figure_6_Top10_State_Profiles\n")
cat("  Figure_7_Ensemble_Affinity_Distribution\n")

cat("\n")
cat("============================================================\n")
