#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2025-12-27
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup & Library Loading
# ------------------------------------------------------------------------------
rm(list = ls())
gc(reset = TRUE, verbose = FALSE)
options(future.globals.maxSize = 1000 * 1024^3) # Set max size to ~1TB
script_start_time <- Sys.time()

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(eoffice)
  library(schard)
  library(tidyverse)
})

# ------------------------------------------------------------------------------
# 2. Global Color 
# ------------------------------------------------------------------------------
group_colors <- c(
  "Baseline" = "#3C5488", "LC" = "#8491B4", "Day0" = "#E64B35", "Day3" = "#ff5B99",
  "Day5" = "#9467BD", "Day7" = "#00A087", "Day14" = "#4DBBD5", "Day21" = "#91D1C2",
  "Day60" = "#7E6148", "CAR" = "#B09C85", "Day75" = "#B09C85",
  "Before treatment" = "#6baed6", "After treatment" = "#fd8d3c"
)

# ------------------------------------------------------------------------------
# 3. Helper Functions
# ------------------------------------------------------------------------------
# 3.1 Custom DotPlot function for consistent aesthetic across different objects
generate_custom_dotplot <- function(seurat_obj, features, group_by_col = 'subcluster_voting', title_text = "Marker Genes Expression") {
  p <- DotPlot(object = seurat_obj, features = features, group.by = group_by_col, col.min = -1, col.max = 1) + 
    scale_color_gradientn(colours = c("#007FB8", "#FFFFFF", "#FF0000")) +
    scale_size(limits = c(0, 100), range = c(1, 6)) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(colour = "grey90", linetype = "dashed"), 
      panel.grid.minor = element_blank(), 
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    guides(
      color = guide_colorbar(title = "Avg Exp", frame.colour = "black", ticks.colour = "black"),
      size = guide_legend(title = "Pct Exp")
    ) +
    labs(title = title_text, x = "", y = "Identity")
  
  return(p)
}

if(!file.exists(qs_file)) {
  seurat_obj <- schard::h5ad2seurat(h5ad_file)
  qsave(seurat_obj, qs_file, nthreads = 512)
} else {
  seurat_obj <- qread(qs_file, nthreads = 512)
}

meta.tumor <- seurat_obj@meta.data
meta.tumor <- meta.tumor[!is.na(meta.tumor$cdr3s_aa) & meta.tumor$majority_voting == 'CD8T', ]


tcr_after_expansion <- meta.tumor %>%
  filter(!is.na(cdr3s_aa) & cdr3s_aa != "") %>%
  filter(group == "After treatment") %>%
  group_by(patient, cdr3s_aa) %>%
  summarise(After_Count = n(), .groups = "drop") %>%
  filter(After_Count > 2) %>%
  arrange(patient, desc(After_Count))

tumor_expanded_clones <- meta.tumor %>%
  filter(!is.na(cdr3s_aa) & cdr3s_aa != "") %>%
  group_by(patient, cdr3s_aa) %>%
  summarise(Tumor_After_Count = n(), .groups = "drop") %>%
  filter(Tumor_After_Count >= 1)


for(task in dotplot_tasks) {
  if(file.exists(task$file)) {
    tmp_obj <- schard::h5ad2seurat(task$file)
    p_dp <- generate_custom_dotplot(tmp_obj, task$features)
    topptx(figure = p_dp, filename = paste0(dir_name, task$out_name), width = task$width, height = 4)
    rm(tmp_obj); gc(verbose = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 6. Cell Type Enrichment / Proportional Distribution (RoE Analysis)
# ------------------------------------------------------------------------------
# Note: Assuming 'mysplit' is defined in 'script/Functions/Functions.R'
# For safety in environments without it, we construct arrays manually.
condition_lvls <- c("M.R. Baseline", "M.R. Day0", "M.R. Day3", "M.R. Day7", "M.R. Day14", "M.R. Day21",
                    "S.R. Baseline", "S.R. Day0", "S.R. Day3", "S.R. Day7", "S.R. Day14", "S.R. Day21")

df_meta <- seurat_obj@meta.data
target_groups <- c("Day0", "Day7", "Baseline", "Day3", "Day14", "Day21")
df_meta <- df_meta[df_meta$group %in% target_groups, ]
df_meta$condition <- paste(df_meta$treatment, df_meta$group)

# Generate overall enrichment plot
p_roe <- distribution_Patient_Mean_Prop(
  meta_data = df_meta,
  celltype_column = "subcluster_voting",
  patient_column = "patient",
  condition_column = "condition",
  celltype_level = sort(unique(df_meta$subcluster_voting)),
  condition_level = condition_lvls
)

# ------------------------------------------------------------------------------
# 7. Clinical Correlation Dynamics (PSA_FC vs Cell Proportions)
# ------------------------------------------------------------------------------
target_days <- c("Baseline", "Day3", "Day7", "Day14", "Day21")
target_cell_type <- "CD8T_Tem_CMC1"
plot_list <- list()

# Note: Depends on 'sort_t' being preloaded from earlier scripts.
if (exists("sort_t")) {
  for (day in target_days) {
    day_prop <- sort_t %>%
      dplyr::filter(group == day & majority_voting == "CD8T") %>%
      dplyr::group_by(patient) %>%
      dplyr::summarise(
        proportion = sum(subcluster_voting == target_cell_type) / n(),
        total_cells = n(),
        .groups = "drop"
      ) %>%
      dplyr::filter(total_cells > 10)
    
    day_plot_data <- day_prop %>%
      dplyr::inner_join(clinical_df, by = "patient") %>%
      dplyr::filter(PSA_FC < 1)
    
    if (nrow(day_plot_data) >= 3) {
      p_cor <- ggplot(day_plot_data, aes(x = proportion, y = PSA_FC)) +
        geom_smooth(method = "lm", color = "black", alpha = 0.2, fill = "gray") +
        geom_point(size = 3, color = "#073248") +
        stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
        theme_bw() +
        labs(title = paste0("Timepoint: ", day), x = "Proportion in CD8T", y = "PSA_FC") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      plot_list[[day]] <- p_cor
    }
  }
  
  if (length(plot_list) > 0) {
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = 3) +
      patchwork::plot_annotation(
        title = paste0("Correlation between ", target_cell_type, " Proportion and PSA_FC"),
        subtitle = "Denominator: Total CD8T cells per patient per day",
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )

    message(">>> Success: Combined correlation plot exported.")
  }
}