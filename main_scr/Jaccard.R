#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2026-03-26
# ==============================================================================

library(pheatmap)
library(eoffice)
library(ggplotify)

# ==============================================================================
# 1. Define global grouping information (lists)
# ==============================================================================
g1 <- c('NK0', 'NK1', 'NK1A', 'NK1B', 'NK1C', 'NK2', 'NK3', 'NKint')

g2 <- c('c59_NK-FCGR3A', 'c60_NK-SELL', 'c61_NK-FCGR3A|IFNG', 'c62_NK-DUSP4', 
        'c63_NK-CD160', 'c64_NK-CD160|IFNG', 'c65_NK-CD160|ICAM1', 'c66_NK-MKI67', 'c67_NK-IL7R')

g3 <- c('CD56dimCD16hi-c1-IL32', 'CD56dimCD16hi-c2-CX3CR1', 'CD56dimCD16hi-c3-ZNF90', 
        'CD56dimCD16hi-c4-NFKBIA', 'CD56dimCD16hi-c5-MKI67', 'CD56dimCD16hi-c6-DNAJB1', 
        'CD56dimCD16hi-c7-NR4A3', 'CD56dimCD16hi-c8-KLRC2', 'CD56brightCD16lo-c1-GZMH', 
        'CD56brightCD16lo-c2-IL7R-RGS1lo', 'CD56brightCD16lo-c3-CCL3', 
        'CD56brightCD16lo-c4-IL7R', 'CD56brightCD16lo-c5-CREM')

g4 <- c('NK_CD16_GZMB', 'NK_CD16_HSPA1A', 'NK_CD16_ISG15', 'NK_CD16_KLRC1', 
        'NK_CD16_MKI67', 'NK_CD16_TXNIP', 'NK_CD56_CD160', 'NK_CD56_CD74', 
        'NK_CD56_FOS', 'NK_CD56_RGS1', 'NK_CD56_SELL')

ordered_cols <- c(g1, g2, g3, g4)

# Define colors (globally applicable)
ann_colors <- list(
  Study_Source = c(
    "Nature_2024" = "#8DD3C7", 
    "CancerCell_2022" = "#FFFFB3", 
    "Cell_2023" = "#BEBADA", 
    "CancerCell_2025" = "#FB8072"
  )
)

# ==============================================================================
# 2. Define a generic plotting function
# ==============================================================================
draw_custom_heatmap <- function(file_path, title_text) {
  
  # A. Read data
  mat <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  # B. Data alignment (handle missing columns)
  valid_cols <- intersect(ordered_cols, colnames(mat))
  mat <- mat[, valid_cols]
  
  # C. Dynamically build annotation bars
  # Calculate 'times' based on the actual number of valid_cols
  len_g1 <- length(intersect(g1, valid_cols))
  len_g2 <- length(intersect(g2, valid_cols))
  len_g3 <- length(intersect(g3, valid_cols))
  len_g4 <- length(intersect(g4, valid_cols))
  
  annotation_col <- data.frame(
    Study_Source = factor(rep(
      c("Nature_2024", "CancerCell_2022", "Cell_2023", "CancerCell_2025"),
      times = c(len_g1, len_g2, len_g3, len_g4)
    ))
  )
  rownames(annotation_col) <- valid_cols
  
  # D. Calculate gap locations
  gap_locations <- cumsum(c(len_g1, len_g2, len_g3))
  
  # E. Generate pheatmap object (note: silent = TRUE)
  p <- pheatmap(mat,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                annotation_col = annotation_col,
                annotation_colors = ann_colors,
                gaps_col = gap_locations,
                display_numbers = TRUE,
                number_format = "%.2f",
                fontsize_number = 6,
                color = colorRampPalette(c("#007FB8", "white", "#FF0000"))(100),
                border_color = "white",
                cellwidth = 15,
                cellheight = 15,
                fontsize_row = 10,
                fontsize_col = 8,
                angle_col = 90,
                main = title_text,
                silent = TRUE # Do not draw the plot directly, return the object instead
  )
  
  # F. Convert to ggplot object for export
  return(as.ggplot(p))
}

# ==============================================================================
# 3. Execute plotting and save
# ==============================================================================

# Generate objects for both plots
p_all <- draw_custom_heatmap(file_all, "All CarNK: Consistency Index")
p_infusion <- draw_custom_heatmap(file_infusion, "Infusion CarNK: Consistency Index")

# 1. Save the first plot (overwrite/create mode)
topptx(figure = p_all, 
       filename = output_ppt, 
       width = 20, 
       height = 8)

# 2. Append the second plot (append mode)
topptx(figure = p_infusion, 
       filename = output_ppt, 
       width = 20, 
       height = 8, 
       append = TRUE)

print(paste("Successfully saved PPT to:", output_ppt))