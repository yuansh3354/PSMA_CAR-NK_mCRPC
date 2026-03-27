#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2025-11-27
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment Setup
# ------------------------------------------------------------------------------
# Clean the workspace
rm(list = ls())
# Free up memory silently
gc(reset = TRUE, verbose = FALSE)

# Create output directory for figures if it does not exist
dir_name <- "figure/"
dir.exists(dir_name) || dir.create(dir_name)

# ------------------------------------------------------------------------------
# 1. Utility Functions for Visualization
# ------------------------------------------------------------------------------

# Function to extract and save a color legend to a PowerPoint file
save_legend_pptx <- function(vec_colors, value_vector, outfile) {
  df <- as.data.frame(table(value_vector))
  colnames(df) <- c("value", "Freq")   # Force consistent column names
  
  p_legend <- ggplot(df, aes(value, Freq, color = value)) +
    geom_point(size = 3) +
    scale_color_manual(values = vec_colors) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  topptx(p_legend, outfile)
}

# Function to clean up the theme of a DimPlot (e.g., UMAP/tSNE)
clean_dimplot <- function(p, colors) {
  p +
    scale_color_manual(values = colors) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank()
    ) +
    NoLegend() +
    labs(title = NULL)
}

# Set default plot dimensions
w <- 6
h <- 6

# ------------------------------------------------------------------------------
# 2. Data Loading and Subsetting
# ------------------------------------------------------------------------------
library(qs)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(eoffice)
library(ggrepel)
library(tidyr)
library(tibble)
library(pheatmap)

# Load infusion CAR-NK data
qs_file = 'infusion.qs'
infusion = qread(qs_file, nthreads = 512)
obs = qread(qs_file, nthreads = 64)

# Ensure patient ID is a factor
obs$patientid = factor(obs$patientid)

# Subset data based on sorting strategy
blood = obs[obs$sort == 'PBMC',]
sort_cat = obs[obs$sort == 'Sorted CAR-NK cells',]
sort_t = obs[obs$sort == 'Sorted T cells',]

# ------------------------------------------------------------------------------
# 3. Cell Type Composition Analysis (Majority Voting)
# ------------------------------------------------------------------------------
{
  # Convert 'obs' to data.table for efficient large-scale data manipulation
  dt_obs <- as.data.table(obs)
  
  # Remove 'Myeloid' cells for this specific visualization
  dt_obs <- dt_obs[dt_obs$majority_voting != 'Myeloid',]
  
  # Define the order of levels for sorting and grouping
  sort_levels <- c("PBMC", "Sorted CAR-NK cells", "Sorted T cells")
  group_levels <- c("Baseline", "LD", "Day0", "Day3", "Day5", "Day7", 
                    "Day14", "Day21", "Day60", "Day75", "Before treatment", "After treatment")
  
  # Apply factor levels to the data
  dt_obs[, sort := factor(sort, levels = sort_levels)]
  dt_obs[, group := factor(group, levels = group_levels)]
  
  # Calculate cell counts for each cell type within each sample
  # Must retain metadata columns used for sorting (sort, patientid, group)
  sample_composition <- dt_obs[, .(cell_count = .N), 
                               by = .(sample, sort, patientid, group, majority_voting)]
  
  # Calculate proportions relative to total cells in each sample
  sample_composition[, proportion := cell_count / sum(cell_count), by = sample]
  
  # Create a sorting logic for the X-axis: sort -> patientid -> group
  sample_composition <- sample_composition[order(sort, patientid, group)]
  
  # Convert 'sample' to an ordered factor to prevent ggplot from sorting alphabetically
  sample_levels <- unique(sample_composition$sample)
  sample_composition[, sample := factor(sample, levels = sample_levels)]
}

# Define facet strip colors for nested plotting
# First level (sort) uses different colors; Second level (patientid) uses a uniform grey
sort_strip_colors <- c(
  "PBMC" = "#FFB6C1",        
  "Whole Blood" = "#ADD8E6", 
  "Bone Marrow" = "#98FB98"  
)
patient_strip_color <- "grey90"

# Generate stacked bar plot for cell composition
p_composition <- ggplot(sample_composition, aes(x = sample, y = proportion, fill = majority_voting)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = major_colors) + # Note: major_colors must be defined in your environment
  
  # Use facet_nested to create grouped strip labels
  facet_nested(
    . ~ sort + patientid,
    scales = "free_x",
    space = "free_x",
    # Define custom background colors for the nested strips
    strip = strip_nested(
      background_x = elem_list_rect(
        fill = list(
          sort_colors,    # Colors for the first level (sort)
          patient_colors  # Colors for the second level (patientid)
        ),
        color = "white"   # Border color for the strip
      )
    )
  ) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.title.x = element_blank(),
    strip.background = element_blank(), # Remove default strip background to show custom colors
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Cell Type Composition Across All Samples",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  )

# Export composition plot to PowerPoint
output_pptx <- paste0(dir_name, "Cell_Composition_Analysis.pptx")
{
  topptx(
    p_composition, 
    filename = output_pptx, 
    width = 20,    # Large width to accommodate many facets
    height = 6, 
    append = FALSE # Create a new file
  )
}

# ------------------------------------------------------------------------------
# 4. Cell Subcluster Composition Analysis (Subcluster Voting)
# ------------------------------------------------------------------------------
{
  dt_obs <- as.data.table(obs)
  
  sort_levels <- c("PBMC", "Sorted CAR-NK cells", "Sorted T cells")
  group_levels <- c("Baseline", "LD", "Day0", "Day3", "Day5", "Day7", 
                    "Day14", "Day21", "Day60", "Day75", "Before treatment", "After treatment")
  
  dt_obs[, sort := factor(sort, levels = sort_levels)]
  dt_obs[, group := factor(group, levels = group_levels)]
  
  # Calculate cell counts and proportions for subclusters
  sample_composition <- dt_obs[, .(cell_count = .N), 
                               by = .(sample, sort, patientid, group, subcluster_voting)]
  
  sample_composition[, proportion := cell_count / sum(cell_count), by = sample]
  
  # Order the data for plotting
  sample_composition <- sample_composition[order(sort, patientid, group)]
  sample_levels <- unique(sample_composition$sample)
  sample_composition[, sample := factor(sample, levels = sample_levels)]
}

# Generate stacked bar plot for subcluster composition
p_composition_sub <- ggplot(sample_composition, aes(x = sample, y = proportion, fill = subcluster_voting)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = subCluster_voting_colors) + # Note: subCluster_voting_colors must be defined
  
  facet_nested(
    . ~ sort + patientid,
    scales = "free_x",
    space = "free_x",
    strip = strip_nested(
      background_x = elem_list_rect(
        fill = list(
          sort_colors,    
          patient_colors  
        ),
        color = "white" 
      )
    )
  ) +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Cell Subcluster Composition Across All Samples",
    y = "Proportion of Total Cells",
    fill = "Cell Subcluster"
  )

# Append subcluster composition plot to the existing PowerPoint file
{
  topptx(
    p_composition_sub, 
    filename = output_pptx, 
    width = 20, 
    height = 10, 
    append = TRUE # Append to the previously created file
  )
}

print(p_composition_sub)

# ------------------------------------------------------------------------------
# 5. Global Cell Type Proportions (Pie Chart)
# ------------------------------------------------------------------------------
{
  # Summarize global cell type counts and calculate proportions
  pie_data_sorted <- obs %>%
    group_by(majority_voting) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>%
    arrange(desc(proportion)) # Sort in descending order
  
  # Lock the factor levels to maintain sorted order in the plot
  pie_data_sorted$majority_voting <- factor(pie_data_sorted$majority_voting, 
                                            levels = pie_data_sorted$majority_voting)
  
  # Format labels with percentages
  pie_data_sorted <- pie_data_sorted %>%
    mutate(percentage = paste0(round(proportion * 100, 2), "%"),
           label = paste0(majority_voting, ": ", percentage))
  
  # Calculate midpoint coordinates for label placement
  # Note: cumulative sum logic matches the direction = -1 in coord_polar
  pie_data_sorted <- pie_data_sorted %>%
    mutate(cumulative = cumsum(proportion),
           midpoint = cumulative - proportion / 2)
  
  # Create pie chart (counter-clockwise layout)
  p_pie_final <- ggplot(pie_data_sorted, aes(x = "", y = proportion, fill = majority_voting)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y", start = 0, direction = -1) + # Convert to polar coordinates
    scale_fill_manual(values = major_colors) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    # Use ggrepel to draw labels for specific cell types to avoid overlapping
    geom_text_repel(
      data = subset(pie_data_sorted, majority_voting %in% c("CD4T", 'Myeloid', 'NK', "CD8T")),
      aes(y = midpoint, label = label),
      size = 4,
      show.legend = FALSE,
      nudge_x = 0.9,       # Push labels outside the pie
      segment.size = 0.5,  # Line thickness
      segment.color = "grey50",
      direction = "y",
      hjust = 0.5,
      force = 3            # Increase repulsion force
    ) +
    labs(fill = "Cell Type")
  
  # Export pie chart
  topptx(p_pie_final, filename = paste0(dir_name, "Blood_Sorted_CCW_Pie.pptx"), width = 11, height = 8)
  print(p_pie_final)
}

# ------------------------------------------------------------------------------
# 6. Myeloid Subcluster Proportions (Pie Chart)
# ------------------------------------------------------------------------------
{
  # Subset to Myeloid cells only
  df = obs[obs$major_subcluster == 'Myeloid',]

  # Calculate proportions for Myeloid subclusters
  pie_data_sorted <- df %>%
    group_by(subcluster_voting) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>%
    arrange(desc(proportion))
  
  # Lock factor levels
  pie_data_sorted$subcluster_voting <- factor(pie_data_sorted$subcluster_voting, 
                                              levels = pie_data_sorted$subcluster_voting)
  
  # Format labels and calculate midpoints
  pie_data_sorted <- pie_data_sorted %>%
    mutate(percentage = paste0(round(proportion * 100, 2), "%"),
           label = paste0(subcluster_voting, ": ", percentage),
           cumulative = cumsum(proportion),
           midpoint = cumulative - proportion / 2)
  
  # Create Myeloid subcluster pie chart
  p_pie_myeloid <- ggplot(pie_data_sorted, aes(x = "", y = proportion, fill = subcluster_voting)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y", start = 0, direction = -1) +
    # Add an extra color for Neutrophils dynamically
    scale_fill_manual(values = c(subCluster_voting_colors, Neutrophil = "#eeea39")) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    labs(fill = "Myeloid Subclusters")
  
  print(p_pie_myeloid)
}

# ------------------------------------------------------------------------------
# 7. Subcluster Correlation Heatmap Across Timepoints
# ------------------------------------------------------------------------------
{
  # Calculate subcluster frequencies grouped by specified timepoints (batches)
  # Note: mysplit() is a custom function assumed to return a character vector
  props_batch <- blood[blood$group %in% mysplit("Day0,Day7,Baseline,Day3,Day14,Day21"),] %>%
    group_by(group, subcluster_voting) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(freq = n / sum(n))
}

{
  # Reshape long data to wide format: rows = subclusters, columns = groups (timepoints)
  cor_matrix_data <- props_batch %>%
    select(group, subcluster_voting, freq) %>%
    pivot_wider(names_from = group, values_from = freq, values_fill = 0) %>%
    column_to_rownames("subcluster_voting")
  
  # Calculate Pearson correlation coefficient matrix between different timepoints
  batch_cor <- cor(cor_matrix_data) 
  
  # Plot correlation heatmap to evaluate compositional stability across batches
  p_heat <- pheatmap(
    batch_cor, 
    display_numbers = TRUE,                       # Display correlation values in cells
    color = colorRampPalette(Strength.cls)(100),  # Use custom color palette defined in environment
    main = "Correlation of Cell Composition Between Timepoints",
    angle_col = 45,                               # Rotate x-axis labels to prevent overlap
    fontsize = 10,
    border_color = "white"                        # White borders for cells
  )
  
  print(p_heat)
}
