#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2026-03-26
# ==============================================================================

{
  library(dplyr)
  
  # ==============================================================================
  # 1. Define Mapping Relationships
  # ==============================================================================
  # Classification mapping based on the blood system.
  # 'Myeloid' includes all myeloid cells; 'B' includes B and plasma cells.
  mapping <- c(
    "B cells"           = "B",
    "Plasma cells"      = "B",
    "CD4T"              = "CD4T",
    "CD8T"              = "CD8T",
    "NK cells"          = "NK",
    "Macrophages"       = "Myeloid",
    "Monocytes"         = "Myeloid",
    "Neutrophils"       = "Myeloid",
    "DC cells"          = "Myeloid",
    # Tumor-specific components retain original names or are grouped as 'Other'
    "Endothelial cells" = "Endothelial",
    "Fibroblasts"       = "Fibroblast",
    "Tumor"             = "Tumor"
  )
  
  # ==============================================================================
  # 2. Execute Classification Conversion
  # ==============================================================================
  # Use the recode function for rapid replacement and create a new column 'blood_style_voting'
  meta.tumor$blood_style_voting <- dplyr::recode(meta.tumor$majority_voting, !!!mapping)
  
  # ==============================================================================
  # 3. Result Validation
  # ==============================================================================
  # Check the quantity distribution after conversion
  print(table(meta.tumor$blood_style_voting))
  
  # Verify presence of Platelet or Erythrocyte (usually rare in tumor data).
  # To force alignment with 7 blood categories, unlisted types can be set to NA.
  target_cats <- c("Myeloid", "CD8T", "Platelet", "CD4T", "Erythrocyte", "B", "NK")
  meta.tumor$majority_voting = meta.tumor$blood_style_voting
  
  # ==============================================================================
  # 1. Build Comprehensive ID Mapping Table
  # ==============================================================================
  # Extract existing mapping relationships (patient -> patientid) from clinical_df
  patient_lookup <- setNames(as.character(clinical_df$patientid), 
                             as.character(clinical_df$patient))
  
  # ==============================================================================
  # 2. Execute Replacement Operation
  # ==============================================================================
  # Perform direct replacement using the lookup vector.
  # If names in meta.tumor$patient are not in the lookup table, NA is returned.
  meta.tumor$patient <- patient_lookup[as.character(meta.tumor$patient)]
  
  # ==============================================================================
  # 3. Check and Validate Results
  # ==============================================================================
  # Check for patients that failed conversion (i.e., resulting in NA)
  missing_after_fix <- unique(meta.tumor$patient[is.na(meta.tumor$patient)])
  
  if(length(missing_after_fix) > 0) {
    message("Warning: The following patients still lack corresponding PatientIDs:")
    print(missing_after_fix)
  } else {
    message("All patient IDs converted successfully!")
  }
  
  # View final unique IDs after conversion
  print(unique(meta.tumor$patient))
  meta.tumor$sample = paste(meta.tumor$patient, gsub(' treatment', '', meta.tumor$group))
}

p1 = visualizeCellTypeDistribution(meta.tumor, 'sample', 'majority_voting') +
  scale_fill_manual(values = major_colors) + 
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray95"),
    panel.grid.major = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

{
  # Convert obs to data.table to process large-scale data efficiently
  dt_obs <- as.data.table(obs)
  dt_obs = dt_obs[dt_obs$majority_voting != 'Myeloid',]
  
  # Define sorting order
  sort_levels <- c("PBMC", "Sorted CAR-NK cells", "Sorted T cells")
  group_levels <- c("Baseline", "LD", "Day0", "Day3", "Day5", "Day7", 
                    "Day14", "Day21", "Day60", "Day75", "Before treatment", "After treatment")
  
  # Convert to factors and set levels
  dt_obs[, sort := factor(sort, levels = sort_levels)]
  dt_obs[, group := factor(group, levels = group_levels)]
  # Note: patientid is already a factor, keeping its original or alphabetical order
  
  # Calculate proportion of majority_voting within each sample.
  # Must retain metadata columns used for sorting (sort, patientid, group)
  sample_composition <- dt_obs[, .(cell_count = .N), 
                               by = .(sample, sort, patientid, group, majority_voting)]
  
  # Calculate percentages
  sample_composition[, proportion := cell_count / sum(cell_count), by = sample]
  
  # --- 2. Determine X-axis sorting logic ---
  # Create a sorting auxiliary column to ensure plotting order follows: sort -> patientid -> group
  sample_composition <- sample_composition[order(sort, patientid, group)]
  
  # Convert sample to an ordered factor to prevent ggplot from applying default alphabetical sorting
  sample_levels <- unique(sample_composition$sample)
  sample_composition[, sample := factor(sample, levels = sample_levels)]
  }

# ==============================================================================
# 2. Plotting Preparation
# ==============================================================================
group_map <- setNames(as.character(sample_composition$group), sample_composition$sample)
group_map[group_map == "Day0"] <- "2Hours"

# ==============================================================================
# 2. Plotting (Applying Label Mappings)
# ==============================================================================
p_composition <- ggplot(sample_composition[sample_composition$sort == 'PBMC',], 
                        aes(x = sample, y = proportion, fill = majority_voting)) +
  # Draw stacked bar plot
  geom_bar(stat = "identity", width = 0.8) +
  # Apply previously defined colors
  scale_fill_manual(values = major_colors) +
  
  # --- Core modification: Replace X-axis display labels ---
  scale_x_discrete(labels = group_map) + 
  
  # --- Facet settings ---
  facet_nested(
    . ~ sort + patientid,
    scales = "free_x",
    space = "free_x",
    strip = strip_nested(
      background_x = elem_list_rect(
        fill = list(sort_colors = sort_colors, patientid = patient_colors)
      )
    )
  ) +
  
  theme_classic() +
  theme(
    # Rotate X-axis labels to prevent overlap
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.title.x = element_blank(),
    
    # Clean up styles
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Cell Type Composition Across PBMC Samples",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  )

{
  # Convert obs to data.table to process large-scale data efficiently
  dt_obs <- as.data.table(obs)
  dt_obs = dt_obs[dt_obs$majority_voting != 'Myeloid',]
  
  # Define sorting order
  sort_levels <- c("PBMC", "Sorted CAR-NK cells", "Sorted T cells")
  group_levels <- c("Baseline", "LD", "Day0", "Day3", "Day5", "Day7", 
                    "Day14", "Day21", "Day60", "Day75", "Before treatment", "After treatment")
  
  # Convert to factors and set levels
  dt_obs[, sort := factor(sort, levels = sort_levels)]
  dt_obs[, group := factor(group, levels = group_levels)]
  
  # Calculate proportion of subcluster_voting within each sample.
  sample_composition <- dt_obs[, .(cell_count = .N), 
                               by = .(sample, sort, patientid, group, subcluster_voting)]
  
  # Calculate percentages
  sample_composition[, proportion := cell_count / sum(cell_count), by = sample]
  
  # --- 2. Determine X-axis sorting logic ---
  sample_composition <- sample_composition[order(sort, patientid, group)]
  
  # Convert sample to ordered factor
  sample_levels <- unique(sample_composition$sample)
  sample_composition[, sample := factor(sample, levels = sample_levels)]
  }

p_composition_sortt <- ggplot(sample_composition[sample_composition$sort == 'Sorted T cells',], 
                              aes(x = sample, y = proportion, fill = subcluster_voting)) +
  # Draw stacked bar plot
  geom_bar(stat = "identity", width = 0.8) +
  # Apply previously defined colors
  scale_fill_manual(values = subCluster_voting_colors) +
  
  # --- Core modification: Replace X-axis display labels ---
  scale_x_discrete(labels = group_map) + 
  
  # --- Facet settings ---
  facet_nested(
    . ~ sort + patientid,
    scales = "free_x",
    space = "free_x",
    strip = strip_nested(
      background_x = elem_list_rect(
        fill = list(sort_colors = sort_colors, patientid = patient_colors)
      )
    )
  ) +
  
  theme_classic() +
  theme(
    # Rotate X-axis labels
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.title.x = element_blank(),
    
    # Clean up styles
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Cell Type Composition Across PBMC Samples",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  )

# Display plot
print(p_composition_sortt)

# Note: Overwriting variable p_composition_sortt for CAR-NK cells
p_composition_sortt <- ggplot(sample_composition[sample_composition$sort == 'Sorted CAR-NK cells',], 
                              aes(x = sample, y = proportion, fill = subcluster_voting)) +
  # Draw stacked bar plot
  geom_bar(stat = "identity", width = 0.8) +
  # Apply previously defined colors
  scale_fill_manual(values = subCluster_voting_colors) +
  
  # --- Core modification: Replace X-axis display labels ---
  scale_x_discrete(labels = group_map) + 
  
  # --- Facet settings ---
  facet_nested(
    . ~ sort + patientid,
    scales = "free_x",
    space = "free_x",
    strip = strip_nested(
      background_x = elem_list_rect(
        fill = list(sort_colors = sort_colors, patientid = patient_colors)
      )
    )
  ) +
  
  theme_classic() +
  theme(
    # Rotate X-axis labels
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.title.x = element_blank(),
    
    # Clean up styles
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Cell Type Composition Across PBMC Samples",
    y = "Proportion of Total Cells",
    fill = "Cell Type"
  )

# Display plot
print(p_composition_sortt)

output_pptx <- paste0(dir_name, "Cell_Composition_Analysis.pptx")
topptx(
  p_composition, 
  filename = output_pptx, 
  width = 20,    # Recommended wider width due to numerous facets
  height = 6, 
  append = FALSE # FALSE means create a new file; TRUE means append to existing
)
topptx(
  p_composition_sortt, 
  filename = output_pptx, 
  width = 20, 
  height = 6, 
  append = FALSE 
)
topptx(
  p1, 
  filename = output_pptx, 
  width = 20, 
  height = 6, 
  append = TRUE 
)

####################################################################################
{
  # ==============================================================================
  # 2. Basic Configuration 
  # ==============================================================================
  # Color mappings
  group_colors <- c(
    "Baseline" = "#3C5488", "LD" = "#8491B4", "Day0" = "#E64B35", 
    "Day3" = "#ff5B99", "Day5" = "#9467BD", "Day7" = "#00A087", 
    "Day14" = "#4DBBD5", "Day21" = "#91D1C2", "Day60" = "#7E6148", 
    "Day75" = "#B09C85", "Before" = "#6baed6", "After" = "#fd8d3c"
  )
  
  # Sorting logic for timepoints
  blood_time_levels <- c("Baseline", "LD", "Day0", "Day3", "Day5", "Day7", "Day14", "Day21", "Day60", "Day75")
  tumor_time_levels <- c("Before", "After")
  
  # ==============================================================================
  # 3. Data Loading and Alignment 
  # ==============================================================================
  
  blood    <- obs[obs$sort == 'PBMC', ]
  sort_t   <- obs[obs$sort == 'Sorted T cells', ]
  sort_cat <- obs[obs$sort == 'Sorted CAR-NK cells', ]
  
  # Tumor data
  meta.tumor <- seurat_obj@meta.data
  
  # Unify mappings
  meta.tumor$patient <- patient_lookup[as.character(meta.tumor$patient)]
  meta.tumor$group_clean <- gsub(" treatment", "", meta.tumor$group)
  
  # ==============================================================================
  # 4. Plotting Function Definition 
  # ==============================================================================
  # Function to generate temporal bar plots, adding geom_text to display specific counts
  plot_temporal <- function(df, title_str, levels_vec) {
    # Preprocess statistical data
    plot_data <- df %>%
      group_by(group) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(group = factor(group, levels = levels_vec)) %>%
      filter(!is.na(group))
    
    ggplot(plot_data, aes(x = group, y = Count, fill = group)) +
      geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
      # Add count annotations: vjust = -0.3 places text above bar, check_overlap prevents clutter
      geom_text(aes(label = scales::comma(Count)), vjust = -0.3, size = 2.5) + 
      scale_fill_manual(values = group_colors) +
      theme_classic() +
      theme(
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(size = 9, face = "bold"),
        # Adjust margins to accommodate text
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      ) +
      # Slightly increase Y-axis upper limit to prevent text cutoff
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
      labs(title = title_str)
  }
  
  # ==============================================================================
  # 5. Generate Subplots
  # ==============================================================================
  # P1: Patient distribution overlap matrix heatmap
  p_blood_ids <- unique(as.character(blood$patientid))
  p_t_ids     <- unique(as.character(sort_t$patientid))
  p_nk_ids    <- unique(as.character(sort_cat$patientid))
  p_tumor_ids <- unique(as.character(meta.tumor$patient))
  all_patients <- sort(unique(c(p_blood_ids, p_t_ids, p_nk_ids, p_tumor_ids)))
  
  p1 <- data.frame(Patient = all_patients) %>%
    mutate(`Blood (PBMC)` = Patient %in% p_blood_ids,
           `Sorted T cells` = Patient %in% p_t_ids,
           `Sorted CAR-NK` = Patient %in% p_nk_ids,
           `Tumor` = Patient %in% p_tumor_ids) %>%
    pivot_longer(cols = -Patient, names_to = "Dataset", values_to = "Present") %>%
    ggplot(aes(x = Dataset, y = Patient, fill = Present)) +
    geom_tile(color = "white", lwd = 0.8) +
    scale_fill_manual(values = c("TRUE" = "#3C5488", "FALSE" = "#E9E9E9")) +
    geom_text(aes(label = ifelse(Present, "Yes", "")), color = "white", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.title = element_blank(), legend.position = "none", panel.grid = element_blank()) +
    labs(title = "Patient Overlap Matrix")
  
  # P2-P4: Blood lineage distribution
  p2 <- plot_temporal(blood, "Blood (PBMC) Counts", blood_time_levels)
  p3 <- plot_temporal(sort_t, "Sorted T Cells Counts", blood_time_levels)
  p4 <- plot_temporal(sort_cat, "Sorted CAR-NK Counts", blood_time_levels)
  
  # P5: Tumor distribution (Fix rename error: overwrite group column directly in mutate)
  p5 <- meta.tumor %>% 
    mutate(group = group_clean) %>% 
    plot_temporal("Tumor Sample Counts", tumor_time_levels)
  
  # ==============================================================================
  # 6. Final Assembly and Saving
  # ==============================================================================
  # Combine plots using patchwork logic
  combined_plot <- (p1 | (p2 / p3 / p4 / p5)) + 
    plot_layout(widths = c(1.2, 1)) +
    plot_annotation(
      title = 'Clinical ScRNA-seq Data Distribution with Exact Cell Counts',
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  print(combined_plot)
  
  topptx(combined_plot, paste0(dir_name, "final_data_summary_with_counts.pptx"), 
         width = 13, height = 11)
  
  message("Visualization completed, including specific bar chart values.")
}