################################################################################
# --- 1. Global Environment Setup ---
library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(rstatix)
library(ggpubr)
library(clusterProfiler)
library(msigdbr)
library(AUCell)
library(patchwork)
library(eoffice)

options(future.globals.maxSize = 1000 * 1024^3)
dir_name <- "analysis_output"
if (!dir.exists(dir_name)) dir.create(dir_name)

# --- 2. Data Loading and Preprocessing ---
# Load main infusion object and reference product object
seurat_obj <- qread('infusion.qs', nthreads = 16)
infusion <- qread('infusion_carnk.qs', nthreads = 16)
car_product <- qread('NK.qs', nthreads = 16)
# Anonymize specific patient identifiers
seurat_obj$patientid <- factor(seurat_obj$patientid)

# Filter timepoints for kinetic analysis
infusion_meta <- infusion@meta.data %>%
  filter(group %in% c("Day0", "Day3", "Day7", "Day14")) %>%
  mutate(condition = paste(treatment, group))

# --- 3. Quality Control (QC) Visualization ---
# Extract QC metrics for violin plots
qc_data <- FetchData(seurat_obj, vars = c(
  'doublet_score', 'subcluster_voting', 'n_genes',
  'total_counts', 'pct_counts_mt'
))
qc_data$genes_per_umi <- qc_data$total_counts / qc_data$n_genes

# Define standard QC plotting function
plot_qc_violin <- function(data, metric, title_label) {
  ggplot(data, aes(x = subcluster_voting, y = .data[[metric]])) +
    geom_violin(trim = FALSE, scale = "width", 
                alpha = 0.6, fill = '#1F77B4', color = NA) +
    geom_boxplot(width = 0.15, fill = "white", 
                 outlier.shape = NA, color = "black") +
    theme_bw() +
    labs(title = title_label, x = NULL, y = metric) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate and save QC plots
qc_metrics <- list(
  "Doublet_Score" = "doublet_score",
  "Gene_Count"    = "n_genes",
  "UMI_Count"     = "total_counts",
  "MT_Percent"    = "pct_counts_mt"
)

lapply(names(qc_metrics), function(n) {
  p <- plot_qc_violin(qc_data, qc_metrics[[n]], n)
  ggsave(file.path(dir_name, paste0("QC_", n, ".png")), p, width = 8, height = 6)
})

# --- 4. Dimensionality Reduction (UMAP) Batch Processing ---
# Function to automate DimPlot generation and legend export
auto_umap_save <- function(obj, group_var, color_vec, file_prefix) {
  # Clean UMAP plot
  p <- DimPlot(obj, group.by = group_var, cols = color_vec, order = TRUE) +
    theme_void() + NoLegend() + labs(title = NULL)
  
  ggsave(file.path(dir_name, paste0(file_prefix, ".pdf")),
         p, width = 6, height = 6)
  
  # Export legend to PowerPoint
  actual_ids <- unique(obj[[group_var, drop = TRUE]])
  legend_data <- as.data.frame(table(obj[[group_var, drop = TRUE]]))
  colnames(legend_data) <- c("ID", "Freq")
  
  p_legend <- ggplot(legend_data, aes(ID, Freq, color = ID)) +
    geom_point(size = 3) +
    scale_color_manual(values = color_vec) +
    theme_void()
  
  topptx(p_legend, file.path(dir_name, paste0(file_prefix, "_legend.pptx")))
}

# Run UMAP visualization for key variables
umap_vars <- list(
  list("patientid", "p1_umap_patient"),
  list("group", "p3_umap_group"),
  list("treatment", "p5_umap_treatment"),
  list("subcluster_voting", "p9_umap_subcluster")
)

# Note: Assumes color vectors (patient_colors, etc.) are pre-defined in the environment
# for (v in umap_vars) auto_umap_save(seurat_obj, v[[1]], get(paste0(v[[1]], "_colors")), v[[2]])

# --- 5. Differential Abundance and Kinetic Trends ---
# Calculate sub-cluster proportion dynamics over time
proportion_trends <- infusion@meta.data %>%
  mutate(group = factor(group, levels = c("Day0", "Day3", "Day7", "Day14"))) %>%
  filter(!is.na(group)) %>%
  group_by(patientid, group, subcluster_voting) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patientid, group) %>%
  mutate(proportion = (n / sum(n)) * 100) %>%
  ungroup()

# Calculate global mean trend
mean_trends <- proportion_trends %>%
  group_by(group, subcluster_voting) %>%
  summarise(mean_prop = mean(proportion), .groups = "drop")

# Visualize dynamics with individual trajectories and mean trend
p_trends <- ggplot() +
  geom_line(data = proportion_trends, aes(x = group, y = proportion, group = patientid), 
            color = "grey85", alpha = 0.6) +
  geom_line(data = mean_trends, aes(x = group, y = mean_prop, group = 1), 
            color = "orange", size = 1.2) +
  geom_point(data = mean_trends, aes(x = group, y = mean_prop), color = "orange") +
  facet_wrap(~subcluster_voting, scales = "free_y", ncol = 4) +
  theme_bw() +
  labs(title = "Sub-cluster Proportion Dynamics", x = "Timepoint", y = "Proportion (%)")

ggsave(file.path(dir_name, "kinetic_trends.pdf"), p_trends, width = 12, height = 8)

# --- 6. Clinical Correlation Analysis ---
# Correlate sub-cluster frequencies with PSA Fold Change (Efficacy)
# Extract specific sub-cluster data at baseline (Day 0)
baseline_corr_data <- proportion_trends %>%
  filter(group == "Day0", subcluster_voting == "CD16_NKFBIA") %>%
  inner_join(clinical_df, by = "patientid")

p_corr <- ggplot(baseline_corr_data, aes(x = proportion, y = PSA_FC)) +
  geom_point(color = "#3C5488", size = 3) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black", fill = "lightgrey") +
  stat_cor(method = "spearman", color = "red") +
  labs(title = "Day0 CD16_NKFBIA vs PSA Fold Change", x = "Cell Frequency (%)", y = "PSA FC") +
  theme_bw()

# --- 7. Functional Scoring (AUCell) ---
# Calculate functional gene set activity across sub-clusters
# rankings calculation
expr_matrix <- GetAssayData(infusion, assay = "RNA", layer = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 16, plotStats = FALSE)

# Compute AUC for specific CAR-NK gene lists
# Assumes CARNK_gene_list is defined
pathway_auc <- AUCell_calcAUC(CARNK_gene_list, cells_rankings, nCores = 16)
auc_matrix <- as.data.frame(t(getAUC(pathway_auc)))

# Integrate AUC scores into metadata
infusion <- AddMetaData(infusion, metadata = auc_matrix)

# Visualize AUC activity as bar charts
auc_summary <- infusion@meta.data %>%
  select(subcluster_voting, all_of(names(CARNK_gene_list))) %>%
  group_by(subcluster_voting) %>%
  summarise(across(everything(), median), .groups = "drop")

plot_auc_bar <- function(data, feature) {
  global_med <- median(data[[feature]])
  ggplot(data, aes(x = .data[[feature]], y = subcluster_voting, fill = subcluster_voting)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_vline(xintercept = global_med, linetype = "dashed") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "#F2F2F2", color = NA), legend.position = "none") +
    labs(title = gsub("_", " ", feature), x = "AUC Score", y = NULL)
}

auc_plots <- lapply(names(CARNK_gene_list), function(f) plot_auc_bar(auc_summary, f))
p_auc_final <- wrap_plots(auc_plots, ncol = 4)
ggsave(file.path(dir_name, "auc_activity_grid.pdf"), p_auc_final, width = 16, height = 8)

# --- 8. Gene Set Enrichment Analysis (GSEA) ---
# Perform GSEA on sub-cluster markers using MSigDB Hallmark sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

# Run GSEA for all clusters
gsea_results_list <- list()
clusters <- unique(degs$cluster) # Assumes degs (markers) is available

for (clus in clusters) {
  cluster_genes <- degs %>%
    filter(cluster == clus) %>%
    arrange(desc(avg_log2FC))
  
  gene_list <- cluster_genes$avg_log2FC
  names(gene_list) <- cluster_genes$gene
  
  res <- GSEA(geneList = gene_list, TERM2GENE = hallmark_sets, pvalueCutoff = 0.05)
  if (!is.null(res) && nrow(res@result) > 0) gsea_results_list[[as.character(clus)]] <- res
}

# --- Final Save ---
qsave(infusion, file.path(dir_name, "infusion_final_processed.qs"))
print("Analysis pipeline completed successfully.")