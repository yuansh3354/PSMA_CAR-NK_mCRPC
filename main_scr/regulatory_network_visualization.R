#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2026-03-26
# ==============================================================================

{
  # --- 1. Check and ensure variable matching ---
  # Verify that the previously defined variable for the mapping table exists
  if (!exists("regulons_gene_lists")) {
    stop("Error: 'regulons_gene_lists' not found. Please run the loom processing step first.")
  }
  
  # --- 2. Preprocess transcription factor list ---
  # Extract the intersection of transcription factors present in the SCENIC results
  tf_to_process <- intersect(tf_list_plus, names(regulons_gene_lists))
  
  # --- 3. Extract TF-Target relationships in parallel ---
  # Ensure target gene expression data is extracted from the RNA Assay
  DefaultAssay(seurat_obj) <- "RNA"
  n_cores <- max(1, parallel::detectCores() - 2)
  
  message(paste0(">>> Starting parallel extraction (", n_cores, " cores)..."))
  
  tf_target_results <- parallel::mclapply(tf_to_process, function(tf) {
    # Extract target genes from the correct variable 'regulons_gene_lists'
    targets <- regulons_gene_lists[[tf]]
    valid_targets <- intersect(targets, rownames(seurat_obj))
    
    if(length(valid_targets) > 0) {
      # Extract expression statistics for filtering
      dp_data <- suppressMessages(
        DotPlot(seurat_obj, features = valid_targets, group.by = 'subcluster_voting')$data
      )
      
      # Filter target genes with significant features in this subcluster
      sig_targets <- dp_data %>%
        dplyr::filter(pct.exp > 75 & avg.exp.scaled > 1) %>%
        dplyr::pull(features.plot) %>%
        unique() %>%
        as.character()
      
      if(length(sig_targets) > 0) {
        return(data.frame(
          from = tf,
          to = sig_targets,
          weight = 1,
          type = "TF-Target"
        ))
      }
    }
    return(NULL)
  }, mc.cores = n_cores)
  
  # Merge parallel computation results
  tf_target_edges <- do.call(rbind, tf_target_results)
  
  # --- 4. Construct comprehensive regulatory network graph ---
  # Calculate activity correlation between TFs
  auc_mat <- as.matrix(seurat_obj[["SCENICplus"]]@data)
  
  # Ensure row names contain the (+) identifier for matching
  if(!any(grepl("\\(\\+\\)", rownames(auc_mat)))) {
    rownames(auc_mat) <- paste0(rownames(auc_mat), "(+)")
  }
  
  # Calculate Spearman correlation for TFs
  tf_tf_cor <- cor(t(auc_mat[intersect(tf_list_plus, rownames(auc_mat)), ]), method = "spearman")
  
  # Convert to edge list and filter connections with strength > 0.3
  tf_tf_edges <- reshape2::melt(tf_tf_cor) %>%
    dplyr::filter(Var1 != Var2 & value > 0.3) %>%
    dplyr::rename(from = Var1, to = Var2, weight = value) %>%
    dplyr::mutate(type = "TF-TF")
  
  # Remove duplicated undirected edges
  tf_tf_edges <- tf_tf_edges[!duplicated(t(apply(tf_tf_edges[,1:2], 1, sort))), ]
  
  # Combine all edge data (TF-TF and TF-Target)
  all_edges <- dplyr::bind_rows(tf_target_edges, tf_tf_edges)
  
  # Output result statistics
  message(">>> Regulatory network construction complete.")
  message(paste0(">>> Total Edges: ", nrow(all_edges)))
}
{
  # --- 1. Load necessary libraries ---
  library(Seurat)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(ggplotify)
  library(eoffice)
  
  # --- 2. Basic parameters and mapping table preparation ---
  # 2.1 Define time point order
  time_order <- c("Day0", "Day3", "Day7", "Day14")
  
  tf_list_plus <- paste0(unlist(tf_list), "(+)")
  
  # 2.3 Construct gene-module mapping table and clean names (prevent XML export errors)
  clean_names <- function(x) {
    x <- gsub("[^[:alnum:]_]", "_", x)
    x <- gsub("^([0-9])", "G_\\1", x)
    return(x)
  }
  
  module_map <- stack(tf_list)
  colnames(module_map) <- c("TF", "Module")
  
  gene_annotation_clean <- all_edges %>%
    dplyr::filter(from %in% tf_list_plus) %>%
    dplyr::mutate(TF_clean = gsub("\\(\\+\\)", "", from),
                  gene_clean = gsub("\\(\\+\\)", "", to)) %>%
    dplyr::select(TF_clean, gene_clean) %>%
    dplyr::rename(gene = gene_clean, TF = TF_clean) %>%
    dplyr::left_join(module_map, by = "TF") %>%
    dplyr::filter(gene %in% rownames(seurat_obj)) %>% 
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(Module = clean_names(Module)) # Clean module names
  
  # --- 3. Core matrix generation (Source for plot_mat_complex) ---
  # 3.1 Calculate Minibulk expression (aggregated by subcluster and time point group)
  exp_matrix_raw <- AverageExpression(seurat_obj, 
                                      features = gene_annotation_clean$gene, 
                                      group.by = c("subcluster_voting", "group"), 
                                      slot = "data")$RNA
  
  # 3.2 Z-score standardization
  # The scale function operates on columns, requiring double transposition to standardize genes (rows)
  plot_mat_complex <- t(scale(t(as.matrix(exp_matrix_raw))))
  # Remove NaNs generated during standardization (genes with constant 0 expression across all groups)
  plot_mat_complex <- plot_mat_complex[complete.cases(plot_mat_complex), ]
  
  # --- 4. Axis sorting and clustering logic ---
  # 4.1 Column sorting: group by subcluster, then temporally Day0 -> Day14 within subclusters
  col_info <- data.frame(ID = colnames(plot_mat_complex)) %>%
    dplyr::mutate(Group = clean_names(sub(".*_", "", ID)),
                  Subcluster = clean_names(sub("_[^_]+$", "", ID))) %>%
    dplyr::mutate(Group = factor(Group, levels = clean_names(time_order))) %>%
    dplyr::arrange(Subcluster, Group)
  
  plot_mat_complex <- plot_mat_complex[, col_info$ID]
  
  col_anno_final <- col_info %>%
    dplyr::select(Subcluster, Group) %>%
    `rownames<-`(col_info$ID)
  
  # 4.2 Row sorting: perform clustering within each module
  row_anno <- gene_annotation_clean %>%
    dplyr::select(gene, Module) %>%
    tibble::column_to_rownames("gene")
  
  final_gene_order <- c()
  row_gaps <- c()
  current_idx <- 0
  
  for (mod in unique(row_anno$Module)) {
    genes_in_mod <- rownames(row_anno)[row_anno$Module == mod]
    genes_in_mod <- intersect(genes_in_mod, rownames(plot_mat_complex))
    
    if (length(genes_in_mod) > 1) {
      hc <- hclust(dist(plot_mat_complex[genes_in_mod, ]))
      final_gene_order <- c(final_gene_order, genes_in_mod[hc$order])
    } else {
      final_gene_order <- c(final_gene_order, genes_in_mod)
    }
    current_idx <- current_idx + length(genes_in_mod)
    row_gaps <- c(row_gaps, current_idx)
  }
  
  # --- 5. Color palette and visualization settings ---
  academic_heatmap_cols <- rev(c("#B2182B", "#D6604D", "#F4A582", "#F7F7F7", "#92C5DE", "#4393C3", "#2166AC"))
  
  # Use the provided variable for subcluster colors
  clean_sub_colors <- subCluster_voting_colors
  names(clean_sub_colors) <- clean_names(names(clean_sub_colors))
  
  anno_colors <- list(
    Group = setNames(c("#E64B35", "#ff5B99", "#00A087", "#4DBBD5"), clean_names(time_order)),
    Module = setNames(brewer.pal(length(unique(row_anno$Module)), "Set3"), unique(row_anno$Module)),
    Subcluster = clean_sub_colors[unique(col_info$Subcluster)]
  )
  
  # Calculate column gaps between subclusters
  col_gaps <- which(!duplicated(col_info$Subcluster))[-1] - 1
  
  p <- pheatmap(plot_mat_complex[final_gene_order, ],
                main = "Complex Minibulk TF-Target Regulation",
                color = colorRampPalette(academic_heatmap_cols)(100),
                annotation_row = row_anno,
                annotation_col = col_anno_final,
                annotation_colors = anno_colors,
                cluster_rows = FALSE, # Manually ordered via intra-module clustering
                cluster_cols = FALSE, # Manually ordered by subcluster and time point
                gaps_row = row_gaps,
                gaps_col = col_gaps,
                show_colnames = FALSE,
                show_rownames = FALSE,
                border_color = NA,
                silent = TRUE)
  
  message(">>> Task completed: Heatmap exported to figure_final folder.")
  
}
{
  # --- 1. Load necessary libraries ---
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(RColorBrewer)
  library(eoffice)
  
  # --- 2. Core parameters and regulatory mapping preparation ---
  # Define time point order
  time_order <- c("Day0", "Day3", "Day7", "Day14")
  
  tf_list_plus <- paste0(unlist(tf_list), "(+)")
  
  # Construct mapping table and clean names
  module_map <- stack(tf_list)
  colnames(module_map) <- c("TF", "Module")
  
  gene_annotation_clean <- all_edges %>%
    dplyr::filter(from %in% tf_list_plus) %>%
    dplyr::mutate(TF_clean = gsub("\\(\\+\\)", "", from),
                  gene_clean = gsub("\\(\\+\\)", "", to)) %>%
    dplyr::select(TF_clean, gene_clean) %>%
    dplyr::rename(gene = gene_clean, TF = TF_clean) %>%
    dplyr::left_join(module_map, by = "TF") %>%
    dplyr::filter(gene %in% rownames(seurat_obj)) %>% 
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(Module))
  
  # --- 3. Calculate Minibulk expression matrix ---
  # Obtain all involved transcription factors and target genes
  unique_genes <- unique(c(gene_annotation_clean$TF, gene_annotation_clean$gene))
  
  # Aggregate by subcluster and time point
  exp_minibulk <- AverageExpression(seurat_obj, 
                                    features = unique_genes, 
                                    group.by = c("subcluster_voting", "group"), 
                                    slot = "data")$RNA %>% as.matrix()
  
  # --- 4. Calculate relative expression activity and correlation ---
  # Calculate the maximum average expression of each gene across all groups
  relative_activity <- apply(exp_minibulk, 1, max)
  
  # Prepare correlation edge table
  network_data <- gene_annotation_clean %>%
    dplyr::filter(TF %in% rownames(exp_minibulk) & gene %in% rownames(exp_minibulk))
  
  # Calculate Pearson correlation between TF and Target across Minibulk samples
  cor_values <- sapply(1:nrow(network_data), function(i) {
    cor(exp_minibulk[network_data$TF[i], ], exp_minibulk[network_data$gene[i], ])
  })
  network_data$correlation <- cor_values
  
  # --- 5. Node statistics and correction (for outliers) ---
  # Mark as significant if involved in any edge with |r| > 0.6
  significant_nodes <- network_data %>%
    dplyr::filter(abs(correlation) > 0.6) %>%
    dplyr::select(TF, gene) %>%
    unlist() %>%
    unique()
  
  gene_stats <- data.frame(
    gene = names(relative_activity),
    rel_exp = relative_activity,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      is_sig = gene %in% significant_nodes,
      # Use log1p to correct outliers (e.g., 70+), ensuring uniform node size distribution
      rel_exp_scaled = log1p(rel_exp),
      plot_color_scaled = ifelse(is_sig, log1p(rel_exp), NA),
      display_label = ifelse(is_sig, gene, "")
    )
  
  # --- 6. Construct network graph attributes ---
  network_data <- network_data %>%
    dplyr::mutate(
      abs_cor = abs(correlation),
      # Dark gray for positive correlation, light gray for negative correlation
      edge_type = ifelse(correlation > 0, "Positive", "Negative")
    )
  
  # --- 7. Render network graph ---
  graph_final <- graph_from_data_frame(network_data, directed = TRUE, vertices = gene_stats)
  
  p_net_final <- ggraph(graph_final, layout = "fr") +
    # Edge settings: thickness and transparency mapped to absolute correlation value
    geom_edge_link(aes(edge_width = abs_cor, edge_alpha = abs_cor, color = edge_type), 
                   show.legend = TRUE) +
    scale_edge_color_manual(values = c("Positive" = "#636363", "Negative" = "#D9D9D9"),
                            name = "Correlation Type") +
    scale_edge_width_continuous(range = c(0.2, 1), name = "Abs Correlation") +
    
    # Node settings: size and color mapped to corrected Log expression
    geom_node_point(aes(size = rel_exp_scaled, color = plot_color_scaled)) +
    
    # Label settings: display only highly correlated nodes
    geom_node_text(aes(label = display_label), 
                   repel = TRUE, size = 1.75, fontface = "bold", max.overlaps = 3) +
    
    # Color gradient and node size mapping
    scale_color_gradientn(colors = c("#4393C3", "#F7F7F7", "#B2182B"), 
                          na.value = "#ECECEC", name = "Log1p(Max Exp)") +
    scale_size_continuous(range = c(1, 6), name = "Node Size (Log Exp)") +
    
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 6),
          legend.position = "right") +
    labs(title = "Optimized TF-Target Regulatory PPI Network",
         subtitle = "Node Size/Color: Log1p Max Minibulk Exp; Gray Nodes: Low Corr (<0.6)")
  
  print(p_net_final)
  
  # --- 8. Export to PPTX ---
  if(!dir.exists("figure_final")) dir.create("figure_final")
  topptx(p_net_final, filename = "figure_final/Optimized_TF_Network_Complete.pptx", width = 12, height = 11)
}
