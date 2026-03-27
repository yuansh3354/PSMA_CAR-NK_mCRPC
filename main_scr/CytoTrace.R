#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2024-02-27
# ==============================================================================

library(CytoTRACE2)
library(Seurat)

# ==============================================================================
# Custom cytotrace2 function (Batching removed, defaults updated)
# ==============================================================================
custom_cytotrace2 <- function(input, species = "human", is_seurat = TRUE, slot_type = "counts", 
                              parallelize_models = TRUE, parallelize_smoothing = TRUE, 
                              ncores = NULL, seed = 14) {
  
  set.seed(seed)
  message("cytotrace2: Started loading data")
  
  # --- 1. Data Loading ---
  if (is_seurat) {
    if (is.character(input) && file.exists(input)) {
      seurat <- readRDS(input)
    } else {
      seurat <- copy(input)
    }
    data <- loadData_fromSeurat(seurat, slot_type)
  } else {
    if (is.character(input) && file.exists(input)) {
      data <- loadData(input)
    } else {
      data <- copy(input)
    }
  }
  
  if (!is.data.frame(data)) {
    message("Attempting to convert the provided input to the required format.")
    data <- as.data.frame(data)
  }
  
  message("Dataset contains ", dim(data)[1], " genes and ", dim(data)[2], " cells.")
  
  # --- 2. Core Configuration ---
  if (parallelize_smoothing || parallelize_models) {
    if (Sys.info()["sysname"] == "Windows") {
      ncores <- 1
    } else if (is.null(ncores)) {
      ncores <- max(1, parallel::detectCores(logical = TRUE) %/% 2)
    }
  } else if (is.null(ncores)) {
    ncores <- 1
  }
  nc <- ncores
  
  parameter_dict <- readRDS(system.file("extdata", "parameter_dict_19.rds", package = "CytoTRACE2"))
  
  # --- 3. Preprocessing (Full Matrix) ---
  message("cytotrace2: Started preprocessing.")
  input_data <- preprocessData(data, species)
  
  ranked_data <- input_data[[1]]
  log2_data <- input_data[[2]]
  
  gene_names <- colnames(ranked_data)
  cell_names <- rownames(ranked_data)
  
  # --- 4. Prediction (Full Matrix) ---
  message("cytotrace2: Started prediction.")
  predicted_df <- predictData(parameter_dict, ranked_data, log2_data, parallelize_models, ncores = nc)
  predicted_df <- predicted_df[cell_names, ]
  
  # --- 5. Dispersion Calculation ---
  num_genes <- ncol(log2_data)
  dispersion_index <- sapply(1:num_genes, function(i) disp_fn(log2_data[, i]))
  top_genes <- gene_names[order(dispersion_index, decreasing = TRUE)[1:min(1000, num_genes)]]
  
  # --- 6. Postprocessing & Smoothing (Full Matrix) ---
  message("cytotrace2: Started postprocessing.")
  # Note: smooth_batch_size is forced to ncol(data) to compute everything in a single pass
  smoothScore <- smoothData(log2_data, predicted_df, top_genes, 
                            ncores = nc, smooth_batch_size = ncol(data), 
                            parallelize_smoothing = parallelize_smoothing, seed = seed)
  
  smoothScore <- smoothScore[cell_names]
  predicted_df$preKNN_CytoTRACE2_Score <- smoothScore
  
  # --- 7. KNN Smoothing & Potency Binning ---
  if (nrow(log2_data) <= 10) {
    predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
    predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
  } else {
    predicted_df <- binData(predicted_df)
    predicted_df <- predicted_df[cell_names, ]
    
    if (nrow(log2_data) <= 100 || sd(log2_data) == 0) {
      predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
      predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
    } else {
      predicted_df <- smoothDatakNN(log2_data, predicted_df, seed, nc)
      predicted_df <- predicted_df[cell_names, ]
    }
  }
  
  # --- 8. Final Relative Score Calculation ---
  predicted_df$CytoTRACE2_Relative <- (predicted_df$CytoTRACE2_Score - min(predicted_df$CytoTRACE2_Score)) / 
    (max(predicted_df$CytoTRACE2_Score) - min(predicted_df$CytoTRACE2_Score))
  
  predicted_df <- predicted_df[c("CytoTRACE2_Score", "CytoTRACE2_Potency", 
                                 "CytoTRACE2_Relative", "preKNN_CytoTRACE2_Score", "preKNN_CytoTRACE2_Potency")]
  
  # --- 9. Return Results ---
  if (!is_seurat) {
    message("cytotrace2: Finished")
    return(predicted_df)
  } else {
    predicted_df <- predicted_df[colnames(seurat), ]
    seurat <- AddMetaData(object = seurat, metadata = predicted_df)
    message("cytotrace2: Finished")
    return(seurat)
  }
}

# ==============================================================================
# Execute Custom CytoTRACE2
# ==============================================================================
# Call the simplified function. Default arguments handle species and object type.
cytotrace2_sce <- custom_cytotrace2(seurat_obj, ncores = 64)