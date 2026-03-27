#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2026-03-26
# ==============================================================================

library(Seurat)
library(qs)
library(CellChat)
library(dplyr)
library(parallel)

# --- Hardware Parallelism Parameters ---
io_threads <- 64 
mc_cores <- 8

# ==============================================================================
# 1. Helper Function: Safe Merge (Compatible with Seurat V4/V5)
# ==============================================================================
SafeJoin <- function(obj) {
  if (inherits(obj[["RNA"]], "StdAssay")) {
    if (length(Layers(obj)) > 2) {
      message("JoinLayers for Seurat V5...")
      obj <- JoinLayers(obj)
    }
  }
  return(obj)
}

# ==============================================================================
# 2. Data Loading
# ==============================================================================
message(">>> Loading Datasets...")
obj_names <- c("seurat_car", "seurat_cd4", "seurat_cd8", "seurat_pbmc")
obj_list <- lapply(obj_names, function(x) qread(paste0(base_path, x, ".qs"), nthreads = io_threads))
names(obj_list) <- c("CAR", "CD4", "CD8", "PBMC")

# Filter and format CAR-NK specific metadata
obj_list$CAR <- subset(obj_list$CAR, subset = sort == 'Sorted CAR-NK cells')
obj_list$CAR$majority_voting <- "CarNK"
obj_list$CAR$subcluster_voting <- paste0("Car_", obj_list$CAR$subcluster_voting)

# ==============================================================================
# 3. Immune Landscape Analysis (With Skip Logic)
# ==============================================================================
treatments <- c("M.R.", "S.R.")

mclapply(treatments, function(treat) {
  result <- tryCatch({
    message(paste("\n[Process] Treatment Group:", treat))
    
    # Define output file paths
    file_maj <- paste0(output_path, "Chat_Total_Immune_", treat, "_majority_voting.qs")
    file_sub <- paste0(output_path, "Chat_Total_Immune_", treat, "_subcluster_voting.qs")
    
    # Extract subsets for the current treatment
    sub_objs <- lapply(obj_list, function(x) {
      if(any(x$treatment == treat)) subset(x, subset = treatment == treat) else NULL
    })
    sub_objs <- Filter(Negate(is.null), sub_objs)
    if (length(sub_objs) < 2) return(NULL)
    
    # --- 3.1 Majority Level Analysis ---
    if (!file.exists(file_maj)) {
      message(">> Level: majority_voting (Calculating...)")
      
      # Note: Ensure 'sub_objs_maj' is properly defined in your environment before this step
      immune_maj <- merge(sub_objs_maj[[1]], y = sub_objs_maj[-1], add.cell.ids = names(sub_objs_maj))
      immune_maj <- SafeJoin(immune_maj)
      
      if (length(unique(as.character(immune_maj$majority_voting))) > 1) {
        chat_maj <- KS_cellchat(immune_maj, group.by = "majority_voting", species = "human")
        qsave(chat_maj, file = file_maj, nthreads = io_threads)
      }
      rm(sub_objs_maj, immune_maj); gc()
    } else {
      message(">> Level: majority_voting (File exists, skipping)")
    }
    
    # --- 3.2 & 3.3 Subcluster Level Analysis (Including Pairwise) ---
    # Check if all relevant subcluster files exist for the current treatment
    pair_targets <- if ("CAR" %in% names(sub_objs)) {
      paste0(output_path, "Chat_Pair_CAR_", setdiff(names(sub_objs), "CAR"), "_", treat, "_subcluster.qs")
    } else { c() }
    all_sub_files <- c(file_sub, pair_targets)
    
    if (!all(file.exists(all_sub_files))) {
      message(">> Level: subcluster (Calculating missing files...)")
      
      # Full subcluster interactions
      if (!file.exists(file_sub)) {
        immune_sub <- merge(sub_objs[[1]], y = sub_objs[-1], add.cell.ids = names(sub_objs))
        immune_sub <- SafeJoin(immune_sub)
        
        if (length(unique(as.character(immune_sub$subcluster_voting))) > 1) {
          chat_sub <- KS_cellchat(immune_sub, group.by = "subcluster_voting", species = "human")
          qsave(chat_sub, file = file_sub, nthreads = io_threads)
        }
        rm(immune_sub); gc()
      }
      
      # Pairwise comparison analysis
      if ("CAR" %in% names(sub_objs)) {
        for (n in setdiff(names(sub_objs), "CAR")) {
          f_pair <- paste0(output_path, "Chat_Pair_CAR_", n, "_", treat, "_subcluster.qs")
          if (!file.exists(f_pair)) {
            message(paste(">> Pairwise: CAR vs", n))
            pair_sub <- merge(sub_objs[["CAR"]], y = sub_objs[[n]], add.cell.ids = c("CAR", n))
            pair_sub <- SafeJoin(pair_sub)
            
            if (length(unique(as.character(pair_sub$subcluster_voting))) > 1) {
              chat_pair <- KS_cellchat(pair_sub, group.by = "subcluster_voting", species = "human")
              qsave(chat_pair, file = f_pair, nthreads = io_threads)
            }
            rm(pair_sub); gc()
          }
        }
      }
    } else {
      message(">> Level: subcluster (All files exist, skipping)")
    }
  }, error = function(e) { message(paste("[Error]:", treat, e$message)); return(e$message) })
  return(result)
}, mc.cores = mc_cores)

# Clean up memory before processing the tumor dataset
rm(obj_list); gc()

# ==============================================================================
# 4. Tumor Dataset Processing
# ==============================================================================
message("\n>>> Processing Tumor Dataset...")
seurat_tumor <- qread(paste0(base_path, "seurat_tumor.qs"), nthreads = io_threads)

mclapply(c("Before treatment", "After treatment"), function(grp) {
  f_tumor <- paste0(output_path, "Chat_Tumor_", gsub(" ", "_", grp), "_subcluster.qs")
  
  if (!file.exists(f_tumor)) {
    tryCatch({
      sub_tumor <- subset(seurat_tumor, subset = group == grp)
      sub_tumor <- SafeJoin(sub_tumor)
      
      if (length(unique(as.character(sub_tumor$subcluster_voting))) > 1) {
        chat_tumor <- KS_cellchat(sub_tumor, group.by = "subcluster_voting", species = "human")
        qsave(chat_tumor, file = f_tumor, nthreads = io_threads)
      }
    }, error = function(e) message(paste("[Error Tumor]:", e$message)))
  } else {
    message(paste(">> Tumor:", grp, "(File exists, skipping)"))
  }
}, mc.cores = mc_cores)

message(">>> All tasks completed.")