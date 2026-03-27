#!/usr/bin/env Rscript
# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2026-03-26
# ==============================================================================

yuansh_cellchat <- function(input_obj,
                        assay = NULL,
                        group.by = NULL,
                        species = c('human', 'mouse'),
                        CellChatDB.use = NULL, # A character vector, which is a subset of c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", "Non-protein Signaling")
                        PPIuse = FALSE,
                        type = "triMean",      # c("triMean", "truncatedMean", "thresholdedMean", "median")
                        min.cells = 10) {
  
  # Initialize CellChat object
  cellchat.obj <- createCellChat(input_obj, assay = assay, group.by = group.by)
  
  # Select database based on species
  if (species == 'human') {
    CellChatDB <- CellChatDB.human
    ppi <- PPI.human
  }
  
  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
    ppi <- PPI.mouse
  }
  
  # Set or subset the CellChat database
  if (is.null(CellChatDB.use)) {
    cellchat.obj@DB <- CellChatDB
  } else {
    CellChatDB <- subsetDB(CellChatDB, search = CellChatDB.use, key = "annotation")
    cellchat.obj@DB <- CellChatDB
  }
  
  # Data preprocessing and gene expression analysis
  cellchat.obj <- subsetData(cellchat.obj) 
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  # Compute communication probabilities
  if (PPIuse == FALSE) {
    cellchat.obj <- computeCommunProb(cellchat.obj, type = type)
  } else {
    cellchat.obj <- projectData(cellchat.obj, ppi)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use = FALSE, type = type)
  }
  
  # Filter communications and compute pathway probabilities
  cellchat.obj <- filterCommunication(cellchat.obj, min.cells = min.cells)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  
  # Aggregate the interaction network (using defaults for all sources and targets)
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  # Compute network centrality
  cellchat.obj <- netAnalysis_computeCentrality(cellchat.obj, slot.name = "netP")
  
  return(cellchat.obj)
}
