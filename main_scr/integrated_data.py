import gc
import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import issparse

# ==========================================
# 1. Environment Setup & Data Loading
# ==========================================
# Configure figure parameters
sc.set_figure_params(dpi=80, frameon=False, vector_friendly=True)

# Define file paths
input_path = "results/01.scanpy_multi/integrated_data.h5ad"
anno_file_path = 'results/annotated_result.csv'

# Load AnnData object
print(f"Loading data from {input_path}...")
adata = sc.read_h5ad(input_path, chunk_size=30 * 6000)

# ==========================================
# 2. Memory Usage Analysis Function
# ==========================================
def inspect_adata_memory(adata):
    """Analyzes and prints the memory consumption of the AnnData object."""
    print(f"=== AnnData Memory Analysis (Total Obs: {adata.n_obs}) ===")
    mem_dict = {}
    
    # 1. Analyze .X (Expression matrix)
    if adata.X is not None:
        if issparse(adata.X):
            mem_dict['.X'] = adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes
            print("Note: .X is a Sparse Matrix")
        else:
            mem_dict['.X'] = adata.X.nbytes
            print("Warning: .X is a Dense Array - High memory consumption")

    # 2. Analyze other components
    mem_dict['.obs'] = adata.obs.memory_usage(deep=True).sum()
    mem_dict['.var'] = adata.var.memory_usage(deep=True).sum()
    mem_dict['.obsm'] = sum(val.nbytes for val in adata.obsm.values())
    mem_dict['.obsp'] = sum((v.data.nbytes + v.indices.nbytes + v.indptr.nbytes) if issparse(v) else v.nbytes for v in adata.obsp.values())
    
    print("-" * 30)
    total_mem = 0
    for k, v in sorted(mem_dict.items(), key=lambda x: x[1], reverse=True):
        gb_size = v / (1024**3)
        total_mem += gb_size
        print(f"{k}: {gb_size:.4f} GB")
    
    print(f"=== Total Estimated: {total_mem:.2f} GB ===\n")

inspect_adata_memory(adata)

# ==========================================
# 3. Metadata Synchronization & Cleaning
# ==========================================
print("Synchronizing external annotations...")
df_anno = pd.read_csv(anno_file_path, index_col=0)
adata.obs = adata.obs.join(df_anno[['majority_voting', 'consistency']], how='left')

# Fix spelling errors and set categorical order for 'group'
if 'group' in adata.obs.columns:
    # Rename 'Baeline' to 'Baseline'
    adata.obs['group'] = adata.obs['group'].replace('Baeline', 'Baseline').astype('category')
    adata.obs['group'] = adata.obs['group'].cat.remove_unused_categories()
    
    # Define clinical timeline order
    target_order = ['Baseline', 'LD', 'Day0', 'Day3', 'Day7', 'Day14', 'Day21', 'Day60', 'Day75']
    valid_order = [x for x in target_order if x in adata.obs['group'].cat.categories]
    remaining = list(set(adata.obs['group'].cat.categories) - set(valid_order))
    
    adata.obs['group'] = adata.obs['group'].cat.reorder_categories(valid_order + remaining, ordered=True)
    print(f"Group order set: {list(adata.obs['group'].cat.categories)}")

# ==========================================
# 4. Sample Classification (Sort Type Logic)
# ==========================================
def define_sort_type(row):
    """Categorizes cells based on patient ID and sample string patterns."""
    patient = str(row['patient'])
    sample_lower = str(row['sample']).lower()
    
    if patient == 'CAR':
        return 'CAR-NK Product'
    if '_carnk_' in sample_lower:
        return 'Sorted CAR-NK cells'
    if any(pattern in sample_lower for pattern in ['_lh_', '-lh', '_lh']):
        return 'PBMC'
    if any(pattern in sample_lower for pattern in ['_tc_', '-t']):
        return 'Sorted T cells'
    return 'Other'

print("Adding 'sort' metadata column...")
adata.obs['sort'] = adata.obs.apply(define_sort_type, axis=1).astype('category')
print(adata.obs['sort'].value_counts())

# ==========================================
# 5. Cell Filtering (Contaminant Removal)
# ==========================================
print(f"Cells before filtering: {adata.n_obs}")

# Filter based on sorting source and expected cell types
# Rule: Sorted CAR-NK should primarily contain ILC/NK; Sorted T cells should be CD4/CD8
mask_sorted_carnk = (adata.obs['sort'] == 'Sorted CAR-NK cells') & (adata.obs['majority_voting'] == 'NK')
mask_sorted_t = (adata.obs['sort'] == 'Sorted T cells') & (adata.obs['majority_voting'].isin(['CD4T', 'CD8T']))
mask_others = ~adata.obs['sort'].isin(['CAR-NK Product', 'Sorted CAR-NK cells', 'Sorted T cells'])

adata = adata[mask_sorted_carnk | mask_sorted_t | mask_others].copy()
print(f"Cells after filtering: {adata.n_obs}")

# ==========================================
# 6. Leiden Cluster Refinement & Renumbering
# ==========================================
old_cats = temp_series.cat.categories
new_cats = [str(i) for i in range(len(old_cats))]
adata.obs['leiden_majority_voting'] = temp_series.cat.rename_categories(new_cats)

# ==========================================
# 7. Visualization (UMAP & Dotplot)
# ==========================================
# Define domain-specific marker genes
marker_genes = {
    "Lymphoid": ["CD3D", "CD3E", "CD4", "CD8A", "GNLY", "NKG7", "MS4A1", "JCHAIN"],
    "Myeloid": ["CD14", "LYZ", "FCGR3A", "MS4A7", "FCER1A", "LILRA4"],
    "Others": ["HBA1", "GYPA", "PPBP", "PF4", "CD34", "MKI67"]
}
# Filter markers to keep only those present in the dataset
valid_markers = {k: [g for g in v if g in adata.var_names] for k, v in marker_genes.items()}

# Plot UMAP overview
sc.pl.umap(
    adata, 
    color=['majority_voting', 'leiden_majority_voting'],
    legend_loc='on data', 
    legend_fontsize=8, 
    legend_fontoutline=2,
    show=False,
    save='_final_annotation_umap.png'
)

# Plot Marker Gene Dotplot
sc.pl.dotplot(
    adata, 
    valid_markers, 
    groupby='leiden_majority_voting', 
    standard_scale='var', 
    dendrogram=True,
    save='_final_marker_dotplot.png'
)

# ==========================================
# 8. Data Export
# ==========================================
# Ensure UMAP coordinates exist before saving
if 'X_umap' not in adata.obsm:
    print("Recalculating neighbors and UMAP...")
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.umap(adata)

output_h5ad = "results/01.scanpy_multi/integrated_data_l1_count.h5ad"
print(f"Saving finalized AnnData to {output_h5ad}...")
adata.write_h5ad(output_h5ad)

# Release memory
gc.collect()
print("Pipeline execution completed successfully.")