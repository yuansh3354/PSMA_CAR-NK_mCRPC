import gc
import os
import pickle
import scanpy.external as sce
import bbknn
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import scanpy as sc
import seaborn as sns
from scipy.sparse import csr_matrix
import re
from pybiomart import Server

# Set plotting parameters
sc.set_figure_params(dpi=80, frameon=False, vector_friendly=True)

## Load blood and tumor tissue data
adata_blood = sc.read_h5ad("results/integrated_data_filtered.h5ad", chunk_size=30*6000)
csv_path = "result/metadata.csv"
metadata_df = pd.read_csv(csv_path, index_col=0)

# Define columns to be synced from ref_df to adata (including columns that may not exist initially)
cols_to_sync = ['majority_voting', 'subcluster_voting', 'cdr3s_aa', 'cdr3s_nt']

# 3. Core processing function (Handle missing column sync, clinical mapping, normalization, and colors)
def process_adata_integrated(adata, ref_df, label):
    print(f"\n>>> Processing {label} data ({adata.n_obs} cells)")
    
    # A. Force metadata synchronization: Resolve missing columns and categorical conflicts
    common_barcodes = adata.obs_names.intersection(ref_df.index)
    aligned_ref = ref_df.loc[common_barcodes, cols_to_sync]
    
    for col in cols_to_sync:
        # Core logic: If column does not exist, initialize with NaN
        if col not in adata.obs.columns:
            adata.obs[col] = np.nan
        
        # Convert column to object to bypass Categorical restrictions, sync data, then convert back to category
        adata.obs[col] = adata.obs[col].astype(object) 
        adata.obs.loc[common_barcodes, col] = aligned_ref[col].values
        adata.obs[col] = adata.obs[col].astype('category')
    
    # B. Clinical mapping
    adata.obs['patientid'] = adata.obs['patient'].astype(str).map(id_map).astype('category')
    adata.obs['PSA_FC'] = adata.obs['patient'].astype(str).map(psa_map).astype(float)
    
    if label == "Tumor":
        adata.obs['group'] = adata.obs['group'].astype('category')
        tumor_order = [c for c in ['Before treatment', 'After treatment'] if c in adata.obs['group'].cat.categories]
        adata.obs['group'] = adata.obs['group'].cat.set_categories(tumor_order, ordered=True)
    
    # D. Preprocessing
    if 'log1p' not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        print("    Normalization and Log transformation complete.")
    
    # E. Cleanup and apply color palettes
    for k in list(adata.uns.keys()):
        if k.endswith('_colors'):
            del adata.uns[k]
            
    for col, palette in palettes.items():
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype('category')
            cats = adata.obs[col].cat.categories
            adata.uns[f"{col}_colors"] = np.array([palette.get(c, "#d3d3d3") for c in cats])
    
    gc.collect()
    print(f">>> {label} processing complete.")

# 4. Execute workflow
# process_adata_integrated(adata_tumor, metadata_df, "Tumor")
process_adata_integrated(adata_blood, metadata_df, "Blood")

# 5. Final verification and saving
output_dir = "sp_data"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# print(f"\nFinal Check - Blood: {adata_blood.n_obs} cells, Tumor: {adata_tumor.n_obs} cells")
gc.collect()

# Environment configuration
sc.set_figure_params(dpi=80, frameon=False, vector_friendly=True)

# --- Helper function: BioMart Filtering ---
def filter_protein_coding_genes(adata, host='http://www.ensembl.org'):
    """Retain only protein-coding genes"""
    print(f"Connecting to BioMart to filter protein-coding genes...")
    try:
        server = Server(host=host)
        dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
        annot = dataset.query(attributes=['external_gene_name', 'gene_biotype'])
        annot.columns = ['gene_id', 'gene_biotype']
        gene_map = dict(zip(annot['gene_id'], annot['gene_biotype']))
        
        adata.var['biotype'] = adata.var_names.map(gene_map).fillna('unknown')
        adata = adata[:, adata.var['biotype'] == 'protein_coding'].copy()
        print(f"Filtering complete, remaining genes: {adata.n_vars}")
    except Exception as e:
        print(f"⚠️ BioMart connection failed, skipping this step: {e}")
    return adata

# --- Core function: Subgroup Re-clustering Analysis ---
def analyze_subgroup(adata_full, cell_type, label):
    print(f"\n{'>'*10} Processing: {label} - {cell_type} {'<'*10}")
    
    # 1. Extract subset
    adata_sub = adata_full[adata_full.obs['majority_voting'] == cell_type].copy()
    if adata_sub.n_obs < 50:
        print(f"Skipping {cell_type}: Cell count too low ({adata_sub.n_obs})")
        return None

    # 2. Restore Raw Counts (if available) and filter low-expression genes
    if adata_sub.raw is not None:
        adata_sub = adata_sub.raw.to_adata()
    sc.pp.filter_genes(adata_sub, min_cells=3)

    # 4. Preprocessing
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)

    # 3. Protein-coding gene filtering & interference gene removal (MT, RP, DNAJ)
    adata_sub = filter_protein_coding_genes(adata_sub)
    patterns = [r'^(MT-|MTRNR2L)', r'^(RPL|RPS)', r'^(DNAJ)']
    combined_pat = "|".join(patterns)
    remove_mask = adata_sub.var_names.str.contains(combined_pat, case=False, regex=True)
    adata_sub = adata_sub[:, ~remove_mask].copy()
    
    gene_counts = (adata_sub.X > 0).sum(axis=1)
    adata_sub.obs['post_filter_n_genes'] = np.array(gene_counts).flatten()
    adata_sub = adata_sub[adata_sub.obs['post_filter_n_genes'] >= 200, :].copy()
    print(f"Final cells remaining for HVG analysis: {adata_sub.n_obs}")
    
    # Threshold: Minimum cells required to keep a batch
    min_cells_threshold = 5 
    
    # Calculate cell counts per batch
    batch_counts = adata_sub.obs['batch'].value_counts()
    
    # Identify valid batches (count >= threshold)
    valid_batches = batch_counts[batch_counts >= min_cells_threshold].index
    
    # Identify batches to be removed for warning notification
    removed_batches = batch_counts[batch_counts < min_cells_threshold].index
    
    if len(removed_batches) > 0:
        print(f"⚠️ Warning: The following batches have fewer than {min_cells_threshold} cells and will be removed to avoid BBKNN errors:")
        for b in removed_batches:
            print(f"  - {b}: {batch_counts[b]} cells")
            
        # Execute filtering, keeping only cells from valid_batches
        adata_sub = adata_sub[adata_sub.obs['batch'].isin(valid_batches)].copy()
    
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=1500)
    
    # 5. Regress out and Scale
    regress_cols = ['total_counts']
    if 'pct_counts_mt' in adata_sub.obs.columns:
        regress_cols.append('pct_counts_mt')
        
    n_cells = adata_sub.shape[0]
    if n_cells > 100000:
        print(f"⚠️ Cell count ({n_cells}) > 100k, skipping time-consuming Regress Out step, performing Scale only...")
        sc.pp.scale(adata_sub, max_value=10)
        sc.tl.pca(adata_sub, n_comps=30)
    else:
        print("Performing regression analysis (Regress Out)...")
        # Standard regression for denoising when cell count is smaller
        sc.pp.regress_out(adata_sub, regress_cols, n_jobs=-1)
        sc.pp.scale(adata_sub, max_value=10)
        sc.tl.pca(adata_sub, n_comps=30)
        
        print("Executing Harmony integration...")
        sce.pp.harmony_integrate(adata_sub, key='batch', basis='X_pca', adjusted_basis='X_pca_harmony')
        rep = 'X_pca_harmony'

    # 8. Neighbor calculation and UMAP (BBKNN integration)
    print("Executing BBKNN integration...")
    bbknn.bbknn(adata_sub, batch_key='batch', use_rep=rep if 'rep' in locals() else 'X_pca')
    sc.tl.umap(adata_sub)
    
    # 9. Save results
    file_name = f"adata_{label}_{cell_type.replace('/', '_')}.h5ad"
    save_path = os.path.join(output_dir, file_name)
    adata_sub.write(save_path)
    print(f"✅ Results saved to: {save_path}")
    
    del adata_sub
    gc.collect()
    return None

# Loop through subgroups for analysis
for sg in subgroups:
    analyze_subgroup(adata_blood, sg, "Blood")