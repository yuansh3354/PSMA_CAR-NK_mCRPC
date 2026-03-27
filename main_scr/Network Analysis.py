import gc
import os
import re
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from matplotlib.lines import Line2D

# ==========================================
# 1. Environment Initialization & Configuration
# ==========================================
# Set scanpy plotting parameters
sc.set_figure_params(dpi=80, frameon=False, vector_friendly=True)

# ==========================================
# 2. Data Loading & Preprocessing
# ==========================================
# Load AnnData objects
adata_tumor = sc.read_h5ad('results/01.scanpy_multi/integrated_data_l2_match_blood.h5ad')
adata_blood = sc.read_h5ad("results/01.scanpy_multi/integrated_data_filtered.h5ad")

def process_adata(adata, ref_df, label):
    """
    Synchronize metadata, apply clinical mapping, and perform normalization.
    """
    print(f"--- Processing {label} Data ({adata.n_obs} cells) ---")
    
    # Sync specific metadata columns from reference dataframe
    cols_to_sync = ['majority_voting', 'subcluster_voting', 'cdr3s_aa', 'cdr3s_nt']
    common_barcodes = adata.obs_names.intersection(ref_df.index)
    
    for col in cols_to_sync:
        adata.obs[col] = adata.obs[col].astype(object) # Avoid categorical assignment issues
        adata.obs.loc[common_barcodes, col] = ref_df.loc[common_barcodes, col].values
        adata.obs[col] = adata.obs[col].astype('category')
    
    # Apply clinical info mapping
    adata.obs['patientid'] = adata.obs['patient'].astype(str).map(id_map).astype('category')
    adata.obs['PSA_FC'] = adata.obs['patient'].astype(str).map(psa_map).astype(float)
    
    # Standard normalization and log-transformation
    if 'log1p' not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
    # Apply global color palettes to uns
    for col, palette in palettes.items():
        if col in adata.obs.columns:
            cats = adata.obs[col].cat.categories
            adata.uns[f"{col}_colors"] = np.array([palette.get(c, "#d3d3d3") for c in cats])
    gc.collect()

# Execute processing
process_adata(adata_tumor, metadata_df, "Tumor")
process_adata(adata_blood, metadata_df, "Blood")

# ==========================================
# 3. Lineage Proportion Dynamics
# ==========================================

def get_proportions(obs, group_cols, target_col):
    """
    Calculate subcluster proportions within a specified lineage group.
    """
    counts = obs.groupby(group_cols + [target_col], observed=True).size().reset_index(name='n')
    totals = counts.groupby(group_cols, observed=True)['n'].transform('sum')
    counts['proportion'] = (counts['n'] / totals) * 100
    return counts

pbmc_props = get_proportions(obs_pbmc, ['patientid', 'group', 'majority_voting'], 'subcluster_voting')
major_lineages = ['CD4T', 'CD8T', 'Myeloid', 'NK', 'B']

# Visualize Median Trends (1x5 Grid)
fig, axes = plt.subplots(1, 5, figsize=(18, 4))
for i, lin in enumerate(major_lineages):
    data = pbmc_props[pbmc_props['majority_voting'] == lin]
    if data.empty: continue
    
    # Plot individual patient trajectories (thin gray lines)
    sns.lineplot(data=data, x='group', y='proportion', units='patientid', estimator=None, 
                 color='gray', alpha=0.2, ax=axes[i])
    
    # Plot group median trend (thick colored line)
    sns.lineplot(data=data, x='group', y='proportion', estimator=np.median, 
                 linewidth=4, marker='o', ax=axes[i], color=palettes['majority_voting'].get(lin, 'black'))
    
    axes[i].set_title(lin, fontweight='bold')
    axes[i].set_ylabel("Proportion (%)" if i==0 else "")
    
    # X-axis cleanup: rename Day0 to 2Hours for better clinical context
    current_labels = [l.get_text() for l in axes[i].get_xticklabels()]
    axes[i].set_xticklabels(["2Hours" if l == "Day0" else l for l in current_labels])

# ==========================================
# 4. Lagged Correlation & Network Analysis
# ==========================================

def get_pivot_matrix(obs, target_sort, target_lineages):
    """
    Generate a pivot matrix of proportions (Subcluster x SampleID).
    """
    df = obs[(obs['sort'] == target_sort) & (obs['majority_voting'].isin(target_lineages))].copy()
    df['sample_id'] = df['patientid'].astype(str) + "_" + df['group'].astype(str)
    
    lineage_totals = df.groupby(['patientid', 'group', 'majority_voting'], observed=True).size().reset_index(name='tot')
    counts = df.groupby(['patientid', 'group', 'majority_voting', 'subcluster_voting'], observed=True).size().reset_index(name='n')
    merged = pd.merge(counts, lineage_totals, on=['patientid', 'group', 'majority_voting'])
    merged['prop'] = (merged['n'] / merged['tot']) * 100
    
    pivot_mat = merged.pivot_table(index='subcluster_voting', columns='sample_id', values='prop').fillna(0)
    # Remove subclusters with zero variance
    return pivot_mat[pivot_mat.std(axis=1) > 1e-10]

# Prepare matrices for cross-correlation
matrix_pbmc = get_pivot_matrix(obs_pbmc, 'PBMC', ['CD4T', 'CD8T', 'Myeloid'])
matrix_carnk = get_pivot_matrix(adata_blood.obs, 'Sorted CAR-NK cells', ['NK'])

def calculate_lagged_corr(m_a, m_b, lags=[0, 1, 2]):
    """
    Calculate time-lagged correlations between A and B, keeping the optimal lag for each pair.
    """
    p_seq = ['Day0', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21', 'Day60']
    c_seq = ['Day0', 'Day3', 'Day7', 'Day14']
    all_res = []
    
    for lag in lags:
        patients = list(set(c.split('_')[0] for c in m_a.columns))
        ref_cols, tgt_cols = [], []
        for pt in patients:
            for t in c_seq:
                try:
                    idx = p_seq.index(t)
                    if idx + lag < len(p_seq):
                        s_ref, s_tgt = f"{pt}_{t}", f"{pt}_{p_seq[idx+lag]}"
                        if s_ref in m_a.columns and s_tgt in m_b.columns:
                            ref_cols.append(s_ref); tgt_cols.append(s_tgt)
                except: continue
        if len(ref_cols) < 5: continue
        
        # Calculate Pearson correlation matrix
        c_mat = np.corrcoef(m_a[ref_cols].values, m_b[tgt_cols].values)[:len(m_a), len(m_a):]
        df = pd.DataFrame(c_mat, index=m_a.index, columns=m_b.index).stack().reset_index()
        df.columns = ['A', 'B', 'r']; df['lag'] = lag
        all_res.append(df)
        
    full_corr = pd.concat(all_res)
    # Return highest absolute correlation across lags
    return full_corr.loc[full_corr.groupby(['A', 'B'])['r'].transform(lambda x: x.abs()).idxmax()]

best_lag_corr = calculate_lagged_corr(matrix_carnk, matrix_pbmc)

# Visualize Network (Hubs identified by degree > 3)
def draw_interaction_net(df, threshold, title, edge_color, out_name):
    """
    Draw a circular topological network of subcluster interactions.
    """
    filtered = df[df['r'].abs() > threshold].copy()
    G = nx.Graph()
    for _, row in filtered.iterrows():
        G.add_edge(row['A'], row['B'], weight=row['r'], lag=row['lag'])
    
    # Plot only modules with more than 3 nodes
    components = [G.subgraph(c).copy() for c in nx.connected_components(G) if len(c) > 3]
    if not components: return
    
    fig, axes = plt.subplots(1, len(components), figsize=(len(components)*6, 6), squeeze=False)
    for i, sub_G in enumerate(components):
        pos = nx.circular_layout(sub_G)
        hubs = [n for n, d in sub_G.degree() if d > 3]
        
        # Draw edges; width reflects correlation strength
        ws = [abs(sub_G[u][v]['weight'])*5 for u, v in sub_G.edges()]
        nx.draw_networkx_edges(sub_G, pos, width=ws, edge_color=edge_color, alpha=0.4, ax=axes[0,i])
        
        # Draw nodes (Hubs as stars, others as circles)
        nx.draw_networkx_nodes(sub_G, pos, nodelist=[n for n in sub_G.nodes() if n not in hubs], 
                               node_size=600, node_color="#D3D3D3", ax=axes[0,i])
        nx.draw_networkx_nodes(sub_G, pos, nodelist=hubs, node_shape='*', node_size=1500, 
                               node_color="#FFD700", edgecolors="black", ax=axes[0,i])
        
        nx.draw_networkx_labels(sub_G, pos, font_size=8, ax=axes[0,i])
        axes[0,i].set_title(f"Module {i+1}", fontweight='bold')
        axes[0,i].axis('off')
    
    plt.savefig(f"figure_final/blood/{out_name}.pdf")
    plt.show()

# ==========================================
# 5. Receptor-Ligand (CellChat-based) Visualization
# ==========================================
# Filter and save subset objects to reduce storage footprint
def save_filtered_subset(adata, name):
    all_genes = [g for sub in pathway_genes.values() for g in sub]
    valid = [g for g in all_genes if g in adata.var_names]
    subset = adata[:, valid].copy()
    subset.write(f"sp_data/{name}_cellchat.h5ad")
    return subset

# Extract subsets for specific cell populations
adata_carnk_chat = save_filtered_subset(adata_blood[adata_blood.obs['sort']=='Sorted CAR-NK cells'], "CARNK")
adata_cd8t_chat = save_filtered_subset(adata_blood[adata_blood.obs['majority_voting']=='CD8T'], "CD8T")

# Visualize with Dotplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
all_viz_genes = [g for sub in pathway_genes.values() for g in sub]

sc.pl.dotplot(adata_carnk_chat, [g for g in all_viz_genes if g in adata_carnk_chat.var_names], 
              groupby='subcluster_voting', title='Ligands in CAR-NK', ax=ax1, show=False, standard_scale='var')
sc.pl.dotplot(adata_cd8t_chat, [g for g in all_viz_genes if g in adata_cd8t_chat.var_names], 
              groupby='subcluster_voting', title='Receptors in CD8T', ax=ax2, show=False, standard_scale='var')

plt.tight_layout()
plt.show()

print("--- Analysis Pipeline Execution Completed ---")