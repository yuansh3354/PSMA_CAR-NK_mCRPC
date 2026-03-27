### ---------------
### Create: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
### Date: 2024-07-16
### ---------------

import os
import argparse
import sys
import pandas as pd
import scanpy as sc
import bbknn
import warnings
import re
import time
import json
from datetime import datetime
warnings.filterwarnings('ignore')

# Set Scanpy parameters
sc.settings.verbosity = 3
sc.logging.print_header()

# Create argument parser
parser = argparse.ArgumentParser(description='Single-cell multi-sample integration and pathway enrichment pipeline')

# Required arguments
parser.add_argument('--sample_file', required=True, help='Sample information file, first column is sample name, second is sample path')
parser.add_argument('--output', required=True, help='Output directory for results')

# Optional arguments
parser.add_argument('--format', default='10X', choices=['h5ad', 'csv', '10X', 'matrix'], help='Input data format')
parser.add_argument('--resolution', type=float, default=1, help='Clustering resolution')
parser.add_argument('--n_top_genes', type=int, default=0, help='Number of highly variable genes')
parser.add_argument('--batch_method', default='harmony', choices=['harmony', 'bbknn'], help='Batch effect correction method')
parser.add_argument('--min_genes', type=int, default=300, help='Minimum number of genes per cell')
parser.add_argument('--max_genes', type=int, default=6000, help='Maximum number of genes per cell')
parser.add_argument('--min_counts', type=int, default=500, help='Minimum counts per cell')
parser.add_argument('--max_counts', type=int, default=50000, help='Maximum counts per cell')
parser.add_argument('--min_cells', type=int, default=5, help='Minimum number of cells per gene')
parser.add_argument('--mito_thresh', type=float, default=10.0, help='Mitochondrial gene threshold percentage')
parser.add_argument('--n_pcs', type=int, default=50, help='Number of principal components for PCA')
parser.add_argument('--n_neighbors_pcs', type=int, default=20, help='Number of principal components used in neighborhood graph construction')
parser.add_argument('--n_neighbors', type=int, default=15, help='Number of neighbors')
parser.add_argument('--random_state', type=int, default=42, help='Random seed')
parser.add_argument('--bbknn_n_pcs', type=int, default=50, help='Number of principal components used for BBKNN')
parser.add_argument('--neighbors_within_batch', type=int, default=5, help='Number of neighbors within batch for BBKNN')
parser.add_argument('--batch_col', default='PatientID', help='Column name for batch information')

args = parser.parse_args()

# Record start time
start_time = time.time()
start_datetime = datetime.now()
print(f"Script started at: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")

# Create output subdirectory
subfolder_name = "scanpy_multi"
output_dir = os.path.join(args.output, subfolder_name)
os.makedirs(output_dir, exist_ok=True)

# Update output path
args.output = output_dir
sc.settings.figdir = output_dir

# Save running parameters
params_dict = {
    'script_name': 'Single-cell multi-sample integration and pathway enrichment pipeline',
    'start_time': start_datetime.strftime('%Y-%m-%d %H:%M:%S'),
    'parameters': vars(args),
    'scanpy_version': sc.__version__,
    'python_version': os.sys.version
}

# Save parameters to JSON file
with open(os.path.join(args.output, 'run_parameters.json'), 'w', encoding='utf-8') as f:
    json.dump(params_dict, f, indent=2)

# Save parameters to text file
with open(os.path.join(args.output, 'run_parameters.txt'), 'w', encoding='utf-8') as f:
    f.write(f"Script Name: {params_dict['script_name']}\n")
    f.write(f"Start Time: {params_dict['start_time']}\n")
    f.write(f"Scanpy Version: {params_dict['scanpy_version']}\n")
    f.write(f"Python Version: {params_dict['python_version']}\n")
    f.write("\nRunning Parameters:\n")
    f.write("-" * 50 + "\n")
    for key, value in vars(args).items():
        f.write(f"{key}: {value}\n")

print(f"Running parameters saved to: {os.path.join(args.output, 'run_parameters.json')}")
print(f"Running parameters saved to: {os.path.join(args.output, 'run_parameters.txt')}")

# 1. Read sample information
print("\n1. Reading sample information...")
sample_df = pd.read_csv(args.sample_file)
sample_list = sample_df['SampleID'].tolist()
path_list = sample_df['PATH'].tolist()
group_list = sample_df['Days'].tolist()
patient_list = sample_df['PatientID'].tolist()
Dosage_list = sample_df['Dosage'].tolist()

print(f"Number of samples: {len(sample_list)}")

# 2. Load sample data
print("\n2. Loading sample data...")
adatas = []
for sample, path, group, patient, Dosage in zip(sample_list, path_list, group_list, patient_list, Dosage_list):
    print(f"  Loading: {sample}")
    if args.format == '10X':
        adata = sc.read_10x_h5(path + '/sample_filtered_feature_bc_matrix.h5')
        sc.pp.filter_genes(adata, min_cells=args.min_cells)
        adata.var_names_make_unique()
    elif args.format == 'h5ad':
        adata = sc.read_h5ad(path)
        sc.pp.filter_genes(adata, min_cells=args.min_cells)
        adata.var_names_make_unique()
    elif args.format == 'csv':
        adata = sc.read_csv(path)
        sc.pp.filter_genes(adata, min_cells=args.min_cells)
        adata.var_names_make_unique()
    elif args.format == 'matrix':
        adata = sc.read_10x_mtx(path, cache=False)
        sc.pp.filter_genes(adata, min_cells=args.min_cells)
        adata.var_names_make_unique()
            
    adata.obs.index = [sample + '_' + i for i in adata.obs.index]
    adata.obs['sample'] = sample
    adata.obs['group'] = group
    adata.obs['patient'] = patient
    adata.obs['Dosage'] = Dosage
    adata.obs['batch'] = adata.obs[args.batch_col].copy()
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adatas.append(adata)

# 3. Concatenate data
print("\n3. Concatenating sample data...")
adata = sc.concat(adatas, join='outer', label='batch', keys=sample_list, index_unique=None)
adata.var_names_make_unique()
keep_genes = [gene for gene in adata.var.index if not re.match(r'^ENSG[0-9]+$', gene)]
adata = adata[:, keep_genes].copy()
keep_genes = [gene for gene in adata.var.index if not re.match(r'^A[A-Z][1-9]+\.[1-9]$', gene)]
adata = adata[:, keep_genes].copy()

# 4. Quality control
print("\n4. Performing quality control...")
print(f"  Before filtering - Cells: {adata.n_obs}, Genes: {adata.n_vars}")
sc.pp.filter_cells(adata, min_genes=args.min_genes)
sc.pp.filter_cells(adata, max_genes=args.max_genes)
sc.pp.filter_cells(adata, min_counts=args.min_counts)
sc.pp.filter_cells(adata, max_counts=args.max_counts)
adata = adata[adata.obs['pct_counts_mt'] < args.mito_thresh, :]

# --- Filter out small batches to prevent Scrublet errors ---
print("  Checking and filtering small batches (required for Scrublet)...")
batch_counts = adata.obs['batch'].value_counts()
min_batch_cells = 75  # Scrublet requires at least 30+, recommended 50+

# Identify non-compliant batches
small_batches = batch_counts[batch_counts < min_batch_cells]
    
if len(small_batches) > 0:
    print(f"  Warning: The following batches have less than {min_batch_cells} cells and will be removed: \n{small_batches}")
    valid_batches = batch_counts[batch_counts >= min_batch_cells].index
    adata = adata[adata.obs['batch'].isin(valid_batches)].copy()
    print(f"  Cells remaining after removal: {adata.n_obs}")
# ----------------------------------------------------

# Remove doublets
print("  Running doublet removal pipeline...")
sc.pp.scrublet(adata, batch_key="batch") 
adata = adata[~adata.obs['predicted_doublet']].copy()

adata.raw = adata
adata.X = adata.raw.X.copy()
print(f"  After filtering - Cells: {adata.n_obs}, Genes: {adata.n_vars}")

# 5. Data preprocessing
print("\n5. Data preprocessing...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

if args.n_top_genes == 0:
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='batch')
else:
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, batch_key='batch')
sc.pp.scale(adata, max_value=10)

# 6. Dimensionality reduction and batch effect correction
# Output number of batches
print(f"  Current number of batches: {adata.obs['batch'].nunique()}")
print("\n6. Dimensionality reduction and batch correction...")
sc.tl.pca(adata, svd_solver='arpack', n_comps=args.n_pcs, random_state=args.random_state)

if args.batch_method == 'harmony':
    print("  Using Harmony for batch correction...")
    try:
        # Always use the 'batch' column for correction
        sc.external.pp.harmony_integrate(adata, 'batch', random_state=args.random_state, adjusted_basis='X_pca_harmony')
        use_rep = 'X_pca_harmony'
        sc.pp.neighbors(adata, n_pcs=args.n_neighbors_pcs, n_neighbors=args.n_neighbors, 
                        random_state=args.random_state, use_rep=use_rep)
    except:
        print("  Harmony is not available, falling back to original PCA.")
        use_rep = 'X_pca'
        sc.pp.neighbors(adata, n_pcs=args.n_neighbors_pcs, n_neighbors=args.n_neighbors, 
                        random_state=args.random_state, use_rep=use_rep)
elif args.batch_method == 'bbknn':
    print("  Using BBKNN for batch correction...")
    try:
        bbknn.bbknn(adata, batch_key='batch', n_pcs=args.bbknn_n_pcs, neighbors_within_batch=args.neighbors_within_batch)
        use_rep = 'X_pca'
    except:
        print("  BBKNN is not available, falling back to original PCA.")
        use_rep = 'X_pca'
        sc.pp.neighbors(adata, n_pcs=args.n_neighbors_pcs, n_neighbors=args.n_neighbors, 
                        random_state=args.random_state, use_rep=use_rep)

# 7. Neighborhood graph and UMAP
print("\n7. Constructing neighborhood graph and UMAP...")
sc.tl.umap(adata, random_state=args.random_state)

# 8. Clustering analysis
print("\n8. Clustering analysis...")
leiden_res = args.resolution
leiden_key = f'leiden_{str(leiden_res).replace(".", "_")}'
sc.tl.leiden(adata, resolution=leiden_res, key_added=leiden_key, random_state=args.random_state, n_iterations=2)

# 9. Visualization
print("\n9. Generating visualization results...")
sc.pl.umap(adata, color=[leiden_key], save='_leiden.png', show=False)
sc.pl.umap(adata, color=['sample'], save='_sample.png', show=False)

# 10. Save results
print("\n10. Saving results...")
adata = adata.raw.to_adata()
del adata.raw
adata.write(os.path.join(args.output, 'integrated_data.h5ad'))
adata.obs[['sample', leiden_key]].to_csv(os.path.join(args.output, 'clusters.csv'))

# Record end time and calculate total elapsed time
end_time = time.time()
end_datetime = datetime.now()
total_time = end_time - start_time

# Update parameter file with end time and total runtime
params_dict['end_time'] = end_datetime.strftime('%Y-%m-%d %H:%M:%S')
params_dict['total_runtime_seconds'] = round(total_time, 2)
params_dict['total_runtime_minutes'] = round(total_time/60, 2)

# Update JSON file
with open(os.path.join(args.output, 'run_parameters.json'), 'w', encoding='utf-8') as f:
    json.dump(params_dict, f, indent=2)

# Update text file
with open(os.path.join(args.output, 'run_parameters.txt'), 'a', encoding='utf-8') as f:
    f.write(f"\nEnd Time: {params_dict['end_time']}\n")
    f.write(f"Total Runtime: {params_dict['total_runtime_seconds']} seconds ({params_dict['total_runtime_minutes']} minutes)\n")

print(f"\nAnalysis completed! Results saved in: {args.output}")
print(f"Run logs updated in: {os.path.join(args.output, 'run_parameters.json')}")
print(f"Run logs updated in: {os.path.join(args.output, 'run_parameters.txt')}")

print("\n==================================================================") 
print(f"\nScript ended at: {end_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total runtime: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
