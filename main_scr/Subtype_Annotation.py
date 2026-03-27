#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# Author: Yuan.Sh, MD (ORCID: 0000-0002-6028-0185)
# Date: 2025-03-26
# ==============================================================================
import argparse
import os
import sys

CORES_TO_USE = "64" 

os.environ["OMP_NUM_THREADS"] = CORES_TO_USE
os.environ["OPENBLAS_NUM_THREADS"] = CORES_TO_USE
os.environ["MKL_NUM_THREADS"] = CORES_TO_USE
os.environ["VECLIB_MAXIMUM_THREADS"] = CORES_TO_USE
os.environ["NUMEXPR_NUM_THREADS"] = CORES_TO_USE

import time
import json
from datetime import datetime
import scanpy as sc
import pandas as pd
import numpy as np
import bbknn
import scipy.sparse as sp
import re
import rapids_singlecell as rsc

import cupy as cp

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    pool_allocator=False,
    managed_memory=True,       # Enable managed memory
    devices=0
)
cp.cuda.set_allocator(rmm_cupy_allocator)
from pybiomart import Server

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Single-cell RNA-seq subpopulation Leiden clustering script (supports Harmony/BBKNN/Both)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--work_dir', type=str, default='', help='Working directory path')
    parser.add_argument('--input', '-i', required=True, help='Input h5ad file path')
    parser.add_argument('--output_dir', type=str, default='./output', help='Output directory path')
    parser.add_argument('--batch_key', default='SampleID', help='Column name used for batch correction')
    
    parser.add_argument(
        '--leiden_res_list', 
        type=str, 
        default='{"Stromal":1.0,"Myeloid":1.0,"Endothelial":1.0,"CD8T":1.0,"B":1.0,"CD4T":1.0,"NK":1.0,"ILC":1.0,"Epithelial":1.0}',
        help='Leiden clustering resolution for different cell types (JSON format)'
    )
    
    parser.add_argument(
        '--n_pcs_list', 
        type=str, 
        default='{"Stromal":200,"Myeloid":200,"Endothelial":200,"CD8T":200,"B":200,"CD4T":200,"NK":200,"ILC":200,"Epithelial":200}',
        help='Number of PCA components for different cell types'
    )
    
    parser.add_argument(
        '--bbknn_n_pcs_list', 
        type=str, 
        default='{"Stromal":50,"Myeloid":50,"Endothelial":50,"CD8T":50,"B":50,"CD4T":50,"NK":50,"ILC":50,"Epithelial":50}',
        help='Number of PCs used when building neighbor graph (if Harmony Only, used as n_pcs parameter)'
    )
    
    parser.add_argument(
        '--n_top_genes_list', 
        type=str, 
        default='{"Stromal":2000,"Myeloid":2000,"Endothelial":2000,"CD8T":2000,"B":2000,"CD4T":2000,"NK":2000,"ILC":2000,"Epithelial":2000}',
        help='Number of highly variable genes'
    )

    # === [Modified] Changed use_harmony_list to integration_method_list ===
    # Optional values: "bbknn", "harmony", "both"
    parser.add_argument(
        '--integration_method_list', 
        type=str, 
        default='{"Stromal":"both","Myeloid":"both","Endothelial":"both","CD8T":"both","B":"both","CD4T":"both","NK":"both","ILC":"both","Epithelial":"both"}',
        help='Batch correction method: bbknn (BBKNN only), harmony (Harmony only), both (Harmony+BBKNN)'
    )
    
    parser.add_argument('--random_seed', type=int, default=42, help='Random seed')
    parser.add_argument('--min_cells_per_batch', type=int, default=1, help='Minimum number of cells per batch (BBKNN neighbors_within_batch)')
    
    return parser.parse_args()

def filter_protein_coding_genes(adata, 
                                host='http://www.ensembl.org', 
                                gene_col=None, 
                                query_attr='external_gene_name'):
    """
    Connect to BioMart and filter AnnData to retain only protein-coding genes.

    Args:
        adata: AnnData object
        host: BioMart server address (default http://www.ensembl.org, if it fails try http://useast.ensembl.org)
        gene_col: Column name containing gene names. If gene names are in the index, keep as None.
        query_attr: BioMart query attribute, 'external_gene_name' (gene name) or 'ensembl_gene_id' (ID).
    
    Returns:
        Filtered AnnData object (copy)
    """
    print(f"Connecting to Ensembl BioMart ({host})...")
    try:
        server = Server(host=host)
        dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    except Exception as e:
        print(f"❌ Failed to connect to BioMart: {e}")
        return adata

    # 1. Get annotations
    print("Querying gene biotypes...")
    # query_attr corresponds to the gene format in user data, gene_biotype is our target
    annot = dataset.query(attributes=[query_attr, 'gene_biotype'])
    annot.columns = ['gene_id', 'gene_biotype'] # Rename for unified processing
    
    # Remove duplicates
    annot = annot.drop_duplicates(subset='gene_id')
    gene_map = dict(zip(annot['gene_id'], annot['gene_biotype']))

    # 2. Map to adata
    # Determine the list of genes to map
    if gene_col:
        target_genes = adata.var[gene_col]
    else:
        target_genes = adata.var_names

    # Map biotypes, fill 'unknown' if not found
    adata.var['biotype'] = target_genes.map(gene_map).fillna('unknown')

    # 3. Statistics and filtering
    n_total = adata.n_vars
    n_coding = (adata.var['biotype'] == 'protein_coding').sum()
    
    print(f"Original number of genes: {n_total}")
    print(f"Number of protein-coding genes: {n_coding}")
    
    if n_coding == 0:
        print("⚠️ Warning: No protein-coding genes found. Please check if the gene name format matches BioMart.")
        return adata

    # 4. Execute filtering
    adata_filtered = adata[:, adata.var['biotype'] == 'protein_coding'].copy()
    print(f"✅ Filtering completed. Remaining genes: {adata_filtered.n_vars}")

    return adata_filtered

def preprocess_data(adata, batch_key, leiden_res, n_pcs, bbknn_n_pcs, n_top_genes, min_cells_per_batch, integration_method):
    """
    Data preprocessing pipeline: QC -> HVG -> Remove interference -> PCA -> [Batch Correction] -> UMAP -> Leiden
    integration_method: 'bbknn', 'harmony', 'both'
    """
    print("Starting data preprocessing...")
    
    # --- 1. Dimension check and correction ---
    var_len = len(adata.var_names)
    if adata.raw is None:
        adata.raw = adata
        adata.X = adata.raw.X.copy()
        print("Warning: adata.raw does not exist, using adata.X as raw counts")
    else:
        adata = adata.raw.to_adata()
        adata.raw = adata
        
    # --- 2. Batch filtering ---
    batch_counts = adata.obs[batch_key].value_counts()
    valid_batches = batch_counts[batch_counts >= min_cells_per_batch].index
    adata = adata[adata.obs[batch_key].isin(valid_batches)].copy()
    
    adata = filter_protein_coding_genes(adata)
    print(f"Retained {len(valid_batches)} batches, {adata.n_obs} cells after filtering")
    
    # --- 6. Remove interference genes ---
    print("Filtering and removing interference genes (MT-, RPL/S, DNAJ)...")
    patterns = [
        re.compile(r'^(MT-|MTRNR2L)', re.IGNORECASE), 
        re.compile(r'^(RPL|RPS)', re.IGNORECASE),      
        re.compile(r'^(DNAJ)', re.IGNORECASE)          
    ]
    remove_mask = np.zeros(adata.n_vars, dtype=bool)
    for pat in patterns:
        remove_mask |= adata.var_names.str.match(pat)
        
    num_removed = remove_mask.sum()
    if num_removed > 0:
        print(f"Detected {num_removed} genes to remove.")
        adata = adata[:, ~remove_mask].copy()
        print(f"Remaining genes after removal: {adata.n_vars}")
    else:
        print("No genes detected for removal.")    
        
    # --- 3. Basic filtering ---
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=3)
    rsc.get.anndata_to_GPU(adata)
    # --- 4. Normalization ---
    rsc.pp.normalize_total(adata, target_sum=1e4)
    rsc.pp.log1p(adata)
    
    # --- 5. Highly variable genes ---
    rsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key=batch_key, n_top_genes=n_top_genes)
      
    # --- 7. Scale & PCA ---
    rsc.pp.scale(adata, max_value=10)
    # sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=True, n_comps=n_pcs)
    rsc.tl.pca(adata, svd_solver="auto", use_highly_variable=True, n_comps=n_pcs)
    # --- 8. Integration and neighbor graph construction (logic branches) ---
    method = integration_method.lower()
    print(f"Executing integration strategy: {method}")

    if method == 'both':
        # 1. Run Harmony first
        print(f"Running Harmony (Key: {batch_key})...")
        rsc.pp.harmony_integrate(adata, 
                                         key=batch_key, 
                                         adjusted_basis='X_pca_harmony')
        
        # 2. Run BBKNN (based on Harmony results)
        print(f"Running BBKNN (Input: X_pca_harmony, n_pcs={bbknn_n_pcs})...")
        rsc.get.anndata_to_CPU(adata)
        bbknn.bbknn(
            adata, 
            batch_key=batch_key, 
            use_rep='X_pca_harmony', 
            n_pcs=bbknn_n_pcs,
            neighbors_within_batch=min_cells_per_batch 
        )

    elif method == 'bbknn':
        # Run BBKNN only (based on original PCA)
        rsc.get.anndata_to_CPU(adata)
        print(f"Running BBKNN Only (Input: X_pca, n_pcs={bbknn_n_pcs})...")
        bbknn.bbknn(
            adata, 
            batch_key=batch_key, 
            use_rep='X_pca', 
            n_pcs=bbknn_n_pcs,
            neighbors_within_batch=min_cells_per_batch 
        )

    elif method == 'harmony':
        # Run Harmony only
        print(f"Running Harmony Only (Key: {batch_key})...")
        rsc.pp.harmony_integrate(adata, key=batch_key, adjusted_basis='X_pca_harmony', max_iter_harmony=20)
        
        # Harmony does not build graph automatically, need to run neighbors manually
        print(f"Building Neighbors Graph (Input: X_pca_harmony)...")
        rsc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=30, n_pcs=bbknn_n_pcs)
    
    else:
        raise ValueError(f"Unknown integration method: {method}. Please use bbknn, harmony, or both.")
    
    rsc.get.anndata_to_GPU(adata)
    # --- 9. UMAP & Leiden ---
    rsc.tl.umap(adata)
    
    leiden_key = f'leiden_{str(leiden_res).replace(".", "_")}'
    rsc.tl.leiden(adata, resolution=leiden_res, key_added=leiden_key)
    rsc.get.anndata_to_CPU(adata)
    print(f"Clustering completed, resolution {leiden_res}, total {len(adata.obs[leiden_key].unique())} clusters")
    return adata, leiden_key

def save_run_parameters(args, anno_dir, start_datetime):
    """Save run parameters"""
    params_dict = {
        'script_name': 'Single-cell RNA-seq subpopulation Leiden clustering script',
        'start_time': start_datetime.strftime('%Y-%m-%d %H:%M:%S'),
        'input_parameters': vars(args)
    }
    json_file = os.path.join(anno_dir, 'run_parameters.json')
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(params_dict, f, indent=2)
    return params_dict, json_file

def create_output_directory(output_dir):
    anno_dir = os.path.join(output_dir, '03.Manual_SUBCLUSTER_ANNO')
    os.makedirs(anno_dir, exist_ok=True)
    return anno_dir

def save_individual_results(subset, anno_dir, cell_type, leiden_key):
    output_file = os.path.join(anno_dir, f'{cell_type}_processed_leiden.h5ad')
    subset.write(output_file)
    print(f"{cell_type} processing results saved to: {output_file}")
    return output_file

def main():
    start_time = time.time()
    start_datetime = datetime.now()
    print(f"Script started at: {start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
    
    sc.settings.verbosity = 3
    args = parse_arguments()
    np.random.seed(args.random_seed)
    
    if args.work_dir and os.path.exists(args.work_dir):
        os.chdir(args.work_dir)
        
    try:
        if not os.path.exists(args.input):
            raise FileNotFoundError(f"Input file does not exist: {args.input}")
        
        anno_dir = create_output_directory(args.output_dir)
        save_run_parameters(args, anno_dir, start_datetime)
        
        print("Loading input data...")
        adata = sc.read_h5ad(args.input, chunk_size=200000)
        
        if 'majority_voting' not in adata.obs.columns:
            if 'majorCluster' in adata.obs.columns:
                 adata.obs['majority_voting'] = adata.obs['majorCluster']
            else:
                 print("Warning: Input data is missing 'majority_voting' column")
        
        # 1. Parse parameters
        leiden_res_dict = json.loads(args.leiden_res_list)
        n_pcs_dict = json.loads(args.n_pcs_list)
        bbknn_n_pcs_dict = json.loads(args.bbknn_n_pcs_list) 
        n_top_genes_dict = json.loads(args.n_top_genes_list)
        # [Modified] Parse integration_method parameters
        integration_method_dict = json.loads(args.integration_method_list)

        # 2. Define target types
        target_cell_types = ['CD8','CD4','NK','ILC','CD8T','Tumor','Neutrophil', 'B', 'CD4T','Stromal','Fibro', 'Myeloid', 'Endothelial',  'Epithelial','CD14+ Monocyte','CD16+ Monocyte']
        
        # 3. Filter types existing in data
        existing_types = adata.obs['majority_voting'].unique()
        candidates = [t for t in target_cell_types if t in existing_types]
        
        # 4. Filter types with complete parameters
        final_target_cell_types = []
        for cell_type in candidates:
            missing_params = []
            if cell_type not in leiden_res_dict: missing_params.append('leiden_res')
            if cell_type not in n_pcs_dict: missing_params.append('n_pcs')
            if cell_type not in bbknn_n_pcs_dict: missing_params.append('bbknn_n_pcs')
            if cell_type not in n_top_genes_dict: missing_params.append('n_top_genes')
            # Check integration_method parameter
            if cell_type not in integration_method_dict: missing_params.append('integration_method')
            
            if len(missing_params) == 0:
                final_target_cell_types.append(cell_type)
            else:
                print(f"Warning: Skipping {cell_type} due to missing configurations in parameter dictionary: {missing_params}")

        target_cell_types = final_target_cell_types
        print(f"Final processing list: {target_cell_types}")

        individual_files = []
        for cell_type in target_cell_types:
            subset = adata[adata.obs['majority_voting'] == cell_type].copy()
            if subset.n_obs == 0: continue
            
            # Get integration method
            current_method = integration_method_dict[cell_type]

            print(f"\n>>> Processing: {cell_type} ({subset.n_obs} cells)")
            print(f"    Integration Method: {current_method}")
            
            subset, leiden_key = preprocess_data(
                subset, 
                batch_key=args.batch_key, 
                leiden_res=leiden_res_dict[cell_type], 
                n_pcs=n_pcs_dict[cell_type],
                bbknn_n_pcs=bbknn_n_pcs_dict[cell_type],
                n_top_genes=n_top_genes_dict[cell_type],
                min_cells_per_batch=args.min_cells_per_batch,
                integration_method=current_method # [Modified] Pass method name
            )
            
            output_file = save_individual_results(subset, anno_dir, cell_type, leiden_key)
            individual_files.append(output_file)
        
        # Finish processing
        total_time = time.time() - start_time
        print(f"\n=== All completed, {len(individual_files)} files generated ===")
        print(f"Total time: {total_time:.2f} seconds")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
