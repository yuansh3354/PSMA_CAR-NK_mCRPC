#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: YuanSh
# Date: 2024-01-27
# Description: Convert all CSV files in the counts_data directory to loom format.
# Suitable for gene expression data (rows as genes, columns as cells).
# The converted loom files are saved in the same directory.

import os
import glob
import loompy as lp
import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path

def csv_to_loom(csv_path, output_dir):
    """
    Convert a CSV file to loom format.
    
    Args:
        csv_path: Path to the CSV file.
        output_dir: Path to the output directory.
        
    Returns:
        Path to the created loom file, or None if failed.
    """
    # Get the file name without extension
    file_name = os.path.splitext(os.path.basename(csv_path))[0]
    print(f"Processing file: {file_name}")
    
    # Read the CSV file
    try:
        print("  Reading CSV file...")
        adata = sc.read_csv(csv_path)
        print(f"  Successfully read CSV, shape: {adata.shape}")
    except Exception as e:
        print(f"  Error reading CSV file: {str(e)}")
        return None
    
    # Check data validity
    if adata.X.shape[0] == 0 or adata.X.shape[1] == 0:
        print("  Data is empty, skipping.")
        return None
    
    # Create the loom file path in the output directory
    loom_path = os.path.join(output_dir, f"{file_name}.loom")
    
    try:
        # Create row and column attributes
        row_attrs = {"Gene": np.array(adata.var_names)}
        col_attrs = {"CellID": np.array(adata.obs_names)}
        
        # Create the loom file (transpose X to make rows=genes, cols=cells)
        print("  Creating loom file...")
        lp.create(loom_path, adata.X.transpose(), row_attrs, col_attrs)
        print(f"  Successfully created loom file: {loom_path}")
        
        return loom_path
    except Exception as e:
        print(f"  Error creating loom file: {str(e)}")
        return None

def main():
    """Main program entry"""
    print("Starting conversion of CSV files to loom format...")
    
    # Define the directory containing the data
    counts_data_dir = "counts_data"
    
    # Check if the target directory exists
    if not os.path.exists(counts_data_dir):
        print(f"Error: Directory '{counts_data_dir}' does not exist.")
        return
    
    # Ensure the target path is a directory
    if not os.path.isdir(counts_data_dir):
        print(f"Error: '{counts_data_dir}' is not a directory.")
        return
    
    # Get all CSV files in the directory
    csv_pattern = os.path.join(counts_data_dir, "*.csv")
    csv_files = glob.glob(csv_pattern)
    
    # Exit if no CSV files are found
    if not csv_files:
        print(f"No CSV files found in directory '{counts_data_dir}'.")
        return
    
    print(f"Found {len(csv_files)} CSV files in directory '{counts_data_dir}'.")
    
    # Iterate and convert each CSV file
    success_count = 0
    for csv_file in csv_files:
        try:
            result = csv_to_loom(csv_file, counts_data_dir)
            if result:
                success_count += 1
                print(f"Successfully converted: {os.path.basename(csv_file)} -> {os.path.basename(result)}")
            else:
                print(f"Conversion failed for: {os.path.basename(csv_file)}")
        except Exception as e:
            print(f"Error processing file {os.path.basename(csv_file)}: {str(e)}")
    
    # Print conversion statistics
    print(f"Conversion complete. Successfully converted {success_count}/{len(csv_files)} files.")
    
    # Calculate and display the size of the generated loom files
    loom_pattern = os.path.join(counts_data_dir, "*.loom")
    loom_files = glob.glob(loom_pattern)
    if loom_files:
        print(f"\nLoom files in directory '{counts_data_dir}' ({len(loom_files)}):")
        for loom_file in loom_files:
            # Convert file size to MB
            loom_size = Path(loom_file).stat().st_size / (1024 * 1024)
            print(f"  {os.path.basename(loom_file)} ({loom_size:.2f} MB)")
    else:
        print(f"\nNo loom files found in directory '{counts_data_dir}'.")

if __name__ == "__main__":
    main()
