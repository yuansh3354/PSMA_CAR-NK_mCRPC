python 01.scanpy_multi.py \
  --sample_file all_samples.csv \
  --output ./results \
  --format matrix \
  --resolution 1 \
  --n_top_genes 0 \
  --batch_method bbknn \
  --min_genes 500 \
  --max_genes 10000 \
  --min_counts 750 \
  --max_counts 80000 \
  --min_cells 5 \
  --mito_thresh 25.0 \
  --n_pcs 150 \
  --n_neighbors_pcs 50 \
  --n_neighbors 15 \
  --random_state 42 \
  --bbknn_n_pcs 50 \
  --neighbors_within_batch 1 \
  --batch_col sample

  
