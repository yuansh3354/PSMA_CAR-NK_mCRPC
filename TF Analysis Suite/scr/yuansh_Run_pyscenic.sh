#!/bin/bash

# Author: YuanSh
# Date: 2024-01-27
# Description: pySCENIC single-cell analysis pipeline using Docker
# Reference: https://lishensuo.github.io/posts/bioinfo/%E5%8D%95%E7%BB%86%E8%83%9E%E6%95%B0%E6%8D%AE%E5%88%86%E6%9E%90--%E8%BD%AC%E5%BD%95%E5%9B%A0%E5%AD%90pyscenic/

# Pull the pySCENIC module from Docker Hub
docker pull aertslab/pyscenic_scanpy:0.12.1_1.9.1

# Download required configuration and reference files
wget -c https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
wget -c https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt
wget -c https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# Start the environment and mount the data directory
# docker run --platform linux/amd64 -it aertslab/pyscenic_scanpy:0.12.1_1.9.1 /bin/bash
docker run --platform linux/amd64 -v /home/user/pySCENIC:/data -it aertslab/pyscenic_scanpy:0.12.1_1.9.1 /bin/bash

# ==============================================================================
# Run pySCENIC analysis inside the container
# ==============================================================================

# Step 1: GRN inference
pyscenic grn \
    --sparse \
    --method grnboost2 \
    --output /dataoutput_adj.tsv \
    /data/pySCENIC_count.loom \
    /data/hs_hgnc_tfs.txt

# Step 2: Regulon analysis (CTX)
pyscenic ctx \
    --output /data/output_regulons.csv \
    --expression_mtx_fname /data/pySCENIC_count.loom \
    --all_modules \
    --mode "dask_multiprocessing" \
    --annotations_fname /data/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --mask_dropouts \
    /data/output_adj.tsv \
    /data/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather

# Step 3: AUCell scoring
pyscenic aucell \
    -o /data/output_aucell.csv \
    /data/pySCENIC_count.loom \
    /data/output_regulons.csv
