#!/bin/bash

# Author: YuanSh
# Date: 2024-01-27
# Description: Automated SCENIC pipeline for processing single-cell data

# Set parameters
TF_FILE="hs_hgnc_tfs.txt"
MOTIF_ANNOTATIONS="motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
RANKINGS_FILE="hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
NUM_WORKERS=32
MODE="dask_multiprocessing"
MAX_PARALLEL=10  # Maximum number of parallel samples
PATTERN="*.loom"  # Loom file pattern
INPUT_DIR="counts_data"  # Input directory
OUTPUT_DIR="scenic_result"  # Output directory

# Subdirectories for the three pipeline steps
GRN_DIR="$OUTPUT_DIR/GRN"
CTX_DIR="$OUTPUT_DIR/CTX"
AUCELL_DIR="$OUTPUT_DIR/AUCell"

# Log file
LOG_FILE="scenic_pipeline.log"

# Function: Log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function: Create all necessary directories
create_directories() {
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$GRN_DIR"
    mkdir -p "$CTX_DIR"
    mkdir -p "$AUCELL_DIR"
    
    log "Creating output directory structure:"
    log "  Main directory: $OUTPUT_DIR"
    log "  GRN results: $GRN_DIR"
    log "  CTX results: $CTX_DIR"
    log "  AUCell results: $AUCELL_DIR"
}

# Function: Process a single sample
process_sample() {
    SAMPLE_ID=$1
    LOOM_FILE="$INPUT_DIR/${SAMPLE_ID}.loom"
    
    log "Starting to process sample: $SAMPLE_ID"
    
    # Step 1: GRN inference
    log "[$SAMPLE_ID] Starting GRN inference..."
    pyscenic grn --num_workers $NUM_WORKERS \
        --sparse \
        --method grnboost2 \
        --output "$GRN_DIR/${SAMPLE_ID}_grn.csv" \
        "$LOOM_FILE" \
        "$TF_FILE"
    
    if [ $? -ne 0 ]; then
        log "[$SAMPLE_ID] GRN inference failed!"
        return 1
    fi
    log "[$SAMPLE_ID] GRN inference completed. Results saved to: $GRN_DIR/${SAMPLE_ID}_grn.csv"
    
    # Step 2: Regulon analysis (CTX)
    log "[$SAMPLE_ID] Starting regulon analysis (CTX)..."
    pyscenic ctx --num_workers $NUM_WORKERS \
        --output "$CTX_DIR/${SAMPLE_ID}_regulons.csv" \
        --expression_mtx_fname "$LOOM_FILE" \
        --all_modules \
        --mode "$MODE" \
        --annotations_fname "$MOTIF_ANNOTATIONS" \
        --mask_dropouts \
        "$GRN_DIR/${SAMPLE_ID}_grn.csv" \
        "$RANKINGS_FILE"
    
    if [ $? -ne 0 ]; then
        log "[$SAMPLE_ID] Regulon analysis failed!"
        return 1
    fi
    log "[$SAMPLE_ID] Regulon analysis completed. Results saved to: $CTX_DIR/${SAMPLE_ID}_regulons.csv"
    
    # Step 3: AUCell scoring
    log "[$SAMPLE_ID] Starting AUCell scoring..."
    pyscenic aucell --num_workers $NUM_WORKERS \
        --output "$AUCELL_DIR/${SAMPLE_ID}_aucell.loom" \
        "$LOOM_FILE" \
        "$CTX_DIR/${SAMPLE_ID}_regulons.csv"
    
    if [ $? -ne 0 ]; then
        log "[$SAMPLE_ID] AUCell scoring failed!"
        return 1
    fi
    log "[$SAMPLE_ID] AUCell scoring completed. Results saved to: $AUCELL_DIR/${SAMPLE_ID}_aucell.loom"
    
    log "Sample $SAMPLE_ID processing completed"
    return 0
}

# Main execution block
main() {
    # Initialize log file
    echo "==== SCENIC Pipeline Started at $(date) ====" > "$LOG_FILE"
    
    # Check if input directory exists
    if [ ! -d "$INPUT_DIR" ]; then
        log "Error: Input directory '$INPUT_DIR' does not exist!"
        exit 1
    fi
    
    # Find all matching loom files
    LOOM_FILES=($(ls $INPUT_DIR/$PATTERN 2>/dev/null))
    
    if [ ${#LOOM_FILES[@]} -eq 0 ]; then
        log "Error: No files matching '$PATTERN' found in '$INPUT_DIR'!"
        exit 1
    fi
    
    log "Found ${#LOOM_FILES[@]} loom files to process in $INPUT_DIR"
    
    # Create all output directories
    create_directories
    
    # Iterate through all loom files and extract sample IDs
    SAMPLES=()
    for LOOM_FILE in "${LOOM_FILES[@]}"; do
        # Extract filename from full path and remove extension
        SAMPLE_ID=$(basename "$LOOM_FILE" .loom)
        SAMPLES+=("$SAMPLE_ID")
    done
    
    log "Sample list: ${SAMPLES[*]}"
    
    # Process all samples
    TOTAL=${#SAMPLES[@]}
    RUNNING=0
    COMPLETED=0
    
    # Sequential or parallel processing
    if [ "$1" == "--sequential" ]; then
        log "Processing samples in sequential mode"
        for SAMPLE in "${SAMPLES[@]}"; do
            process_sample "$SAMPLE"
            if [ $? -eq 0 ]; then
                COMPLETED=$((COMPLETED+1))
            fi
        done
    else
        log "Processing samples in parallel mode. Maximum parallel jobs: $MAX_PARALLEL"
        
        # Track PIDs of all background processes
        declare -A PIDS
        
        for SAMPLE in "${SAMPLES[@]}"; do
            # Wait until running processes are less than maximum parallel limit
            while [ $RUNNING -ge $MAX_PARALLEL ]; do
                # Check which processes have completed
                for PID in "${!PIDS[@]}"; do
                    if ! kill -0 $PID 2>/dev/null; then
                        SAMPLE_DONE="${PIDS[$PID]}"
                        wait $PID
                        EXIT_CODE=$?
                        
                        if [ $EXIT_CODE -eq 0 ]; then
                            log "Sample $SAMPLE_DONE completed successfully"
                            COMPLETED=$((COMPLETED+1))
                        else
                            log "Sample $SAMPLE_DONE processing failed with exit code: $EXIT_CODE"
                        fi
                        
                        unset PIDS[$PID]
                        RUNNING=$((RUNNING-1))
                    fi
                done
                
                # Short sleep to avoid high CPU usage
                sleep 5
            done
            
            # Start processing a new sample
            log "Starting process for sample $SAMPLE..."
            process_sample "$SAMPLE" &
            PID=$!
            PIDS[$PID]="$SAMPLE"
            RUNNING=$((RUNNING+1))
            
            log "Currently running: $RUNNING, Completed: $COMPLETED, Total: $TOTAL"
        done
        
        # Wait for all remaining processes to complete
        log "Waiting for all remaining processes to complete..."
        for PID in "${!PIDS[@]}"; do
            SAMPLE="${PIDS[$PID]}"
            wait $PID
            EXIT_CODE=$?
            
            if [ $EXIT_CODE -eq 0 ]; then
                log "Sample $SAMPLE completed successfully"
                COMPLETED=$((COMPLETED+1))
            else
                log "Sample $SAMPLE processing failed with exit code: $EXIT_CODE"
            fi
        done
    fi
    
    # Generate result statistics report
    log "All samples processed. Successfully completed: $COMPLETED/$TOTAL"
    log "Result files are saved in: $OUTPUT_DIR"
    
    # Count files in each directory
    GRN_COUNT=$(ls -1 "$GRN_DIR"/*.csv 2>/dev/null | wc -l)
    CTX_COUNT=$(ls -1 "$CTX_DIR"/*.csv 2>/dev/null | wc -l)
    AUCELL_COUNT=$(ls -1 "$AUCELL_DIR"/*.loom 2>/dev/null | wc -l)
    
    log "Result statistics:"
    log "  GRN files: $GRN_COUNT (saved in $GRN_DIR)"
    log "  CTX files: $CTX_COUNT (saved in $CTX_DIR)"
    log "  AUCell files: $AUCELL_COUNT (saved in $AUCELL_DIR)"
    
    echo "==== SCENIC Pipeline Ended at $(date) ====" >> "$LOG_FILE"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --tf-file)
            TF_FILE="$2"
            shift; shift
            ;;
        --motif-annotations)
            MOTIF_ANNOTATIONS="$2"
            shift; shift
            ;;
        --rankings-file)
            RANKINGS_FILE="$2"
            shift; shift
            ;;
        --num-workers)
            NUM_WORKERS="$2"
            shift; shift
            ;;
        --mode)
            MODE="$2"
            shift; shift
            ;;
        --max-parallel)
            MAX_PARALLEL="$2"
            shift; shift
            ;;
        --input-dir)
            INPUT_DIR="$2"
            shift; shift
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            # Update subdirectory paths
            GRN_DIR="$OUTPUT_DIR/GRN"
            CTX_DIR="$OUTPUT_DIR/CTX"
            AUCELL_DIR="$OUTPUT_DIR/AUCell"
            shift; shift
            ;;
        --sequential)
            SEQUENTIAL="true"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --tf-file <file>             Transcription factor list file"
            echo "  --motif-annotations <file>   Motif annotations file"
            echo "  --rankings-file <file>       Ranking database file"
            echo "  --num-workers <number>       Number of worker processes (Default: 64)"
            echo "  --mode <mode>                Processing mode (Default: dask_multiprocessing)"
            echo "  --max-parallel <number>      Maximum parallel samples (Default: 5)"
            echo "  --input-dir <dir>            Input directory (Default: counts_data)"
            echo "  --output-dir <dir>           Output directory (Default: scenic_result)"
            echo "  --sequential                 Use sequential processing mode"
            exit 1
            ;;
    esac
done

# Run main program
if [ "$SEQUENTIAL" == "true" ]; then
    main --sequential
else
    main
fi
