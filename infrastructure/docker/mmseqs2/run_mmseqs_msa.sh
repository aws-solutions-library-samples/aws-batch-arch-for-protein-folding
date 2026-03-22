#!/bin/bash

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# MMSeqs2 MSA script
# This script runs the MMSeqs2 workflow to generate multiple sequence alignments (MSAs)
# against specified reference databases.

set -e

# Source logging library
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"
source "$SCRIPT_DIR/logging.sh"

# Default values
THREADS=16
EVALUE="0.1"
COV=0.5
COV_MODE=2
MAX_SEQS=10000
SENSITIVITY=7.5
NUM_ITERATIONS=3
OUTPUT_FORMAT="a3m"
WRITE_HITS_TABLE=1
GPU=0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --query_fasta)
            QUERY_FASTA="$2"
            shift 2
            ;;
        --run_dir)
            RUN_DIR="$2"
            shift 2
            ;;
        --databases_paths)
            DATABASES_PATHS="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --evalue)
            EVALUE="$2"
            shift 2
            ;;
        --cov)
            COV="$2"
            shift 2
            ;;
        --cov_mode)
            COV_MODE="$2"
            shift 2
            ;;
        --max_seqs)
            MAX_SEQS="$2"
            shift 2
            ;;
        --s)
            SENSITIVITY="$2"
            shift 2
            ;;
        --num_iterations)
            NUM_ITERATIONS="$2"
            shift 2
            ;;
        --output_format)
            OUTPUT_FORMAT="$2"
            shift 2
            ;;
        --write_hits_table)
            WRITE_HITS_TABLE="$2"
            shift 2
            ;;
        --gpu)
            GPU="$2"
            shift 2
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$QUERY_FASTA" ] || [ -z "$RUN_DIR" ] || [ -z "$DATABASES_PATHS" ]; then
    log_error "--query_fasta, --run_dir, and --databases_paths are required"
    exit 1
fi

# Build resolved parameters JSON
RESOLVED_PARAMS_JSON=$(cat <<EOF
{
  "query_fasta": "$QUERY_FASTA",
  "run_dir": "$RUN_DIR",
  "databases_paths": "$DATABASES_PATHS",
  "threads": $THREADS,
  "evalue": "$EVALUE",
  "cov": $COV,
  "cov_mode": $COV_MODE,
  "max_seqs": $MAX_SEQS,
  "sensitivity": $SENSITIVITY,
  "num_iterations": $NUM_ITERATIONS,
  "output_format": "$OUTPUT_FORMAT",
  "write_hits_table": $WRITE_HITS_TABLE,
  "gpu": $GPU
}
EOF
)

# Get total RAM in GB, use 80% of it
TOTAL_MEM_GB=$(free -g | awk '/^Mem:/{print $2}')
SPLIT_MEM_GB=$(( TOTAL_MEM_GB * 80 / 100 ))
SPLIT_MEM_LIMIT="${SPLIT_MEM_GB}G"

# Create directory structure
mkdir -p "$RUN_DIR/tmp"

# Write resolved parameters to config
echo "$RESOLVED_PARAMS_JSON" > "$RUN_DIR/mmseqs_config.json"

# Setup logging - tee output to log file
LOG_FILE="$RUN_DIR/mmseqs_logs.log"
exec 1> >(tee -a "$LOG_FILE")
exec 2>&1

log_section "MMSeqs2 MSA Generation Job Started"
log_info "Query FASTA: $QUERY_FASTA"
log_info "Run Directory: $RUN_DIR"
log_info "Databases Paths: $DATABASES_PATHS"
log_debug "Threads: $THREADS"
log_debug "E-value: $EVALUE"
log_debug "Coverage: $COV"
log_debug "Coverage Mode: $COV_MODE"
log_debug "Max Sequences: $MAX_SEQS"
log_debug "Sensitivity: $SENSITIVITY"
log_debug "Number of Iterations: $NUM_ITERATIONS"
log_debug "Output Format: $OUTPUT_FORMAT"
log_debug "Write Hits Table: $WRITE_HITS_TABLE"
log_debug "Split Memory Limit: $SPLIT_MEM_LIMIT"
if [ "$GPU" -eq 1 ]; then
    log_info "GPU Mode: ENABLED (will use all available GPUs)"
else
    log_debug "GPU Mode: DISABLED"
fi


# Create query database (only once)
log_info "Creating query database..."
run_quiet mmseqs createdb "$QUERY_FASTA" "$RUN_DIR/queryDB"

# Process each database - search and store results
# Format: dbname:dbpath,dbname:dbpath (e.g., uniref90:uniref90/uniref90,swissprot:swissprot/SwissProt)
IFS=',' read -ra DB_PATH_ARRAY <<< "$DATABASES_PATHS"
RESULT_DBS=()  # Track result databases for later merging
DB_NAMES=()    # Track database names for summary

for db_entry in "${DB_PATH_ARRAY[@]}"; do
    db_entry=$(echo "$db_entry" | xargs)  # trim whitespace
    # Parse dbname:dbpath format
    dbname="${db_entry%%:*}"
    dbpath="${db_entry#*:}"
    # Create safe name for result files (replace / with _)
    safe_dbname=$(echo "$dbname" | tr '/' '_')

    log_section "Processing database: $dbname (path: $dbpath)"

    TARGET_DB="/database/${dbpath}"

    # Check if target database exists
    if [ ! -f "$TARGET_DB.dbtype" ]; then
        log_error "Target database not found at $TARGET_DB"
        exit 1
    fi

    mkdir -p "$RUN_DIR/tmp/${safe_dbname}"

    # Step 1: Search (with filtering parameters)
    # Build GPU flag if enabled (GPU only applies to search step)
    GPU_FLAG=""
    if [ "$GPU" -eq 1 ]; then
        GPU_FLAG="--gpu 1"
    fi

    log_info "Running search against $dbname..."
    run_quiet mmseqs search "$RUN_DIR/queryDB" "$TARGET_DB" \
        "$RUN_DIR/result_${safe_dbname}" "$RUN_DIR/tmp/${safe_dbname}" \
        --threads "$THREADS" \
        -e "$EVALUE" \
        --max-seqs "$MAX_SEQS" \
        -s "$SENSITIVITY" \
        -c "$COV" \
        --cov-mode "$COV_MODE" \
        --num-iterations "$NUM_ITERATIONS" \
        --db-load-mode 0 \
        $GPU_FLAG

    # Step 2: Convert to MSA for this database (temp file)
    log_info "Converting results to MSA format for $dbname..."
    run_quiet mmseqs result2msa "$RUN_DIR/queryDB" "$TARGET_DB" \
        "$RUN_DIR/result_${safe_dbname}" \
        "$RUN_DIR/tmp/${safe_dbname}.${OUTPUT_FORMAT}"

    # Generate hits table for this database if requested
    if [ "$WRITE_HITS_TABLE" -eq 1 ]; then
        log_info "Generating hits table for $dbname..."
        run_quiet mmseqs convertalis "$RUN_DIR/queryDB" "$TARGET_DB" \
            "$RUN_DIR/result_${safe_dbname}" \
            "$RUN_DIR/mmseqs_hits_${safe_dbname}.tsv" \
            --format-mode 4 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits,qcov,tcov,theader \
            --threads "$THREADS"

        # Count and log hits
        HIT_COUNT=$(tail -n +2 "$RUN_DIR/mmseqs_hits_${safe_dbname}.tsv" | wc -l)
        log_info "Found $HIT_COUNT hits in $dbname"
    fi

    RESULT_DBS+=("$safe_dbname")
    DB_NAMES+=("$dbname")

    log_info "Completed search for $dbname"
done

# Combine all MSA files into single combined.a3m
log_section "Combining MSA results"

FIRST=true
for safe_dbname in "${RESULT_DBS[@]}"; do
    MSA_FILE="$RUN_DIR/tmp/${safe_dbname}.${OUTPUT_FORMAT}"
    if [ -f "$MSA_FILE" ]; then
        if [ "$FIRST" = true ]; then
            # First file: copy entire content (includes query header)
            cat "$MSA_FILE" > "$RUN_DIR/mmseqs_output.${OUTPUT_FORMAT}"
            FIRST=false
        else
            # Subsequent files: append sequences (skip first 2 lines which are query header)
            tail -n +3 "$MSA_FILE" >> "$RUN_DIR/mmseqs_output.${OUTPUT_FORMAT}"
        fi
    fi
done

log_info "Combined MSA written to: $RUN_DIR/mmseqs_output.${OUTPUT_FORMAT}"

# Create combined hits table with database source column
if [ "$WRITE_HITS_TABLE" -eq 1 ]; then
    log_section "Creating combined hits table"

    COMBINED_FILE="$RUN_DIR/mmseqs_hits_combined.tsv"
    FIRST_DB=true

    for safe_dbname in "${RESULT_DBS[@]}"; do
        TSV_FILE="$RUN_DIR/mmseqs_hits_${safe_dbname}.tsv"
        if [ -f "$TSV_FILE" ]; then
            if [ "$FIRST_DB" = true ]; then
                # First file: add "database" column to header
                head -1 "$TSV_FILE" | awk '{print $0"\tdatabase"}' > "$COMBINED_FILE"
                # Add data rows with database name
                tail -n +2 "$TSV_FILE" | awk -v db="$safe_dbname" '{print $0"\t"db}' >> "$COMBINED_FILE"
                FIRST_DB=false
            else
                # Subsequent files: skip header, append data with database name
                tail -n +2 "$TSV_FILE" | awk -v db="$safe_dbname" '{print $0"\t"db}' >> "$COMBINED_FILE"
            fi
        fi
    done

    TOTAL_HITS=$(tail -n +2 "$COMBINED_FILE" | wc -l)
    log_info "Combined hits table created: $COMBINED_FILE"
    log_info "Total hits across all databases: $TOTAL_HITS"

    # Generate summary statistics
    log_info "Generating summary statistics..."

    SUMMARY_FILE="$RUN_DIR/mmseqs_summary.json"

    # Start JSON structure
    echo "{" > "$SUMMARY_FILE"
    echo "  \"query_fasta\": \"$QUERY_FASTA\"," >> "$SUMMARY_FILE"
    echo "  \"databases\": [$(printf '\"%s\",' "${DB_NAMES[@]}" | sed 's/,$//')]," >> "$SUMMARY_FILE"
    echo "  \"parameters\": {" >> "$SUMMARY_FILE"
    echo "    \"evalue\": \"$EVALUE\"," >> "$SUMMARY_FILE"
    echo "    \"coverage\": $COV," >> "$SUMMARY_FILE"
    echo "    \"cov_mode\": $COV_MODE," >> "$SUMMARY_FILE"
    echo "    \"max_seqs\": $MAX_SEQS," >> "$SUMMARY_FILE"
    echo "    \"sensitivity\": $SENSITIVITY," >> "$SUMMARY_FILE"
    echo "    \"num_iterations\": $NUM_ITERATIONS" >> "$SUMMARY_FILE"
    echo "  }," >> "$SUMMARY_FILE"
    echo "  \"per_database_stats\": {" >> "$SUMMARY_FILE"

    # Calculate per-database statistics
    FIRST=true
    for safe_dbname in "${RESULT_DBS[@]}"; do
        TSV_FILE="$RUN_DIR/mmseqs_hits_${safe_dbname}.tsv"
        if [ -f "$TSV_FILE" ]; then
            HIT_COUNT=$(tail -n +2 "$TSV_FILE" | wc -l)
            AVG_IDENTITY=$(tail -n +2 "$TSV_FILE" | awk '{sum+=$3; count++} END {if(count>0) printf "%.4f", sum/count; else print "0"}')
            AVG_QCOV=$(tail -n +2 "$TSV_FILE" | awk '{sum+=$15; count++} END {if(count>0) printf "%.4f", sum/count; else print "0"}')

            [ "$FIRST" = false ] && echo "," >> "$SUMMARY_FILE"
            echo "    \"${safe_dbname}\": {" >> "$SUMMARY_FILE"
            echo "      \"hit_count\": $HIT_COUNT," >> "$SUMMARY_FILE"
            echo "      \"avg_identity\": $AVG_IDENTITY," >> "$SUMMARY_FILE"
            echo "      \"avg_query_coverage\": $AVG_QCOV" >> "$SUMMARY_FILE"
            echo -n "    }" >> "$SUMMARY_FILE"
            FIRST=false
        fi
    done

    # Calculate combined statistics
    echo "" >> "$SUMMARY_FILE"
    echo "  }," >> "$SUMMARY_FILE"
    AVG_IDENTITY_COMBINED=$(tail -n +2 "$COMBINED_FILE" | awk '{sum+=$3; count++} END {if(count>0) printf "%.4f", sum/count; else print "0"}')
    echo "  \"combined_stats\": {" >> "$SUMMARY_FILE"
    echo "    \"total_hits\": $TOTAL_HITS," >> "$SUMMARY_FILE"
    echo "    \"avg_identity\": $AVG_IDENTITY_COMBINED" >> "$SUMMARY_FILE"
    echo "  }" >> "$SUMMARY_FILE"
    echo "}" >> "$SUMMARY_FILE"

    log_info "Summary statistics written to: $SUMMARY_FILE"
fi

# Cleanup ALL intermediate files
log_debug "Cleaning up intermediate files..."
rm -f "$RUN_DIR/queryDB"*
rm -rf "$RUN_DIR/tmp"
for safe_dbname in "${RESULT_DBS[@]}"; do
    rm -f "$RUN_DIR/result_${safe_dbname}"*
done

log_section "MMSeqs2 MSA Generation Job Completed"
log_info "Outputs available at:"
log_info "  MSA: $RUN_DIR/mmseqs_output.${OUTPUT_FORMAT}"
if [ "$WRITE_HITS_TABLE" -eq 1 ]; then
    log_info "  Combined Hits: $RUN_DIR/mmseqs_hits_combined.tsv"
    for safe_dbname in "${RESULT_DBS[@]}"; do
        log_info "    ${safe_dbname} Hits: $RUN_DIR/mmseqs_hits_${safe_dbname}.tsv"
    done
    log_info "  Summary Stats: $RUN_DIR/mmseqs_summary.json"
fi
log_info "  Config: $RUN_DIR/mmseqs_config.json"
log_info "  Logs: $LOG_FILE"
