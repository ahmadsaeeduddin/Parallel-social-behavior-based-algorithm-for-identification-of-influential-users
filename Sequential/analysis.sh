#!/bin/bash

# Load Intel VTune environment
source /opt/intel/oneapi/vtune/latest/env/vars.sh

# Define your sequential executable
EXEC=./psaiim_seq

# Output directory for results
OUTPUT_DIR=vtune_seq

echo "[RUNNING VTUNE HOTSPOTS ANALYSIS: SEQUENTIAL]"
vtune -collect hotspots -result-dir $OUTPUT_DIR $EXEC

echo "[DONE] Results saved to $OUTPUT_DIR"
