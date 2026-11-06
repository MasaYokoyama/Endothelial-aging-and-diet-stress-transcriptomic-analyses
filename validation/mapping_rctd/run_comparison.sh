#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -ne 1 ]; then
  echo 'Usage: bash run_comparison.sh config.local.R' >&2
  exit 1
fi
Rscript 01_preflight.R "$1"
Rscript 02_compare_mappings.R "$1"
Rscript 03_rerun_paired_rctd.R "$1"
Rscript 04_summarize_results.R "$1"
