#!/bin/bash

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --path_mdata) path_mdata="$2"; shift ;;
    --out_dir)    out_dir="$2";    shift ;;
    --path_out)   path_out="$2";   shift ;;
    *) echo "Unknown parameter: $1"; exit 1 ;;
  esac
  shift
done


export TMPDIR="/workdir/vangysel/tmp/${SLURM_JOB_ID:-tmp}"
mkdir -p "$TMPDIR"

Rscript "workflow/scripts/mth/deepmaps/deepmaps.R" \
  -m "$path_mdata" \
  -d "$out_dir" \
  -o "$path_out"

Rscript "workflow/scripts/mth/deepmaps/grn.R" \
  -d "$out_dir" \
  -o "$path_out"