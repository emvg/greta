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


Rscript "workflow/scripts/mth/deepmaps/deepmaps.r" \
  -m "$path_mdata" \
  -d "$out_dir" \
  -o "$path_out"