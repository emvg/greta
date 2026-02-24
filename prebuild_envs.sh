#!/bin/bash
# Usage: ./prebuild_envs.sh 
# Builds only specified YAMLs into .snakemake/conda/<hash>_

set -euo pipefail

SNK_CONDA=".snakemake/conda"
mkdir -p "$SNK_CONDA"

# list of config files that require a conda env
CONDA_ENV=("workflow/envs/scdori.yaml" "workflow/envs/dictys.yaml" "workflow/envs/scgpt.yaml")

if [ "$#" -eq 0 ]; then
    echo "Error: Please provide at least one conda YAML file to prebuild."
    exit 1
fi

for yaml in "$@"; do
    if [ ! -f "$yaml" ]; then
        echo "Warning: YAML file '$yaml' not found, skipping."
        continue
    fi

    # Compute hash (same as Snakemake does)
    HASH=$(md5sum "$yaml" | awk '{print $1}')_
    TARGET="$SNK_CONDA/$HASH"

    echo "Prebuilding Conda environment for $yaml -> $TARGET"

    # Remove any leftover folder / yaml
    rm -rf "$TARGET" "$SNK_CONDA/$HASH.yaml"

    # Create environment
    mamba env create --file "$yaml" --prefix "$TARGET"

    echo "Done: $TARGET"
done

echo "All requested environments prebuilt in $SNK_CONDA"