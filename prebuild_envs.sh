#!/bin/bash
SNAKEMAKE=/home/ucl/inma/vangysel/.conda/envs/greta/bin/snakemake

$SNAKEMAKE --ri --conda-create-envs-only \
    --conda-frontend conda \
    --profile config/slurm/ \
    dts/hg38/pbmc10k/cases/all/runs/o_dictys.o_dictys.o_dictys.o_dictys.mdl.csv \
    dts/hg38/pbmc10k/cases/all/runs/o_scdori.o_scdori.o_scdori.o_scdori.mdl.csv \
    dts/hg38/pbmc10k/cases/all/runs/o_scgpt.o_scgpt.o_scgpt.o_scgpt.mdl.csv