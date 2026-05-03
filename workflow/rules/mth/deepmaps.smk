rule mdl_o_deepmaps:
    threads: 12
    singularity: 'workflow/envs/deepmaps.sif'
    input:
        img='workflow/envs/deepmaps.sif',
        mdata=rules.extract_case.output.mdata
    output:
        dir=directory('dts/{org}/{dat}/cases/{case}/runs/deepmaps/'),
        out='dts/{org}/{dat}/cases/{case}/runs/o_deepmaps.o_deepmaps.o_deepmaps.o_deepmaps.mdl.csv'
    params:
        script='workflow/scripts/mth/deepmaps/deepmaps.sh'
    resources:
        mem_mb= lambda wildcards, attempt: restart_mem(wildcards, attempt) * 4,
        runtime=300,
        partition='gpu',
        slurm="--gres=gpu:GeForceRTX2080Ti:1" #"gres=gpu:1",
    shell:
        """
        export HDF5_USE_FILE_LOCKING=FALSE
        mkdir -p {output.dir}
        set -e
        timeout $(({resources.runtime}-20))m \
        bash {params.script} \
        --path_mdata {input.mdata} \
        --out_dir {output.dir} \
        --path_out {output.out}
        """


# snakemake --profile config/slurm dts/hg38/pbmc10k/cases/all/runs/o_deepmaps.o_deepmaps.o_deepmaps.o_deepmaps.mdl.csv
# snakemake -j 16 dts/{org}/{dat}/cases/{case}/runs/o_deepmaps.o_deepmaps.o_deepmaps.o_deepmaps.mdl.csv