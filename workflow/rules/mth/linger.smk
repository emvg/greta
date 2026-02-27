rule mdl_o_linger:
    threads: 4 #32 for dragon2
    singularity: 'workflow/envs/linger.sif'
    input:
        img='workflow/envs/linger.sif',
        linger_GRN=rules.linger_prior.output.dir,
        mdata=rules.extract_case.output.mdata
    output:
        dir=directory('dts/{org}/{dat}/cases/{case}/runs/linger/'),
        out='dts/{org}/{dat}/cases/{case}/runs/o_linger.o_linger.o_linger.o_linger.mdl.csv'
    params:
        version=config['methods']['linger']['version'],
        script='workflow/scripts/mth/linger/linger.sh'
    resources:
        mem_mb=32000,
        runtime=30
    shell:
        """
        mkdir -p {output.dir}
        set -e
        timeout $(({resources.runtime}-20))m \
        bash {params.script} \
        --linger_GRN {input.linger_GRN} \
        --out_dir {output.dir} \
        --path_mdata {input.mdata} \
        --version {params.version} \
        --path_out {output.out}
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """