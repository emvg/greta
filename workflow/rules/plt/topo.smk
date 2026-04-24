localrules: plt_topo

rule plt_topo:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        sims='anl/topo/{org}.{dat}.{case}.sims_mult.csv',
        stats='anl/topo/{org}.{dat}.{case}.stats_mult.csv',
    output: 'plt/topo/{org}.{dat}.{case}.topo.pdf'
    shell:
        """
        python workflow/scripts/plt/topo.py \
        -s {input.sims} \
        -t {input.stats} \
        -o {output}
        """
