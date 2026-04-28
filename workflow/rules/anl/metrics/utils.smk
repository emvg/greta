localrules: aggr_metric, metric_summ


"""
snakemake run_mech_metrics --profile config/slurm/
snakemake run_pred_metrics --profile config/slurm/
snakemake run_prior_metrics --profile config/slurm/
snakemake metric_aggr --profile config/slurm/
""" 

"""
# Run the similarity metrics
snakemake --profile config/slurm/ anl/topo/hg38.pbmc10k.all.sims_mult.csv

# Run the graph
snakemake --profile config/slurm/ plt/topo/hg38.pbmc10k.all.topo.pdf
"""

"""
# Run stability metrics on LINGER

snakemake --profile config/slurm/ anl/stab/pitupair.ovc.csv -n
"""

rule aggr_metric:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        lambda w: make_combs_rules(w=w, rule_name='{typ}_{tsk}'.format(typ=w.type, tsk=w.task), do_decoupling=True)
    output:
        'anl/metrics/{type}/{task}/{db}/{org}.{dat}.{case}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/aggregate.py \
        -i {input} \
        -o {output}
        """

# Datasets to run metrics on
metric_dts = ['pbmc10k']  # , 'brain', 'rpe_choroid'

# Methods to run metrics on   'figr',
metric_mths = [
    'celloracle',
    'collectri',  
    'deepmaps',
    'dictys',
    'dorothea',
    'granie',
    'linger',
    'linger_baseline',
    'pando',
    'random',
    'scenic',
    'scenicplus'
]


def make_mech_rules(dat):
    org = config['dts'][dat]['organism']
    case = 'all'
    return [
        f'anl/metrics/mech/prt/knocktf/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/mech/tfa/knocktf/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/mech/sss/sss/{org}.{dat}.{case}.scores.csv',
    ]


def make_pred_rules(dat):
    org = config['dts'][dat]['organism']
    case = 'all'
    paths = [
        f'anl/metrics/pred/omics/gtf/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/pred/omics/cretf/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/pred/omics/gcre/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/pred/gsets/hall/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/pred/gsets/reac/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/pred/gsets/prog/{org}.{dat}.{case}.scores.csv',
    ]
    if org == 'hg38':
        paths += [f'anl/metrics/pred/gsets/kegg/{org}.{dat}.{case}.scores.csv']
    return paths


def make_prior_rules(dat):
    org = config['dts'][dat]['organism']
    case = 'all'
    paths = [
        f'anl/metrics/prior/grn/collectri/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/genom/tfb/chipatlas/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/genom/tfb/remap2022/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/genom/tfb/unibind/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/genom/cre/blacklist/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/genom/cre/encode/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/genom/cre/phastcons/{org}.{dat}.{case}.scores.csv',
        f'anl/metrics/genom/cre/promoters/{org}.{dat}.{case}.scores.csv',
    ]
    if org == 'hg38':
        paths += [
            f'anl/metrics/prior/tfm/hpa/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfm/tfmdb/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfp/europmc/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfp/intact/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/gwascatalogue/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/zhang21/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/c2g/eqtlcatalogue/{org}.{dat}.{case}.scores.csv',
        ]
    return paths

def make_metric_rules(dat):
    org = config['dts'][dat]['organism']
    case = 'all'
    if org == 'hg38':
        return [
            f'anl/metrics/mech/prt/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/tfa/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/sss/sss/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gtf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/cretf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gcre/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/kegg/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/hall/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/reac/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/prog/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfm/hpa/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfm/tfmdb/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfp/europmc/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/tfp/intact/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/grn/collectri/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/chipatlas/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/remap2022/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/unibind/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/blacklist/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/encode/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/gwascatalogue/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/phastcons/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/zhang21/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/promoters/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/c2g/eqtlcatalogue/{org}.{dat}.{case}.scores.csv',
        ]
    elif org == 'mm10':
        return [
            f'anl/metrics/mech/prt/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/tfa/knocktf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/mech/sss/sss/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gtf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/cretf/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/omics/gcre/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/hall/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/reac/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/pred/gsets/prog/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/prior/grn/collectri/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/chipatlas/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/remap2022/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/tfb/unibind/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/blacklist/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/encode/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/phastcons/{org}.{dat}.{case}.scores.csv',
            f'anl/metrics/genom/cre/promoters/{org}.{dat}.{case}.scores.csv',
        ]

rule run_mech_metrics:
    input: [make_mech_rules(dat) for dat in metric_dts]

rule run_pred_metrics:
    input: [make_pred_rules(dat) for dat in metric_dts]

rule run_prior_metrics:
    input: [make_prior_rules(dat) for dat in metric_dts]

rule metric_aggr:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        [make_metric_rules(dat=dat) for dat in metric_dts]
    output:
        'anl/metrics/summary/metrics.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/aggr_all.py {output}
        """


rule metric_summ:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        aggr=[make_metric_rules(dat=dat) for dat in config['dts'].keys()],
        scale='anl/stab/pitupair.ovc.csv',
        dfs='anl/topo/hg38.pitupair.topo_diff_seed.csv'
    output:
        metrics='anl/metrics/summary/metrics.csv',
        scale='anl/metrics/summary/scalability.csv',
        pair='anl/metrics/summary/pair.csv',
    shell:
        """
        python workflow/scripts/anl/metrics/aggr_all.py {output.metrics} && \
        python workflow/scripts/anl/metrics/scalability.py {input.scale} {input.dfs} {output.scale} && \
        python workflow/scripts/anl/metrics/pair.py {output.metrics} {output.pair}
        """
