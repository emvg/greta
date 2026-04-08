rule download_rpe_choroid:
    threads: 8
    singularity: 'workflow/envs/figr.sif'
    input: 'workflow/envs/figr.sif'
    output:
        frags=expand(
            'dts/hg38/rpe_choroid/{celltype}.frags.tsv.gz',
            celltype=config['dts']['rpe_choroid']['celltypes']
        ),
        tbis=expand(
            'dts/hg38/rpe_choroid/{celltype}.frags.tsv.gz.tbi',
            celltype=config['dts']['rpe_choroid']['celltypes']
        ),
        barcodes=temp(local('dts/hg38/rpe_choroid/barcodes.tsv.gz')),
        features=temp(local('dts/hg38/rpe_choroid/features.tsv.gz')),
        matrix=temp(local('dts/hg38/rpe_choroid/matrix.mtx.gz')),
    params:
        celltypes=config['dts']['rpe_choroid']['celltypes'],
        url_barcodes=config['dts']['rpe_choroid']['url']['barcodes'],
        url_features=config['dts']['rpe_choroid']['url']['features'],
        url_matrix=config['dts']['rpe_choroid']['url']['matrix'],
        url_frags_base=config['dts']['rpe_choroid']['url']['frags_base'],
    resources:
        mem_mb=8000,
    shell:
        """
        data_path=$(dirname {output.barcodes})
        # Download RNA matrix files
        wget --no-verbose '{params.url_barcodes}' -O '{output.barcodes}'
        wget --no-verbose '{params.url_features}' -O '{output.features}'
        wget --no-verbose '{params.url_matrix}' -O '{output.matrix}'
        # Download per-celltype ATAC fragment files
        for celltype in {params.celltypes}; do
            url="{params.url_frags_base}GSE262151%5F${{celltype}}%5Ffrags%2Etsv%2Egz"
            echo "Downloading fragments for $celltype: $url"
            wget --no-verbose "$url" -O "$data_path/${{celltype}}.frags.tsv.gz"
        done
        # Format all fragment files
        ls $data_path/*.frags.tsv.gz | xargs -n 1 -P {threads} bash workflow/scripts/dts/format_frags.sh
        """


rule prc_annot_rpe_choroid:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        img='workflow/envs/gretabench.sif',
        frags=rules.download_rpe_choroid.output.frags,
        barcodes=rules.download_rpe_choroid.output.barcodes
    output:
        annot=temp(local('dts/hg38/rpe_choroid/annot.csv')),
    shell:
        """
        python workflow/scripts/dts/rpe_choroid/prc_annot.py \
        -f {input.frags} \
        -b {input.barcodes} \
        -o {output.annot}
        """


rule callpeaks_rpe_choroid:
    threads: 8
    singularity: 'workflow/envs/gretabench.sif'
    input:
        frags=rules.download_rpe_choroid.output.frags,
        annot=rules.prc_annot_rpe_choroid.output.annot,
    output:
        peaks=temp(local('dts/hg38/rpe_choroid/peaks.h5ad'))
    resources:
        mem_mb=110000,
        runtime=400,
    shell:
        """
        python workflow/scripts/dts/callpeaks.py \
        -f {input.frags} \
        -a {input.annot} \
        -t /workdir/vangysel/tmp/ \
        -n {threads} \
        -o {output.peaks}
        """


rule annotate_rpe_choroid:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        path_peaks=rules.callpeaks_rpe_choroid.output.peaks,
        path_annot=rules.prc_annot_rpe_choroid.output.annot,
        barcodes=rules.download_rpe_choroid.output.barcodes,
        features=rules.download_rpe_choroid.output.features,
        matrix=rules.download_rpe_choroid.output.matrix,
        gid=rules.gen_gid_ensmbl.output.hg38,
    output: out='dts/hg38/rpe_choroid/annotated.h5mu'
    resources: mem_mb=32000
    shell:
        """
        python workflow/scripts/dts/rpe_choroid/rpe_choroid.py \
        -a {input.barcodes} \
        -b {input.features} \
        -c {input.matrix} \
        -d {input.path_peaks} \
        -e {input.path_annot} \
        -f {input.gid} \
        -o {output.out}
        """