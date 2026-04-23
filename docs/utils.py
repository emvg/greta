import sys

def create_output_name(dataset, method):
    return "dts/hg38/{}/cases/all/runs/o_{}.o_{}.o_{}.o_{}.mdl.csv".format(
        dataset, method, method, method, method
    )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python utils.py <dataset> <method>")
        sys.exit(1)

    dataset = sys.argv[1]
    method = sys.argv[2]

    output_name = create_output_name(dataset, method)
    print(output_name)

# export HDF5_USE_FILE_LOCKING=FALSE

"""
dts/hg38/pbmc10k/cases/all/runs/o_celloracle.o_celloracle.o_celloracle.o_celloracle.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_collectri.o_collectri.o_collectri.o_collectri.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_crema.o_crema.o_crema.o_crema.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_dictys.o_dictys.o_dictys.o_dictys.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_directnet.o_directnet.o_directnet.o_directnet.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_dorothea.o_dorothea.o_dorothea.o_dorothea.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_figr.o_figr.o_figr.o_figr.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_granie.o_granie.o_granie.o_granie.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_grnboost.o_grnboost.o_grnboost.o_grnboost.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_hummus.o_hummus.o_hummus.o_hummus.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_linger.o_linger.o_linger.o_linger.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_pando.o_pando.o_pando.o_pando.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_pearson.o_pearson.o_pearson.o_pearson.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_random.o_random.o_random.o_random.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_scenic.o_scenic.o_scenic.o_scenic.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_scenicplus.o_scenicplus.o_scenicplus.o_scenicplus.mdl.csv
dts/hg38/pbmc10k/cases/all/runs/o_spearman.o_spearman.o_spearman.o_spearman.mdl.csv
"""