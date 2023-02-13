import shutil
shutil.copyfile(snakemake.input.PretextMap, snakemake.output.pretextCP2keyResults)
with open(snakemake.output.IndividualPretextMD, "w") as out:
	print("", file=out)
	print("###",snakemake.params.assemblyName," HiC Heatmap", file=out)
	print("![](", snakemake.input.PretextMap, "){ height=30% }", file=out)
