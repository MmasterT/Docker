import pandas as pd
from tabulate import tabulate
samples=pd.read_csv(snakemake.input.results, dtype=str, index_col=0, delim_whitespace=True, skip_blank_lines=True)
# samples=samples.reset_index
turn2FloatAndMb=['Bases_Mb','Est_Size_Mb','Longest_Scaff_Mb','Scaff_N50_Mb','Scaff_NG50_Mb','Scaff_N95_Mb','Scaff_NG95_Mb','Longest_Cont_Mb','Cont_N50_Mb','Cont_NG50_Mb','Cont_N95_Mb','Cont_NG95_Mb']
roundDecimals=['Comp_BUSCOs_%','Comp_Single_BUSCOs_%','Het_%','GC_%','QV','Completeness','Bases_Mb','Est_Size_Mb','Longest_Scaff_Mb','Scaff_N50_Mb','Scaff_NG50_Mb','Scaff_N95_Mb','Scaff_NG95_Mb','Longest_Cont_Mb','Cont_N50_Mb','Cont_NG50_Mb','Cont_N95_Mb','Cont_NG95_Mb']
print('this is samples',samples)
sampleTransposed=samples.T
print('this is sampleTransposed',sampleTransposed)
for header in turn2FloatAndMb:
	sampleTransposed[header]=sampleTransposed[header].astype(float).div(1000000)
for i in range(0,4):
	sampleTransposed[roundDecimals[i]]=sampleTransposed[roundDecimals[i]].str.replace('%','')
for roundHeader in roundDecimals:
	sampleTransposed[roundHeader]=sampleTransposed[roundHeader].astype(float).round(2)
with open(snakemake.output[0], "w") as out_markdown:
	print("", file=out_markdown)
	print("\\blandscape", file=out_markdown)
	print("", file=out_markdown)
	print("\\tiny", file=out_markdown)
	print(tabulate(sampleTransposed.rename_axis('ASM_ID'), headers='keys',tablefmt="pipe", showindex=True), file=out_markdown)
	print("\\elandscape", file=out_markdown)
	print("", file=out_markdown)

with open(snakemake.output[1], "w") as out_plain:
	print(tabulate(sampleTransposed.rename_axis('ASM_ID'), headers='keys',tablefmt="plain", showindex=True), file=out_plain)
