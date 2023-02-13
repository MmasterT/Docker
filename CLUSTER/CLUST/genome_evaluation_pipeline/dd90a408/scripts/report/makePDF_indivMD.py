import os
import numpy as np
import pandas as pd
from tabulate import tabulate

key_stats=snakemake.input[0]
#key_stats="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/L_sceleratus_aggregatedResults.tsv"
key_stats_table = pd.read_csv(key_stats, dtype=str, delim_whitespace=True, skip_blank_lines=True)
key_stats_table = key_stats_table.iloc[:-4 , :]

genomescope_linear=snakemake.input[1]
#genomescope_linear="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/L_sceleratus_gScope_linear_plot.png"
genomescope_log=snakemake.input[2]
#genomescope_log="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/L_sceleratus_gScope_log_plot.png"

QV_score=snakemake.input[3]
#QV_score="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/L_sceleratus_merq.qv"
QV_score_table = pd.read_csv(QV_score, dtype=str, delim_whitespace=True, skip_blank_lines=True, names=["assembly_name","Assembly_Only", "Total_Kmers", "QV_Score", "Error"])
QV_score_table	= QV_score_table.drop(['Error'], axis=1).set_index(['assembly_name'])

kmer_completeness=snakemake.input[4]
##kmer_completeness="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/L_sceleratus_merq.completeness.stats"
kmer_completeness_table = pd.read_csv(kmer_completeness, dtype=str, delim_whitespace=True, skip_blank_lines=True, names=["assembly_name","all","Kmers_Assembly", "Kmers_Reads", "%"])
kmer_completeness_table	= kmer_completeness_table.drop(["all"], axis=1).set_index(['assembly_name'])


busco=snakemake.input[5]
#busco="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/only_busco_scores.txt"

with open(busco) as scoresOnlyFile:
    lines = scoresOnlyFile.readlines()
    lines = [line.rstrip() for line in lines]

kmer_flat=snakemake.input[6]
#kmer_flat="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/L_sceleratus_merq.L_sceleratus.spectra-cn.fl.png"

kmer_stacked=snakemake.input[7]
#kmer_stacked="/srv/public/users/james94/data/Lsceleratus/GEP_results/L_sceleratus/5_Key_Results/L_sceleratus_merq.L_sceleratus.spectra-cn.st.png"


params_assemblyID=snakemake.params[0]
#params_assemblyID="Lagocephalus sceleratus"

params_merylKmer=snakemake.params[1]
# params_merylKmer="21"

params_buscoDB=snakemake.params[2]
# params_buscoDB="Vertebrata"



with open(snakemake.output[0], 'w') as outFile:
# print("---", "title: 'GEP Quick-View'", "author: James Sullivan", "date: 11/09/2021", "geometry: margin=2cm", "classoption: table", "documentclass: extarticle", "urlcolor: blue", "colorlinks: true", "header-includes: |", "  \\rowcolors{2}{gray!10}{gray!25}" , "  ```{=latex}", "---", sep="\n")
	print("", file=outFile)
	print("\\twocolumn", file=outFile)
	print("\\Large", file=outFile)
	print("# Assembly:",params_assemblyID, file=outFile)
	print("", file=outFile)
	print(" - db built with kmer size ", params_merylKmer, "bp", file=outFile)
	print("", file=outFile)
	print("### Genome Statistics", file=outFile)
	print("", file=outFile)
	print(tabulate(key_stats_table, headers='keys',tablefmt="pipe", showindex=False), file=outFile)
	print("", file=outFile)
	print("### Genomescope2 Profile (Linear)", file=outFile)
	print("![](", genomescope_linear, "){ width=37% }", file=outFile)
	print("", file=outFile)
	print("### Genomescope2 Profile (Log)", file=outFile)
	print("![](", genomescope_log, "){ width=37% }", file=outFile)
	print("\\", file=outFile)
	print("\\pagebreak", file=outFile)
	print("\\", file=outFile)
	print("", file=outFile)
	print("\\large", file=outFile)
	print("### QV Score", file=outFile)
	print(tabulate(QV_score_table, headers='keys',tablefmt="pipe", showindex=True), file=outFile)
	print("", file=outFile)
	print("\\large", file=outFile)
	print("### Kmer Completeness", file=outFile)
	print(tabulate(kmer_completeness_table, headers='keys',tablefmt="pipe", showindex=True), file=outFile)
	print("", file=outFile)
	print("\\Large", file=outFile)
	print("### BUSCOv5 (database: ", params_buscoDB, ")", file=outFile)
	print("```", file=outFile)
	for line in lines:
		print(line, file=outFile)
	print("```", file=outFile)
	print("\\", file=outFile)
	print("\\large", file=outFile)
	print("", file=outFile)
	print("### K-mer Multiplicity PRI Only (Stacked)", file=outFile)
	print("![](", kmer_flat, "){ width=44% }", file=outFile)
	print("\\", file=outFile)
	print("", file=outFile)
	print("### K-mer Multiplicity PRI + ALT (Stacked)", file=outFile)
	print("![](", kmer_stacked, "){ width=44% }", file=outFile)
	print("\\", file=outFile)
	print("\\pagebreak", file=outFile)


# ---
# title: "GEP Quick-View"
# date: 11/09/2021
# geometry: a4paper
# classoption: table
# documentclass: extarticle
# urlcolor: blue
# #pdf_document: null
# colorlinks: true
# header-includes: |
#   \rowcolors{2}{gray!10}{gray!25}
#   ```{=latex}
#   \usepackage{pdflscape}
#   \usepackage{longtable}
#   \newcommand{\blandscape}{\begin{landscape}}
#   \newcommand{\elandscape}{\end{landscape}}
#     ```
# output: pdf_document
# ---
#
# print("---", "title: 'GEP Quick-View'", "author: James Sullivan", "date: 11/09/2021", "geometry: margin=2cm", "classoption: table", "documentclass: extarticle", "urlcolor: blue", "colorlinks: true", "header-includes: |", "  \\rowcolors{2}{gray!10}{gray!25}" , "  ```{=latex}", "" sep="\n")
