


def merylDB(wildcards):
	return samples.loc[(wildcards.asmID), "merylDB"]


localrules: symlink_UnzippedFastq_R1_HiC, 					\
 			symlink_UnzippedFastq_R2_HiC, 					\
			symlink_UnzippedFasta_PRI, 						\
			symlink_UnzippedFasta_ALT, 						\
			KeyResults_GenomescopeProfiles, 				\
			KeyResults, 									\
			Tables_TSV, 									\
			IndividualResults_md, 							\
			PretextMaps_md, 								\
			Table_md, 										\
			Reports_md, 									\
			Reports_pdf,									\
			ColouredTable_html, 							\
			HeatmapTable_html,			\
			BothTables_pdf,						\
			ConcatAll_pdfs

def HiC_R1_gzipped(wildcards):
	return yesGzip_HiC_R1.loc[(wildcards.asmID), "HiC_R1"]

rule unzipFastq_R1_HiC:
	input:
		assembly=HiC_R1_gzipped,
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R1.fastq")),
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/pigzUnzip.HiC.R1.{asmID}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	threads:
		resource['unzipFastq_R1_HiC']['threads']
	resources:
		mem_mb=resource['unzipFastq_R1_HiC']['mem_mb'],
		time=resource['unzipFastq_R1_HiC']['time'],
	shell:
		"""
		pigz -p {threads} -c -d -k {input.assembly} > {output} 2> {log}
		"""

def HiC_R1_unzipped(wildcards):
	return noGzip_HiC_R1.loc[(wildcards.asmID), "HiC_R1"]

rule symlink_UnzippedFastq_R1_HiC:
	input:
		assembly=HiC_R1_unzipped,
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R1.fastq")),
	container:
		None
	shell:
		"""
		ln -sf {input} {output}
		"""


def HiC_R2_gzipped(wildcards):
	return yesGzip_HiC_R2.loc[(wildcards.asmID), "HiC_R2"]

rule unzipFastq_R2_HiC:
	input:
		assembly=HiC_R2_gzipped,
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R2.fastq")),
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/pigzUnzip..HiC.R2.{asmID}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	threads:
		resource['unzipFastq_R2_HiC']['threads']
	resources:
		mem_mb=resource['unzipFastq_R2_HiC']['mem_mb'],
		time=resource['unzipFastq_R2_HiC']['time']
	shell:
		"""
		pigz -p {threads} -c -d -k {input.assembly} > {output} 2> {log}
		"""

def HiC_R2_unzipped(wildcards):
	return noGzip_HiC_R2.loc[(wildcards.asmID), "HiC_R2"]

rule symlink_UnzippedFastq_R2_HiC:
	input:
		assembly=HiC_R2_unzipped,
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R2.fastq")),
	container:
		None
	shell:
		"""
		ln -sf {input} {output}
		"""





rule indexFasta_PRI:
	input:
		assemblyPRI=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
	output:
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.0123"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.amb"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.ann"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.bwt.2bit.64"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.pac"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.fai")
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/indexASM.PRI.{asmID}.log")
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['indexFasta_PRI']['threads']
	resources:
		mem_mb=resource['indexFasta_PRI']['mem_mb'],
		time=resource['indexFasta_PRI']['time']
	shell:
		"""
		(bwa-mem2 index {input.assemblyPRI}) &> {log}
		(samtools faidx {input.assemblyPRI}) &>> {log}
		"""

rule convertFastqTObam_R1_HiC:
	input:
		HiC_R1=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R1.fastq"),
		assembly=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
		indexes= [ os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.0123"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.amb"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.ann"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.bwt.2bit.64"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.pac"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.fai") ]
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R1.bam"))
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['convertFastqTObam_R1_HiC']['threads']
	resources:
		mem_mb=resource['convertFastqTObam_R1_HiC']['mem_mb'],
		time=resource['convertFastqTObam_R1_HiC']['time']
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/fastq2bam.HiC.R1.{asmID}.log")
	shell:
		"""
		(bwa-mem2 mem -t {threads} -B8 {input.assembly} {input.HiC_R1} | samtools view -Sb - > {output})  &> {log}
		"""

rule convertFastqTObam_R2_HiC:
	input:
		HiC_R2=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R2.fastq"),
		assembly=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
		indexes= [ os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.0123"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.amb"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.ann"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.bwt.2bit.64"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.pac"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.fai") ]
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R2.bam"))
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['convertFastqTObam_R2_HiC']['threads']
	resources:
		mem_mb=resource['convertFastqTObam_R2_HiC']['mem_mb'],
		time=resource['convertFastqTObam_R2_HiC']['time']
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/fastq2bam.HiC.R2.{asmID}.log")
	shell:
		"""
		(bwa-mem2 mem -t {threads} -B8 {input.assembly} {input.HiC_R2} | samtools view -Sb - > {output})  &> {log}
		"""

rule filter5PrimeEnd_R1_HiC:
	input:
		HiC_R1_bam=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R1.bam"),
		assembly=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
		indexes= [ os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.0123"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.amb"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.ann"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.bwt.2bit.64"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.pac"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.fai") ]
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R1.FILTERED.bam"))
	params:
		script=os.path.join(workflow.basedir, "scripts/process_HiC/filter_five_end.pl")
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['filter5PrimeEnd_R1_HiC']['threads']
	resources:
		mem_mb=resource['filter5PrimeEnd_R1_HiC']['mem_mb'],
		time=resource['filter5PrimeEnd_R1_HiC']['time']
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/filtered.HiC.R1.{asmID}.log")
	shell:
		"""
 		(samtools view -h {input.HiC_R1_bam}| perl {params.script} | samtools view -@{threads} -Sb - > {output}) &> {log}
		"""

rule filter5PrimeEnd_R2_HiC:
	input:
		HiC_R2_bam=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R2.bam"),
		assembly=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
		indexes= [ os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.0123"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.amb"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.ann"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.bwt.2bit.64"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.pac"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta.fai") ]
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R2.FILTERED.bam"))
	params:
		script=os.path.join(workflow.basedir, "scripts/process_HiC/filter_five_end.pl")
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['filter5PrimeEnd_R2_HiC']['threads']
	resources:
		mem_mb=resource['filter5PrimeEnd_R2_HiC']['mem_mb'],
		time=resource['filter5PrimeEnd_R2_HiC']['time']
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/filtered.HiC.R1.{asmID}.log")
	shell:
		"""
 		(samtools view -h {input.HiC_R2_bam}| perl {params.script} | samtools view -@{threads} -Sb - > {output}) &> {log}
		"""


rule pairAndCombineFiltered_HiC:
	input:
		R1=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R1.FILTERED.bam"),
		R2=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.R2.FILTERED.bam")
	output:
		temp(os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.COMBINED.FILTERED.bam"))
	params:
		script=os.path.join(workflow.basedir, "scripts/process_HiC/two_read_bam_combiner.pl")
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['pairAndCombineFiltered_HiC']['threads']
	resources:
		mem_mb=resource['pairAndCombineFiltered_HiC']['mem_mb'],
		time=resource['pairAndCombineFiltered_HiC']['time']
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/combine.filtered.HiC.{asmID}.log")
	shell:
		"""
		(perl {params.script} {input.R1} {input.R2} | samtools view -@{threads} -Sb > {output}) &>{log}
		"""

rule pretextMap:
	input:
		HiC_alignment=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.COMBINED.FILTERED.bam")
	output:
		pretextFile=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.COMBINED.FILTERED.pretext")
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['pretextMap']['threads']
	resources:
		mem_mb=resource['pretextMap']['mem_mb'],
		time=resource['pretextMap']['time']
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/PretextMap.{asmID}.log")
	shell:
		"""
		(samtools view -h {input.HiC_alignment} | PretextMap -o {output.pretextFile} --sortby length --mapq 10) &> {log}
		"""

rule pretextSnapshot:
	input:
		pretextFile=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.COMBINED.FILTERED.pretext")
	output:
		pretextSnapshotFULL=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.COMBINED.FILTERED_FullMap.png"),
	params:
		outDirectory=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/")
	conda:
		os.path.join(workflow.basedir, "envs/HiC_CONTACT_MAPS.yaml")
	threads:
		resource['pretextSnapshot']['threads']
	resources:
		mem_mb=resource['pretextSnapshot']['mem_mb'],
		time=resource['pretextSnapshot']['time']
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/PretextSnapshot.{asmID}.log")
	shell:
		"""
		(PretextSnapshot -m {input.pretextFile} -o {params.outDirectory}) &> {log}
		"""

#######################################################################################################################################

def PRI_asm_gzipped(wildcards):
	return yesGzip_PRI.loc[(wildcards.asmID), "PRI_asm"]


rule unzipFasta_PRI:
	input:
		assembly=PRI_asm_gzipped,
	output:
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/pigzUnzip.PRI.{asmID}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	threads:
		resource['unzipFasta_PRI']['threads']
	resources:
		mem_mb=resource['unzipFasta_PRI']['mem_mb'],
		time=resource['unzipFasta_PRI']['time']
	shell:
		"""
		pigz -p {threads} -c -d -k {input.assembly} > {output} 2> {log}
		"""

def PRI_asm_unzipped(wildcards):
	return noGzip_PRI.loc[(wildcards.asmID), "PRI_asm"]

rule symlink_UnzippedFasta_PRI:
	input:
		assembly=PRI_asm_unzipped,
	output:
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
	container:
		None
	shell:
		"""
		ln -fs {input} {output}
		"""

################

def ALT_asm_gzipped(wildcards):
	return yesGzip_ALT.loc[(wildcards.asmID), "ALT_asm"]

rule unzipFasta_ALT:
	input:
		assembly=ALT_asm_gzipped,
	output:
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.ALT.fasta"),
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/pigzUnzip.ALT.{asmID}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	threads:
		resource['unzipFasta_ALT']['threads']
	resources:
		mem_mb=resource['unzipFasta_ALT']['mem_mb'],
		time=resource['unzipFasta_ALT']['time']
	shell:
		"""
		pigz -p {threads} -c -d -k {input.assembly} > {output} 2> {log}
		"""


def ALT_asm_unzipped(wildcards):
	return noGzip_ALT.loc[(wildcards.asmID), "ALT_asm"]

rule symlink_UnzippedFasta_ALT:
	input:
		assembly=ALT_asm_unzipped,
	output:
		os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.ALT.fasta"),
	container:
		None
	shell:
		"""
		ln -s {input} {output}
		"""


####################


def altFile(wildcards):
	if samples.loc[(wildcards.asmID), "ALT_present"] == True:
		return os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.ALT.fasta")
	else:
		return os.path.join(workflow.basedir, "scripts/standIN_files/ALT_missing.fasta")

# rule symlinkMerylDB:
# 	input:
# 		merylDB_provided=merylDB
# 	output:
# 		symlink_merylDB=directory(os.path.join(config['Results'], "1_evaluation/{asmID}/04_merquryQVandKAT/merylDB_providedFor_{asmID}.meryl"))
# 	container:
# 		None
# 	shell:
# 		"""
# 		ln -s {input.merylDB_provided} {output.symlink_merylDB}
# 		"""

rule merqury:
	input:
		assemblyPRI=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
		assemblyALT=altFile,
		merylDB_provided=merylDB,
		# merylDB=os.path.join(config['Results'], "1_evaluation/{asmID}/04_merquryQVandKAT/merylDB_providedFor_{asmID}.meryl")
	params:
		outFile= "{asmID}" + "_merqOutput",
		outPath= os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA"),
		symlink_merylDB=directory(os.path.join(config['Results'], "1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/merylDB_providedFor_{asmID}.meryl"))
#		symlink_merylDB=directory(os.path.join(config['Results'], "1_evaluation/{asmID}/04_merquryQVandKAT/merylDB_providedFor_{asmID}.meryl"))
	threads:
		resource['merqury']['threads']
	resources:
		mem_mb=resource['merqury']['mem_mb'],
		time=resource['merqury']['time']
	output:
		os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.qv"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.completeness.stats"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.{asmID}.PRI.spectra-cn.st.png"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.{asmID}.PRI.spectra-cn.fl.png"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.spectra-cn.st.png"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.spectra-cn.fl.png"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/merylDB_providedFor_{asmID}.hist")
	log:
		os.path.join(config['Results'],"1_evaluation/{asmID}/logs/merqury.{asmID}.log")
	priority:
		3
	conda:
		os.path.join(workflow.basedir, "envs/MERYL_MERQURY.yaml")
	shell:
		"""
		ln -fs {input.merylDB_provided} {params.symlink_merylDB}
		cd {params.outPath}
		export OMP_NUM_THREADS={threads}
		(merqury.sh {params.symlink_merylDB} {input.assemblyPRI} {input.assemblyALT} {params.outFile}) &> {log}
		"""



#######################################################################################################################################



rule busco5:
	input:
		assembly=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
		lineage = os.path.join(workflow.basedir, "buscoLineage/" + buscoDataBaseName + "_odb10"),
	params:
		assemblyName = "{asmID}",
		chngDir = os.path.join(config['Results'], "1_evaluation/{asmID}/BUSCOs")
	threads:
		resource['busco5']['threads']
	resources:
		mem_mb=resource['busco5']['mem_mb'],
		time=resource['busco5']['time']
	output:
		os.path.join(config['Results'], "1_evaluation/{asmID}/BUSCOs/{asmID}/short_summary.specific." + buscoDataBaseName + "_odb10.{asmID}.txt"),
	conda:
		os.path.join(workflow.basedir, "envs/BUSCO.yaml")
	log:
		os.path.join(config['Results'], "1_evaluation/{asmID}/logs/busco5.{asmID}.log")
	priority:
		20
	shell:
		"""
		cd {params.chngDir}
		(busco -m genome --offline --in {input.assembly} -o {params.assemblyName} -l {input.lineage} -c {threads} -f --limit 5) &> {log}
		"""





rule GenomeScope2Profiles:
	input:
		hist=os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/merylDB_providedFor_{asmID}.hist"),
	params:
		outFolder=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/"),
		name="{asmID}_k{kmer}",
		kmer="{kmer}",
		cpHist=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/merylDB_providedFor_{asmID}_10000.hist")
	output:
		summary=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_summary.txt"),
		logPlot=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_log_plot.png"),
		linearPlot=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_linear_plot.png"),
		estimatedSize=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_sizeEst.txt")
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_R_SCRIPTS.yaml")
	log:
		os.path.join(config['Results'],"1_evaluation/{asmID}/logs/genomescopeProfiles.{asmID}.{kmer}.log")
	threads:
		resource['GenomeScope2Profiles']['threads']
	resources:
		mem_mb=resource['GenomeScope2Profiles']['mem_mb'],
		time=resource['GenomeScope2Profiles']['time']
	shell:
		"""
		head -n 10000 {input.hist} > {params.cpHist}
		(genomescope2 -i {params.cpHist} -o {params.outFolder} -k {params.kmer} -n {params.name}) &> {log}
		grep "Genome Haploid Length" {output.summary} | awk {{'print $6'}} | sed 's/,//g'> {output.estimatedSize}
		"""

rule assemblyStats:
	input:
		assembly=os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_FASTAS/{asmID}.PRI.fasta"),
		estGenome=lambda wildcards: expand(os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_sizeEst.txt"), asmID=wildcards.asmID, kmer=dictSamples[wildcards.asmID][4]),
	output:
		scaffStats=os.path.join(config['Results'],"1_evaluation/{asmID}/ASSEMBLY_STATISTICS/{asmID}_scaffold_stats.tsv"),
		contStats=os.path.join(config['Results'],"1_evaluation/{asmID}/ASSEMBLY_STATISTICS/{asmID}_contig_stats.tsv")
	params:
		script = os.path.join(workflow.basedir, "scripts/assembly_stats/stats.py"),
		path = os.path.join(config['Results'], "1_evaluation/{asmID}/ASSEMBLY_STATISTICS/"),
		filename="{asmID}",
		given_size=lambda wildcards: expand("{genomeSize}", genomeSize=dictSamples[wildcards.asmID][5])
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_PYTHON_SCRIPTS.yaml")
	log:
		os.path.join(config['Results'],"1_evaluation/{asmID}/logs/assemblyStats.{asmID}.log")
	threads:
		resource['assemblyStats']['threads']
	resources:
		mem_mb=resource['assemblyStats']['mem_mb'],
		time=resource['assemblyStats']['time']
	shell:
		"""
		(python {params.script} {input.assembly} {input.estGenome} {params.filename} {params.given_size} {output.scaffStats} {output.contStats}) &> {log}
		"""




rule KeyResults_GenomescopeProfiles:
	input:
		gscopeSum=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_summary.txt"),
		gscopeLog=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_log_plot.png"),
		gscopeLin=os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_linear_plot.png")
	output:
		gscopeSum=os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/{asmID}_k{kmer}_summary.txt"),
		gscopeLog=os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/{asmID}_k{kmer}_log_plot.png"),
		gscopeLin=os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/{asmID}_k{kmer}_linear_plot.png")
	container:
	 	None
	shell:
		"""
		cp {input.gscopeSum} {output.gscopeSum}
		cp {input.gscopeLog} {output.gscopeLog}
		cp {input.gscopeLin} {output.gscopeLin}
		"""
#



rule KeyResults:
	input:
		gscopeSum=lambda wildcards: expand(os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_summary.txt"), asmID=wildcards.asmID, kmer=dictSamples[wildcards.asmID][4]),
		# gscopeLog=lambda wildcards: expand(os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_log_plot.png"), asmID=wildcards.asmID, kmer=dictSamples[wildcards.asmID][4]),
		# gscopeLin=lambda wildcards: expand(os.path.join(config['Results'], "1_evaluation/{asmID}/GENOMESCOPE_PROFILES/{asmID}_k{kmer}_linear_plot.png"), asmID=wildcards.asmID, kmer=dictSamples[wildcards.asmID][4]),
		busco=os.path.join(config['Results'], "1_evaluation/{asmID}/BUSCOs/{asmID}/short_summary.specific." + buscoDataBaseName + "_odb10.{asmID}.txt"),
		qv=os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.qv"),
		completeness=os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.completeness.stats"),
		spectraStacked_PRI=os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.{asmID}.PRI.spectra-cn.st.png"),
		spectraFlat_PRI=os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.{asmID}.PRI.spectra-cn.fl.png"),
		spectraStacked_both=os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.spectra-cn.st.png"),
		spectraFlat_both=os.path.join(config['Results'],"1_evaluation/{asmID}/QV.KMER-COMPLETENESS.CN-SPECTRA/{asmID}_merqOutput.spectra-cn.fl.png"),
		scaffStats=os.path.join(config['Results'],"1_evaluation/{asmID}/ASSEMBLY_STATISTICS/{asmID}_scaffold_stats.tsv"),
		contStats=os.path.join(config['Results'],"1_evaluation/{asmID}/ASSEMBLY_STATISTICS/{asmID}_contig_stats.tsv"),
	params:
		asm_level=lambda wildcards: expand("{ASM_LEVEL}", ASM_LEVEL=dictSamples[wildcards.asmID][0]),
		resultsPath=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS"),
		keyValues2=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/results_withoutRowNames_2.tsv"),
		keyValues=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/results_withoutRowNames.tsv"),
		rowNames=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/rowNames.tsv")
	output:
		keyValuesWithRowNames=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_aggregatedSTATS.tsv"),
		busco=os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/short_summary.specific." + buscoDataBaseName + "_odb10.{asmID}.txt"),
		buscoScores=os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/only_buscoScores_{asmID}.txt"),
		qv=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.qv"),
		completeness=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.completeness.stats"),
		spectraStacked_PRI=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.{asmID}.PRI.spectra-cn.st.png"),
		spectraFlat_PRI=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.{asmID}.PRI.spectra-cn.fl.png"),
		spectraStacked_both=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.spectra-cn.st.png"),
		spectraFlat_both=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.spectra-cn.fl.png"),
		scaffStats=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_scaffold_stats.tsv"),
		contStats=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_contig_stats.tsv"),
	container:
		None
	shell:
		"""
		head -n 15 {input.busco} | tail -n 7 | sed "s/^['\t']*//" > {output.buscoScores}
		cp {input.busco} {output.busco}
		cp {input.completeness} {output.completeness}
		cp {input.qv} {output.qv}
		cp {input.spectraStacked_PRI} {output.spectraStacked_PRI}
		cp {input.spectraFlat_PRI} {output.spectraFlat_PRI}
		cp {input.spectraFlat_both} {output.spectraFlat_both}
		cp {input.spectraStacked_both} {output.spectraStacked_both}
		cp {input.scaffStats} {output.scaffStats}
		cp {input.contStats} {output.contStats}
		echo "{wildcards.asmID}" >> {params.keyValues}
		echo "{params.asm_level}" >> {params.keyValues}
		echo "$(grep 'total_bps' {input.scaffStats} | awk {{'print $3'}})" >> {params.keyValues}
		echo "$(grep 'Genome Haploid Length' {input.gscopeSum} | awk {{'print $6'}})" >> {params.keyValues}
		echo "$(grep 'Heterozygous (ab)' {input.gscopeSum} | awk {{'print $4'}})" >> {params.keyValues}
		echo "$(grep 'gc_content' {input.scaffStats} | awk {{'print $3'}})" >> {params.keyValues}
		echo "$(grep 'sequence_count' {input.scaffStats} | awk {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'sequence_count' {input.contStats} | awk {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'number_of_gaps' {input.contStats} | awk {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'longest (bp)' {input.scaffStats} | awk {{'print $3'}})" >> {params.keyValues}
		echo "$(grep 'N50' {input.scaffStats} | awk {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'NG50' {input.scaffStats} | awk {{'print $4'}})" >> {params.keyValues}
		echo "$(grep 'N95' {input.scaffStats} | awk {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'NG95' {input.scaffStats} | awk {{'print $4'}})" >> {params.keyValues}
		echo "$(grep 'longest (bp)' {input.contStats} | awk {{'print $3'}})" >> {params.keyValues}
		echo "$(grep 'N50' {input.contStats} | awk {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'NG50' {input.contStats} | awk {{'print $4'}})" >> {params.keyValues}
		echo "$(grep 'N95' {input.contStats} | awk {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'NG95' {input.contStats} | awk {{'print $4'}})" >> {params.keyValues}
		echo "$(awk {{'print $4'}} {input.qv})" | head -n 1 >> {params.keyValues}
		echo "$(awk {{'print $5'}} {input.completeness})" | head -n 1 >> {params.keyValues}
		echo "$(grep 'C:' {input.busco} | awk -F'[:\\[,]' {{'print $2'}})" >> {params.keyValues}
		echo "$(grep 'C:' {input.busco} | awk -F'[:\\[,]' {{'print $4'}})" >> {params.keyValues}
		echo -e "ASM_ID\nASM_LEVEL\nBases\nEst_Size\nHet_%\nGC_%\nScaff\nCont\nGaps\nLongest_Scaff\nScaff_N50\nScaff_NG50\nScaff_N95\nScaff_NG95\nLongest_Cont\nCont_N50\nCont_NG50\nCont_N95\nCont_NG95\nQV\nCompleteness\nComp_BUSCOs_%\nComp_Single_BUSCOs_%" > {params.rowNames}
		paste -d'\t' {params.rowNames} {params.keyValues} > {output.keyValuesWithRowNames}
		rm {params.rowNames} {params.keyValues}
		"""

rule Tables_TSV:
	input:
		allResults=expand(os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_aggregatedSTATS.tsv"), asmID=list(dictSamples.keys())),
		sampleSheet= config['samplesTSV'],
		config=os.path.join(workflow.basedir, "configuration/config.yaml")
	output:
		results=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS.tsv"),
		newSampleSheet=os.path.join(config['Results'],"1_evaluation/finalResults/saved_configuration/savedSampleSheet.tsv"),
		newConfigFile=os.path.join(config['Results'],"1_evaluation/finalResults/saved_configuration/savedConfig.yaml")
	container:
		None
	shell:
		"""
		cp {input.sampleSheet} {output.newSampleSheet}
		cp {input.config} {output.newConfigFile}
		echo -e "ASM_ID\nASM_LEVEL\nBases_Mb\nEst_Size_Mb\nHet_%\nGC_%\nScaff\nCont\nGaps\nLongest_Scaff_Mb\nScaff_N50_Mb\nScaff_NG50_Mb\nScaff_N95_Mb\nScaff_NG95_Mb\nLongest_Cont_Mb\nCont_N50_Mb\nCont_NG50_Mb\nCont_N95_Mb\nCont_NG95_Mb\nQV\nCompleteness\nComp_BUSCOs_%\nComp_Single_BUSCOs_%" | \
		paste -d'\t' - {input.allResults} | \
		awk -F'\t' '{{s="";for (i=1;i<=NF;i+=2) {{s=s?s FS $i:$i}} print s}}' | \
		column -t > {output.results}
		sed -i 's/,//g' {output.results}
		"""

rule IndividualResults_md:
	input:
		os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_aggregatedSTATS.tsv"),
		lambda wildcards: expand(os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/{asmID}_k{kmer}_log_plot.png"),asmID=wildcards.asmID, kmer=dictSamples[wildcards.asmID][4]),
		lambda wildcards: expand(os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/{asmID}_k{kmer}_linear_plot.png"),asmID=wildcards.asmID, kmer=dictSamples[wildcards.asmID][4]),
		os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.qv"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.completeness.stats"),
		os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/only_buscoScores_{asmID}.txt"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.{asmID}.PRI.spectra-cn.st.png"),
		os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_merqOutput.spectra-cn.st.png")
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_PYTHON_SCRIPTS.yaml")
	params:
		"{asmID}",
		"{kmer}",
		config['busco5Lineage']
	output:
		os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_k{kmer}_markdownForReport.md")
	script:
		os.path.join(workflow.basedir, "scripts/report/makePDF_indivMD.py")



rule PretextMaps_md:
	input:
		PretextMap=os.path.join(config['Results'], "1_evaluation/{asmID}/HiC_MAPS/{asmID}.HiC.COMBINED.FILTERED_FullMap.png")
	output:
		pretextCP2keyResults=os.path.join(config['Results'], "1_evaluation/{asmID}/KEY_RESULTS/{asmID}.pretext_hiC_FullMap.png"),
		IndividualPretextMD=os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_PRETEXT.md")
	params:
		assemblyName='{asmID}'
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_PYTHON_SCRIPTS.yaml")
	script:
		os.path.join(workflow.basedir, "scripts/report/pretextMapsToMarkdown.py")

rule Table_md:
	input:
		results=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS.tsv")
	output:
		os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS.md"),
		os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_roundedMB.tsv")
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_PYTHON_SCRIPTS.yaml")
	script:
		os.path.join(workflow.basedir, "scripts/report/addFullTableForReport.py")

rule ColouredTable_html:
	input:
		os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_roundedMB.tsv")
	params:
		os.path.join(workflow.basedir, "scripts/compare_results/colouredHeatmap_legend.tsv")
	output:
		os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_COLOURED.html"),
	# threads:
	# 	resource['fullTable_heatmap_external_create']['threads']
	# resources:
	# 	mem_mb=resource['fullTable_heatmap_external_create']['mem_mb'],
	# 	time=resource['fullTable_heatmap_external_create']['time']
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_R_SCRIPTS.yaml")
	script:
		os.path.join(workflow.basedir, "scripts/compare_results/fullTable_heatmap_external.R")

rule HeatmapTable_html:
	input:
		os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_roundedMB.tsv")
	params:
		os.path.join(workflow.basedir, "scripts/compare_results/internalComparison_legend.tsv")
	output:
		os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_GRADIENT.html")
	# threads:
	# 	resource['fullTable_heatmap_external_create']['threads']
	# resources:
	# 	mem_mb=resource['fullTable_heatmap_external_create']['mem_mb'],
	# 	time=resource['fullTable_heatmap_external_create']['time']
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_R_SCRIPTS.yaml")
	script:
		os.path.join(workflow.basedir, "scripts/compare_results/fullTable_heatmap_internalComparison.R")

rule BothTables_pdf:
	input:
		coloured=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_COLOURED.html"),
		gradient=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_GRADIENT.html")
	params:
		css=os.path.join(workflow.basedir, "scripts/report/tableOnSamePage.css")
	log:
		os.path.join(config['Results'], "1_evaluation/logs/ComparisonTables_createPDF.log")
	output:
		coloured=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_COLOURED.pdf"),
		gradient=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_GRADIENT.pdf")
	# threads:
	# 	resource['fullTable_heatmap_external_createPDF']['threads']
	# resources:
	# 	mem_mb=resource['fullTable_heatmap_external_createPDF']['mem_mb'],
	# 	time=resource['fullTable_heatmap_external_createPDF']['time']
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_PYTHON_SCRIPTS.yaml")
	shell:
		"""
		(pandoc -V papersize:a3 -V margin-top=1.5cm -V margin-left=1.5cm -V margin-right=0 -V margin-bottom=0 -c {params.css} -o {output.coloured} --pdf-engine=wkhtmltopdf {input.coloured}) &>> {log}
		(pandoc -V papersize:b3 -V margin-top=1.5cm -V margin-left=1.5cm -V margin-right=0 -V margin-bottom=0 -c {params.css} \
		-o {output.gradient} --pdf-engine=wkhtmltopdf --pdf-engine-opt="-O" --pdf-engine-opt="Landscape" {input.gradient}) &>> {log}
		"""

rule Reports_md:
	input:
		indivMD=[expand(os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_k{kmer}_markdownForReport.md"), asmID=key, kmer=value5) for key, [value1, value2, value3, value4, value5, value6, value7, value8, value9, value10, value11, value12, value13, value14, value15, value16] in dictSamples.items()],
		landingPage=os.path.join(workflow.basedir, "scripts/report/reportLandingPage.md"),
		# landingPageTABLE=os.path.join(workflow.basedir, "scripts/reportTABLELandingPage.md"),
		endTableMD=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS.md"),
		IndividualPretextMD=[expand(os.path.join(config['Results'],"1_evaluation/{asmID}/KEY_RESULTS/{asmID}_PRETEXT.md"), asmID=list(testDictPRETEXT.keys()))]
	output:
		FullMarkdown=os.path.join(config['Results'],"1_evaluation/finalResults/individual_reports/ALL_individual_REPORTS.md"),
		# endTableMD=os.path.join(config['Results'],"1_evaluation/finalResults/FullTableMarkdown_wPreamble.md")
	container:
		None
	shell:
		"""
		cat {input.landingPage} {input.indivMD} {input.IndividualPretextMD} >> {output.FullMarkdown}
		"""
#		cat {input.landingPageTABLE} {input.endTableMD} >> {output.endTableMD}


rule Reports_pdf:
	input:
		md_report=os.path.join(config['Results'],"1_evaluation/finalResults/individual_reports/ALL_individual_REPORTS.md"),
#		md_comparison_table=os.path.join(config['Results'],"1_evaluation/finalResults/FullTableMarkdown_wPreamble.md")
	output:
		pdf_report=os.path.join(config['Results'],"1_evaluation/finalResults/individual_reports/ALL_individual_REPORTS.pdf"),
#		pdf_comparison_table=os.path.join(config['Results'],"1_evaluation/finalResults/FULL_TABLE_PDF.pdf")
	log:
		os.path.join(config['Results'], "1_evaluation/logs/ReportWithoutComparisonTables_createPDF.log")
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_PYTHON_SCRIPTS.yaml")
	# threads:
	# 	resource['ReportWithoutComparisonTables_createPDF']['threads']
	# resources:
	# 	mem_mb=resource['ReportWithoutComparisonTables_createPDF']['mem_mb'],
	# 	time=resource['ReportWithoutComparisonTables_createPDF']['time']
	shell:
		"""
		(pandoc -o {output.pdf_report} {input.md_report} --pdf-engine=tectonic) &>> {log}
		"""
	#	(pandoc -o {output.pdf_comparison_table} {input.md_comparison_table} --pdf-engine=tectonic) &>> {log}

rule ConcatAll_pdfs:
	input:
		pdf_report=os.path.join(config['Results'],"1_evaluation/finalResults/individual_reports/ALL_individual_REPORTS.pdf"),
		coloured=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_COLOURED.pdf"),
		gradient=os.path.join(config['Results'],"1_evaluation/finalResults/tables/TABLE_OF_RESULTS_GRADIENT.pdf")
#		md_comparison_table=os.path.join(config['Results'],"1_evaluation/finalResults/FullTableMarkdown_wPreamble.md")
	output:
		pdf_report=os.path.join(config['Results'],"1_evaluation/finalResults/GEP_FINAL_REPORT.pdf"),
#		pdf_comparison_table=os.path.join(config['Results'],"1_evaluation/finalResults/FULL_TABLE_PDF.pdf")
	log:
		os.path.join(config['Results'], "1_evaluation/logs/ConcatAll_pdfs.log")
	conda:
		os.path.join(workflow.basedir, "envs/AUXILIARY_PYTHON_SCRIPTS.yaml")
	# threads:
	# 	resource['makePDF']['threads']
	# resources:
	# 	mem_mb=resource['makePDF']['mem_mb'],
	# 	time=resource['makePDF']['time']
	shell:
		"""
		(gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output.pdf_report} {input.pdf_report} {input.gradient} {input.coloured}) &>> {log}
		"""
