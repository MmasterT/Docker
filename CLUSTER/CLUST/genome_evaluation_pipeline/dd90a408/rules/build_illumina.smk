



localrules: symlink_UnzippedFastq_R1_illumina,\
 			symlink_UnzippedFastq_R2_illumina, \
			symLink_Trim10xBarcodes_noSequencingAdaptTrim_illumina, \
			symlink_No10xWithSequencingAdaptTrim_illumina, \
			symlink_No10xOrSequencingAdaptTrim_illumina, \
			symlink_Trim10xBarcodes_R2_illumina, \
			multiQC_illumina




def R1_gzipped(wildcards):
	return yesGzip_R1.loc[(wildcards.sample, wildcards.readCounter), "Library_R1"]

def R1_notgzipped(wildcards):
	return noGzip_R1.loc[(wildcards.sample, wildcards.readCounter), "Library_R1"]

def R2_gzipped(wildcards):
	return yesGzip_R2.loc[(wildcards.sample, wildcards.readCounter), "Library_R2"]

def R2_notgzipped(wildcards):
	return noGzip_R2.loc[(wildcards.sample, wildcards.readCounter), "Library_R2"]

rule unzipFastq_R1_illumina:
	input:
		assembly=R1_gzipped,
	output:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.{trim10x}.{trimAdapters}_R1.fastq"),
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{readCounter}.{trim10x}.{trimAdapters}_R1_pigzUnzip.log"),
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	threads:
		resource['unzipFastq_R1_illumina']['threads']
	resources:
		mem_mb=resource['unzipFastq_R1_illumina']['mem_mb'],
		time=resource['unzipFastq_R1_illumina']['time']
	shell:
		"""
		pigz -p {threads} -c -d -k {input.assembly} > {output} 2> {log}
		"""

rule symlink_UnzippedFastq_R1_illumina:
	input:
		assembly=R1_notgzipped,
	output:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.{trim10x}.{trimAdapters}_R1.fastq"),
	container:
		None
	shell:
		"""
		ln -s {input} {output}
		"""

rule unzipFastq_R2_illumina:
	input:
		assembly=R2_gzipped,
	output:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.{trim10x}.{trimAdapters}_R2.fastq"),
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{readCounter}.{trim10x}.{trimAdapters}_R2_pigzUnzip.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	threads:
		resource['unzipFastq_R2_illumina']['threads']
	resources:
		mem_mb=resource['unzipFastq_R2_illumina']['mem_mb'],
		time=resource['unzipFastq_R2_illumina']['time']
	shell:
		"""
		pigz -p {threads} -c -d -k {input.assembly} > {output} 2> {log}
		"""

rule symlink_UnzippedFastq_R2_illumina:
	input:
		assembly=R2_notgzipped,
	output:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.{trim10x}.{trimAdapters}_R2.fastq"),
	container:
		None
	shell:
		"""
		ln -s {input} {output}
		"""


rule trim10xBarcodes_illumina:
	input:
		read1=os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.10xTrimmed.{trimAdapters}_R1.fastq"),
	output:
		read1=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.10xTrimmed.{trimAdapters}_Read1.fastq"),
	threads:
		resource['trim10xBarcodes_illumina']['threads']
	resources:
		mem_mb=resource['trim10xBarcodes_illumina']['mem_mb'],
		time=resource['trim10xBarcodes_illumina']['time']
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{readCounter}.10xTrimmed.10BarcodeRemoval_Trimmomatic.{trimAdapters}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	shell:
		"""
		(trimmomatic SE -threads {threads} {input.read1} {output.read1} HEADCROP:23) &> {log}
		"""


rule symlink_Trim10xBarcodes_R2_illumina:
	input:
		read2=os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.10xTrimmed.{trimAdapters}_R2.fastq")
	output:
		read2=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.10xTrimmed.{trimAdapters}_Read2.fastq")
	container:
		None
	shell:
		"""
		ln -s {input.read2} {output.read2}
		"""



rule symLink_Trim10xBarcodes_noSequencingAdaptTrim_illumina:
	input:
		read1=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.10xTrimmed.notAdaptTrimmed_Read1.fastq"),
		read2=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.10xTrimmed.notAdaptTrimmed_Read2.fastq"),
	output:
		read1=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.10xTrimmed.notAdaptTrimmed_val_1.fq"),
		read2=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.10xTrimmed.notAdaptTrimmed_val_2.fq")
	container:
		None
	shell:
		"""
		ln -s {input.read1} {output.read1}
		ln -s {input.read2} {output.read2}
		"""

rule symlink_No10xOrSequencingAdaptTrim_illumina:
	input:
		read1=os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.not10xTrimmed.notAdaptTrimmed_R1.fastq"),
		read2=os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.not10xTrimmed.notAdaptTrimmed_R2.fastq")
	output:
		read1=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.not10xTrimmed.notAdaptTrimmed_val_1.fq"),
		read2=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.not10xTrimmed.notAdaptTrimmed_val_2.fq")
	container:
		None
	shell:
		"""
		ln -s {input.read1} {output.read1}
		ln -s {input.read2} {output.read2}
		"""

rule symlink_No10xWithSequencingAdaptTrim_illumina:
	input:
		read1=os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.not10xTrimmed.AdaptTrimmed_R1.fastq"),
		read2=os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/{readCounter}.not10xTrimmed.AdaptTrimmed_R2.fastq")
	output:
		read1=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.not10xTrimmed.AdaptTrimmed_Read1.fastq"),
		read2=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.not10xTrimmed.AdaptTrimmed_Read2.fastq")
	container:
		None
	shell:
		"""
		ln -s {input.read1} {output.read1}
		ln -s {input.read2} {output.read2}
		"""

rule trimSequencingAdapters_illumina:
	input:
		read1= os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.AdaptTrimmed_Read1.fastq"),
		read2= os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.AdaptTrimmed_Read2.fastq"),
	params:
		outputDir=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/"),
		r1_prefix="{readCounter}.{trim10x}.AdaptTrimmed",
	output:
		read1=temp(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.AdaptTrimmed_val_1.fq")),
		read2=temp(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.AdaptTrimmed_val_2.fq"))
	threads:
		resource['trimSequencingAdapters_illumina']['threads']
	resources:
		mem_mb=resource['trimSequencingAdapters_illumina']['mem_mb'],
		time=resource['trimSequencingAdapters_illumina']['time'],
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{readCounter}.{trim10x}.AdaptTrimmed_tGalore.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	shell:
		"""
		(trim_galore -j {threads} --basename {params.r1_prefix} --dont_gzip --length 65 -o {params.outputDir} --paired {input.read1} {input.read2}) &> {log}
		"""

rule fastQC_illumina:
	input:
		read1=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.{trimAdapters}_val_1.fq"),
		read2=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.{trimAdapters}_val_2.fq")
	params:
		folder2out=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/fastqc")
	output:
		os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/fastqc/{readCounter}.{trim10x}.{trimAdapters}_val_1_fastqc.html"),
		os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/fastqc/{readCounter}.{trim10x}.{trimAdapters}_val_2_fastqc.html")
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{readCounter}.{trim10x}.{trimAdapters}_fastqc.log")
	threads:
		resource['fastQC_illumina']['threads']
	resources:
		mem_mb=resource['fastQC_illumina']['mem_mb'],
		time=resource['fastQC_illumina']['time'],
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	shell:
		"(fastqc {input} -o {params.folder2out} -t {threads}) &> {log}"

rule multiQC_illumina:
	input:
		read1=lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/fastqc/{readCounter}.{trim10x}.{trimAdapters}_val_1_fastqc.html"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], trim10x=dictSamples[wildcards.sample][1], trimAdapters=dictSamples[wildcards.sample][2]),
		read2=lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/fastqc/{readCounter}.{trim10x}.{trimAdapters}_val_2_fastqc.html"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], trim10x=dictSamples[wildcards.sample][1], trimAdapters=dictSamples[wildcards.sample][2])
	params:
		folder2qc=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/fastqc/"),
		folder2out=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/multiqc/"),
		filename="{sample}.multiqcReport.html"
	output:
		os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/QC/multiqc/{sample}.multiqcReport.html")
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{sample}.multiqc.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	shell:
		"(multiqc {params.folder2qc} -o {params.folder2out} -n {params.filename}) &> {log}"



rule merylCount_R1_illumina:
	input:
		read1= os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.{trimAdapters}_val_1.fq")
	params:
		kmer = "{kmer}",
	threads:
		resource['merylCount_R1_illumina']['threads']
	resources:
		mem_mb=resource['merylCount_R1_illumina']['mem_mb'],
		time=resource['merylCount_R1_illumina']['time'],
	output:
		temp(directory(os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/merylDb/{readCounter}.{trim10x}.{trimAdapters}_R1.{kmer}.meryl")))
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{readCounter}.{trim10x}.{trimAdapters}_meryl_R1.{kmer}.log")
	priority:
		10
	conda:
		os.path.join(workflow.basedir, "envs/MERYL_MERQURY.yaml")
	shell:
		"""
		(meryl count k={params.kmer} threads={threads} {input.read1} output {output}) &> {log}
		"""

rule merylCount_R2_illumina:
	input:
		read2= os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.{trimAdapters}_val_2.fq")
	params:
		kmer = "{kmer}",
	threads:
		resource['merylCount_R2_illumina']['threads']
	resources:
		mem_mb=resource['merylCount_R2_illumina']['mem_mb'],
		time=resource['merylCount_R2_illumina']['time'],
	output:
		temp(directory(os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/merylDb/{readCounter}.{trim10x}.{trimAdapters}_R2.{kmer}.meryl")))
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{readCounter}.{trim10x}.{trimAdapters}_meryl_R2.{kmer}.log")
	priority:
		10
	conda:
		os.path.join(workflow.basedir, "envs/MERYL_MERQURY.yaml")
	shell:
		"""
		(meryl count k={params.kmer} threads={threads} {input.read2} output {output}) &> {log}
		"""


rule merylUnion_illumina:
	input:
		# removeReads1=lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.{trimAdapters}_val_1.fq"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], trim10x=dictSamples[wildcards.sample][1], trimAdapters=dictSamples[wildcards.sample][2]),
		# removeReads2=lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/{readCounter}.{trim10x}.{trimAdapters}_val_2.fq"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], trim10x=dictSamples[wildcards.sample][1], trimAdapters=dictSamples[wildcards.sample][2]),
		read1=lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/merylDb/{readCounter}.{trim10x}.{trimAdapters}_R1.{kmer}.meryl"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], trim10x=dictSamples[wildcards.sample][1], kmer=dictSamples[wildcards.sample][0], trimAdapters=dictSamples[wildcards.sample][2]),
		read2=lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/merylDb/{readCounter}.{trim10x}.{trimAdapters}_R2.{kmer}.meryl"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], trim10x=dictSamples[wildcards.sample][1], kmer=dictSamples[wildcards.sample][0], trimAdapters=dictSamples[wildcards.sample][2])
	params:
		removeReadDIR_trimmed=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_trimReads/"),
		removeReadDIR_unzipped=os.path.join(config['Results'],"0_buildDatabases/{sample}/illuminaReads/temp_unzipFastqs/"),
		kmer = "{kmer}",
		path= os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/merylDb/")
	threads:
		resource['merylUnion_illumina']['threads']
	resources:
		mem_mb=resource['merylUnion_illumina']['mem_mb'],
		time=resource['merylUnion_illumina']['time'],
	output:
		directory(os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/merylDb/complete_illumina.{sample}.{kmer}.meryl")),
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/illuminaReads/logs/{sample}.illumina_meryl.{kmer}.log")
	priority:
		10
	conda:
		os.path.join(workflow.basedir, "envs/MERYL_MERQURY.yaml")
	shell:
		"""
		(meryl union-sum {input.read1} {input.read2} output {output}) &> {log}
		rm -r {params.removeReadDIR_trimmed}
		rm -r {params.removeReadDIR_unzipped}
		"""
