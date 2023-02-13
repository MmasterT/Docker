
localrules: symlink_UnzippedFastq_hifi, \
			symlink_noSMRTBellAdaptTrim_hifi, \
			multiQC_hifi



# def fq_to_trimSMRTbell(wildcards):
#   return trimSMRTbell.loc[(wildcards.sample, wildcards.readCounter), "hifi_reads"]
#
#
# def fq_to_notTrimSMRTbell(wildcards):
#   return notrimSMRTbell.loc[(wildcards.sample, wildcards.readCounter), "hifi_reads"]


def hifi_gzipped(wildcards):
	return yesGzip_hifi.loc[(wildcards.sample, wildcards.readCounter), "hifi_reads"]


def hifi_notgzipped(wildcards):
	return noGzip_hifi.loc[(wildcards.sample, wildcards.readCounter), "hifi_reads"]



rule unzipFastq_hifi:
	input:
		fastq=hifi_gzipped,
	output:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/temp_unzipFastqs/{readCounter}.{smrtornot}.fastq"),
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/logs/pigzUnzip.{readCounter}.{smrtornot}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	threads:
		resource['unzipFastq_hifi']['threads']
	resources:
		mem_mb=resource['unzipFastq_hifi']['mem_mb'],
		time=resource['unzipFastq_hifi']['time'],
	shell:
		"""
		pigz -p {threads} -c -d -k {input.fastq} > {output} 2> {log}
		"""

rule symlink_UnzippedFastq_hifi:
	input:
		fastq=hifi_notgzipped,
	output:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/temp_unzipFastqs/{readCounter}.{smrtornot}.fastq"),
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/logs/pigzUnzip.{readCounter}.{smrtornot}.log")
	container:
		None
	shell:
		"""
		ln -s {input.fastq} {output}
		echo "{input.fastq} not gzipped. Symlink created in place of expected decompressed file." > {log}
		"""

rule trimSMRTBellAdapters_hifi:
	input:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/temp_unzipFastqs/{readCounter}.smrtTrimmed.fastq"),
	output:
		outputFile=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/temp_trimReads/{readCounter}.smrtTrimmed.fastq")
	threads:
		resource['trimSMRTBellAdapters_hifi']['threads']
	resources:
		mem_mb=resource['trimSMRTBellAdapters_hifi']['mem_mb'],
		time=resource['trimSMRTBellAdapters_hifi']['time'],
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/logs/trimSMRTbell.{readCounter}.log")
	priority:
		15
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	shell:
		"""
		(cutadapt -j {threads} -o {output.outputFile} {input} -b AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA --overlap 35 -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT --overlap 45 --revcomp -e 0.1 --discard-trimmed) &> {log}
		"""

rule symlink_noSMRTBellAdaptTrim_hifi:
	input:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/temp_unzipFastqs/{readCounter}.notsmrtTrimmed.fastq"),
	output:
		outputFile=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/temp_trimReads/{readCounter}.notsmrtTrimmed.fastq")
	container:
		None
	shell:
		"""
		ln -s {input} {output.outputFile}
		"""

rule fastQC_hifi:
	input:
		os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/temp_trimReads/{readCounter}.{smrtornot}.fastq")
	params:
		folder2out=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/QC/fastqc")
	output:
		os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/QC/fastqc/{readCounter}.{smrtornot}_fastqc.html")
	threads:
		resource['fastQC_hifi']['threads']
	resources:
		mem_mb=resource['fastQC_hifi']['mem_mb'],
		time=resource['fastQC_hifi']['time'],
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/logs/fastQC.hifi.{readCounter}.{smrtornot}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	shell:
		"""
		(fastqc {input} -o {params.folder2out} -t {threads}) &> {log}
		"""

rule multiQC_hifi:
	input:
		lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/QC/fastqc/{readCounter}.{smrtornot}_fastqc.html"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], smrtornot=dictSamples[wildcards.sample][1])
	params:
		folder2qc=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/QC/fastqc/"),
		folder2OUT=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/QC/multiqc/"),
		filename="{sample}.multiqcReport.html"
	output:
		os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/QC/multiqc/{sample}.multiqcReport.html")
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/logs/multiQC.{sample}.log")
	conda:
		os.path.join(workflow.basedir, "envs/UNZIP_and_QC.yaml")
	shell:
		"(multiqc {params.folder2qc} -o {params.folder2OUT} -n {params.filename}) &> {log}"


rule merylCount_hifi:
	input:
		reads=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/temp_trimReads/{readCounter}.{smrtornot}.fastq")
	params:
		kmer = "{kmer}"
	threads:
		resource['merylCount_hifi']['threads']
	resources:
		mem_mb=resource['merylCount_hifi']['mem_mb'],
		time=resource['merylCount_hifi']['time'],
	output:
		temp(directory(os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/merylDb/" + "{readCounter}" + "_hifi_dB.{smrtornot}.{kmer}.meryl"))),
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/logs/meryl_hifi_count.{readCounter}.{kmer}.{smrtornot}.log")
	priority:
		10
	conda:
		os.path.join(workflow.basedir, "envs/MERYL_MERQURY.yaml")
	shell:
		"""
		(meryl count k={params.kmer} threads={threads} {input.reads} output {output}) &> {log}
		"""

rule merylUnion_hifi:
	input:
		lambda wildcards: expand(os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/merylDb/{readCounter}_hifi_dB.{smrtornot}.{kmer}.meryl/"), sample=wildcards.sample, readCounter=dictReadCounter[wildcards.sample], kmer=dictSamples[wildcards.sample][0], smrtornot=dictSamples[wildcards.sample][1])
	params:
		kmer = "{kmer}",
		removeReadDIR_trimmed=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/temp_trimReads/"),
		removeReadDIR_unzipped=os.path.join(config['Results'],"0_buildDatabases/{sample}/hifiReads/temp_unzipFastqs/")
	threads:
		resource['merylUnion_hifi']['threads']
	resources:
		mem_mb=resource['merylUnion_hifi']['mem_mb'],
		time=resource['merylUnion_hifi']['time'],
	output:
		directory(os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/merylDb/complete_hifi.{sample}.{kmer}.meryl")),
	log:
		os.path.join(config['Results'], "0_buildDatabases/{sample}/hifiReads/logs/meryl_hifi_combine.{sample}.{kmer}.log")
	priority:
		10
	conda:
		os.path.join(workflow.basedir, "envs/MERYL_MERQURY.yaml")
	shell:
		"""
		(meryl union-sum {input} output {output}) &> {log}
		rm -r {params.removeReadDIR_trimmed}
		rm -r {params.removeReadDIR_unzipped}
		"""
