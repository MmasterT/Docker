FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="49e58cca4906e028f46b0ae5f6c1dac3f5dee26aa638b69ea037d04804b01d82"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/R_scripts.yaml
#   prefix: /conda-envs/264f76b767c6fb4787a76b073e3c1550
#   channels:
#    - bioconda
#    - conda-forge
#   dependencies:
#     - python >=3.10.4
#     - r-argparse=2.1.5=r41hc72bb7e_0
#     - r-base >=4.1.3
#     - r-minpack.lm=1.2_2=r41h1463581_0
#     - genomescope2=2.0=py310r41hdfd78af_5
#     - r-dplyr==1.0.9
#     - r-formattable==0.2.1

RUN mkdir -p /conda-envs/264f76b767c6fb4787a76b073e3c1550
COPY envs/AUXILIARY_R_SCRIPTS.yaml /conda-envs/264f76b767c6fb4787a76b073e3c1550/environment.yaml

# Conda environment:
#   source: envs/busco_and_assembly.yaml
#   prefix: /conda-envs/4e57ed29df8b6f849000ab15b5c719f2
#   channels:
#    - conda-forge
#    - bioconda
#    - anaconda
#    - defaults
#   dependencies:
#     - augustus >=3.3
#     - biopython
#     - blast >=2.10.1
#     - fonts-conda-ecosystem
#     - hmmer >=3.1b2
#     - metaeuk
#     - pandas
#     - prodigal
#     - python >=3.3
#     - r-base
#     - r-ggplot2 >=2.2.1
#     - sepp >=4.3.10
#     - wget
#     - busco == 5.3.0
RUN mkdir -p /conda-envs/4e57ed29df8b6f849000ab15b5c719f2
COPY envs/BUSCO.yaml /conda-envs/4e57ed29df8b6f849000ab15b5c719f2/environment.yaml

# Conda environment:
#   source: envs/merylMerq_2.yaml
#   prefix: /conda-envs/bac9eedd8f673dc515402b61559f1582
#   channels:
#    - bioconda
#    - defaults
#    - conda-forge
#   dependencies:
#     - merqury == 1.3
#     - bedtools == 2.29.2
#     - meryl == 1.3
#     - cairo == 1.16
#     - openjdk == 11.0.13
#     - r-argparse == 2.1.5
#     - r-base >=4.0,<4.1.0a0
#     - r-ggplot2=3.3.2=r40hc72bb7e_1
#     - r-scales == 1.2.0
#     - samtools == 1.15.1
#     - xorg-libxrender == 0.9.10
#     - xorg-libxext == 1.3.4
#     - xorg-libxau == 1.0.9

RUN mkdir -p /conda-envs/bac9eedd8f673dc515402b61559f1582
COPY envs/MERYL_MERQURY.yaml /conda-envs/bac9eedd8f673dc515402b61559f1582/environment.yaml

# Conda environment:
#   source: envs/pigz.yaml
#   prefix: /conda-envs/9f8a8ffad4001cd6453433a4ececbe5b
#   channels:
#    - anaconda
#    - bioconda
#    - conda-forge
#   dependencies:
#    - pigz=2.4
#    - trim-galore=0.6.6
#    - pandoc=2.11
#    - multiqc=1.9
#    - python=3.8.5
#    - trimmomatic=0.39

RUN mkdir -p /conda-envs/9f8a8ffad4001cd6453433a4ececbe5b
COPY envs/UNZIP_and_QC.yaml /conda-envs/9f8a8ffad4001cd6453433a4ececbe5b/environment.yaml

# Conda environment:
#   source: envs/pretext.yaml
#   prefix: /conda-envs/d0744b07096b67261c96e38bff814c84
#   channels:
#    - anaconda
#    - bioconda
#    - conda-forge
#   dependencies:
#    - pretextmap == 0.1.9
#    - pretextsnapshot == 0.0.4
#    - bwa-mem2 == 2.2.1
#    - samtools == 1.14

RUN mkdir -p /conda-envs/d0744b07096b67261c96e38bff814c84
COPY envs/HiC_CONTACT_MAPS.yaml /conda-envs/d0744b07096b67261c96e38bff814c84/environment.yaml

# Conda environment:
#   source: envs/python_scripts.yaml
#   prefix: /conda-envs/4d7a7861b78c18f5b15454f62f5a7f70
#   channels:
#    - conda-forge
#    - bioconda
#    - anaconda
#   dependencies:
#    - python>=3.9.10
#    - tabulate=0.8.7
#    - beautifulsoup4=4.9
#    - mamba=0.15.2
#    - pandoc=2.15.*
#    - tectonic
#    - wkhtmltopdf
#    - pandas
#    - numpy
#    - ghostscript

RUN mkdir -p /conda-envs/4d7a7861b78c18f5b15454f62f5a7f70
COPY envs/AUXILIARY_PYTHON_SCRIPTS.yaml /conda-envs/4d7a7861b78c18f5b15454f62f5a7f70/environment.yaml

# Step 2: Generate conda environments
RUN conda config --set remote_read_timeout_secs 600
RUN conda env create --prefix /conda-envs/264f76b767c6fb4787a76b073e3c1550 --file /conda-envs/264f76b767c6fb4787a76b073e3c1550/environment.yaml && \
    conda env create --prefix /conda-envs/4e57ed29df8b6f849000ab15b5c719f2 --file /conda-envs/4e57ed29df8b6f849000ab15b5c719f2/environment.yaml && \
    conda env create --prefix /conda-envs/bac9eedd8f673dc515402b61559f1582 --file /conda-envs/bac9eedd8f673dc515402b61559f1582/environment.yaml && \
    conda env create --prefix /conda-envs/9f8a8ffad4001cd6453433a4ececbe5b --file /conda-envs/9f8a8ffad4001cd6453433a4ececbe5b/environment.yaml && \
    conda env create --prefix /conda-envs/d0744b07096b67261c96e38bff814c84 --file /conda-envs/d0744b07096b67261c96e38bff814c84/environment.yaml && \
    conda env create --prefix /conda-envs/4d7a7861b78c18f5b15454f62f5a7f70 --file /conda-envs/4d7a7861b78c18f5b15454f62f5a7f70/environment.yaml && \
    conda clean --all -y
