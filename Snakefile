# LongCallR - phase long-read RNA-seq for allele-specific tests

# Jack Humphrey 2026
import pandas as pd
import os

lcr_threads = 16

# for working on minerva
#shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge;")

samtools = config["samtools_path"]
longcallr_path = config["longcallr_path"]
longcallr = os.path.join(longcallr_path, "target/release/longcallR")
asj_to_bed = os.path.join(longcallr_path, "allele_specific/asj_to_bed.py")

rediportal = "/sc/arion/projects/ad-omics/data/references/editing/TABLE1_hg38_v3.txt.gz"
region_string = ""

ref_genome = config["ref_genome"] + ".fa"
ref_gtf = config["ref_gtf"]
metadata = config["metadata"]
data_code = config["data_code"]
outFolder = config["out_folder"]
prep = config["prep"]
if "phased_vcf" in config.keys():
    phased_vcf = config["phased_vcf"] # make optional later
    phased_vcf_string = "--input-vcf " + phased_vcf + " --direct-haplotag"
else:
    phased_vcf_string = ""
#gwas = "/sc/arion/projects/ad-omics/data/references/GWAS/Bellenguez_AD/Bellenguez_2021.processed.tsv.gz"
# file either TSV or XLSX
if ".tsv" in metadata:
    meta_df = pd.read_csv(metadata, sep = '\t')
if ".xlsx" in metadata:
    meta_df = pd.read_excel(metadata)

# turn off asediting for now
allele_tests = ["ase","asj"]

samples = meta_df['sample']

# for testing
#samples = "16-078_2_MFG"
#region_string = "-r chr16:74821372-89210978"
#outFolder = "test_lcr_editing"

metadata_dict = meta_df.set_index("sample").T.to_dict()

if prep == "pacbio":
    longcallr_string = "hifi-masseq"

if prep == "nanopore_direct":
    longcallr_string = "ont-drna"

if prep == "nanopore_cdna":
    longcallr_string = "ont-cdna"

rule all:
    input:
        expand(outFolder + "/{sample}/phased/{sample}.{test}.tsv", sample = samples, test = allele_tests ),
        expand(outFolder + "/{sample}/phased/{sample}.phasing_rate.tsv", sample = samples)

rule longcallR:
    output: "{outFolder}/{sample}/phased/{sample}.phased.bam"
    params:
        prefix = "{outFolder}/{sample}/phased/{sample}"
    run:
        input_bam = metadata_dict[wildcards.sample]["bam_path"]
        shell("{longcallr} -b {input_bam} -f {ref_genome} {phased_vcf_string} -p {longcallr_string} -t {lcr_threads} -o {params.prefix} {region_string}")
        shell("{samtools} index {output}")

# for each phased bam output table on how many reads were able to be phased
rule get_phasing_rate:
    input:
        "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        "{outFolder}/{sample}/phased/{sample}.phasing_rate.tsv"
    params:
        script = "scripts/get_phasing_rate.py"
    run:
        shell("python {params.script} --results-dir {outFolder} --samples {wildcards.sample} --out {output} -t {lcr_threads}")

rule split_bam:
    input:
        bam =  "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        h1 = "{outFolder}/{sample}/phased/{sample}.phased.hap1.bam",
        h2 = "{outFolder}/{sample}/phased/{sample}.phased.hap2.bam"
    run:
        shell("{samtools} view -h -b -d HP:1 {input.bam} > {output.h1}\
            {samtools} index {output.h1}\
            {samtools} view -h -b -d HP:2 {input.bam} > {output.h2}\
            {samtools} index {output.h2}")

rule allele_specific_splicing:
    input:
        bam = "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        bed = "{outFolder}/{sample}/phased/{sample}.asj.0.05.bed",
        tsv = "{outFolder}/{sample}/phased/{sample}.asj.tsv"
    params:
        prefix = "{outFolder}/{sample}/phased/{sample}"
    run:
        shell("{longcallr} asj -a {ref_gtf} -b {input.bam} -f {ref_genome} -o {params.prefix} -t {lcr_threads}")
        shell("{asj_to_bed} {output.tsv} 0.05 > {output.bed}")

rule allele_specific_expression:
    input:
        bam = "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        ase = "{outFolder}/{sample}/phased/{sample}.ase.tsv"
    params:
        prefix = "{outFolder}/{sample}/phased/{sample}"
    run:
        shell("{longcallr} ase -a {ref_gtf} -b {input.bam} -o {params.prefix} -t {lcr_threads}")

rule allele_specific_editing:
    input:
        bam = "{outFolder}/{sample}/phased/{sample}.phased.bam"
    output:
        ased = "{outFolder}/{sample}/phased/{sample}.asediting.tsv"
    params:
        script = "scripts/longcallR-asediting_v13.py",
        prefix = "{outFolder}/{sample}/phased/{sample}"
    run:
        shell("python {params.script} \
        -b {input.bam} \
        -r {rediportal} \
        -o {params.prefix} \
        -a {ref_gtf} \
        -t {lcr_threads} \
        --min_coverage 10 \
        --min_editing_rate 0.01")

#rule assemble_{outFolder}:
# in R
#x <- list.files(pattern = "TSPAN14.*asj.tsv", recursive = TRUE); names(x) <- dirname(x)
#d <- map_df(x, read_tsv, .id = "sample", col_types = "cccnnnnnnnllc") %>% arrange(P_value)

#x <- list.files(pattern = "TSPAN14.*asediting.tsv", recursive = TRUE); names(x) <- dirname(x)
#d2 <- map_df(x, read_tsv, .id = "sample", col_types = "cnccnnnnnnnnnnnl") %>% arrange(P_value)
