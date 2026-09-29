# Integrate allele-specific features with GWAS
# Jack Humphrey 2026
import pandas as pd
import os

lcr_threads = 16

shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate isoseq-pipeline;")

#longcallr = "/sc/arion/projects/ad-omics/data/software/longcallR/target/release/longcallR"

longcallr = config["longcallr_path"]
asj_to_bed = longcallr + "/allele_specific/asj_to_bed.py"

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

samples = meta_df['sample']


# sanity check - can the GWAS and QTL data be found in the database?
gwas_data = config["gwas"]

gwas_df = pd.read_excel("/sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx", sheet_name = "GWAS")

gwas_check = all(i in gwas_df["dataset"].tolist() for i in gwas_data )

if not gwas_check:
    print(" * GWAS cannot be found in database")
    sys.exit()

rule all:
    input:
        expand(outFolder + "/GWAS/{GWAS}.{out_type}", GWAS = gwas_data, out_type = ["ase_gwas_joint.tsv", "asj_gwas_joint.tsv", "asediting_gwas_joint.tsv"] )

# GWAS testing
# for a GWAS, extract genome-wide significant SNPs
# overlap with SNPs used to phase RNA
# for each ASE, ASJ, ASED, find features (genes, junctions, edSites) within those haplotypes
# get all AS features the right way around, with respect to the GWAS SNP
# compare allelic FCs and do binomial test to get P-value of consistent direct of allelic bias
# do inverse-variance weighting meta-analysis when present in 2 or more samples
# scripts]$ python longcallR-gwas-joint.py -h
#

rule integrate_gwas:
    input:
        ase = expand(outFolder + "/{sample}/phased/{sample}.ase.tsv", sample = samples),
        asj = expand(outFolder + "/{sample}/phased/{sample}.asj.tsv", sample = samples),
        ased = expand(outFolder + "/{sample}/phased/{sample}.asediting.tsv", sample = samples)
    output:
        ase = outFolder + "/GWAS/{GWAS}.ase_gwas_joint.tsv",
        asj = outFolder + "/GWAS/{GWAS}.asj_gwas_joint.tsv",
        ased = outFolder + "/GWAS/{GWAS}.asediting_gwas_joint.tsv"
    params:
        script = "scripts/longcallR-gwas-joint.py"
    run:
        gwas = gwas_df.loc[ gwas_df["dataset"] == wildcards.GWAS, "full_processed_path"].values[0]
        gwas_build = gwas_df.loc[ gwas_df["dataset"] == wildcards.GWAS, "build"].values[0]
        if gwas_build == "hg19":
            liftover_string = "--liftover hg19"
        else:
            liftover_string = ""
        gwas_out = outFolder + "/GWAS/" + wildcards.GWAS
        shell("python {params.script} \
            --results-dir {outFolder} \
            --mode ase \
            --gwas {gwas} {liftover_string}\
            --min-samples 1 \
            --out {gwas_out}" )
        shell("python {params.script} \
            --results-dir {outFolder} \
            --mode asj \
            --min-samples 1 \
            --gwas {gwas} {liftover_string}\
            --out {gwas_out}" )
        shell("python {params.script} \
            --results-dir {outFolder} \
            --mode asediting \
            --min-samples 1 \
            --gwas {gwas} {liftover_string}\
            --out {gwas_out}" )


