import pandas as pd

# Read sample sheet
df = pd.read_csv(config["samples"], dtype=str)
samples = list(df["sample name"].unique())
# Filter samles:
if len(config["include"]) > 0:
    samples = config["include"]
    
def get_patients(wildcards):
    inds = df[config['pooled_sample_col']]==wildcards.sample
    s = df.loc[inds, config['individual_sample_col']]
    return s.str.extract(config["regex"])[0].values


def get_patient_gvcf(wildcards):
    return [f"{config['datadir']}/{patient}.gvcf.gz" \
        for patient in get_patients(wildcards)]

def gvcfs_wrapper(wildcards):
    return [f"-V {f}" for f in get_patient_gvcf(wildcards)]

# Merge gvcf
rule dbimport:
    input:
        gvcfs=get_patient_gvcf
    output:
        "{outdir}/{sample}_db"
    threads: config['threads']
    params:
        regions=[f" -L {r}" for r in config['regions'].split(',')]
    shell:
        "gatk GenomicsDBImport "
        "$(echo {input.gvcfs} | sed 's/^/ -V /') " 
        "--genomicsdb-workspace-path {output} "
        "{input.regions} "
        "--max-num-intervals-to-import-in-parallel {threads}"

# Joint calling on merged database 
rule genotype:
    input:
        db="{outdir}/{sample}_db",
        fa=expand("{refdir}/hg19.fa", refdir=config['refdir'])
    output:
        "{outdir}/{sample}.vcf.gz"
    shell:
        "gatk GenotypeGVCFs "
        "-R {input.fa} "
        "-V gendb://{input.db} "
        "-O {output}"

# Extract variants, remove indels
rule extract_variants:
    input:
        "{outdir}/{sample}.vcf.gz"
    output:
        "{outdir}/{sample}.var.vcf.gz"
    conda:
        "../envs/crossmap.yaml"
    shell:
        "vcftools --gzvcf {input} --remove-indels --recode --recode-INFO-all " 
        "--stdout | bgzip > {output} && "
        "bcftools index {output}"

# Keep variants in exons
rule filter_exons:
    input:
        vcf="{outdir}/{sample}.var.vcf.gz",
        bed=expand("{refdir}/hg19_exonp250.bed", refdir=config['refdir'])
    output:
        "{outdir}/{sample}.var.exonp250_hg19.vcf.gz"
    conda:
        "../envs/crossmap.yaml"
    shell:
        "bcftools view {input.vcf} -R {input.bed} -O u | "
        "bcftools sort -O z -o {output} && "
        "bcftools index {output}"

# Liftover from hg19 to hg38
rule liftover:
    input:
        chain=expand("{refdir}/hg19ToHg38.over.chain.gz", refdir=config['refdir']),
        fa=expand("{refdir}/hg38.fa", refdir=config['refdir']),
        vcf="{outdir}/{sample}.var.exomep250_hg19.vcf.gz"
    output:
        tvcf=temp("{outdir}/{sample}.var.exomep250_hg38_temp.vcf"),
    conda:
        "../envs/crossmap.yaml"
    shell:
        "CrossMap.py vcf {input.chain} {input.vcf} {input.fa} {output.tvcf}"
        
# Remove troublesome KI27* contigs. Filter by regions.
rule filter_regions:
    input:
        tvcf="{outdir}/{sample}.var.exomep250_hg38_temp.vcf"
    output:
        tbcf=temp("{outdir}/{sample}.var.exomep250_hg38_temp.bcf"),
        tbcfi=temp("{outdir}/{sample}.var.exomep250_hg38_temp.bcf.csi"),
        vcf="{outdir}/{sample}.var.exomep250_hg38.vcf.gz",
        index="{outdir}/{sample}.var.exomep250_hg38.vcf.gz.csi"
    params:
        reg=config['regions']
    conda:
        "../envs/crossmap.yaml"
    shell:
        "sed -i '/_KI27/d' {input.tvcf} && "
        "bcftools sort {input.tvcf} -O b -o {output.tbcf} && "
        "bcftools index {output.tbcf} && "
        "bcftools view {output.tbcf} --regions {params.reg} -O z -o {output.vcf} && "
        "bcftools index {output.vcf}"








