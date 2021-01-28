import pandas as pd

# Read sample sheet
df = pd.read_csv(config["samples"], dtype=str)
samples = list(df["sample name"].unique())

def get_patients(wildcards):
    inds = df[config['pooled_sample_col']]==wildcards.sample
    s = df.loc[inds, config['individual_sample_col']]
    return s.str.extract(config["regex"])[0].values

# Convert gvcf to vcf
rule gvcf_to_vcf:
    input:
        expand("{datadir}/{patient}.gvcf.gz", datadir=config['datadir'], 
            allow_missing=True)
    output:
        vcf=temp("{outdir}/{patient}.var.vcf.gz"),
        index=temp("{outdir}/{patient}.var.vcf.gz.csi")
    conda:
        "../envs/crossmap.yaml"
    shell:
        "zcat {input} | extract_variants | bgzip > {output.vcf} && "
        "bcftools index {output.vcf} -o {output.index}"

## Filter variants 
rule filter_snps:
    input:
        bed=expand("{refdir}/hg19_exonp250.bed", refdir=config['refdir']),
        vcf="{outdir}/{patient}.var.vcf.gz",
        index="{outdir}/{patient}.var.vcf.gz.csi"
    output:
        vcf="{outdir}/{patient}.exome.vcf.gz",
        index="{outdir}/{patient}.exome.vcf.gz.csi"
    # params:
    #     flt='MIN(FMT/DP)>20 & MIN(FMT/GQ)>20'
    conda:
        "../envs/crossmap.yaml"
    shell:
        "FLT='MIN(FMT/DP)>20 & MIN(FMT/GQ)>20' && "
        "bcftools view {input.vcf} -R {input.bed} -i $FLT -O u | "
        "bcftools sort -O z -o {output.vcf} && "
        "bcftools index {output.vcf} -o {output.index}"

def get_patient_vcf(wildcards):
    return [f"{wildcards.outdir}/{patient}.exome.vcf.gz" for patient in get_patients(wildcards)]

## Merge. Remove small chr. Remove small number of records with <NON_REF>.
rule merge:
    input:
        get_patient_vcf
    output:
        "{outdir}/{sample}_exomep250_hg19.vcf.gz"
    params:
        reg='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,'+\
        'chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY'
    conda:
        "../envs/crossmap.yaml"
    shell:
        "bcftools merge -m all -0 {input} --regions {params.reg} -O u | "
        "bcftools view -a -c 1 -O v | "
        "sed '/<NON_REF>/d' | "  
        "bcftools norm -d all -O z -o {output}"     

## Convert VCF to hg38.  Remove troublesome KI27* contigs.    
rule liftover:
    input:
        chain=expand("{refdir}/hg19ToHg38.over.chain.gz", refdir=config['refdir']),
        fa=expand("{refdir}/hg38.fa", refdir=config['refdir']),
        vcf="{outdir}/{sample}_exomep250_hg19.vcf.gz"
    output:
        tvcf=temp("{outdir}/{sample}_exomep250_hg38_temp.vcf"),
        tbcf=temp("{outdir}/{sample}_exomep250_hg38_temp.bcf"),
        tbcfi=temp("{outdir}/{sample}_exomep250_hg38_temp.bcf.csi"),
        vcf="{outdir}/{sample}_exomep250_hg38.vcf.gz",
        index="{outdir}/{sample}_exomep250_hg38.vcf.gz.csi"
    params:
        reg='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,'+\
        'chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY'
    conda:
        "../envs/crossmap.yaml"
    shell:
        "CrossMap.py vcf {input.chain} {input.vcf} {input.fa} {output.tvcf} && "
        "sed -i '/_KI27/d' {output.tvcf} && "
        "bcftools sort {output.tvcf} -O b -o {output.tbcf} && "
        "bcftools index {output.tbcf} && "
        "bcftools view {output.tbcf} --regions {params.reg} -O z -o {output.vcf} && "
        "bcftools index {output.vcf}"


