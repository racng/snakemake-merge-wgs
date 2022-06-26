import pandas as pd

# Read sample sheet
df = pd.read_csv(config["samples"], dtype=str)
samples = list(df[config["pooled_sample_col"]].unique())
#samples = list(df["sample name"].unique())
print("Samples in sample sheet: ", *samples)
# Filter samples
if len(config["include"]) > 0:
    samples = config["include"]
    
def get_patients(wildcards):
    inds = df[config['pooled_sample_col']]==wildcards.sample
    s = df.loc[inds, config['individual_sample_col']]
    return s.str.extract(config["regex"])[0].unique()


def get_patient_gvcf(wildcards):
    return [f"{config['datadir']}/{patient}.gvcf.gz" \
        for patient in get_patients(wildcards)]

# Merge gvcf
rule dbimport:
    input:
        gvcfs=get_patient_gvcf
    output:
        temp(directory("{outdir}/{sample}_db"))
    log: 
        "{outdir}/log/{sample}_dbimport.log"
    threads: config['threads']
    params:
        #threads=config['threads'],
        regions=[f" -L {r}" for r in config['regions'].split(',')],
        gvcfs=lambda wildcards, input: [f" -V {v}" for v in input["gvcfs"]]
    shell:
        "gatk GenomicsDBImport "
        "{params.gvcfs} " 
        "--genomicsdb-workspace-path {output} "
        "{params.regions} "
        "--max-num-intervals-to-import-in-parallel {threads} &> {log}"

# Joint calling on merged database 
rule genotype:
    input:
        db="{outdir}/{sample}_db",
        fa=expand("{refdir}/hg19.fa", refdir=config['refdir'])
    output:
        vcf=temp("{outdir}/{sample}.genotype.vcf.gz"),
        idx=temp("{outdir}/{sample}.genotype.vcf.gz.tbi")
    log: 
        "{outdir}/log/{sample}_genotype.log"
    shell:
        "gatk GenotypeGVCFs "
        "-R {input.fa} "
        "-V gendb://{input.db} "
        "-O {output.vcf} &> {log}"

# Extract variants, remove indels
rule extract_variants:
    input:
        vcf="{outdir}/{sample}.genotype.vcf.gz",
        index="{outdir}/{sample}.genotype.vcf.gz.tbi"
    output:
        vcf=temp("{outdir}/{sample}.var.vcf.gz"),
        idx=temp("{outdir}/{sample}.var.vcf.gz.csi")
    log: 
        "{outdir}/log/{sample}_extract_var.log"
    conda:
        "../envs/crossmap.yaml"
    shell:
        "vcftools --gzvcf {input.vcf} --remove-indels --recode --recode-INFO-all " 
        "--stdout | bgzip > {output.vcf} 2> {log} && "
        "bcftools index {output.vcf}"

# Keep variants in exons
rule filter_exons:
    input:
        vcf="{outdir}/{sample}.var.vcf.gz",
        idx="{outdir}/{sample}.var.vcf.gz.csi",
        bed=expand("{refdir}/hg19_exonp250.bed", refdir=config['refdir'])
    output:
        vcf=temp("{outdir}/{sample}.var.exonp250_hg19.vcf.gz"),
        idx=temp("{outdir}/{sample}.var.exonp250_hg19.vcf.gz.csi")
    conda:
        "../envs/crossmap.yaml"
    shell:
        "bcftools view {input.vcf} -R {input.bed} -O u | "
        "bcftools sort -O z -o {output.vcf} && "
        "bcftools index {output.vcf}"

# Liftover from hg19 to hg38
rule liftover:
    input:
        chain=expand("{refdir}/hg19ToHg38.over.chain.gz", refdir=config['refdir']),
        fa=expand("{refdir}/hg38.fa", refdir=config['refdir']),
        vcf="{outdir}/{sample}.var.exonp250_hg19.vcf.gz",
        idx="{outdir}/{sample}.var.exonp250_hg19.vcf.gz.csi"
    output:
        tvcf=temp("{outdir}/{sample}.var.exonp250_hg38_temp.vcf"),
        unmap=temp("{outdir}/{sample}.var.exonp250_hg38_temp.vcf.unmap")
    log: 
        "{outdir}/log/{sample}_liftover.log"
    conda:
        "../envs/crossmap.yaml"
    shell:
        "CrossMap.py vcf {input.chain} {input.vcf} {input.fa} {output.tvcf} &> {log}"
        
# Remove troublesome KI27* contigs. Filter by regions.
rule filter_regions:
    input:
        tvcf="{outdir}/{sample}.var.exonp250_hg38_temp.vcf"
    output:
        tbcf=temp("{outdir}/{sample}.var.exonp250_hg38_temp.bcf"),
        tbcfi=temp("{outdir}/{sample}.var.exonp250_hg38_temp.bcf.csi"),
        vcf="{outdir}/{sample}.var.exonp250_hg38.vcf.gz",
        index="{outdir}/{sample}.var.exonp250_hg38.vcf.gz.csi"
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


