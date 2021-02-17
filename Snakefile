configfile: "config.yaml"

include: "rules/merge.smk"
include: "rules/download.smk"

rule all:
    input:
        # Download references if necessary
        ancient(expand("{refdir}/hg19_exonp250.bed", refdir=config['refdir'])),
	    ancient(expand("{refdir}/hg38.fa", refdir=config['refdir'])),
        ancient(expand("{refdir}/hg38.fa.fai", refdir=config['refdir'])),
        ancient(expand("{refdir}/hg19ToHg38.over.chain.gz", refdir=config['refdir'])),

        # Merged output per sample
        expand("{outdir}/{sample}.var.exonp250_hg38.vcf.gz", 
            outdir=config['outdir'], sample=samples),



