rule download_bed:
    output:
        "{refdir}/hg19_exonp250.bed"
    shell:
        """
        URL='https://docs.google.com/uc?export=download&id=1_ALrh2bxyNsizYKhNdSP_gfIlq7CGphw'
        wget --no-check-certificate $URL -O {output}
        """

rule download_fa:
    output:
        "{refdir}/hg38.fa",
        "{refdir}/hg38.fa.fai"
    shell:
        "URL='https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz'"
        "wget $URL -O {refdir}/hg38.fa.gz"
        "gunzip {refdir}/hg38.fa.gz"
        "samtools faidx {refdir}/hg38.fa -o {refdir}/hg38.fa.fai"

rule download_chain:
    output:
        "{refdir}/hg19ToHg38.over.chain.gz"
    shell:
        "URL='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'"
        "wget $URL -O {refdir}/hg38.fa.gz"
        