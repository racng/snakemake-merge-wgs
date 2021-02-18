# Snakemake workflow: Merging GVCFs variants
Snakemake pipeline for processing [GVCF files](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format) to generate multi-sample VCF of exon variants. 
The workflow is built for GVCFs based on hg19 alignment and outputs variants that are lifted over to hg38. 

Current workflow is based on [GATK tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035889971--How-to-Consolidate-GVCFs-for-joint-calling-with-GenotypeGVCFs) for joint calling on consolidated GVCFs:
1. Merge GVCFs using GenomicsDBImport for specfified intervals.
2. Joint calling on merged GenomicsDB using GenotypeGVCFs.
3. Extract variants (remove indels) using vcftools.
4. Filter for exon variants using bcftools and bed file.
5. Liftover to hg38 using CrossMap.
6. Filter VCF by region and remove Ki27 contig.

Previous workflow extracted variants from each GVCF, lifted over to hg38, filtered by region and quality, and then finally merged into one vcf (old workflow in `variants.smk`). This is deprecated due to loss of variants from independent filtering per sample. Assuming reference genotype for missing variants might have introduced artificial genotypes during the merging step. 

## Requirements 
Dependencies are managed by snakemake via conda environment. See `envs/crossmap.yaml` for full list of dependencies installed via conda. Missing linux libraries may cause error during installing of CrossMap (ex. libcurl4-openssl-dev). 
Install them with `sudo apt-get install` command before continuing running snakemake pipeline.

[GATK4 installation](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4) is required. 

## Usage
In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository.

### Step 1: Install workflow

Clone this repository:
```
git clone https://github.com/racng/snakemake-merge-wgs.git
```

### Step 2: Configure workflow

Prepare a samplesheet csv file with the individual sample ids and their pooled sample name. Column names are flexible but must be consistent with your configuration file. 

Here is an example sample sheet, where `PBMC sample origin` contains the names of individual samples, and `sample name` contains the pooled sample label to group by. Note that multiple samples can come from the same patient, so we provide additional string parsing argument for extracting patient id (see regex below).

| PBMC sample origin | sample name |
| :---: | :---: |
| INCOV001-1 | B1 |
| INCOV002-1 | B1 |
| INCOV003-2 | B2 |
| INCOV004-1 | B2 |
| 1234BW | B2 |

Configure the workflow according to your needs via editing `config.yaml`.
- `datadir`: Directory containing GVCFs formated as "{patient_id}.gvcf.gz".
- `refdir`: Directory containing references
- `outdir`: Output directory 
- `samples`: Path to sample sheet csv file
- `pooled_sample_col`: Column name for pooled sample labels (ex. PBMC sample origin)
- `individual_sample_col`:  Column name for individual sample ids (ex. sample name)
- `include`: Optional, specified list of pooled samples to process. Put empty list `[]` if processing all pooled samples.
- `regex`: Regex pattern for parsing patient id from the individual_sample_col in samplesheet. 
ex. `"(INCOV\\d{3}|\\d{4}BW)"` would parse out `INCOV001, INCOV002, INCOV003, INCOV004, 1234BW` from the example sample sheet. If no string parsing is necessary, put `"(.+)"`. Note that backslash are special characters for YAML files.
- `regions`: regions/intervals to include when running GenomicsDBImport and bcftools filtering. 
ex. `"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"`
- `threads`: Number of parallel threads allowed (used by GenomicsDBImport)

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [mamba](https://github.com/mamba-org/mamba):

    conda install -c conda-forge mamba
    mamba create -c conda-forge -c bioconda -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally using `$N` cores (or 'all' cores)

    snakemake --use-conda --cores $N

Use `--configfile` option to select a different config file.
To download packages to a dedicated folder for multiple runs, use the `--conda-prefix` option.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

To run snakemake in the background and pipe the terminal output to file, wrap the above commands 

    nohup [INSERT SNAKEMAKE COMMAND] > log.txt & echo $! > pid.txt

Process ID is saved in `pid.txt`.
Check background job and its jobid with `jobs`. Bring back to foreground with `fg [jobid]`. 

## Authors

- Rachel Ng (@racng) - Heath Lab, Institute for Systems Biology
- Chris Lausted - Hood Lab, Institute for Systems Biology
