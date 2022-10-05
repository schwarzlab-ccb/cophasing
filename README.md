# Pipeline-Cophasing (PCP) #

Nextflow-based pipeline for cophasing of GAM reads.

## Pipeline

### Input:
1. GAM experiemt samples in BAM format
2. Reference genome in FASTA format
3. Known SNPs in VCF format. **NOTE**: The VCF file must be in plaintext (not compressed) and end with `.vcf`.

### Output
A set of tables of the form `{genome_name}.{bin_size}.table` for different bin sizes with the number of reads per bin. The columns are:
1. `chrom`: chromosome name
2. `start`: start position of the bin
3. `stop`: the stop position of the bin
4. `{sample}`: the number of reads in the bin
5. `{sample_hap1}`: the number but for the first haplotype
6. `{sample_hap2}`: the number but for the second haplotype

### Process

#### Align reads to the reference haplotype
1. filter out positions without phasing information
2. convert VCF to BED format
3. find the closest snip for each read and merge

#### Create bins
1. bin the genome into fixed sized, non-overlapping windows of desired resolution, eg 50kb, 100kb, 200kb 

#### Calculat coverage 
1. for comparison, calculate the coverage of each window using the original unsplit GAM samples,
2. calculate coverage files of all split GAM samples for each haplotype,
3. combine coverage files of all samples into one coverage table, per resolution

![Pipeline-Cophasing](./doc/pipeline_cophasing_chart.png)

## Requirements:
You can skip installing the tools if you use conda with the provided environment (see below).

Requirements:
* java-jre
* nextflow
* bedtools
* bfctools
* vcftools
* samtools
* bedops
* gatk4

## Execution

To execture run 
`nextflow run pipeline-cophasing.nf [parameters]`

### Output

By default the results are written to the `./out` folder.

### Test run

Random testing data are provided as a part of the package. The default test execution can be done as following:

* Download the reference data using `sh ./DownloadRefData.sh`.
* Execute based on the environment you are using:
    * [Conda] `nextflow -C test_data.config run pipeline-cophasing.nf -with-conda pcp-env.yaml`
    * [System] `nextflow -C test_data.config run pipeline-cophasing.nf`

**NOTE:** The parameters for the execution are stored in the Nextflow configuration file `test_data.config`.

### Parameters

#### Mandatory
* `--fa path` reference file either as `.fa`  or `.fa.gz`.
* `--bam path` alignment files either as `.bam` or `.sam`. This can be a glob pattern (e.g. `sample*.bam`). All files matching the pattern are used then.

**NOTE:** Read and variants are aligned by name. For example `sample5.bam` is matched to `sample5.vcf`. If you have single vcf for all samples, use the `vcf` parameter detailed below.

#### Default 
* `--bins [int]` the bin sizes to be used, `default=[50000, 100000, 200000]`,
* `--vcf path` variant call files either as `.vcf` or `.vcf.gz`. This can be a glob pattern (e.g. `sample*.vcf`), `default={bam_filè}.vcf` for each bam file provided.
* `--out path` a path to a folder where the output is stored, `default=./out`.

## Input data 

Note that GATK requires `.gz` files to be compressed with `bgzip`, not `gzip`.

#### Tested HG38 version
Downloaded from: `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`

## Contact
Email questions, feature requests and bug reports to **Adam Streck, adam.streck@mdc-berlin.de**.

## License
Pipeline-CoPhasing is available under the MIT License. 
