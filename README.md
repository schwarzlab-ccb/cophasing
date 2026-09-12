# Cophasing (v1.0.0)

Nextflow-based pipeline for cophasing of GAM reads.

## Cophasing Pipeline

### Input:
1. GAM experiment samples in the BAM format,
2. Reference genome in the FASTA format,
3. **Phased** SNPs in the VCF format.

### Output
A set of tables of the form `{genome_name}.{bin_size}.{hap}.coverage.tsv` and `{genome_name}.{bin_size}.{hap}.segregation.tsv` for different bin sizes with the number of reads per bin per `sample`. Haplotype `hap` is either `hap1`, `hap2`, or `both`.  The columns are:
1. `chrom`: chromosome name
2. `start`: start position of the bin
3. `stop`: the stop position of the bin
4. `[sample]`: The number of reads in the bin for `coverage.tsv`, and 1 or 0 for `segregation.tsv`, with 1 meaning the bin has been determined to be covered after filtering, 0 the opposite.

**Note**: For segregation, only autosomes are considered. The filter requires that the chromosomes are named in either the UCSC (chr1-chr22) or Ensembl format (1-22).

### Process

#### Align reads to the reference haplotype
1. filter out positions without phasing information
2. convert VCF to BED format
3. find the closest SNP for each read and merge

#### Create bins
1. bin the genome into fixed sized, non-overlapping windows of desired resolution, eg 50kb, 100kb, 200kb 

#### Calculate coverage 
1. for comparison, calculate the coverage of each window using the original unsplit GAM samples,
2. calculate coverage files of all split GAM samples for each haplotype,
3. combine coverage files of all samples into one coverage table, per resolution

The coverage table describes the number of **bases** covered by reads in each window; therefore, it is possible that the sum of coverage from both haplotypes is higher than the total number of bases covered in a window, as some bases may be covered by both haplotypes.

![Pipeline-Cophasing](./doc/pipeline_cophasing_chart.png)

#### Segregation
The segregation algorithm removes spurious bins in particular the following two steps:
1. Finds a separation threshold for bin count to remove noise reads.
2. Remove orphan bins (bins with no neighbours).

## Requirements:
Use the provided conda environment file `pcp-env.yaml` to install all the required software.

## Execution

To execute, run:

`nextflow run main.nf [parameters]`

### Output

By default the results are written to the `./out` folder.

### Test run

Random testing data are provided as a part of the package. Run the following command to test the pipeline:

`nextflow run main.nf -c test_data.config` 

**NOTE:** The parameters for the execution are stored in the Nextflow configuration file `test_data.config`.

### Parameters

#### Mandatory
* `--bam path` alignment files either as `.bam` or `.sam`. This can be a glob pattern (e.g. `sample*.bam`). All files matching the pattern are used then.
* `--fa path` reference file either as `.fa`  or `.fa.gz`.
* `--vcf path` a variant call file either as `.vcf` or `.vcf.gz`. Must contain phased GT information.

**Note**: GATK requires `.gz` files to be compressed with `bgzip`, not `gzip`.

#### Default 
* `--name string` the prefix that will be given to the samples, `default=<the name of the reference file>`,
* `--bins [int]` the bin sizes to be used, `default=[50000, 100000, 200000]`,
* `--out path` a path to a folder where the output is stored, `default=./out`,
* `--min_depth int` a minimum read depth per SNP to be included, `default=1`,
* `--min_ratio int` tests that the dominant base for a SNP is at least `min_ratio` times more often present than the remaining observed bases [BCFTools](https://samtools.github.io/bcftools/bcftools.html#expressions), `default=5`
* `--cutoff int` maximum distance from a read to a closest so that the read is still matched to the SNP, `default=0`.

### Visualization

A sample visualization method is shown in the notebook `NMPI_matrix_vis.ipynb`. This notebook loads a singular segregation table and plots a contact map for one chromosome in a specified region. 

## Permutation Test Pipeline

A downstream Nextflow-based pipeline that takes the segregation tables produced by the Cophasing pipeline and identifies differential contacts between haplotypes using permutation testing.


### Input
The output directory of the main co-phasing pipeline (`--out`).

### Output
Results are written to subdirectories inside `output_dir`:
* `01_curated_segregation_tables/` — segregation tables after WDF-based dat curation
* `02_permutation_test/` — per-chromosome permutation test results (`.pkl`)
* `03_thresholds/` — knee-point threshold tables (`.tsv`) per chromosome
* `04_cool_files/` — contact matrices in `.cool` format for differential contacts

### Process
1. **Curation** — filters segregation tables based on windowed deviation from median (WDF), removing low-quality bins using `high_factor` and `low_factor` thresholds
2. **Permutation test** — runs a chromosome-level permutation test (`num_perm` permutations) to calculate empirical p-values for each contact
3. **Threshold identification** — identifies a significance threshold per chromosome
4. **Cool file generation** — converts significant permutation test results into `.cool` contact matrices

```mermaid
---
config:
  flowchart:
    nodeSpacing: 2
    rankSpacing: 5
    diagramPadding: 1
    padding: 1
  themeVariables:
    fontSize: 8px
---
flowchart TD
    IN[/"CoPhasing output\nhap1 · hap2 · both\nsegregation tables"/]

    IN ==> S1

    subgraph S1["`① **Curation**`"]
        P1["Window detection frequency-based bin filtering\ncurates segregation tables and generates masking file"]
    end

    S1 ==> CUR[/"curated_*_cutoff.*.hap1.segregation.tsv\ncurated_*_cutoff.*.hap2.segregation.tsv\ncurated_*_cutoff.*.both.segregation.tsv"/]
    S1 ==> BINS[/"bins_rm_hap1_hap2.tsv"/]

    CUR ==> S2
    BINS ==> S2

    subgraph S2["`② **Permutation Test**`"]
        P2["per chromosome<br/>• calculate contact frequencies (CF) per haplotype<br/>• calculate ∆CF (hap1-hap2)<br/>• build null distribution by randomly reshuffling segregation tables across haplotypes"]
    end

    S2 ==> PKL[/"permutation_test_results\n_chr_multiprocessing.pkl"/]

    PKL ==> S3
    BINS ==> S3

    subgraph S3["`③ **Threshold Identification**`"]
        P3["• p-value calculation<br/>• normalized Gaussian smoothing<br/>• knee-point detection per chromosome"]
    end

    S3 ==> KNEE[/"knee_point_table_threshold_steps.tsv"/]
    S3 ==> RATIO[/"ratios_table_threshold_steps.tsv\nratio_differential_contacts_info.tsv"/]

    CUR ==> S4
    PKL ==> S4
    KNEE ==> S4
    BINS ==> S4

    subgraph S4["`④ **Cool File Generation**`"]
        P4["per chromosome<br/>• Identifiction of differential contacts based on thresholding<br/>• save results as long format .tsv and .cool files"]
    end

    S4 ==> COOL[/"permutation_test_results_*.cool\npermutation_test_results_*.tsv.gz"/]

    classDef input fill:#b5c7eb,stroke:#0d47a1,color:#000,stroke-width:2px;
    classDef process fill:#f5f5f5,stroke:#555555,color:#000,stroke-width:2px;
    classDef file fill:#99c5c4,stroke:#008080,color:#000,stroke-width:2px;
    classDef output fill:#ffe2a6,stroke:#ffcc80,color:#000,stroke-width:2px;


    class IN input;
    class P1,P2,P3,P4 process;
    class CUR,BINS,PKL,KNEE,RATIO file;
    class COOL output;

    classDef leftAlign text-align:left;
    class P1,P2,P3,P4 leftAlign;
```

### Execution
`nextflow run pipeline-permutation-test.nf -c permutation_test.config`


## Authors
This pipeline has been developed at Max Delbrück Center for Molecular Medicine, Berlin. Authors:

* Dr. Adam Streck: cophasing pipeline development,
* Dr. Julia Markowski: creator of the co-phasing method,
* Dr. Alexander Kukalev: creator of the separation algorithm,
* Claudia Robens: permutation test pipeline development.

## Contact
Email questions, feature requests and bug reports to **Adam Streck, adam.streck@iccb-cologne.org**.

## License
This repository is available under the MIT License. 
