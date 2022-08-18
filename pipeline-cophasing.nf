nextflow.enable.dsl=2

params.out_dir = "out"
params.debug_out = "out"
scripts_folder = "${projectDir}/scripts"


process filterUnphased 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf) 

    output:
    tuple val(name), path("${name}.recode.vcf") 

    script:
    """
    vcftools --vcf $vcf --out $name --phased --recode
    """
}

process convertBamToBed 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(bam) 
    
    output:
    tuple val(name), path("${name}.bed")

    script:
    """
    bedtools bamtobed -i $bam > ${name}.bed
    """
}

process bcftoolsPileup 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(ref_genome)
    tuple val(name), path(bam), path(vcf)
    
    output:
    tuple val(name), path("${name}.pileup")

    script:
    """
    bcftools mpileup -f $ref_genome -T $vcf -a FORMAT/AD,INFO/AD -O v $bam > ${name}.pileup 
    """
}

process SNPsPerSample 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(ref_genome)
    tuple val(name), path(bam), path(vcf)
    
    output:
    tuple val(name), path("${name}.pileup")

    script:
    """
    $scripts_folder/co-phasing/create_sample_spec_SNP_seeds_to_find_closest_obs_SNP.R 
    """
}

process getGenomeSizes 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(ref_genome)
    
    output:
    path("${ref_genome.baseName}.sizes")

    script:
    """
    samtools faidx ${ref_genome}
    cut -f1,2 ${ref_genome}.fai > ${ref_genome.baseName}.sizes
    """
}

process binGenome 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(genome_sizes)
    val(bin)
    
    output:
    tuple val(bin), path("${genome_sizes.baseName}.${bin}.bed")

    script:
    """
    bedtools makewindows -g $genome_sizes -w ${bin} > ${genome_sizes.baseName}.${bin}.bed
    """
}

process calcCoverage
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    tuple val(bin), path(bin_bed), val(name), path(read_bed)
    
    output:
    tuple val(bin), val(name), path("${name}.${bin}.cov")
    
    script:
    """
    bedtools coverage -a $bin_bed -b $read_bed | awk \"{print \\\$5}\" > ${name}.${bin}.cov
    """
}


process createCoverageTables
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    tuple val(bin), path(genome_bin), val(names), path(covs)
    
    output:
    path("${bin}.table")
    
    script:
    header = "chrom\tstart\tstop\t" + names.join("\t")
    """
    # find all coverage files and paste them by column with the resolution genomic bins
    echo -e \"$header\" > ${bin}.table
    paste $genome_bin $covs >> ${bin}.table
    """
}

def getVcfFiles(bam_names) 
{
    if (params.vcf) 
    {
        vcf_files = Channel.fromPath(params.vcf)
        vcf_files.ifEmpty{error "No VCF files found matching the naming of the bam files."}
        vcf = bam_names.combine(vcf_files)
        vcf.ifEmpty{error "Failed to match VCF files onto BAM files."}
    }
    else 
    {
        vcf_filename = params.bam - ".bam" + ".vcf"
        vcf = Channel.fromFilePairs(vcf_filename, size: 1)
    }
    return vcf
}

workflow 
{
    bam = Channel.fromFilePairs(params.bam, size: 1)
    bam_names = bam.map{name, file -> name}
    vcf = getVcfFiles(bam_names)
    ref_fa = file(params.fa)
    bin_sizes = Channel.from(50000, 100000, 200000)
    
    filtered_vcf = filterUnphased(vcf)
    sample_beds = convertBamToBed(bam)
    genome_size = getGenomeSizes(ref_fa)
    genome_bins = binGenome(genome_size, bin_sizes)
    samples_with_bins = genome_bins.combine(sample_beds)
    coverages = calcCoverage(samples_with_bins)
    covs_by_bin = genome_bins.join(coverages.groupTuple())
    cov_tables = createCoverageTables(covs_by_bin)

    pileup = bcftoolsPileup(ref_fa, bam.join(filtered_vcf))
    // seeds = seedsPerSample(vcf.join(pileup))
}