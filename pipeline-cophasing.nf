nextflow.enable.dsl=2

// Default parameter values
params.bins = [50000, 100000, 200000]
params.out = "out"
params.debug_out = ""

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
    tuple val(name), path("${name}.bam.bed")

    script:
    """
    bedtools bamtobed -i $bam |  sort -k1,1V -k2,2n -k3,3n > ${name}.bam.bed
    """
}

process convertVcfToBed 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf) 
    
    output:
    tuple val(name), path("${name}.vcf.bed")

    script:
    """
    cat $vcf | vcf2bed |  sort -k1,1V -k2,2n -k3,3n > ${name}.vcf.bed
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

process findClosesBed 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(bam_bed), path(vcf_bed) 
    
    output:
    tuple val(name), path("${name}.closest.bed")
    
    script:
    """
    bedtools closest -d -t all -k 2 -a $bam_bed -b $vcf_bed > ${name}.closest.bed
    """
}

process splitBamFilesToHaps 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(bam), path(closest_bed), val(hap)
    
    output:
    tuple val("${name}_${hap}"), path("${name}_${hap}.bam.bed")
    
    script:
    // TODO: Check if this is correct splitting!   
    filter = ("${hap}" == "hap1") ?  "1|0" : "0|1"
    """
    grep -F '$filter' $closest_bed | cut -f4 > ${name}.${hap}.bed.list
    gatk FilterSamReads -I $bam -O ${name}_${hap}.bam -READ_LIST_FILE ${name}.${hap}.bed.list -FILTER includeReadList    
    bedtools bamtobed -i ${name}_${hap}.bam > "${name}_${hap}.bam.bed"
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

// find all coverage files and paste them by column with the resolution genomic bins
process createCoverageTables
{
    publishDir params.out, mode: "copy"
    
    input:
    val(genome_name)
    tuple val(bin), path(genome_bin), val(names), path(covs)
    
    output:
    path("${genome_name}.${bin}.table")
    
    script:
    // Names of samples and files are lexicographically sorted first
    names.sort()
    cov_files = "$covs".split(" ");
    cov_files.sort()
    header = "chrom\tstart\tstop\t" + names.join("\t")
    """
    echo -e \"$header\" > ${genome_name}.${bin}.table
    paste $genome_bin ${cov_files.join(" ")} >> ${genome_name}.${bin}.table
    """
}

def getVcfFiles(bam_names) 
{
    if (params.containsKey("vcf")) 
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
    bin_sizes = Channel.from(params.bins)
    
    genome_size = getGenomeSizes(ref_fa)
    genome_bins = binGenome(genome_size, bin_sizes)
    
    filtered_vcf = filterUnphased(vcf)
    
    sample_beds = convertBamToBed(bam)

    // pileup = bcftoolsPileup(ref_fa, bam.join(filtered_vcf))
    vcf_beds = convertVcfToBed(filtered_vcf)
    closest_beds = findClosesBed(sample_beds.join(vcf_beds))
    hap_beds = splitBamFilesToHaps(bam.join(closest_beds).combine(Channel.from("hap1", "hap2")))

    coverages = calcCoverage(genome_bins.combine(sample_beds.mix(hap_beds)))
    covs_by_bin = genome_bins.join(coverages.groupTuple())
    cov_tables = createCoverageTables(ref_fa.baseName, covs_by_bin)
}