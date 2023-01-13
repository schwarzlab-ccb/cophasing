nextflow.enable.dsl=2

// IMPROVEMENT: consider monoallelic if ratio more than 10:1 (less strict)
// Default parameter values
params.bins = [50000, 100000, 200000]
params.out = "out"
params.debug_out = ""
params.cutoff = 0
params.min_depth = 10

process filterUnphased 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple path(vcf) 

    output:
    tuple path("${vcf.simpleName}.filtered.vcf.gz"), path("${vcf.simpleName}.filtered.vcf.gz.tbi") 

    script:
    vcftool_type = "${vcf.name}".endsWith(".vcf.gz") ? "gzvcf" : "vcf"
    """
    vcftools --$vcftool_type $vcf --phased --recode --stdout | bgzip -c > ${vcf.simpleName}.filtered.vcf.gz
    tabix ${vcf.simpleName}.filtered.vcf.gz
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
    bedtools bamtobed -i $bam | sort -k1,1 -k2,2n > ${name}.bam.bed
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
    obs_ref = ":[[:digit:]]+,0" // OBSERVED REFERENCE
    obs_alt = ":0,[[:digit:]]+,0" // OBSERVED ALTERNATIVE
    """
    cat $vcf | vcf2bed | cut -f1-3,11 | sort -k1,1 -k2,2nn > ${name}.vcf.bed    
    sed -E -i "s/0\\|1.*${obs_ref}/hap1/g" ${name}.vcf.bed 
    sed -E -i "s/1\\|0.*${obs_alt}/hap1/g" ${name}.vcf.bed 
    sed -E -i "s/0\\|1.*${obs_alt}/hap2/g" ${name}.vcf.bed 
    sed -E -i "s/1\\|0.*${obs_ref}/hap2/g" ${name}.vcf.bed 
    """
}

process bcftoolsPileup 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(ref_genome)
    tuple val(name), path(bam), path(vcf), path(tb)
    
    output:
    tuple val(name), path("${name}.pileup.vcf")

    script:
    """
    bcftools mpileup -f $ref_genome -T $vcf -a FORMAT/DP,FORMAT/AD -O v $bam | vcf-sort | bgzip > ${name}.temp.vcf.gz 
    tabix ${name}.temp.vcf.gz
    echo `bcftools query -l $vcf` `bcftools query -l ${name}.temp.vcf.gz` > samples.txt
    if [ `wc -l < samples.txt` != 1 ]; then echo "there must be exactly one sample in the VCF ${vcf}"; exit 1; fi;
    bcftools annotate -a $vcf -c ALT,FORMAT/GT ${name}.temp.vcf.gz -S samples.txt > ${name}.pileup.vcf
    """
}

// Keep only sites where are reads are allocated to one allele
process filterSites {    
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf)
    
    output:
    tuple val(name), path("${name}.mono")
    
    script:
    filter_mono = "(FORMAT/AD[0:0] > 0 && FORMAT/AD[0:1] == 0) || (FORMAT/AD[0:0] == 0 && FORMAT/AD[0:2] == 0)"
    filter_depth = "(FORMAT/DP[0:0] > $params.min_depth)"
    """
    bcftools filter -i "$filter_mono && $filter_depth" $vcf > ${name}.mono
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
    bedtools closest -d -t all -k 1 -a $bam_bed -b $vcf_bed > ${name}.closest.bed
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
    """
    grep -F '$hap' $closest_bed | cut -f4 > ${name}.${hap}.bed.list
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


workflow 
{
    bam = Channel.fromFilePairs(params.bam, size: 1)
    ref_fa = file(params.fa)
    vcf = file(params.vcf)

    bin_sizes = Channel.from(params.bins)    

    filtered_vcf = filterUnphased(vcf)
    combined = bam.combine(filtered_vcf)    
    pileup = bcftoolsPileup(ref_fa, combined)
    filtered_pileup = filterSites(pileup)
    vcf_beds = convertVcfToBed(filtered_pileup)

    sample_beds = convertBamToBed(bam)
    genome_size = getGenomeSizes(ref_fa)
    genome_bins = binGenome(genome_size, bin_sizes)
    closest_beds = findClosesBed(sample_beds.join(vcf_beds))
    hap_beds = splitBamFilesToHaps(bam.join(closest_beds).combine(Channel.from("hap1", "hap2")))

    coverages = calcCoverage(genome_bins.combine(sample_beds.mix(hap_beds)))
    covs_by_bin = genome_bins.join(coverages.groupTuple())
    cov_tables = createCoverageTables(ref_fa.baseName, covs_by_bin)
}
