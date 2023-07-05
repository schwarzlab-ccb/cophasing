nextflow.enable.dsl=2

// IMPROVEMENT: consider monoallelic if ratio more than 10:1 (less strict)

// Default parameter values
params.bins = [50000, 100000, 200000] // Bin sizes to use for analysis, each bin size will be analyzed separately
params.out = "out" // Output directory containing the results
params.debug_out = "" // If set, will output intermediate files to this directory
params.cutoff = 0 // Maximum distance between a read and a variant to be considered for analysis
params.min_depth = 1 // Minimum required read depth per variant to be considered for analysis
params.name = "" // Will default to the name of the FA file if not set

process filterUnphased 
{
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(vcf) 

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
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

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
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf) 
    
    output:
    tuple val(name), path("${name}.vcf.bed")

    script:
    obs_ref = ":(([[:digit:]]+,0\$)|([[:digit:]]+,[[:digit:]],0))\$" // OBSERVED REFERENCE
    obs_alt = ":((0,[[:digit:]]+,0\$)|(0,[[:digit:]]+,[[:digit:]],0))\$" // OBSERVED ALTERNATIVE
    """
    cat $vcf | vcf2bed | cut -f1-3,11 | sort -k1,1 -k2,2nn > ${name}.vcf.bed    
    sed -E -i "s/1\\|0.*${obs_alt}/hap1/g" ${name}.vcf.bed 
    sed -E -i "s/0\\|1.*${obs_alt}/hap2/g" ${name}.vcf.bed 
    sed -E -i "s/0\\|1.*${obs_ref}/hap1/g" ${name}.vcf.bed 
    sed -E -i "s/1\\|0.*${obs_ref}/hap2/g" ${name}.vcf.bed 
    """
}

process bcftoolsPileup 
{
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

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
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf)
    
    output:
    tuple val(name), path("${name}.mono")
    
    // Retain if either monoallelic, or biallelic with no more than one read on the second allele and more than one on first allele
    script:
    filter_mono = "(FORMAT/AD[0:0] > 1 && FORMAT/AD[0:1] < 2 && FORMAT/AD[0:2] == 0) || (FORMAT/AD[0:0] > 0 && FORMAT/AD[0:1] == 0) || (FORMAT/AD[0:0] == 0 && FORMAT/AD[0:2] < 2)"
    filter_depth = "(FORMAT/DP[0:0] >= $params.min_depth)"
    """
    bcftools filter -i "$filter_mono && $filter_depth" $vcf > ${name}.mono
    """
}


process findClosesBed 
{
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(bam_bed), path(vcf_bed) 
    
    output:
    tuple val(name), path("${name}.closest.bed")
    
    script:
    """
    bedtools closest -d -t all -k 1 -a $bam_bed -b $vcf_bed | awk '\$NF <= $params.cutoff' > ${name}.closest.bed
    """
}

// TODO: This should be improved - currently the original reads are filtered by info in the closest.bed, but the file itself could be used
process splitBamFilesToHaps 
{
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(bam), path(closest_bed), val(hap)
    
    output:
    tuple val("${name}_${hap}"), path("${name}_${hap}.bam.bed")
    
    script:
    """
    grep -F '$hap' $closest_bed | cut -f1-4,10 > ${name}_${hap}.bam.bed
    """
}

process getGenomeSizes 
{
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

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
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

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
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    tuple val(bin), path(bin_bed), val(name), path(read_bed)
    
    output:
    tuple val("${bin}.${name[-4..-1]}"), val(name), path("${name}.${bin}.cov")
    
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
    path("${genome_name}.${bin}.coverage.tsv")
    
    script:
    // Names of samples and files are lexicographically sorted first
    names.sort()
    cov_files = "$covs".split(" ");
    cov_files.sort()
    header = "chrom\tstart\tstop\t" + names.collect{it[0..-6]}.join("\t") // remove the haplotype from the sample name
    table_file = "${genome_name}.${bin}.coverage.tsv"
    """
    echo -e \"$header\" > ${table_file}
    paste $genome_bin ${cov_files.join(" ")} >> ${table_file}
    """
}

def groupFilesBySizeAndType(cov_tables) 
{
	return cov_tables.map 
    { file -> 
        def matcher = file.name =~ /(.*)\.(hap1|hap2|both)\.coverage\.tsv/
        if (matcher.matches()) {
            def size = matcher[0][1]
            def type = matcher[0][2]
            return tuple(size, type, file)
        } else {
            println("File ${file} did not match pattern")
            return null
        }
    }
    .filter { it != null }
    .groupTuple()
    .map 
    { size, types, files ->
        def fileMap = ['hap1': null, 'hap2': null, 'both': null]
        types.eachWithIndex { type, i -> 
            fileMap[type] = files[i]
        }
        return tuple(size, fileMap['hap1'], fileMap['hap2'], fileMap['both'])
    }
    .filter { it[1] != null && it[2] != null && it[3] != null }
}

process createSegregationTables
{
    publishDir params.out, mode: "copy"

    input:
    tuple val(group_name), path(hap1_cov), path(hap2_cov), path(both_cov)

    output:
    path("${group_name}.*.tsv")

    script:
    """
    python $projectDir/scripts/get_segregation_tables.py $group_name $hap1_cov $hap2_cov $both_cov
    """
}

workflow 
{
    // Inputs
    bam = Channel.fromFilePairs(params.bam, size: 1)
    ref_fa = file(params.fa)
    vcf = file(params.vcf)
    bin_sizes = Channel.from(params.bins)    

    // Create and filter pileup to obtain phased variant sites observed in the reads
    filtered_vcf = filterUnphased(vcf)
    combined = bam.combine(filtered_vcf)    
    pileup = bcftoolsPileup(ref_fa, combined)
    filtered_pileup = filterSites(pileup)
    vcf_beds = convertVcfToBed(filtered_pileup)

    // Calculate bins
    genome_size = getGenomeSizes(ref_fa)
    genome_bins = binGenome(genome_size, bin_sizes)

    // Split reads into haplotypes 
    sample_beds = convertBamToBed(bam)
    closest_beds = findClosesBed(sample_beds.join(vcf_beds))
    named_beds = sample_beds.map { it -> [it[0] + "_both", it[1]]}
    hap_beds = splitBamFilesToHaps(bam.join(closest_beds).combine(Channel.from("hap1", "hap2")))
    all_beds = named_beds.mix(hap_beds)

    // Calculate coverage for each sample and bin
    window_sample_pairs = genome_bins.combine(all_beds)
    coverages = calcCoverage(window_sample_pairs)
    tables = genome_bins.combine(Channel.from(["both", "hap1", "hap2"])).map { it -> [it[0] + "." + it[2], it[1]] }
    covs_by_bin = tables.join(coverages.groupTuple())
    output_name = params.name != "" ? params.name : ref_fa.baseName
    cov_tables = createCoverageTables(output_name, covs_by_bin)

	// Create the segregation table
	grouped_tables = groupFilesBySizeAndType(cov_tables)    
    segregation_tables = createSegregationTables(grouped_tables)
}

