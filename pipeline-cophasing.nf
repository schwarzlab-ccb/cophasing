nextflow.enable.dsl=2

// Default parameter values
params.bins = [50000, 100000, 200000] // Bin sizes to use for analysis, each bin size will be analyzed separately
params.out = "out" // Output directory containing the results
params.debug_out = "" // If set, will output intermediate files to this directory
params.cutoff = 0 // Maximum distance between a read and a variant to be considered for analysis
params.min_depth = 1 // Minimum required read depth per variant to be considered for analysis
params.min_ratio = 5 // Main base must be at least {min_ratio} times more often represented than the remaining bases (or the only one represented)
params.name = "" // Output file name, will default to the name of the FA file if not set


process bgzip 
{    
    input:
    path(fa)

    output:
    path("*.gz")
    
    """
    bgzip $fa 
    """  
}


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
    obs_ref = "FORMAT/AD[0:0]>=FORMAT/AD[0:1]" // OBSERVED REFERENCE, note that equal should not ever happen, even though its permitted here
    obs_alt = "FORMAT/AD[0:0]<FORMAT/AD[0:1]" // OBSERVED ALTERNATIVE
    // Below we print out the position twice, awk is used to increment the second position by 1 (end)
    """
    touch ${name}.vcf.bed    
    bcftools query $vcf -i '(GT="1|0" && $obs_ref) || (GT="0|1" && $obs_alt) ' -f '%CHROM %POS %POS hap1\n' | awk '{ \$3 = \$3 + 1 } 1' | sed 's/ /\t/g' >> ${name}.vcf.bed   
    bcftools query $vcf -i '(GT="1|0" && $obs_alt) || (GT="0|1" && $obs_ref) ' -f '%CHROM %POS %POS hap2\n' | awk '{ \$3 = \$3 + 1 } 1' | sed 's/ /\t/g' >> ${name}.vcf.bed      
    sort ${name}.vcf.bed -k1,1 -k2,2nn -o ${name}.vcf.bed
    """
}

process bcftoolsPileup 
{
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(ref_genome)
    tuple path(vcf), path(tb)
    tuple val(name), path(bam)
    
    output:
    tuple val(name), path("${name}.pileup.vcf")

    script:
    """
    bcftools mpileup -f $ref_genome -T $vcf -a FORMAT/DP,FORMAT/AD -O v $bam | bcftools sort > ${name}.pileup.vcf
    """
}

// If ALT is not known in pileup, take it from VCF. Annotate only where ALT matches between pileup and VCF
process bcftoolsAnnotate
{
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple path(vcf), path(tb)
    tuple val(name), path(pileup)
    
    output:
    tuple val(name), path("${name}.genotyped.vcf")

    script:
    // The sample name is not the same between the Pileup and the reference, needs to be matched, hence the samples.txt file
    """
    echo `bcftools query -l $vcf` `bcftools query -l $pileup` > samples.txt
    if [ `wc -l < samples.txt` != 1 ]; then echo "there must be exactly one sample in the VCF $vcf"; exit 1; fi;
    awk -F'\t' 'BEGIN {OFS="\t"} {split(\$5, a, ","); \$5 = a[1]; print \$0}' $pileup > ${name}.cut.vcf
    bgzip ${name}.cut.vcf
    tabix ${name}.cut.vcf.gz
    bcftools annotate -a $vcf -c ALT ${name}.cut.vcf.gz -i "FORMAT/AD[0:1]<=0" -k -S samples.txt > ${name}.fill.vcf
    bgzip ${name}.fill.vcf
    tabix ${name}.fill.vcf.gz    
    bcftools annotate -a $vcf -c FORMAT/GT ${name}.fill.vcf.gz -S samples.txt > ${name}.annotated.vcf
    bcftools view -i 'FORMAT/GT!="."' ${name}.annotated.vcf > ${name}.genotyped.vcf
    """
}

// Keep only sites where are reads are allocated to one allele
process filterSites {    
    publishDir "${params.debug_out}/${task.process}", mode: "copy", enabled: params.debug_out != ""

    input:
    tuple val(name), path(vcf)
    
    output:
    tuple val(name), path("${name}.mono")
    
    // Filter only if the read depth is above the minimum and the most represented base occurs at least params.max_ratio more often than the rest of the bases
    script:
    filter_four = "(FORMAT/AD[0:3] > 0)"
    filter_three = "(FORMAT/AD[0:2] > 0) && (FORMAT/AD[0:0] > 0)"
    filter_two = "(FORMAT/AD[0:0] > 0) && (FORMAT/AD[0:1] > 0)"
    filter_depth = "FORMAT/DP[0:0] < $params.min_depth"
    """
    bcftools filter -e "($filter_depth) || ($filter_two) || ($filter_three) || ($filter_four)" $vcf > ${name}.mono
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

// Function to check if a variable is an integer
def isInteger(value) {
    if (value == null) {
        return false
    }
    if (value instanceof Integer) {
        return value >= 0
    }
    else {
        return false
    }
}

workflow 
{
    // Asserts
    if (params.debug_out != "") {
        def outputFolder = file(params.debug_out)
        if (!outputFolder.exists()) {
            outputFolder.mkdirs()
        }
    }

    if (!isInteger(params.min_ratio)) {
        error "Provided min_ratio must be a non-negative integer, is ${params.min_ratio}."
        System.exit(1)
    }
    if (!isInteger(params.min_depth)) {
        error "Provided min_depth must be a non-negative integer, is ${params.min_depth}."
        System.exit(2)
    }
    if (!isInteger(params.cutoff)) {
        error "Provided cutoff must be a non-negative integer, is ${params.cutoff}."
        System.exit(3)
    }

    // Inputs
    bam = Channel.fromFilePairs(params.bam, size: 1)
    ref_fa = file(params.fa)
    fa_is_zipped = params.fa.endsWith(".gz")
    fasta = fa_is_zipped ? ref_fa : bgzip(ref_fa)
    vcf = file(params.vcf)
    bin_sizes = Channel.from(params.bins)    

    // Create and filter pileup to obtain phased variant sites observed in the reads
    filtered_vcf = filterUnphased(vcf)
    pileup = bcftoolsPileup(fasta, filtered_vcf, bam)
    annotated_pileup = bcftoolsAnnotate(filtered_vcf, pileup)
    filtered_pileup = filterSites(annotated_pileup)
    vcf_beds = convertVcfToBed(filtered_pileup)

    // Calculate bins
    genome_size = getGenomeSizes(fasta)
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

