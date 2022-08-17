nextflow.enable.dsl=2

params.out_dir = "out"
params.debug_out = "out"
scripts_folder = "${projectDir}/scripts"

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
    
    output:
    tuple path("*.50000.bed"), path("*.100000.bed"), path("*.200000.bed")

    script:
    """
    bedtools makewindows -g $genome_sizes -w 50000 > ${genome_sizes.baseName}.50000.bed
    bedtools makewindows -g $genome_sizes -w 100000 > ${genome_sizes.baseName}.100000.bed
    bedtools makewindows -g $genome_sizes -w 200000 > ${genome_sizes.baseName}.200000.bed
    """
}

process calcCoverage
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    tuple val(name), path(bed)
    tuple path(binned_genome_50000), path(binned_genome_100000), path(binned_genome_200000)
    
    output:
    tuple path("*.50000.cov"), path("*.100000.cov"), path("*.200000.cov")
    
    script:
    """
    bedtools coverage -a ${binned_genome_50000} -b $bed | awk \"{print \\\$5}\" > ${name}.50000.cov
    bedtools coverage -a ${binned_genome_100000} -b $bed | awk \"{print \\\$5}\" > ${name}.100000.cov
    bedtools coverage -a ${binned_genome_200000} -b $bed | awk \"{print \\\$5}\" > ${name}.200000.cov
    """
}


process createCoverageTables
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""
    
    input:
    tuple val(name), path(bed)
    tuple path(coverage_50000), path(coverage_100000), path(coverage_200000)
    
    output:
    path("${name}.table")
    
    script:
    """
    echo "$coverage_50000 $coverage_100000 $coverage_200000" > ${name}.table
    """
    // """
    // # find all coverage files and paste them by column with the resolution genomic bins
    // covlist=\`find -type f -name '*.cov'\`
    // paste /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/F123_binned_genomes/F123_sorted_autosomes.window_size.$res.bed \$covlist > coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table

    // # then make and add a header including the sample names 
    // echo -e \"chrom\\tstart\\tstop\" \$covlist > header.txt
    // tr ' ' \\\t < header.txt | sponge header.txt
    // cat header.txt coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table | sponge coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table
    // cp coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table $output_dir
    // """
}

process pythonTest
{
    output:
    stdout

    script:
    """
    $scripts_folder/test.py
    """
}

process RTest
{
    output:
    stdout

    script:
    """
    $scripts_folder/test.py
    """
}

workflow 
{
    bam = Channel.fromFilePairs(params.bam, size: 1)
    vcf = Channel.fromFilePairs(params.vcf, size: 1)   
    ref_fa = file(params.fa)
    
    pythonTest().view()
    bed = convertBamToBed(bam)
    genome_size = getGenomeSizes(ref_fa)
    genome_bins = binGenome(genome_size)
    coverage = calcCoverage(bed, genome_bins).view()
    cov_table = createCoverageTables(bed, coverage)

    // pileup = bcftoolsPileup(ref_fa, bam.join(vcf))
    // seeds = seedsPerSample(vcf.join(pileup))
}