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
    ref_bed = file(params.bed)
    pythonTest().view()

    // bed = convertBamToBed(bam)
    // pileup = bcftoolsPileup(ref_fa, bam.join(vcf))
    // seeds = seedsPerSample(vcf.join(pileup))
}