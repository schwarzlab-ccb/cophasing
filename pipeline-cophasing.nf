nextflow.enable.dsl=2

params.out_dir = "out"
params.debug_out = "out"
scripts_folder = "./scripts/"

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

process bcftoolsMpileup 
{
    publishDir "${params.debug_out}", mode: "copy", enabled: params.debug_out != ""

    input:
    path(ref_genome)
    tuple val(name), path(bam), path(vcf)
    
    output:
    tuple val(name), path("${name}.pileup.vcf")

    script:
    """
    bcftools mpileup -f $ref_genome -T $vcf -a FORMAT/AD,INFO/AD -O v $bam > ${name}.pileup.vcf 
    """
}

workflow 
{
    bam = Channel.fromFilePairs(params.bam, size: 1)
    vcf = Channel.fromFilePairs(params.vcf, size: 1)   
    ref_fa = file(params.fa)
    ref_bed = file(params.bed)

    bed = convertBamToBed(bam)
    pileup = bcftoolsMpileup(ref_fa, bam.join(vcf))

}