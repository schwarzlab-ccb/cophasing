nextflow.enable.dsl=2

scripts_folder = "${projectDir}/scripts"

process pythonTest
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
    bam.map{name, file -> name}.view()
    if (params.vcf) 
    {
        vcf_files = Channel.fromPath(params.vcf)
        comb = bam.combine(vcf_files)
    }
    else 
    {
        vcf_filename = params.bam - ".bam" + ".vcf"
        vcf_pairs = Channel.fromFilePairs(vcf_filename, size: 1)
        comb = bam.join(vcf_pairs)
    }
    bam_vcf_pairs = comb.flatten().collate(3)
    // vcf_files = vcf.map{a, b -> b}.toList()
    // comb = bam.combine(vcf_files)
    // bam_vcf = bam.join(vcf).ifEmpty(comb)
    // bam_vcf.view()
    // pythonTest().view()
}