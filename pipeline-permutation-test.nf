nextflow.enable.dsl=2

process CurateSegregationTables {
    tag "curation"

    output:
    path "01_curation.log"

    script:
    """
    python ${projectDir}/permutation_test_scripts/01_segregation_table_curation_by_mean_WDF.py \\
        --input_dir '${params.input_dir}' \\
        --output_dir '${params.output_dir}' \\
        --cutoff ${params.cutoff} \\
        --resolution ${params.resolution} \\
        --high_factor ${params.high_factor} \\
        --low_factor ${params.low_factor} \\
        > 01_curation.log 2>&1
    """
}

process PermutationTestChromosomeLevel {
    tag { chr }

    publishDir "${params.output_dir}/02_permutation_test", mode: 'copy'

    input:
    val chr
    val _ready

    output:
    tuple val(chr), path("permutation_test_results_${chr}_multiprocessing.pkl"), path("${chr}.log")


    script:
    """
    python ${projectDir}/permutation_test_scripts/02_permutation_test_chromosome_level.py \\
        --output_dir '${params.output_dir}' \\
        --cutoff ${params.cutoff} \\
        --resolution ${params.resolution} \\
        --chr ${chr} \\
        --num_perm ${params.num_perm} \\
        --pseudocount ${params.pseudocount} \\
        --out permutation_test_results_${chr}_multiprocessing.pkl \\
        > ${chr}.log 2>&1
    """
}

process IdentifyThresholds {
    tag "threshold_identification"

    input:
    val n_results

    output:
    path "03_permutation_test_threshold_identification_contact_ratio.log"

    script:
    """
    python ${projectDir}/permutation_test_scripts/03_permutation_test_threshold_identification_contact_ratio.py \\
        --output_dir '${params.output_dir}' \\
        --chromosomes '${params.chromosomes}' \\
        --gaussian_kernel_size ${params.gaussian_kernel_size} \\
        > 03_permutation_test_threshold_identification_contact_ratio.log 2>&1
    """
}

process GenerateCoolFromPerm {
    tag { chr }

    input:
    val(chr)
    val _ready

    output:
    path "04_cool_file_${chr}.log"

    script:
    """
    python ${projectDir}/permutation_test_scripts/04_generate_cool_file_permutation_test_results.py \\
        --chr ${chr} \\
        --output_dir '${params.output_dir}' \\
        --cutoff ${params.cutoff} \\
        --resolution ${params.resolution} \\
        --gaussian_kernel_size ${params.gaussian_kernel_size} \\
        > 04_cool_file_${chr}.log 2>&1
    """
}

workflow {
    def chrs = params.chromosomes.tokenize(',')
    def chr_channel = Channel.fromList(chrs)

    curation_done = CurateSegregationTables()

    perm_chr = PermutationTestChromosomeLevel(chr_channel, curation_done.map { true })

    trigger_perm_complete = perm_chr.count()

    thresh_done_log = IdentifyThresholds(trigger_perm_complete)

    thresh_ready = thresh_done_log.map { true }
    GenerateCoolFromPerm(chr_channel, thresh_ready)
}
