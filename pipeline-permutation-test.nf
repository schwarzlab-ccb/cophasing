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
    cpus params.num_workers    

    publishDir "${params.output_dir}/02_permutation_test", mode: 'copy'

    input:
    val chr
    path _ready  // curation log path — cache invalidated when step 1 reruns

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
        --num_workers ${params.num_workers} \\
        --pseudocount ${params.pseudocount} \\
        --out permutation_test_results_${chr}_multiprocessing.pkl \\
        > ${chr}.log 2>&1
    """
}

process IdentifyThresholds {
    tag "threshold_identification"

    input:
    path perm_results  // actual pkl files — content hash changes when step 2 reruns

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
    tuple val(chr), path(_ready)  // threshold log path — cache invalidated when step 3 reruns

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

    perm_chr = PermutationTestChromosomeLevel(chr_channel, curation_done)

    // Collect actual pkl files so their content hashes drive IdentifyThresholds cache
    all_pkl_files = perm_chr.map { _chr, pkl, _log -> pkl }.collect()
    thresh_done_log = IdentifyThresholds(all_pkl_files)

    // combine() pairs each chromosome with the threshold log — works for any number of chromosomes
    GenerateCoolFromPerm(chr_channel.combine(thresh_done_log))
}
