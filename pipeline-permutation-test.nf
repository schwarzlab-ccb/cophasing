nextflow.enable.dsl=2

// Every process stages its inputs and declares the files it produces, so Nextflow
// hashes real content and `-resume` reflects the state of the results rather than
// of the work directory alone. Results reach params.output_dir through publishDir.
// Each script is handed `--output_dir .`, i.e. its own task directory, and the
// subdirectory layout it writes there is what gets published.

process CurateSegregationTables {
    tag "curation"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path segregation_tables

    output:
    path "01_curated_segregation_tables", emit: curated
    path "01_curation.log"

    script:
    """
    python ${projectDir}/permutation_test_scripts/01_segregation_table_curation_by_mean_WDF.py \\
        --input_dir . \\
        --output_dir . \\
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
    path curated

    output:
    tuple val(chr), path("permutation_test_results_${chr}_multiprocessing.pkl"), emit: pkl
    path "${chr}.log"

    script:
    """
    python ${projectDir}/permutation_test_scripts/02_permutation_test_chromosome_level.py \\
        --output_dir . \\
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

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path curated
    path perm_results, stageAs: '02_permutation_test/*'

    output:
    path "03_thresholds", emit: thresholds
    path "03_permutation_test_threshold_identification_contact_ratio.log"

    script:
    """
    python ${projectDir}/permutation_test_scripts/03_permutation_test_threshold_identification_contact_ratio.py \\
        --output_dir . \\
        --chromosomes '${params.chromosomes}' \\
        --gaussian_kernel_size ${params.gaussian_kernel_size} \\
        > 03_permutation_test_threshold_identification_contact_ratio.log 2>&1
    """
}

process GenomicWindowsBed {
    tag "genomic_windows"

    publishDir "${params.output_dir}/04_cool_files", mode: 'copy'

    output:
    path "res_${params.resolution}_genomic_windows.bed"

    script:
    """
    python ${projectDir}/permutation_test_scripts/make_genomic_windows.py \\
        --resolution ${params.resolution} \\
        --out res_${params.resolution}_genomic_windows.bed
    """
}

process GenerateCoolFromPerm {
    tag { chr }

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(chr), path(perm_result, stageAs: '02_permutation_test/*')
    path curated
    path thresholds
    path genomic_windows

    output:
    path "04_cool_files/cool_files_npmi_permutation_test_${chr}"
    path "04_cool_file_${chr}.log"

    script:
    """
    python ${projectDir}/permutation_test_scripts/04_generate_cool_file_permutation_test_results.py \\
        --chr ${chr} \\
        --output_dir . \\
        --cutoff ${params.cutoff} \\
        --resolution ${params.resolution} \\
        --gaussian_kernel_size ${params.gaussian_kernel_size} \\
        --pseudocount ${params.pseudocount} \\
        --bed ${genomic_windows} \\
        > 04_cool_file_${chr}.log 2>&1
    """
}

workflow {
    if (!file(params.input_dir).isAbsolute())
        error "input_dir must be an absolute path, got: ${params.input_dir}"
    if (!file(params.output_dir).isAbsolute())
        error "output_dir must be an absolute path, got: ${params.output_dir}"

    def chrs = params.chromosomes.tokenize(',')
    def chr_channel = Channel.fromList(chrs)

    // Staged by content, so editing or replacing a segregation table reruns curation
    seg_tables = Channel.fromPath(
        "${params.input_dir}/*.${params.resolution}.*.segregation*", checkIfExists: true).collect()

    curated = CurateSegregationTables(seg_tables).curated

    perm_chr = PermutationTestChromosomeLevel(chr_channel, curated).pkl

    all_pkl_files = perm_chr.map { _chr, pkl -> pkl }.collect()
    thresholds = IdentifyThresholds(curated, all_pkl_files).thresholds

    genomic_windows = GenomicWindowsBed()

    GenerateCoolFromPerm(perm_chr, curated, thresholds, genomic_windows)
}
