#!/bin/bash

# this script is an adaptation from Sashas code on creating coverage tables. Here we will only cover the first part, using bedtools coverage to create a coverage file per bam file with a given binned genome. In the following script, all individual coverage files are gonna be joined into one cov table.
# we need to create coverage tables from the original, non allele-specific bam files first, to ascertain covered and not covered bins of the genome, in multiple resolutions / bin or window sizes. 
# next we create coverage tables for each of the parental assignments of reads (p1, p2 / CAST,J129) and each of the different approaches, again in different resolutions. Together with the coverage file of the orig, we can figure out in subsequent steps if the allele-specific coverage is high enough to create meaningful allele-specific chromatin contact maps. 

# Sasha used bedtools makewindows to create equal bins along the ref genome. I did this separately, outside of the script, in order to avoid producing them over and over again for the different approaches.
# next bedtools coverage is used to calculate the coverage with options -a, where the previously produced binned genome is given and -b, with the bam files as input.  


# include USAGE message, if no arguments are given to the script
if [ "$#" -ne 3 ]; then
        echo "USAGE: $0 <input_dir> <parent> <output_dir>"
        exit
fi

input_dir=$1 #eg /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123_2kb_SNP_dist
parent=$2
output_dir=$3 #eg /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/coverage_files/subsampled_2kb_SNP_dist

# exit if no input is given
if [ ! -d "$input_dir" ]; then
    echo "directory $input_dir not found. Exiting."
    exit 1
fi

# if outputdirectory is not existing yet, make it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

if [ ! -d "$output_dir/log" ]; then
    mkdir -p "$output_dir/log"
fi

dist_cutoffs=`find $input_dir/bam_files_haplotype_assigned_reads_1/ -mindepth 1 -maxdepth 1 -type d -printf '%f '`

files=`find $input_dir/bam_files_haplotype_assigned_reads_1/100000.bp/$parent -name '*.bam'`
#files="/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123_2kb_SNP_dist/bam_files_haplotype_assigned_reads_3/100000.bp/CAST/f123_1NP_GAM180914_H09_1_S72_R1_001.fastq.gz.rmdup.trimmed_last_base.assigned_CAST.bam"
#files=/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/bam_files_haplotype_assigned_reads_1/100000.bp/J129/f123_1NP_GAM180207_H02_1_S16_R1_001.fastq.gz.rmdup.trimmed_last_base.assigned_J129.bam

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $input_dir/*.bam, writing to $output_dir at $(date +%x_%r)" >> $output_dir/run_log.txt

# we create small job scripts which are send off with sbatch:

job_directory=/fast/users/jmarkow_m/scratch/job

k=1
for f in $files; do
	j=$k$parent
	
        echo "${j}"
        job_file="${job_directory}/${j}_cov.job"

        b=`basename $f .bam`

        echo "#!/bin/bash
#SBATCH --job-name=${j}_cov.job
#SBATCH --output=$output_dir/log/${j}_cov.out
#SBATCH --error=$output_dir/log/${j}_cov.err
#SBATCH --mem=10G

source /fast/users/jmarkow_m/work/miniconda/etc/profile.d/conda.sh
conda activate PhasingGAM #H1

for n in {1..10};do
        echo \$n
        for dist in $dist_cutoffs; do
                echo \"dist \$dist\"
                for res in 50000 100000 200000; do
                        echo \"res \$res\"
                        bedtools coverage -a /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/F123_binned_genomes/F123_sorted_autosomes.window_size.\$res.bed -b $input_dir/bam_files_haplotype_assigned_reads_\$n/\$dist/$parent/$b.bam | awk '{print \$5}' > $output_dir/subsampled_F123_\$n/dist_cutoff_\$dist/res_\$res/$parent/$b.cov
                done
        done
done
" > $job_file
        sbatch $job_file
        ((k++))
done
