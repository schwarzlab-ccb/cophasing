#!/bin/bash

input=$1
input_dir=/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/coverage_files/subsampled_${input}kb_SNP_dist/subsampled_F123_1
output_dir=/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/coverage_files/subsampled_${input}kb_SNP_dist/cov_tables_iter1_full
job_directory=/fast/users/jmarkow_m/scratch/job

# if outputdirectory is not existing yet, make it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

if [ ! -d "$output_dir/log" ]; then
    mkdir -p "$output_dir/log"
fi

k=1
for dist in 0.bp 50000.bp 100000.bp 200000.bp basic; do
	for res in 50000 100000 200000; do
		for parent in CAST J129;do
			j=$k.$dist.$res.$parent
			echo "${j}"
			job_file="${job_directory}/${j}.job"
			echo "#!/bin/bash
#SBATCH --job-name=${j}.job
#SBATCH --output=${output_dir}/log/${j}.out
#SBATCH --error=${output_dir}/log/${j}.err
#SBATCH --mem=10G

source /fast/users/jmarkow_m/work/miniconda/etc/profile.d/conda.sh
conda activate H1

# change into the directory where coverage files are stored 
cd $input_dir/dist_cutoff_$dist/res_$res/$parent/

# find all coverage files and paste them by column with the resolution genomic bins
covlist=\`find -type f -name '*.cov'\`
paste /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/F123_binned_genomes/F123_sorted_autosomes.window_size.$res.bed \$covlist > coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table

# then make and add a header including the sample names 
echo -e \"chrom\\tstart\\tstop\" \$covlist > header.txt
tr ' ' \\\t < header.txt | sponge header.txt
cat header.txt coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table | sponge coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table
cp coverage_table_cophased_subsampled_F123_dist_cutoff_$dist.res_$res.$parent.iteration_1.210521.table $output_dir
" > $job_file
			sbatch $job_file
			((k++))
		done
	done
done
