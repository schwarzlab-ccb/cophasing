#!/bin/bash
# this script is an adaptation from Sashas code on creating coverage tables. This is the second part, combining the previously generated coverage files per sample per parent into coverage tables of the whole dataset. 
# IMPORTANT: we need to create coverage tables from the original, non allele-specific bam files first, to ascertain covered and not covered bins of the genome, in multiple resolutions / bin or window sizes. 
# here we create coverage tables for each of the parental assignments of reads (CAST,J129) and each of the different approaches, again in different resolutions. Together with the coverage file of the orig, we can figure out in subsequent steps if the allele-specific coverage is high enough to create meaningful allele-specific chromatin contact maps. 

input_dir=$1 # /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/coverage_files/subsampled_2kb_SNP_dist/
if [ ! -d "$input_dir/log" ]; then
    mkdir -p "$input_dir/log"
fi
if [ ! -d "$input_dir/cov_tables" ]; then
    mkdir -p "$input_dir/cov_tables"
fi

job_directory=/fast/users/jmarkow_m/scratch/job

k=1
for i in {1..10};do
        for dist in 0.bp 50000.bp 100000.bp 200000.bp basic; do
		for res in 50000 100000 200000; do
			for parent in CAST J129;do
				j=$k.$i.$dist.$res.$parent
				echo "${j}"
				job_file="${job_directory}/${j}.job"
				echo "#!/bin/bash
#SBATCH --job-name=${j}.job
#SBATCH --output=${input_dir}/log/${j}_covtab.out
#SBATCH --error=${input_dir}/log/${j}_covtab.err
#SBATCH --mem=10G

source /fast/users/jmarkow_m/work/miniconda/etc/profile.d/conda.sh
conda activate H1

# change into the directory where coverage files are stored 
cd $input_dir/subsampled_F123_$i/dist_cutoff_$dist/res_$res/$parent/
#rm header.txt
find *.table | xargs rm

# find all coverage files and paste them by column with the resolution genomic bins
covlist=\`find -type f -name '*.cov'\`
paste /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/F123_binned_genomes/F123_sorted_autosomes.window_size.$res.bed \$covlist > coverage_table_cophased_subsampled_F123_3NP_dist_cutoff_$dist.res_$res.$parent.iteration_$i.table

# then make and add a header including the sample names 
echo -e \"chrom\\tstart\\tstop\" \$covlist > header.txt
tr ' ' \\\t < header.txt | sponge header.txt
cat header.txt coverage_table_cophased_subsampled_F123_3NP_dist_cutoff_$dist.res_$res.$parent.iteration_$i.table | sponge coverage_table_cophased_subsampled_F123_3NP_dist_cutoff_$dist.res_$res.$parent.iteration_$i.table
cp coverage_table_cophased_subsampled_F123_3NP_dist_cutoff_$dist.res_$res.$parent.iteration_$i.table $input_dir/cov_tables/
" > $job_file
				sbatch $job_file
				((k++))
			done
		done
	done
done
