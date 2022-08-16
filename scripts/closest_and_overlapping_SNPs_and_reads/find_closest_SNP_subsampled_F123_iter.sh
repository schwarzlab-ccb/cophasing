#!/bin/bash

# this script uses bedtools version "bedtools v2.26.0" to find the closest (incl overlapping) SNP for each read in a GAM sample

# include USAGE message, if no arguments are given to the script
if [ "$#" -ne 3 ]; then
        echo "USAGE: $0 <input_dir_bed_files_GAM_reads> <input_dir_bed_files_sample_SNP_seeds> <output_dir>"
        exit
fi

indir_reads=$1 ## directory with input bed files of GAM reads
indir_SNPs=$2 ## directory with input bed files of seed SNPs, which have been observed in the respective GAM sample
outdir=$3 ## directory to store output

# exit if no input is given
if [ ! -d "$indir_reads" ]; then
    echo "directory $indir_reads not found. Exiting."
    exit 1
fi

if [ ! -d "$indir_SNPs" ]; then
    echo "directory $indir_SNPs not found. Exiting."
    exit 1
fi

# if outputdirectory is not existing yet, make it
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

if [ ! -d "$outdir/log" ]; then
    mkdir -p "$outdir/log"
fi

for m in {1..10};do
        mkdir $outdir/closest_SNP_subsampled_$m
done

read_files=`find $indir_reads -name '*.autosomes.bed'` # or which ever extension of files it is
#read_files=`find $indir_reads -name '*.autosomes.bed.sorted'` # or which ever extension of files it is
#read_files=/fast/groups/ag_schwarz/Projects/project-gam/data/F123_bam_files_1123_1NP_mm10_Sept_2019_last_base_trimmed_bed/f123_1NP_GAM171024_C01_1_S3_R1_001.fastq.gz.rmdup.trimmed_last_base.autosomes.bed.sorted
#read_files=/fast/groups/ag_schwarz/Projects/project-gam/data/F123_bam_files_1123_1NP_mm10_Sept_2019_last_base_trimmed_bed/f123_1NP_GAM180207_H02_1_S16_R1_001.fastq.gz.rmdup.trimmed_last_base.autosomes.bed.sorted

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir_reads/*.autosomes.bed.sorted and $indir_SNPs/*.autosomes.pileup.closest_SNP_seed.bed, writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

# we create small job scripts which are send off with sbatch:
job_directory=/fast/users/jmarkow_m/scratch/job

i=1
for f in $read_files; do

        echo "${i}"
        job_file="${job_directory}/${i}_closest_SNP.job"
        b=`basename $f .autosomes.bed`
#        b=`basename $f .autosomes.bed.sorted`

        echo "#!/bin/bash
#SBATCH --job-name=${i}_closest_SNP.job
#SBATCH --output=$outdir/log/${i}_closest_SNP.out
#SBATCH --error=$outdir/log/${i}_closest_SNP.err
#SBATCH --partition=medium
#SBATCH --mem=32G

source /fast/users/jmarkow_m/work/miniconda/etc/profile.d/conda.sh
conda activate H1 # PhasingGAM

sort -V -k1,1 -k2,2 $f > $f.sorted

for n in {1..10};do
	echo \$n
	indir_SNPs_dir=${indir_SNPs}/closest_SNP_seed_files_subsampled_\$n
	echo \$indir_SNPs_dir
	snp=\`find \$indir_SNPs_dir -name $b*\`
	echo \$snp
	echo bedtools closest -d -t all -g /fast/groups/ag_schwarz/Projects/project-gam/data/mm10_sorted_autosomes.chrom.sizes -a $f.sorted -b \$snp > $outdir/closest_SNP_subsampled_\$n/$b.closest_SNP.bed
#       echo bedtools closest -d -t all -g /fast/groups/ag_schwarz/Projects/project-gam/data/mm10_sorted_autosomes.chrom.sizes -a $f -b \$snp > $outdir/closest_SNP_subsampled_\$n/$b.closest_SNP.bed
        bedtools closest -d -t all -g /fast/groups/ag_schwarz/Projects/project-gam/data/mm10_sorted_autosomes.chrom.sizes -a $f.sorted -b \$snp > $outdir/closest_SNP_subsampled_\$n/$b.closest_SNP.bed
#	bedtools closest -d -t all -g /fast/groups/ag_schwarz/Projects/project-gam/data/mm10_sorted_autosomes.chrom.sizes -a $f -b \$snp > $outdir/closest_SNP_subsampled_\$n/$b.closest_SNP.bed
done
" > $job_file
        sbatch $job_file
        ((i++))
done
