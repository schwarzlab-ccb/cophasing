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


read_files=`find $indir_reads -name '*.autosomes.bed'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/GAM_beds/Repl1/H1hESC_3NP_R1_GAM_171214_H06_2_S112_R1_001.rmdup.autosomes.bed

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir_reads/*.autosomes.bed and $indir_SNPs/*.autosomes.pileup.closest_SNP_seed.bed, writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:


job_directory=/fast/groups/ag_schwarz/Projects/project-gam/src/.job

i=1
for f in $read_files; do

        echo "${i}"
        job_file="${job_directory}/${i}.job"
        b=`basename $f .bed`
	snp=`find $indir_SNPs -name $b*`

        echo "#!/bin/bash
#SBATCH --job-name=${i}.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err
#SBATCH --partition=medium
#SBATCH --mem=16G

module load /fast/users/jmarkow_m/work/miniconda
source activate PhasingGAM

bedtools closest -d -t all -a $f -b $snp > $outdir/$b.closest_SNP.bed" > $job_file
        sbatch $job_file
        ((i++))
done
