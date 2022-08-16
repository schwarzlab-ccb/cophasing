#!/bin/bash

# this script uses samtools version 1.9 to filter the GAM bams for autosomes, aka remove chromosomes X and Y

# include USAGE message, if no arguments are given to the script
if [ "$#" -ne 2 ]; then
        echo "USAGE: $0 <input_dir> <output_dir>"
        exit
fi

indir=$1 ## directory with input files
outdir=$2 ## directory to store output

# exit if no input is given
if [ ! -d "$indir" ]; then
    echo "directory $indir not found. Exiting."
    exit 1
fi

# if outputdirectory is not existing yet, make it
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

if [ ! -d "$outdir/log" ]; then
    mkdir -p "$outdir/log"
fi

files=`find $indir -name '*.bam'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/results/trimmed_F123_3NP/f123_3NP_GAM161008_B7_S50_R1_001.fastq.gz.rmdup.trimmed_last_base.bam

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir/*.rmdup.autosomes.bed, writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:


job_directory=/fast/users/jmarkow_m/scratch/job

i=1
for f in $files; do

        echo "${i}"
        job_file="${job_directory}/${i}_filter_autosomes.job"

        b=`basename $f .bam`
        echo "#!/bin/bash
#SBATCH --job-name=${i}_filter_autosomes.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err
#SBATCH --partition=medium
#SBATCH --mem=16G

source /fast/users/jmarkow_m/work/miniconda/etc/profile.d/conda.sh
conda activate PhasingGAM
samtools view -b -h -L /fast/groups/ag_schwarz/Projects/project-gam/data/mm10_chromsize.bed $f > $outdir/$b.autosomes.bam" > $job_file
        sbatch $job_file
        ((i++))
done
