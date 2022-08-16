#!/bin/bash

# this script uses awk to curate the intermediate files of finding the distance between reads and their closest SNP, also identifying overlapping SNPs
# here we will remove all ties, so there is only one closest SNP per read (and thus only one shortest distance per read)

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

files=`find $indir -name '*.closest_SNP.bed'` # or which ever extension of files it is


# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir/*.closest_SNP.bed, writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:

job_directory=/fast/groups/ag_schwarz/Projects/project-gam/src/.job

i=1
for f in $files; do

        echo "${i}"
        job_file="${job_directory}/${i}.job"

	b=`basename $f .closest_SNP.bed`
        echo "#!/bin/bash
#SBATCH --job-name=${i}.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err
#SBATCH --partition=medium
#SBATCH --mem=16G

module load /fast/users/jmarkow_m/work/miniconda
#conda activate PhasingGAM
source activate PhasingGAM
# unique the first ID column to remove ties, so we are left with one closest SNP per input read 
awk '!_[\$4]++' $f > $outdir/$b.closest_SNP.dedup.bed" > $job_file
        sbatch $job_file
        ((i++))
done

