#!/bin/bash

# this script uses bedtools version "bedtools v2.26.0" to find the closest (incl overlapping) SNP for each read in a GAM sample

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

files=`find $indir -name '*.autosomes.bed'` # or which ever extension of files it is
#files=`find $indir -name '*.rmdup.autosomes.bed'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/GAM_beds/Repl1/H1hESC_3NP_R1_GAM_171214_H06_2_S112_R1_001.rmdup.autosomes.bed

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir/*.autosomes.bed, writing to $outdir at $(date +%x_%r)" > $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:


job_directory=/fast/groups/ag_schwarz/Projects/project-gam/src/.job

i=1
for f in $files; do

        echo "${i}"
        job_file="${job_directory}/${i}.job"

        b=`basename $f .bed`
        echo "#!/bin/bash
#SBATCH --job-name=${i}.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err
#SBATCH --partition=medium
#SBATCH --mem=16G

module load /fast/users/jmarkow_m/work/miniconda
#conda activate PhasingGAM
source activate PhasingGAM
# was thinking that "closest" against itself would show the closest read, but it will only return the identical read -.-
# use -k 2 option, to report the 2 closest reads. afterwards, go over them with awk to delete those lines were the identical read was reported
bedtools closest -d -t all -k 2 -a $f -b $f > $outdir/$b.closest_read_incl_same.bed" > $job_file
        sbatch $job_file
        ((i++))
done
