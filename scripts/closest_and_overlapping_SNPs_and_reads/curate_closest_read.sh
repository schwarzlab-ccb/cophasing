#!/bin/bash

# this script uses awk to curate the intermediate files of finding the distance between closest reads, also identifying overlapping reads
# in the output of the previous step, we find the read itself as being closest, then the second closest read and then all ties
# for now, the closest read is always the reads itself, we will remove this (and safe the intermediate results)
# then we will remove all ties, so there is only one closest read per read (and thus only one shortest distance per read)
# and we remove all closest readpairs in which the coordinates of the second read are smaller than those of the first, so we do not double count the distances. (e.g read 1 is closest to read 2 with 3 bp difference, and read 2 is closest to read 1 with , again, 3bp difference)

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

files=`find $indir -name '*.closest_read_incl_same.bed'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/closest_reads_intermediate_files/Repl1/H1hESC_3NP_R1_GAM_171121_A01_1_S1_R1_001.rmdup.autosomes.closest_read_incl_same.bed


# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir/*.closest_read_incl_same.bed, writing to $outdir at $(date +%x_%r)" > $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:

job_directory=/fast/groups/ag_schwarz/Projects/project-gam/src/.job

i=1
for f in $files; do

        echo "${i}"
        job_file="${job_directory}/${i}.job"

	b=`basename $f .closest_read_incl_same.bed`
        echo "#!/bin/bash
#SBATCH --job-name=${i}.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err
#SBATCH --partition=medium
#SBATCH --mem=16G

module load /fast/users/jmarkow_m/work/miniconda
#conda activate PhasingGAM
source activate PhasingGAM
# remove lines where the read IDs are identical, then unique the first ID column to remove ties, so we are left with one closest read per input read 
awk '\$4 != \$10' $f | awk '!_[\$4]++' > $outdir/$b.closest_read_excl_same.dedup.bed" > $job_file
#awk '\$4 != \$10' $f | awk '\$2 <= \$8' > $outdir/$b.closest_read_excl_same_smaller_pos.bed" > $job_file
        sbatch $job_file
        ((i++))
done

# too slow:
#i=1
#for f in $files; do

#	echo "${i}"

#        b=`basename $f .closest_read_incl_same.bed`
	
#	awk '$4 != $10' $f | awk '$2 <= $8' > $outdir/$b.closest_read_excl_same_smaller_pos.bed
#	((i++))
#done
