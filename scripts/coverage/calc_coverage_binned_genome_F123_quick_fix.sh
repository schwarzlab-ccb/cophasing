#
#!/bin/bash

# this script is an adaptation from Sashas code on creating coverage tables. Here we will only cover the first part, using bedtools coverage to create a coverage file per bam file with a given binned genome. In the following script, all individual coverage files are gonna be joined into one cov table.
# IMPORTANT: at this point I will no longer divide between Repl1 and Repl2, but just handle them together, thus I need two input directories
# we need to create coverage tables from the original, non allele-specific bam files first, to ascertain covered and not covered bins of the genome, in multiple resolutions / bin or window sizes. 
# next we create coverage tables for each of the parental assignments of reads (p1, p2) and each of the different approaches, again in different resolutions. Together with the coverage file of the orig, we can figure out in subsequent steps if the allele-specific coverage is high enough to create meaningful allele-specific chromatin contact maps. 

# Sasha used bedtools makewindows to create equal bins along the ref genome. I did this separately, outside of the script, in order to avoid producing them over and over again for the different approaches.
# next bedtools coverage is used to calculate the coverage with options -a, where the previously produced binned genome is given and -b, with the bam files as input.  


# include USAGE message, if no arguments are given to the script
#if [ "$#" -ne 3 ]; then
#        echo "USAGE: $0 <input_bam_dir> <binned_genome> <output_dir>"
#        exit
#fi

indir_bam=$1 ## directory with input bam files of GAM experiments
#binned_genome=$2 ## binned genome with window / bin size showing the resolution of final chromatin contact map
outdir=$2 ## directory to store output

# exit if no input is given
if [ ! -d "$indir_bam" ]; then
    echo "directory $indir_bam not found. Exiting."
    exit 1
fi

# if outputdirectory is not existing yet, make it
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

if [ ! -d "$outdir/log" ]; then
    mkdir -p "$outdir/log"
fi

files=`find $indir_bam -name '*autosomes.bam'` 
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60_autosomes/Repl1/H1hESC_3NP_R1_GAM_171121_A03_2_S113_R1_001.rmdup.autosomes.bam

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir_bam/*.autosomes.bam, with $binned_genome writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:


job_directory=/fast/groups/ag_schwarz/Projects/project-gam/src/.job

i=1
for f in $files; do

        echo "${i}"
        job_file="${job_directory}/${i}.job"

        b=`basename $f .bam`

        echo "#!/bin/bash
#SBATCH --job-name=${i}.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err
#SBATCH --mem=10G

source /fast/users/jmarkow_m/work/miniconda/etc/profile.d/conda.sh
conda activate H1

bedtools coverage -a /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/F123_binned_genomes/F123_sorted_autosomes.window_size.100000.bed -b $f | awk '{print \$5}' > /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/coverage_files/orig/res_100kb/$b.cov
bedtools coverage -a /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/F123_binned_genomes/F123_sorted_autosomes.window_size.200000.bed -b $f | awk '{print \$5}' > /fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/coverage_files/orig/res_200kb/$b.cov" > $job_file
        sbatch $job_file
        ((i++))
done
