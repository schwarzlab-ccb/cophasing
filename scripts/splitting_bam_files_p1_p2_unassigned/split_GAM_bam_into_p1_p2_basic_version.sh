
#!/bin/bash

# this script uses grep to find and picards FilterSamReads to extract the reads assigned to either parental chromosome p1 or parental chromosome p2 based on the observed haplotype of their closest SNP. 
# the base of this assignment is (1) the pileup of the reads in the GAM sample at known H1 SNP positions (liftover from dixon). From this pileup, we know which SNPS are observed in this sample and which haplotype / parental chromosome was covered. We used these observed SNPs per sample to create (2) a file reporting the distance of each read in the sample to its closest SNP
# now we use the information about the observed haplotype of the covered SNPs (sample specific) and the information about the closest SNP to each read to assign the reads to the haplotype of their closest SNP
# in the first version, as a base to improve from, each read will be assigned to the haplotype of the closest SNP, no reads will be unassigned (except in the unlikely case that there are reads from the chromosome, but no single SNP covered) 
# in the next versions we will introduce distance cutoffs and more complex strategies where more factors are included in the decision making


# include USAGE message, if no arguments are given to the script
if [ "$#" -ne 3 ]; then
        echo "USAGE: $0 <input_bam_dir> <input_closestSNP_dir> <output_dir>"
        exit
fi

indir_bam=$1 ## directory with input bam files of GAM experiments
indir_closestSNP=$2 ## directory with input bed files of reads and their closest SNP
outdir=$3 ## directory to store output

tempdir=/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/splitting_bams_read_IDs_temps

# exit if no input is given
if [ ! -d "$indir_bam" ]; then
    echo "directory $indir_bam not found. Exiting."
    exit 1
fi

if [ ! -d "$indir_closestSNP" ]; then
    echo "directory $indir_closestSNP not found. Exiting."
    exit 1
fi

# if outputdirectory is not existing yet, make it
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

if [ ! -d "$outdir/log" ]; then
    mkdir -p "$outdir/log"
fi

files=`find $indir_bam -name '*.autosomes.bam'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60_autosomes/Repl1/H1hESC_3NP_R1_GAM_171121_A03_2_S113_R1_001.rmdup.autosomes.bam

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir_bam/*.rmdup.autosomes.bam, with $indir_closestSNP/*.autosomes.closest_SNP.bed writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:


job_directory=/fast/groups/ag_schwarz/Projects/project-gam/src/.job

i=1
for f in $files; do

        echo "${i}"
        job_file="${job_directory}/${i}.job"

	b=`basename $f .bam`
	closest_snp=`find $indir_closestSNP -name $b*`

        echo "#!/bin/bash
#SBATCH --job-name=${i}.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err
#SBATCH --partition=medium
#SBATCH --mem=64G

module load /fast/users/jmarkow_m/work/miniconda
source activate H1
grep -F 'p1' $closest_snp | cut -f4 > $tempdir/$b.p1
grep -F 'p2' $closest_snp | cut -f4 > $tempdir/$b.p2
java -jar /fast/users/jmarkow_m/work/miniconda/envs/H1/share/picard-2.23.8-0/picard.jar FilterSamReads I=$f O=$outdir/p1/$b.assigned_p1.bam READ_LIST_FILE=$tempdir/$b.p1 FILTER=includeReadList
java -jar /fast/users/jmarkow_m/work/miniconda/envs/H1/share/picard-2.23.8-0/picard.jar FilterSamReads I=$f O=$outdir/p2/$b.assigned_p2.bam READ_LIST_FILE=$tempdir/$b.p2 FILTER=includeReadList
rm $tempdir/$b.p1
rm $tempdir/$b.p2" > $job_file
        sbatch $job_file
        ((i++))
done

