#!/bin/bash

# this script uses multiple instances of grep to find and extract the reads assigned to either parental chromosome p1 or parental chromosome p2 or can not be assigned. 
# the base of this assignment is (1) the pileup of the reads in the GAM sample at known H1 SNP positions (liftover from dixon). From this pileup, we know which SNPS are observed in this sample and which haplotype / parental chromosome was covered. We used these observed SNPs per sample to create (2) a file reporting the distance of each read in the sample to its closest SNP
# now we use the information about the observed haplotype of the covered SNPs (sample specific) and the information about the closest SNP to each read to assign the reads to the haplotype of their closest SNP
# in the first version, as a base to improve from, each read was assigned to the haplotype of the closest SNP, no reads were unassigned (except in the unlikely case that there are reads from the chromosome, but no single SNP covered) 
# in this next versions I introduce hard distance cutoffs, which can be given to the script as an additional argument
# more complex strategies where more factors are included in the decision making will be introduced in following more advanced and complex scripts. 


# include USAGE message, if no arguments are given to the script
if [ "$#" -ne 4 ]; then
        echo "USAGE: $0 <input_bam_dir> <input_closestSNP_dir> <distance_cutoff> <output_dir>"
        exit
fi

indir_bam=$1 ## directory with input bam files of GAM experiments
indir_closestSNP=$2 ## directory with input bed files of reads and their closest SNP
dist_cutoff=$3 ## distance cut off in base pairs. reads with a larger distance to the nearest SNP will not be phased / assigned to a parental haplotype but remain unassigned
outdir=$4 ## directory to store output

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

# check that the distance cutoff is given as an integer
re='^[0-9]+$'
if ! [[ $dist_cutoff =~ $re ]] ; then
   echo "given distance cutoff must be an integer, please provide the distance cutoff in bp. Exiting." >&2; exit 1
fi

# if outputdirectory is not existing yet, make it
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

if [ ! -d "$outdir/log" ]; then
    mkdir -p "$outdir/log"
fi

mkdir -p "$outdir/$dist_cutoff.bp/p1"
mkdir -p "$outdir/$dist_cutoff.bp/p2"

files=`find $indir_bam -name '*.autosomes.bam'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60_autosomes/Repl1/H1hESC_3NP_R1_GAM_171121_A03_2_S113_R1_001.rmdup.autosomes.bam


# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir_bam/*.autosomes.bam, with $indir_closestSNP/*.autosomes.closest_SNP.bed writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

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

#module load /fast/users/jmarkow_m/work/miniconda
source activate H1
# first take all lines where the distance between read and closest SNP is smaller or equal to the given distance cutoff

# then proceede to grep the read IDs according to the observed haplotype of the closest SNP
awk -v dist="$dist_cutoff" 'BEGIN {FS=\"\t\"}; {if ( \$14 <= dist ) {print \$0} }' $closest_snp | grep -F 'p1' | cut -f4  > $tempdir/$b.$dist_cutoff.p1
awk -v dist="$dist_cutoff" 'BEGIN {FS=\"\t\"}; {if ( \$14 <= dist ) {print \$0} }' $closest_snp | grep -F 'p2' | cut -f4  > $tempdir/$b.$dist_cutoff.p2
java -jar /fast/users/jmarkow_m/work/miniconda/envs/H1/share/picard-2.23.8-0/picard.jar FilterSamReads I=$f O=$outdir/$dist_cutoff.bp/p1/$b.assigned_p1.bam READ_LIST_FILE=$tempdir/$b.$dist_cutoff.p1 FILTER=includeReadList
java -jar /fast/users/jmarkow_m/work/miniconda/envs/H1/share/picard-2.23.8-0/picard.jar FilterSamReads I=$f O=$outdir/$dist_cutoff.bp/p2/$b.assigned_p2.bam READ_LIST_FILE=$tempdir/$b.$dist_cutoff.p2 FILTER=includeReadList
rm $tempdir/$b.$dist_cutoff.p1
rm $tempdir/$b.$dist_cutoff.p2" > $job_file
        sbatch $job_file
        ((i++))
done

