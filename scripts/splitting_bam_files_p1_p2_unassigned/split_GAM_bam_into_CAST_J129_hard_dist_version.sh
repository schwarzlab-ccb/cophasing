#!/bin/bash

# this script uses multiple instances of grep to find and extract the reads assigned to either parental chromosome p1 or parental chromosome p2 or can not be assigned. 
# the base of this assignment is (1) the pileup of the reads in the GAM sample at known H1 SNP positions (liftover from dixon). From this pileup, we know which SNPS are observed in this sample and which haplotype / parental chromosome was covered. We used these observed SNPs per sample to create (2) a file reporting the distance of each read in the sample to its closest SNP
# now we use the information about the observed haplotype of the covered SNPs (sample specific) and the information about the closest SNP to each read to assign the reads to the haplotype of their closest SNP
# in the first version, as a base to improve from, each read was assigned to the haplotype of the closest SNP, no reads were unassigned (except in the unlikely case that there are reads from the chromosome, but no single SNP covered) 
# in this next versions I introduce hard distance cutoffs, which can be given to the script as an additional argument
# more complex strategies where more factors are included in the decision making will be introduced in following more advanced and complex scripts. 


# include USAGE message, if no arguments are given to the script
if [ "$#" -ne 3 ]; then
        echo "USAGE: $0 <input_bam_dir> <input_closestSNP_dir> <output_dir>"
        exit
fi

indir_bam=$1 ## directory with input bam files of GAM experiments
indir_closestSNP=$2 ## directory with input bed files of reads and their closest SNP
#dist_cutoff=$3 ## distance cut off in base pairs. reads with a larger distance to the nearest SNP will not be phased / assigned to a parental haplotype but remain unassigned
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

# check that the distance cutoff is given as an integer
#re='^[0-9]+$'
#if ! [[ $dist_cutoff =~ $re ]] ; then
#   echo "given distance cutoff must be an integer, please provide the distance cutoff in bp. Exiting." >&2; exit 1
#fi

# if outputdirectory is not existing yet, make it
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

if [ ! -d "$outdir/log" ]; then
    mkdir -p "$outdir/log"
fi

for n in {1..10};do
	mkdir -p "$outdir/bam_files_haplotype_assigned_reads_$n/"
	for dist in 0 50000 100000 200000; do
		mkdir -p "$outdir/bam_files_haplotype_assigned_reads_$n/$dist.bp/CAST"
		mkdir -p "$outdir/bam_files_haplotype_assigned_reads_$n/$dist.bp/J129"
	done
done

files=`find $indir_bam -name '*.autosomes.bam'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/data/F123_3NP/F123_3NPs_last_base_trimmed_autosomes/f123_3NP_GAM161008_A7_S49_R1_001.fastq.gz.rmdup.trimmed_last_base.autosomes.bam
#files=/fast/groups/ag_schwarz/Projects/project-gam/data/F123_bam_files_1123_1NP_mm10_Sept_2019_last_base_trimmed/f123_1NP_GAM171024_C01_1_S3_R1_001.fastq.gz.rmdup.trimmed_last_base.autosomes.bam
#files=/fast/groups/ag_schwarz/Projects/project-gam/data/F123_bam_files_1123_1NP_mm10_Sept_2019_last_base_trimmed/f123_1NP_GAM180207_H02_1_S16_R1_001.fastq.gz.rmdup.trimmed_last_base.autosomes.bam

# write a little log file, that states with which input the script as running and when, and store it in the outputdirectory
echo "Running analysis with $0 on $indir_bam/*.autosomes.bam, with $indir_closestSNP/closest_SNP_subsampled_\$n/*.closest_SNP.bed writing to $outdir at $(date +%x_%r)" >> $outdir/run.log

# the cluster job scheduler is now slurm, therefore the commands need to be changed from qsub
# we create small job scripts which are send off with sbatch:

job_directory=/fast/users/jmarkow_m/scratch/job

i=1
for f in $files; do

        echo "${i}"
        job_file="${job_directory}/${i}_split.job"

	b=`basename $f .autosomes.bam`

        echo "#!/bin/bash
#SBATCH --job-name=${i}_split.job
#SBATCH --output=$outdir/log/${i}_split.out
#SBATCH --error=$outdir/log/${i}_split.err

source /fast/users/jmarkow_m/work/miniconda/etc/profile.d/conda.sh
conda activate H1

for n in {1..10};do
        echo \$n
        indir_closestSNP_dir=$indir_closestSNP/closest_SNP_subsampled_\$n
        echo \$indir_closestSNP_dir
        closest_snp=\`find \$indir_closestSNP_dir -name $b*\`
        echo \$closest_snp

        for dist in 0 50000 100000 200000; do
		# first take all lines where the distance between read and closest SNP is smaller or equal to the given distance cutoff
		# then proceede to grep the read IDs according to the observed haplotype of the closest SNP
		awk -v dist="\$dist" 'BEGIN {FS=\"\t\"}; {if ( \$13 <= dist && \$13 >= 0 ) {print \$0} }' \$closest_snp | grep -F 'CAST' | cut -f4  > $tempdir/\$n.$b.\$dist.CAST
		awk -v dist="\$dist" 'BEGIN {FS=\"\t\"}; {if ( \$13 <= dist && \$13 >= 0 ) {print \$0} }' \$closest_snp | grep -F 'J129' | cut -f4  > $tempdir/\$n.$b.\$dist.J129
		java -jar /fast/users/jmarkow_m/work/miniconda/envs/H1/share/picard-2.23.8-0/picard.jar FilterSamReads I=$f O=$outdir/bam_files_haplotype_assigned_reads_\$n/\$dist.bp/CAST/$b.assigned_CAST.bam READ_LIST_FILE=$tempdir/\$n.$b.\$dist.CAST FILTER=includeReadList
		java -jar /fast/users/jmarkow_m/work/miniconda/envs/H1/share/picard-2.23.8-0/picard.jar FilterSamReads I=$f O=$outdir/bam_files_haplotype_assigned_reads_\$n/\$dist.bp/J129/$b.assigned_J129.bam READ_LIST_FILE=$tempdir/\$n.$b.\$dist.J129 FILTER=includeReadList
		rm $tempdir/\$n.$b.\$dist.CAST
		rm $tempdir/\$n.$b.\$dist.J129
        done
done
" > $job_file
        sbatch $job_file
        ((i++))
done

