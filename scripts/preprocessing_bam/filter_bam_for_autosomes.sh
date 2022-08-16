if [ "$#" -ne 2 ]; then
	echo "USAGE: $0 <input_dir> <output_dir>"
	exit
fi

indir=$1 ## directory with input files 
outdir=$2 ## directory to store output 
#autosomes=/fast/groups/ag_schwarz/Projects/project-gam/H1/data/hg38_sorted_autosomes.bed
autosomes=/fast/groups/ag_schwarz/Projects/project-gam/data/mm10_chromsize.bed

if [ ! -d "$indir" ]; then
    echo "directory $indir not found. Exiting."
    exit 1
fi 

if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
    mkdir "$outdir/log"
fi


#files=`find $indir -name '*.bam'` # or which ever extension of files it is 
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60/Repl1/H1hESC_3NP_R1_GAM_171121_A03_2_S113_R1_001.rmdup.bam
files=/fast/groups/ag_schwarz/Projects/project-gam/results/trimmed_F123_3NP/f123_3NP_GAM161008_B7_S50_R1_001.fastq.gz.rmdup.trimmed_last_base.bam

echo "Running analysis with $0 on $indir/*.bam, writing to $outdir at $(date +%x_%r)" > $outdir/run.log


if [ ! -d "$outdir/log" ]; then
    mkdir -p "$outdir/log"
fi

job_directory=/fast/users/jmarkow_m/scratch/job
i=1
for f in $files; do
        b=`basename $f .bam`
        echo "${i}"
        job_file="${job_directory}/${i}_filter_autosomes.job"

        echo "#!/bin/bash
#SBATCH --job-name=${i}_filter_autosomes.job
#SBATCH --output=$outdir/log/${i}.out
#SBATCH --error=$outdir/log/${i}.err

        samtools view -b -h -L "$autosomes" -o "$outdir/$b.autosomes.bam" "$f"
" > $job_file
#        sbatch $job_file
        ((i++))
done

