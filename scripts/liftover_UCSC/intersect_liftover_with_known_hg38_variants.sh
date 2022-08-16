if [ "$#" -ne 2 ]; then
	echo "USAGE: $0 <input_dir> <output_dir>"
	exit
fi

indir=$1 ## directory with input files 
outdir=$2 ## directory to store output 

if [ ! -d "$indir" ]; then
    echo "directory $indir not found. Exiting."
    exit 1
fi 

if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
    mkdir "$outdir/log"
fi

files=`find $indir -name '*[0-9].rs_sorted.bed'`

echo "Running analysis with $0 on $indir/*.bam, writing to $outdir at $(date +%x_%r)" > $outdir/run.log

for f in $files; do
    	b=`basename $f .rs_sorted.bed`
	qsub_cmd="join -1 4 -2 4 -o 0,1.1,1.3,2.1,2.3 /fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/GSE52457_H1_haps_onlySNPs_liftover_hg18_to_hg38_sorted_rs_only_rs.bed "$f""
	echo $qsub_cmd '>' "$outdir/$b.liftover_overlap.bed" | qsub -cwd -V -l h_vmem=120g -o "$outdir/log/$b.log" -e "$outdir/log/$b.log"
	echo
done
