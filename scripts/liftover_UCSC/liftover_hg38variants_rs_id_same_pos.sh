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

files=`find $indir -name '*.bed'`


echo "Running analysis with $0 on $indir/*.bam, writing to $outdir at $(date +%x_%r)" > $outdir/run.log

for f in $files; do
    	b=`basename $f .bed`
#	qsub -cwd -V -l h_vmem=120g -o "$outdir/log/$b.log" -e "$outdir/log/$b.log" -- awk '$2 == $4 && $3 == $5' "$f" > "$outdir/$b.liftover_hg38variants_identical.bed"
#	qsub -cwd -V -l h_vmem=120g -o "$outdir/log/$b.log" -e "$outdir/log/$b.log" <<EOF
#	awk '$2 == $4 && $3 == $5' "$f" > "$outdir/$b.liftover_hg38variants_identical.bed"
#EOF
#	awk '$2 == $4 && $3 == $5' "$f" > "$outdir/$b.liftover_hg38variants_identical.bed"
	awk '$2 == $4 && $3 == $5' "$f" | awk '{print $2, $3, $1}' | sort -n -k2 > "$outdir/$b.liftover_hg38variants_identical.bed"
	echo
done

