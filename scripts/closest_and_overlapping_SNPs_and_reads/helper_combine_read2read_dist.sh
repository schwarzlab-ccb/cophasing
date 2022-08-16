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

files=`find $indir -name '*.closest_read_excl_same.dedup.bed'` # or which ever extension of files it is
#files=/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/closest_reads_intermediate_files/Repl1/H1hESC_3NP_R1_GAM_171121_A01_1_S1_R1_001.rmdup.autosomes.closest_read_excl_same_smaller_pos.bed

for f in $files; do
	awk '{print $13}' $f >> $outdir/combined_distances.txt
done


