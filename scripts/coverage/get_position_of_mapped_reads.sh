# read in list of files,
# open a new output file,
# loop over list of files,
# into that new  output file, paste the filename
# then paste the chr, position and CIGAR string (length) of mapped reads in that file, comma separated

files=`cat /fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60_autosomes/list_of_filenames.txt`

outfile=/fast/groups/ag_schwarz/Projects/project-gam/H1/results/coverage/positions_of_mapped_reads.csv

touch $outfile

cd /fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60_autosomes/

for f in $files; do
	samtools view $f | awk '{print $3,$4,$6,"'$f'"}' >> $outfile
done

