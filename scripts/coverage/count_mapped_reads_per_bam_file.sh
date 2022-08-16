# read in list of files,
# open a new output file,
# loop over list of files,
# into that new  output file, paste the filename
# then paste the number of mapped reads in that file, comma separated

files=`cat /fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60_autosomes/list_of_filenames.txt`

outfile=/fast/groups/ag_schwarz/Projects/project-gam/H1/results/coverage/mapped_reads_per_bam.csv

touch $outfile

# the lidt of files has relative paths stored, so in order tof ind them, we need to change to the directory
cd /fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_3NP_GAM_13032020/H1hESC_hg38_PMR20_OW60_autosomes/

for f in $files; do
	counts=`samtools view -c -F 260 $f`
	echo $f,$counts >> $outfile
done

