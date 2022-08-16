# create SNP seeds specifically per GAM sample

# load the pileup results of each sample again, to be 100% sure the resulting SNP seed file is named correctly
# again, clean / filter combined SNP set for matching REF and ALT bases, as well as phased SNPs and Mappability of 1
# create bed file of filtered, observed SNPs for each GAM sample separately

################################################################################
################################################################################
################################################################################
# libraries
library(GAMIBHEAR)
library(stringr)
library(data.table)

# load pileup results
vcf_dir_1 = '/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/results/pileup/Repl1'
vcf_dir_2 = '/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/results/pileup/Repl2'
files_1 = list.files(path = vcf_dir_1, pattern = '*.pileup.vcf$')
files_2 = list.files(path = vcf_dir_2, pattern = '*.pileup.vcf$')

VCF_list_1 = load_sample_VCF(files_1,vcf_dir_1)
VCF_list_2 = load_sample_VCF(files_2,vcf_dir_2)

processed_SNVs = process_sample_VCF(c(VCF_list_1,VCF_list_2))
# dim(processed_SNVs) [1] 1911075       7
processed_SNVs$repl=1
processed_SNVs$repl[which(as.numeric(str_sub(processed_SNVs$SAMPLE,7))>length(files_1))]=2

# load truth, merge
H1_truth=readRDS("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/H1_dixon_haplotypes_liftover_hg38.rds")
processed_SNVs_merged_truth=merge(processed_SNVs,H1_truth,by=c('CHROM','POS'),all.x=T,all.y=F)
# dim(processed_SNVs_merged_truth) [1] 1911075      15

# load mappability, merge
mapp=fread('/fast/groups/ag_schwarz/Projects/project-gam/H1/Genotyping/intermediate/mappability/mapinter_v3')
mappdb=as.data.frame(mapp)
mappdb=mappdb[,c(1,2,7)]
colnames(mappdb)=c('CHROM','POS','MAPPABILITY')
processed_SNVs_merged_truth_mapp=merge(processed_SNVs_merged_truth,mappdb,by=c('CHROM','POS'),all.x=T,all.y=F)
# dim(processed_SNVs_merged_truth_mapp) [1] 1911075      16

# remove unphased SNPs
processed_SNVs_merged_truth_mapp_phased = processed_SNVs_merged_truth_mapp[substr(processed_SNVs_merged_truth_mapp$GT,2,2)=='|',]
# dim(processed_SNVs_merged_truth_mapp_phased) [1] 1843340      16

# remove SNPs with mappability lower 1 
processed_SNVs_merged_truth_mapp_phased_minmapp1 = processed_SNVs_merged_truth_mapp_phased[which(processed_SNVs_merged_truth_mapp_phased$MAPPABILITY==1),]
# dim(processed_SNVs_merged_truth_mapp_phased_minmapp1) [1] 1703137      16

# remove SNPs where hg38 REF allele from the pileup does not match the given Dixon SNP REF allele
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF = processed_SNVs_merged_truth_mapp_phased_minmapp1[processed_SNVs_merged_truth_mapp_phased_minmapp1$REF.x==processed_SNVs_merged_truth_mapp_phased_minmapp1$REF.y,]
# dim(processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF) [1] 1699689      16

# remove SNPs where observed ALT allele from the pileup does not match the given Dixon SNP ALT allele / or N, in case it was not observed
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT = processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF$ALT.x==processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF$ALT.y|processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF$ALT.x=="N",]
# dim(processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT) [1] 1695744      16


# that finishes the filtering / cleaning
# now delete irrelevant columns, add some more information like total counts and observed haplotypes 
# delete START and TYPE column
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT=processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT[, !names(processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT) %in% c("START","TYPE")]
# add total counts
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$TOTAL_COUNT = processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$REF_COUNT + processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$ALT_COUNT
# check if both alleles are observed 
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$obs_both_alleles=processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$REF_COUNT>0 & processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$ALT_COUNT>0
# assign the observed haplotype 
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$p1=0
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$p2=0
# p1 = 1|0 
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$p1[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="0|1"]=processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$REF_COUNT[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="0|1"]
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$p1[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="1|0"]=processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$ALT_COUNT[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="1|0"]
# p2 = 0|1 
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$p2[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="0|1"]=processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$ALT_COUNT[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="0|1"]
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$p2[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="1|0"]=processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$REF_COUNT[processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$GT=="1|0"]
# chrom, pos sort
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT = processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT[order(processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$POS),]
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT = processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT[order(as.numeric(substring(processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$CHROM, 4))),]
# save
saveRDS(processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT, '/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/H1_dixon_haplotypes_liftover_hg38_filter_phased_obs_ref_alt_mapp.rds')
#processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT=readRDS('/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/H1_dixon_haplotypes_liftover_hg38_filter_phased_obs_ref_alt_mapp.rds')
################################################################################
# write out the file names
# select each sample, write out the SNP seed set in bed format with the correct file name, no deduplication needed

#substr file name wo get rid off the file ending 
#H1hESC_3NP_R1_GAM_171121_A03_2_S113_R1_001.rmdup.autosomes.pileup.vcf > H1hESC_3NP_R1_GAM_171121_A01_1_S1_R1_001.rmdup.autosomes.pileup.closest_SNP_seed.bed


processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$START = processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$POS-1
processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT$PARENT = c("p1","p2")[apply(processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT[,c(17,18)],1,which.max)]
SNP_seed_beds=processed_SNVs_merged_truth_mapp_phased_minmapp1_sameREF_sameALT[,c('CHROM','START','POS','ID','GT','PARENT','SAMPLE')]

files=c(files_1,files_2)
outfilenames=sub('\\.vcf$', '.closest_SNP_seed.bed', files)
samples_sort = unique(SNP_seed_beds$SAMPLE)[order(as.numeric(substring(unique(SNP_seed_beds$SAMPLE), 7)))] # e.g. 'sample7'

for (i in 1:length(files)){
    print(paste(i,' / ', length(files)))
    sample_seed_SNPs=SNP_seed_beds[which(SNP_seed_beds$SAMPLE==samples_sort[i]),]
    if(i <= length(files_1)){
        write.table(sample_seed_SNPs,paste('/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/closest_SNP_seed_files/Repl1/',outfilenames[i],sep = ''),
                    sep = '\t',
                    quote = F,
                    row.names = F,
                    col.names = F)
    }else{
        write.table(sample_seed_SNPs,paste('/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/closest_SNP_seed_files/Repl2/',outfilenames[i],sep = ''),
                    sep = '\t',
                    quote = F,
                    row.names = F,
                    col.names = F)
    }
}

################################################################################
