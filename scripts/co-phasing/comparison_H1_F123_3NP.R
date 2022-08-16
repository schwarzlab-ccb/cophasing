# comparison H1 and F123 3NP data
################################################################################
### load in data
#/Users/jmarkow/Desktop/GAM/scratch/cophasing/
    
number_reads_F123_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_reads_F123_3NP_per_chr.rds")
number_reads_F123_3NP = read.csv("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_reads_F123_3NP.txt",header = F,sep = " ")
number_reads_H1_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_reads_H1_3NP_per_chr.rds")
number_reads_H1_3NP = read.csv("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_reads_H1_3NP.txt",header = F,sep = " ")
number_SNPs_F123_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_SNPs_F123_3NP_per_chr.rds")
number_SNPs_F123_3NP = read.csv("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_SNPs_F123_3NP.txt",header = F,sep = " ")
number_SNPs_H1_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_SNPs_H1_3NP_per_chr.rds")
number_SNPs_H1_3NP = read.csv("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_SNPs_H1_3NP.txt",header = F,sep = " ")

full_SNPs_mm10 = readRDS('/Users/jmarkow/Desktop/GAM/data/clean_mm10_gt_merged_mappability_minmap1.rds')
nrow(full_SNPs_mm10) #18150228
nr_snps_F132_all = nrow(full_SNPs_mm10) 
full_SNPs_mm10$CHROM = factor(full_SNPs_mm10$CHROM, levels=unique(full_SNPs_mm10$CHROM))
full_SNPs_mm10_chr = split(full_SNPs_mm10,full_SNPs_mm10$CHROM)
nr_snps_F132_all_chr = unlist(lapply(full_SNPs_mm10_chr,nrow))

SNPs_mm10=readRDS('/Users/jmarkow/Desktop/GAM/data/clean_mm10_gt_merged_mappability_minmap1.rds')
nrow(SNPs_mm10)

#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/results/mara-scripts/H1_truths_mapp_filtered_H1_GT.rds /Users/jmarkow/Desktop/GAM/data/
SNPs_H1 = readRDS("/Users/jmarkow/Desktop/GAM/data/H1_truths_mapp_filtered_H1_GT.rds")
nrow(SNPs_H1) #1506858
nr_snps_H1_all = nrow(SNPs_H1) 
SNPs_H1_possort = SNPs_H1[order(SNPs_H1$POS),]
SNPs_H1 = SNPs_H1_possort[order(as.numeric(substring(SNPs_H1_possort$CHROM, 4))),]
SNPs_H1$CHROM = factor(SNPs_H1$CHROM, levels=unique(SNPs_H1$CHROM))
SNPs_H10_chr = split(SNPs_H1,SNPs_H1$CHROM)
nr_snps_H1_all_chr = unlist(lapply(SNPs_H10_chr,nrow))

### number of reads per 3NP sample
hist(number_reads_F123_3NP$V1,
     breaks=seq(-0.5,(max(number_reads_F123_3NP$V1)+50000.5),50000),
     col = adjustcolor("grey",alpha.f = 0.5), 
     freq = T,
     main = "number of reads per 3NP sample",
     xlab = "number of reads")
hist(number_reads_H1_3NP$V1,breaks=seq(-0.5,(max(number_reads_F123_3NP$V1)+50000.5),50000),freq = T,add=T,col=adjustcolor("gold",0.5))
abline(v=mean(number_reads_F123_3NP$V1),col='grey',lwd=3)
abline(v=mean(number_reads_H1_3NP$V1),col='gold',lwd=3)
legend("topright",fill=c(adjustcolor("grey",alpha.f = 0.5),adjustcolor("gold",0.5)),legend=c("F123 (1019)","H1 (381)"),bty='n')


hist(number_reads_F123_3NP$V1,
     breaks=seq(-0.5,(max(number_reads_F123_3NP$V1)+50000.5),50000),
     col = adjustcolor("grey",alpha.f = 0.5), 
     freq = F,
     ylim = c(0,12e-7),
     main = "number of reads per 3NP sample",
     xlab = "number of reads")
hist(number_reads_H1_3NP$V1,breaks=seq(-0.5,(max(number_reads_F123_3NP$V1)+50000.5),50000),freq = F,add=T,col=adjustcolor("gold",0.5))
abline(v=mean(number_reads_F123_3NP$V1),col='grey',lwd=3)
abline(v=mean(number_reads_H1_3NP$V1),col='gold',lwd=3)
legend("topright",fill=c(adjustcolor("grey",alpha.f = 0.5),adjustcolor("gold",0.5)),legend=c("F123 (1019)","H1 (381)"),bty='n')


### number of chromosomes observed (reads)



### number of reads per chromosome

### number / percent of SNPS observed per 3NP sample
nr_snps_F132_all
nr_snps_H1_all
nr_snps_F132_all_chr
nr_snps_H1_all_chr

breaksize = 1000
hist(number_SNPs_F123_3NP$V1,
     breaks=seq(-0.5,(max(number_SNPs_F123_3NP$V1)+breaksize+0.5),breaksize),
     col = adjustcolor("grey",alpha.f = 0.5), 
     freq = T,
     main = "number of obseved SNPs per 3NP sample",
     xlab = "number of SNPs")
hist(number_SNPs_H1_3NP$V1,
     breaks=seq(-0.5,(max(number_SNPs_H1_3NP$V1)+breaksize+0.5),breaksize),
     col = adjustcolor("gold",alpha.f = 0.5), 
     freq = T,
     main = "number of obseved SNPs per 3NP sample",
     xlab = "number of SNPs")


breaksize = 0.001
hist(number_SNPs_F123_3NP$V1/nr_snps_F132_all,
     breaks=seq(0,(max(number_SNPs_F123_3NP$V1/nr_snps_F132_all)+breaksize),breaksize),
     xlim = c(0,max(max(number_SNPs_F123_3NP$V1/nr_snps_F132_all),max(number_SNPs_H1_3NP$V1/nr_snps_H1_all))),
     col = adjustcolor("grey",alpha.f = 0.5), 
     freq = T,
     main = "percent of obseved SNPs per 3NP sample",
     xlab = "percent of SNPs")
hist(number_SNPs_H1_3NP$V1/nr_snps_H1_all,
     breaks=seq(0,(max(number_SNPs_H1_3NP$V1/nr_snps_H1_all)+breaksize),breaksize),
     freq = T,
     add = T,
     col=adjustcolor("gold",0.5))
legend("topright",fill=c(adjustcolor("grey",alpha.f = 0.5),adjustcolor("gold",0.5)),legend=c("F123 (1019)","H1 (381)"),bty='n')

breaksize = 0.001
hist(number_SNPs_F123_3NP$V1/nr_snps_F132_all,
     breaks=seq(0,(max(number_SNPs_F123_3NP$V1/nr_snps_F132_all)+breaksize),breaksize),
     xlim = c(0,max(max(number_SNPs_F123_3NP$V1/nr_snps_F132_all),max(number_SNPs_H1_3NP$V1/nr_snps_H1_all))),
     ylim = c(0,350),
     col = adjustcolor("grey",alpha.f = 0.5), 
     freq = F,
     main = "percent of obseved SNPs per 3NP sample",
     xlab = "percent of SNPs")
hist(number_SNPs_H1_3NP$V1/nr_snps_H1_all,
     breaks=seq(0,(max(number_SNPs_H1_3NP$V1/nr_snps_H1_all)+breaksize),breaksize),
     freq = F,
     add = T,
     col=adjustcolor("gold",0.5))
legend("topright",fill=c(adjustcolor("grey",alpha.f = 0.5),adjustcolor("gold",0.5)),legend=c("F123 (1019)","H1 (381)"),bty='n')




### number / percent of SNPS observed per chromosome per 3NP sample

### number of chromosomes observed (SNPs)

################################################################################
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/processed_F123_SNPs.rds /Users/jmarkow/Desktop/GAM/results/processed_SNPs_latest/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_3NP_mm10_all_files_and_processed_observed_SNPs.rds /Users/jmarkow/Desktop/GAM/results/processed_SNPs_latest/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/H1_dixon_haplotypes_liftover_hg38_filter_phased_obs_ref_alt_mapp.rds /Users/jmarkow/Desktop/GAM/results/processed_SNPs_latest/

processed_SNVs_F123_1NP=readRDS("/Users/jmarkow/Desktop/GAM/results/processed_SNPs_latest/processed_F123_SNPs.rds") # 35879529        7
processed_SNVs_F123_3NP=readRDS('/Users/jmarkow/Desktop/GAM/results/processed_SNPs_latest/F123_3NP_mm10_all_files_and_processed_observed_SNPs.rds')
processed_SNVs_F123_3NP = processed_SNVs_F123_3NP[[2]] # 60919164        7
processed_SNVs_H1_3NP = readRDS('/Users/jmarkow/Desktop/GAM/results/processed_SNPs_latest/H1_dixon_haplotypes_liftover_hg38_filter_phased_obs_ref_alt_mapp.rds')
#[1] 1695744      18
    
sample_contribution_F123_1NP = as.data.frame(sort(table(processed_SNVs_F123_1NP$SAMPLE),decreasing = T))
percentage_SNPs_per_sample_F123_1NP = (sample_contribution_F123_1NP$Freq/nr_snps_F132_all)*100

sample_contribution_F123_3NP = as.data.frame(sort(table(processed_SNVs_F123_3NP$SAMPLE),decreasing = T))
percentage_SNPs_per_sample_F123_3NP = (sample_contribution_F123_3NP$Freq/nr_snps_F132_all)*100

sample_contribution_H1_3NP = as.data.frame(sort(table(processed_SNVs_H1_3NP$SAMPLE),decreasing = T))
percentage_SNPs_per_sample_H1_3NP = (sample_contribution_H1_3NP$Freq/nr_snps_H1_all)*100

cex_plot=1
cex_rel=1.2#0.8
par(cex=cex_plot)

# sum(number_SNPs_F123_3NP%in%sample_contribution$Freq)
hist(sample_contribution_F123_3NP$Freq,
     breaks = seq(0,max(sample_contribution_F123_3NP$Freq)+2000,2000),
     xlab='Number of SNVs per NuP',
     #xlab='Number of SNVs per sample',
     ylab='Number of NuPs',
     ylim=c(0,130),
     #ylab='Number of samples',
     main='',
     #,main='Number of observed SNPs \n  per sample'
     col=adjustcolor('red',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
axis(3,labels = paste(round((seq(0,300000,50000)/nr_snps_F132_all)*100,2),'%',sep=''),at=seq(0,300000,50000),tcl=0.5,mgp = c(0, -2, 0),cex.axis=cex_rel)
abline(v=mean(sample_contribution_F123_3NP$Freq),col='red',lwd=3) # (mean(sample_contribution_F123_3NP$Freq)/nrow(SNPs_mm10))*100 # [1] 0.3293803
# hist(number_SNPs_F123_3NP$V1,
#      breaks = seq(0,max(number_SNPs_F123_3NP$V1)+2000,2000),
#      add=T,col="red")

hist(sample_contribution_F123_1NP$Freq,
     breaks = seq(0,max(sample_contribution_F123_1NP$Freq)+2000,2000),
     xlab='Number of SNVs per NuP',
     #xlab='Number of SNVs per sample',
     ylab='Number of NuPs',
     #ylab='Number of samples',
     main='',
     #,main='Number of observed SNPs \n  per sample'
     col=adjustcolor('grey80',alpha.f = 0.3),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1,
     add=T)
abline(v=mean(sample_contribution_F123_1NP$Freq),col='grey80',lwd=3)

hist(sample_contribution_H1_3NP$Freq,
     breaks = seq(0,max(sample_contribution_H1_3NP$Freq)+2000,2000),
     xlab='Number of SNVs per NuP',
     #xlab='Number of SNVs per sample',
     ylab='Number of NuPs',
     #ylab='Number of samples',
     main='',
     #,main='Number of observed SNPs \n  per sample'
     col=adjustcolor("gold",0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1,
     add=T)
axis(3,line=3,labels = paste(round((seq(0,300000,50000)/nr_snps_H1_all)*100,2),'%',sep=''),at=seq(0,300000,50000),tcl=0.5,mgp = c(0, -2, 0),cex.axis=cex_rel)
abline(v=mean(sample_contribution_H1_3NP$Freq),col='gold',lwd=3) #(mean(sample_contribution$Freq)/nr_snps_H1_all)*100 #0.2953677


mean(sample_contribution_H1_3NP$Freq)/nr_snps_H1_all*100 #0.2953677 # comparable to F123 3NP
mean(sample_contribution_F123_1NP$Freq)/nr_snps_F132_all*100 #0.1760293 
mean(sample_contribution_F123_3NP$Freq)/nr_snps_F132_all*100 #0.3293803 #1.87 times the 1NP


#### percent

binsize=0.02
hist(sample_contribution_F123_3NP$Freq/nr_snps_F132_all*100,
     breaks = seq(0,max(sample_contribution_F123_3NP$Freq/nr_snps_F132_all*100)+binsize,binsize),
     xlab='Percent of SNVs observed per sample',
     #xlab='Number of SNVs per sample',
    # ylab='Number of NuPs',
    ylim=c(0,130),
     #ylab='Number of samples',
     main='',
     #,main='Number of observed SNPs \n  per sample'
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
#axis(3,labels = paste(round((seq(0,300000,50000)/nr_snps_F132_all)*100,2),'%',sep=''),at=seq(0,300000,50000),tcl=0.5,mgp = c(0, -2, 0),cex.axis=cex_rel)
abline(v=mean(sample_contribution_F123_3NP$Freq)/nr_snps_F132_all*100,col='grey',lwd=3) # (mean(sample_contribution_F123_3NP$Freq)/nrow(SNPs_mm10))*100 # [1] 0.3293803
# hist(number_SNPs_F123_3NP$V1,
#      breaks = seq(0,max(number_SNPs_F123_3NP$V1)+2000,2000),
#      add=T,col="red")

hist(sample_contribution_F123_1NP$Freq/nr_snps_F132_all*100,
     breaks = seq(0,max(sample_contribution_F123_1NP$Freq/nr_snps_F132_all*100)+binsize,binsize),
     col=adjustcolor('red',alpha.f = 0.3),
     add=T)
abline(v=mean(sample_contribution_F123_1NP$Freq)/nr_snps_F132_all*100,col='red',lwd=3)

hist(sample_contribution_H1_3NP$Freq/nr_snps_H1_all*100,
     breaks = seq(0,max(sample_contribution_H1_3NP$Freq/nr_snps_H1_all*100)+binsize,binsize),
     col=adjustcolor("gold",0.5),
     add=T)
abline(v=mean(sample_contribution_H1_3NP$Freq)/nr_snps_H1_all*100,col='gold',lwd=3) 
legend("topright",fill=c(adjustcolor("red",alpha.f = 0.5),adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 1NP (1123)", "F123 3NP (1019)","H1 3NP (381)"),bty='n')



# density

hist(sample_contribution_F123_3NP$Freq/nr_snps_F132_all*100,
     breaks = seq(0,max(sample_contribution_F123_3NP$Freq/nr_snps_F132_all*100)+binsize,binsize),
     xlab='Percent of SNVs observed per sample',
     ylim=c(0,6),
     freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(sample_contribution_F123_3NP$Freq)/nr_snps_F132_all*100,col='grey',lwd=3) # (mean(sample_contribution_F123_3NP$Freq)/nrow(SNPs_mm10))*100 # [1] 0.3293803

hist(sample_contribution_F123_1NP$Freq/nr_snps_F132_all*100,
     breaks = seq(0,max(sample_contribution_F123_1NP$Freq/nr_snps_F132_all*100)+binsize,binsize),
     col=adjustcolor('red',alpha.f = 0.3),
     freq = F,
     add=T)
abline(v=mean(sample_contribution_F123_1NP$Freq)/nr_snps_F132_all*100,col='red',lwd=3)

hist(sample_contribution_H1_3NP$Freq/nr_snps_H1_all*100,
     breaks = seq(0,max(sample_contribution_H1_3NP$Freq/nr_snps_H1_all*100)+binsize,binsize),
     col=adjustcolor("gold",0.5),
     freq = F,
     add=T)
abline(v=mean(sample_contribution_H1_3NP$Freq)/nr_snps_H1_all*100,col='gold',lwd=3) 
legend("topright",fill=c(adjustcolor("red",alpha.f = 0.5),adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 1NP (1123)", "F123 3NP (1019)","H1 3NP (381)"),bty='n')

### try with the data I got from the pileups

# number_SNPs_F123_3NP
# number_SNPs_H1_3NP

hist(number_SNPs_F123_3NP$V1/nr_snps_F132_all*100,
     breaks = seq(0,max(number_SNPs_F123_3NP$V1/nr_snps_F132_all*100)+binsize,binsize),
     xlab='Percent of SNVs observed per sample',
     ylim=c(0,4),
     freq = F,
     main='',
     col=adjustcolor('red',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(number_SNPs_F123_3NP$V1)/nr_snps_F132_all*100,col='red',lwd=3) # 0.3475509

hist(number_SNPs_H1_3NP$V1/nr_snps_H1_all*100,
     breaks = seq(0,max(number_SNPs_H1_3NP$V1/nr_snps_H1_all*100)+binsize,binsize),
     col=adjustcolor('gold',alpha.f = 0.3),
     freq = F,
     add=T)
abline(v=mean(number_SNPs_H1_3NP$V1)/nr_snps_H1_all*100,col='gold',lwd=3) #0.350017

# from this data the percetnages are almost identical, so we have filtered out quite some with the cleaning of the h1 positions.
# still a difference, there are more samples with lower percetnage of SNPs in H1
legend("topright",fill=c(adjustcolor("red",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

################################################################################
# per chromosome
# number_reads_F123_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_reads_F123_3NP_per_chr.rds")
# number_reads_H1_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_reads_H1_3NP_per_chr.rds")
# number_SNPs_F123_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_SNPs_F123_3NP_per_chr.rds")
# number_SNPs_H1_3NP_per_chr = readRDS("/Users/jmarkow/Desktop/GAM/scratch/cophasing/number_SNPs_H1_3NP_per_chr.rds")

hist(number_reads_F123_3NP_per_chr[,2:ncol(number_reads_F123_3NP_per_chr)])
number_reads_F123_3NP_per_chr_df = number_reads_F123_3NP_per_chr[,2:ncol(number_reads_F123_3NP_per_chr)]
sum(is.na(number_reads_F123_3NP_per_chr_df)) #0
table(number_reads_F123_3NP_per_chr_df>0) # all TRUE 19361
hist(number_reads_F123_3NP_per_chr_df)
hist(as.data.frame(number_reads_F123_3NP_per_chr_df[,1:7]))
is.list(number_reads_F123_3NP_per_chr_df)
is.data.frame(number_reads_F123_3NP_per_chr_df)
hist(number_reads_F123_3NP_per_chr_df)
boxplot(t(number_reads_F123_3NP_per_chr_df))
hist(t(number_reads_F123_3NP_per_chr_df))
table(is.na(t(number_reads_F123_3NP_per_chr_df))) 


hist(t(number_reads_F123_3NP_per_chr_df))




number_reads_F123_3NP_per_chr_df = number_reads_F123_3NP_per_chr[,2:ncol(number_reads_F123_3NP_per_chr)]
rownames(number_reads_F123_3NP_per_chr_df) = number_reads_F123_3NP_per_chr[,1]

boxplot(t(number_reads_F123_3NP_per_chr_df),ylim=c(0,max(number_reads_F123_3NP_per_chr_df,number_reads_H1_3NP_per_chr_df)))
mean(number_reads_F123_3NP_per_chr_df)
number_reads_F123_3NP_per_chr_dm = data.matrix(number_reads_F123_3NP_per_chr_df)
mean(number_reads_F123_3NP_per_chr_dm) # 45122.67
binsize = 10000
hist(number_reads_F123_3NP_per_chr_dm,
     breaks = seq(0,max(number_reads_F123_3NP_per_chr_dm)+binsize, binsize),
     xlab='Reads per chromosome per sample',
     #ylim=c(0,4),
     #freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(number_reads_F123_3NP_per_chr_dm),col='grey',lwd=3) # 45122.67


hist(log10(number_reads_F123_3NP_per_chr_dm+1),
     breaks = seq(0,7,0.1),
     xlab='log10(Reads per chromosome per sample+1)',
     #ylim=c(0,4),
     #freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(log10(number_reads_F123_3NP_per_chr_dm+1)),col='grey',lwd=3) 

### H1
number_reads_H1_3NP_per_chr_df = number_reads_H1_3NP_per_chr[,2:ncol(number_reads_H1_3NP_per_chr)]
rownames(number_reads_H1_3NP_per_chr_df) = number_reads_H1_3NP_per_chr[,1]

table(is.na(number_reads_H1_3NP_per_chr_df))
table(number_reads_H1_3NP_per_chr_df==0)
number_reads_H1_3NP_per_chr_df[,370:381]
boxplot(t(number_reads_H1_3NP_per_chr_df),ylim=c(0,max(number_reads_F123_3NP_per_chr_df,number_reads_H1_3NP_per_chr_df)))
number_reads_H1_3NP_per_chr_dm = data.matrix(number_reads_H1_3NP_per_chr_df)
mean(number_reads_H1_3NP_per_chr_dm) # 41253.38

binsize = 10000
hist(number_reads_H1_3NP_per_chr_dm,
     breaks = seq(0,max(number_reads_H1_3NP_per_chr_dm)+binsize, binsize),
     xlab='Reads per chromosome per sample',
     #ylim=c(0,4),
     #freq = F,
     main='',
     col=adjustcolor('gold',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(number_reads_H1_3NP_per_chr_dm),col='gold',lwd=3) # 45122.67


hist(log10(number_reads_H1_3NP_per_chr_dm+1),
     breaks = seq(0,7,0.1),
     xlab='log10(Reads per chromosome per sample+1)',
     #ylim=c(0,4),
     #freq = F,
     main='',
     col=adjustcolor('gold',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(log10(number_reads_H1_3NP_per_chr_dm+1)),col='gold',lwd=3) 



# combined
hist(log10(number_reads_F123_3NP_per_chr_dm+1),
     breaks = seq(0,7,0.1),
     xlab='log10(Reads per chromosome per sample+1)',
     #freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(log10(number_reads_F123_3NP_per_chr_dm+1)),col='grey',lwd=3) 
hist(log10(number_reads_H1_3NP_per_chr_dm+1),
     breaks = seq(0,7,0.1),
     #freq = F,
     main='',
     col=adjustcolor('gold',alpha.f = 0.5),
     add = T)
abline(v=mean(log10(number_reads_H1_3NP_per_chr_dm+1)),col='gold',lwd=3) 
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

hist(log10(number_reads_F123_3NP_per_chr_dm+1),
     breaks = seq(0,7,0.1),
     xlab='log10(Reads per chromosome per sample+1)',
     ylim = c(0,1),
     freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(log10(number_reads_F123_3NP_per_chr_dm+1)),col='grey',lwd=3) 
hist(log10(number_reads_H1_3NP_per_chr_dm+1),
     breaks = seq(0,7,0.1),
     freq = F,
     main='',
     col=adjustcolor('gold',alpha.f = 0.5),
     add = T)
abline(v=mean(log10(number_reads_H1_3NP_per_chr_dm+1)),col='gold',lwd=3) 
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

##### number SNPs per chr
# number_SNPs_F123_3NP_per_chr 
# number_SNPs_H1_3NP_per_chr


number_SNPs_F123_3NP_per_chr_df = number_SNPs_F123_3NP_per_chr[,2:ncol(number_SNPs_F123_3NP_per_chr)]
rownames(number_SNPs_F123_3NP_per_chr_df) = number_SNPs_F123_3NP_per_chr[,1]
number_SNPs_F123_3NP_per_chr_df_dm = data.matrix(number_SNPs_F123_3NP_per_chr_df)
table(is.na(number_SNPs_F123_3NP_per_chr_df_dm))
# FALSE  TRUE 
# 19319    42 

number_SNPs_H1_3NP_per_chr_df = number_SNPs_H1_3NP_per_chr[,2:ncol(number_SNPs_H1_3NP_per_chr)]
rownames(number_SNPs_H1_3NP_per_chr_df) = number_SNPs_H1_3NP_per_chr[,1]
number_SNPs_H1_3NP_per_chr_dm = data.matrix(number_SNPs_H1_3NP_per_chr_df)
#number_SNPs_H1_3NP_per_chr_dm[,365:381]
table(is.na(number_SNPs_H1_3NP_per_chr_dm))
# FALSE  TRUE 
# 8350    32
table(number_SNPs_H1_3NP_per_chr_dm==0) # all 8350 FALSE 

number_SNPs_F123_3NP_per_chr_df_dm_NA0 = number_SNPs_F123_3NP_per_chr_df_dm
number_SNPs_F123_3NP_per_chr_df_dm_NA0[which(is.na(number_SNPs_F123_3NP_per_chr_df_dm_NA0))]=0
number_SNPs_H1_3NP_per_chr_dm_NA0 = number_SNPs_H1_3NP_per_chr_dm
number_SNPs_H1_3NP_per_chr_dm_NA0[which(is.na(number_SNPs_H1_3NP_per_chr_dm_NA0))]=0
binsize = 1000
hist(number_SNPs_F123_3NP_per_chr_df_dm_NA0,
     breaks = seq(0,max(number_SNPs_F123_3NP_per_chr_df_dm_NA0)+binsize,binsize),
     xlab='observed SNPs per chromosome per sample',
     #freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(number_SNPs_F123_3NP_per_chr_df_dm_NA0),col='grey',lwd=3) # 3320.067

binsize = 100
hist(number_SNPs_H1_3NP_per_chr_dm_NA0,
     breaks = seq(0,max(number_SNPs_H1_3NP_per_chr_dm_NA0)+binsize,binsize),
     xlab='observed SNPs per chromosome per sample',
     #freq = F,
     main='',
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1,
     col=adjustcolor('gold',alpha.f = 0.5))
abline(v=mean(number_SNPs_H1_3NP_per_chr_dm_NA0),col='gold',lwd=3)  # 239.7391 per chr per sample
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

###

hist(log10(number_SNPs_F123_3NP_per_chr_df_dm_NA0+1),
     breaks = seq(0,5,0.1),
     xlab='log10(SNPs per chromosome per sample+1)',
     #freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(log10(number_SNPs_F123_3NP_per_chr_df_dm_NA0+1)),col='grey',lwd=3) 

hist(log10(number_SNPs_H1_3NP_per_chr_dm_NA0+1),
     breaks = seq(0,5,0.1),
     #freq = F,
     main='',
     col=adjustcolor('gold',alpha.f = 0.5),
     add = T)
abline(v=mean(log10(number_SNPs_H1_3NP_per_chr_dm_NA0+1)),col='gold',lwd=3) 
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

### percent
percent_SNPs_F123_3NP_per_chr_df_dm_NA0 = apply(number_SNPs_F123_3NP_per_chr_df_dm_NA0,2,function(col){return(col/nr_snps_F132_all_chr)})
percent_SNPs_H1_3NP_per_chr_df_dm_NA0 = apply(number_SNPs_H1_3NP_per_chr_dm_NA0,2,function(col){return(col/nr_snps_H1_all_chr)})
binsize=0.001
hist(percent_SNPs_F123_3NP_per_chr_df_dm_NA0,
     breaks = seq(0,max(percent_SNPs_F123_3NP_per_chr_df_dm_NA0)+binsize,binsize),
     xlab=' % observed SNPs per chromosome per sample',
     #freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(percent_SNPs_F123_3NP_per_chr_df_dm_NA0),col='grey',lwd=3) # 0.003495725


hist(percent_SNPs_H1_3NP_per_chr_df_dm_NA0,
     breaks = seq(0,max(percent_SNPs_H1_3NP_per_chr_df_dm_NA0)+binsize,binsize),
     xlab='observed SNPs per chromosome per sample',
     #freq = F,
     main='',
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1,
     col=adjustcolor('gold',alpha.f = 0.5),
     add = T)
abline(v=mean(percent_SNPs_H1_3NP_per_chr_df_dm_NA0),col='gold',lwd=3)  # 0.003611573
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

###

hist(percent_SNPs_F123_3NP_per_chr_df_dm_NA0,
     breaks = seq(0,max(percent_SNPs_F123_3NP_per_chr_df_dm_NA0)+binsize,binsize),
     xlab=' % observed SNPs per chromosome per sample',
     freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
abline(v=mean(percent_SNPs_F123_3NP_per_chr_df_dm_NA0),col='grey',lwd=3) # 0.003495725


hist(percent_SNPs_H1_3NP_per_chr_df_dm_NA0,
     breaks = seq(0,max(percent_SNPs_H1_3NP_per_chr_df_dm_NA0)+binsize,binsize),
     xlab='observed SNPs per chromosome per sample',
     freq = F,
     main='',
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1,
     col=adjustcolor('gold',alpha.f = 0.5),
     add = T)
abline(v=mean(percent_SNPs_H1_3NP_per_chr_df_dm_NA0),col='gold',lwd=3)  # 0.003611573
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

### log10 percent
hist(log10(percent_SNPs_F123_3NP_per_chr_df_dm_NA0),
     breaks = seq(-7,0,0.1),
     xlab='log10( % observed SNPs per chromosome per sample)',
     #freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
#abline(v=mean(log10(percent_SNPs_F123_3NP_per_chr_df_dm_NA0),na.rm=T),col='grey',lwd=3) 

hist(log10(percent_SNPs_H1_3NP_per_chr_df_dm_NA0),
     breaks = seq(-7,0,0.1),
     #freq = F,
     main='',
     col=adjustcolor('gold',alpha.f = 0.5),
     add = T)
#abline(v=mean(log10(percent_SNPs_H1_3NP_per_chr_df_dm_NA0)),col='gold',lwd=3) 
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')

hist(log10(percent_SNPs_F123_3NP_per_chr_df_dm_NA0),
     breaks = seq(-7,0,0.1),
     ylim = c(0,0.9),
     xlab='log10( % observed SNPs per chromosome per sample)',
     freq = F,
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
#abline(v=mean(log10(percent_SNPs_F123_3NP_per_chr_df_dm_NA0),na.rm=T),col='grey',lwd=3) 

hist(log10(percent_SNPs_H1_3NP_per_chr_df_dm_NA0),
     breaks = seq(-7,0,0.1),
     freq = F,
     main='',
     col=adjustcolor('gold',alpha.f = 0.5),
     add = T)
#abline(v=mean(log10(percent_SNPs_H1_3NP_per_chr_df_dm_NA0)),col='gold',lwd=3) 
legend("topright",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')




#################################################################################
# nr chrom observed per sample
nr_obs_chr_per_sample_F123 = unlist(apply(number_SNPs_F123_3NP_per_chr_df_dm_NA0,2,function(col){return(sum(sign(col)))}))
nr_obs_chr_per_sample_H1 = unlist(apply(number_SNPs_H1_3NP_per_chr_dm_NA0,2,function(col){return(sum(sign(col)))}))

# table(nr_obs_chr_per_sample_F123==19)
# FALSE  TRUE 
# 33   986
# table(nr_obs_chr_per_sample_H1==22)
# FALSE  TRUE 
# 24   357 

# hist(nr_obs_chr_per_sample_F123,
#      breaks = seq(0,max(nr_obs_chr_per_sample_F123),1),
#      main='',
#      col=adjustcolor('grey',alpha.f = 0.5),
#      cex=cex_plot,
#      cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
#      las=1)
# hist(nr_obs_chr_per_sample_H1,
#      breaks = seq(0,max(nr_obs_chr_per_sample_H1),1),
#      main='',
#      col=adjustcolor('gold',alpha.f = 0.5),
#      cex=cex_plot,
#      cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
#      las=1)


hist(nr_obs_chr_per_sample_H1,
     breaks = seq(0,max(nr_obs_chr_per_sample_H1),1),
     ylim = c(0,1),
     main='',
     xlab = 'number of chr with obs SNPs per sample',
     freq = F,
     col=adjustcolor('gold',alpha.f = 0.5),
     cex=cex_plot,
     cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
     las=1)
hist(nr_obs_chr_per_sample_F123,
     breaks = seq(0,max(nr_obs_chr_per_sample_F123),1),
     main='',
     col=adjustcolor('grey',alpha.f = 0.5),
     freq = F,
     add = T)
legend("topleft",fill=c(adjustcolor("grey",0.5),adjustcolor("gold",0.5)),legend=c("F123 3NP (1019)","H1 3NP (381)"),bty='n')





  
################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################
# look at some fischgrätenplots aka up and down facing read counts per sample per chr

# for that I need to merge the processed SNPs with the truth (F132 at least, for H1 it is already done)
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/prelim/clean_mm10_gt_merged_mappability_minmap1.rds /Users/jmarkow/Desktop/GAM/data/
truth_f123 = readRDS('/Users/jmarkow/Desktop/GAM/data/clean_mm10_gt_merged_mappability_minmap1.rds')
processed_SNVs_F123_3NP_merged_truth = merge(x = processed_SNVs_F123_3NP,y = truth_f123[,c('CHROM','POS','F123_GT')],by=c('CHROM','POS'),all.x=T,all.y=F)
processed_SNVs_F123_3NP_merged_truth = processed_SNVs_F123_3NP_merged_truth[order(processed_SNVs_F123_3NP_merged_truth$POS),]
processed_SNVs_F123_3NP_merged_truth = processed_SNVs_F123_3NP_merged_truth[order(as.numeric(substring(processed_SNVs_F123_3NP_merged_truth$CHROM, 4))),]

processed_SNVs_F123_3NP_merged_truth$CAST_COUNT=0
processed_SNVs_F123_3NP_merged_truth$J129_COUNT=0
# CAST = 1|0 
processed_SNVs_F123_3NP_merged_truth$CAST_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='0|1']=processed_SNVs_F123_3NP_merged_truth$REF_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='0|1']
processed_SNVs_F123_3NP_merged_truth$CAST_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='1|0']=processed_SNVs_F123_3NP_merged_truth$ALT_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='1|0']
# J129 = 0|1 
processed_SNVs_F123_3NP_merged_truth$J129_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='0|1']=processed_SNVs_F123_3NP_merged_truth$ALT_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='0|1']
processed_SNVs_F123_3NP_merged_truth$J129_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='1|0']=processed_SNVs_F123_3NP_merged_truth$REF_COUNT[processed_SNVs_F123_3NP_merged_truth$F123_GT=='1|0']

processed_SNVs_F123_3NP_merged_truth$TOTAL_COUNT=processed_SNVs_F123_3NP_merged_truth$REF_COUNT + processed_SNVs_F123_3NP_merged_truth$ALT_COUNT

processed_SNVs_F123_3NP_merged_truth$sign=sign(processed_SNVs_F123_3NP_merged_truth$CAST_COUNT)
processed_SNVs_F123_3NP_merged_truth$sign[processed_SNVs_F123_3NP_merged_truth$CAST_COUNT==0]=-1

###
chromsize_mm10=read.table('/Users/jmarkow/Desktop/RNAseq/data/chromsize.tsv',header = T)
chromsize_hg38=read.table('/Users/jmarkow/Desktop/GAM/data/hg38_sorted_autosomes.bed',header = F)
chromsize_hg38 = chromsize_hg38[,c(1,3)]
colnames(chromsize_hg38) = c("chr","size")
###    



#pdf("/Users/jmarkow/Desktop/GAM/figures/F123_3NP_CAST_J129_reacounts.pdf")
colours_allele = c(rgb(251,176,59,maxColorValue = 255),rgb(158,0,93,maxColorValue = 255)) # CAST, J129
cex_plot=1
cex_rel=1.2#0.8
par(cex=cex_plot)
par(las=1)
barwidth=0.1

processed_SNVs_F123_3NP_merged_truth_chr = split(processed_SNVs_F123_3NP_merged_truth,processed_SNVs_F123_3NP_merged_truth$CHROM)
processed_SNVs_F123_3NP_merged_truth_chr = processed_SNVs_F123_3NP_merged_truth_chr[c(paste('chr',1:19,sep=''))]
par(mfrow=c(4,4))
i=1
for (chr in processed_SNVs_F123_3NP_merged_truth_chr){
    samples_raw_chr = unique(chr$SAMPLE)
    samples_sort_chr = samples_raw_chr[order(as.numeric(substring(samples_raw_chr, 7)))] 
    for (sample in samples_sort_chr){
        plot_sample=chr[chr$SAMPLE==sample,]
        max_y = max(abs(plot_sample$TOTAL_COUNT))
        
        print(paste(i,sample,sep=" : "))
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_mm10[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        # plot the bars (first ref_count, then alt_count) at the genomic position width the height corresponding to the amount of readcounts at that SNP
        #axis(1,at=c(0,50000000,100000000,150000000,200000000),labels=c(0,50,100,150,200))
        axis(1,at=c(0,1e7,2e7,3e7,4e7,5e7,6e7),labels=seq(0,60,10),cex.axis=cex_rel)
        axis(2,at=c(-100,-50,0,50,100),labels=c(100,50,0,50,100),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$CAST_COUNT),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$J129_COUNT, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        
        # label the haplotypes 
        text(x=(chromsize_mm10[i,2]/2),y=-(max_y),labels = 'CAST',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_mm10[i,2]/2),y=max_y,labels = 'J129',cex=1.1,col = colours_allele[2])
        
    }
    i=i+1
}
#dev.off()

###
# compare this to the H1 3NPs
inch=0.393701
processed_SNVs_H1_3NP_basic = processed_SNVs_H1_3NP[,c("CHROM","POS","SAMPLE","TOTAL_COUNT","p1","p2")]
#pdf("/Users/jmarkow/Desktop/GAM/figures/H1_P1_P2_reacounts.pdf",width=15*inch,,height=15*inch)
#cex_value=1.75
#cex=1.75
par(mar=c(3,3,0.1,0.1)) #bltr # margins
par(mgp=c(1.8,0.5,0)) #titel,label,line # axis placement
cex_plot=1
cex_rel=1.2#0.8
par(cex=cex_plot)
par(las=1)
barwidth=0.1

colours_allele = c(rgb(251,176,59,maxColorValue = 255),rgb(158,0,93,maxColorValue = 255)) # P1, P2


#processed_SNVs_H1_3NP_basic_chr = split(processed_SNVs_H1_3NP_basic,processed_SNVs_H1_3NP_basic$CHROM)
#processed_SNVs_H1_3NP_basic_chr = processed_SNVs_H1_3NP_basic_chr[c(paste('chr',1:22,sep=''))]
par(mfrow=c(2,2))
i=1
for (chr in processed_SNVs_H1_3NP_basic_chr){
    samples_raw_chr = unique(chr$SAMPLE)
    samples_sort_chr = samples_raw_chr[order(as.numeric(substring(samples_raw_chr, 7)))] 
    samples_sort_chr = samples_sort_chr[1:20]
    for (sample in samples_sort_chr){
        plot_sample=chr[chr$SAMPLE==sample,]
        max_y = max(abs(plot_sample$TOTAL_COUNT))
        
        print(paste(i,sample,sep=" : "))
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_hg38[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             
             
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        # plot the bars (first ref_count, then alt_count) at the genomic position width the height corresponding to the amount of readcounts at that SNP
        #axis(1,at=c(0,50000000,100000000,150000000,200000000),labels=c(0,50,100,150,200))
        # axis(1,at=c(0,1e7,2e7,3e7,4e7,5e7,6e7),labels=seq(0,60,10),cex.axis=cex_rel)
        # axis(2,at=c(-100,-50,0,50,100),labels=c(100,50,0,50,100),cex.axis=cex_rel)
        # max(processed_SNVs_H1_3NP_basic$TOTAL_COUNT) # 138
        axis(1,at=seq(0,chromsize_hg38[i,2]+1e7,1e7),labels=seq(0,round((chromsize_hg38[i,2]+1e7)/1e6),10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$p1),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$p2, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        
        # label the haplotypes 
        text(x=(chromsize_hg38[i,2]/2),y=-(max_y),labels = 'P1',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_hg38[i,2]/2),y=max_y,labels = 'P2',cex=1.1,col = colours_allele[2])
        
    }
    i=i+1
}
#dev.off()
sample="sample4"
par(mfrow=c(1,1))



################################################################################
################################################################################
# merge the processed SNP data with the known and calculated haplotypes for H1 and F123

#GAMIBHEAR phasing results
H1_GAMIBHEAR_PHASE=readRDS("/fast/groups/ag_schwarz/Projects/project-gam/results/mara-scripts/H1_results/H1_compiled_graph_haplotypes_new_parameter_lower_bound.rds")
F123_3NP_GAMIBHEAR_PHASE=readRDS("/fast/groups/ag_schwarz/Projects/project-gam/results/mara-scripts/dataframes/large_3NP_graph_haplotypes_new_parameter_lower_bound.rds")
F123_1NP_GAMIBHEAR_PHASE=readRDS("/fast/groups/ag_schwarz/Projects/project-gam/results/phasing_1123_gd_t.rds")

# processed GAM SNP data
processed_SNVs_F123_1NP=readRDS("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/processed_F123_SNPs.rds") # 35879529        7
processed_SNVs_F123_3NP=readRDS('/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_3NP_mm10_all_files_and_processed_observed_SNPs.rds')
processed_SNVs_F123_3NP = processed_SNVs_F123_3NP[[2]] # 60919164        7
processed_SNVs_H1_3NP = readRDS('/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/H1_dixon_haplotypes_liftover_hg38_filter_phased_obs_ref_alt_mapp.rds')
processed_SNVs_H1_3NP_basic = processed_SNVs_H1_3NP[,c("CHROM","POS","REF_COUNT","ALT_COUNT","TOTAL_COUNT","SAMPLE","GT","p1","p2")]

# merge in GAMIBHEAR results for H1
processed_SNVs_H1_3NP_basic_chr = split(processed_SNVs_H1_3NP_basic,processed_SNVs_H1_3NP_basic$CHROM)
processed_SNVs_H1_3NP_basic_chr = processed_SNVs_H1_3NP_basic_chr[c(paste('chr',1:22,sep=''))]

# names(processed_SNVs_H1_3NP_basic_chr)
# [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22"
# > names(H1_GAMIBHEAR_PHASE)
# [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22"

processed_SNVs_H1_3NP_basic_chr_GAMIBHEAR = lapply(c(1:22),function(chr){
    merged = merge(processed_SNVs_H1_3NP_basic_chr[[chr]],H1_GAMIBHEAR_PHASE[[chr]][,c("CHROM","POS","GT")],by=c("CHROM","POS"))
    colnames(merged) = c("CHROM","POS","REF_COUNT","ALT_COUNT","TOTAL_COUNT","SAMPLE","GT_dix","p1_dix","p2_dix","GT_gam")
    merged$p1_gam=0
    merged$p2_gam=0
    merged$p1_gam[merged$GT_gam=='0|1']=merged$REF_COUNT[merged$GT_gam=='0|1']
    merged$p1_gam[merged$GT_gam=='1|0']=merged$ALT_COUNT[merged$GT_gam=='1|0']
    merged$p2_gam[merged$GT_gam=='0|1']=merged$ALT_COUNT[merged$GT_gam=='0|1']
    merged$p2_gam[merged$GT_gam=='1|0']=merged$REF_COUNT[merged$GT_gam=='1|0']
    merged = merged[order(merged$POS),]

    return(merged)
})

# merge in truth and GAMIBHEAR results for F123
truth_f123 = readRDS('/fast/groups/ag_schwarz/Projects/project-gam/prelim/clean_mm10_gt_merged_mappability_minmap1.rds')
truth_f123_chr = split(truth_f123,truth_f123$CHROM)
truth_f123_chr = truth_f123_chr[c(paste('chr',1:19,sep=''))]

processed_SNVs_F123_3NP_chr =  split(processed_SNVs_F123_3NP,processed_SNVs_F123_3NP$CHROM)
processed_SNVs_F123_3NP_chr = processed_SNVs_F123_3NP_chr[c(paste('chr',1:19,sep=''))]

# names(processed_SNVs_F123_3NP_chr)
# [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" 
# names(truth_f123_chr)
# [1] "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" 

processed_SNVs_F123_3NP_chr_GAMIBHEAR = lapply(c(1:19),function(chr){
    merged = merge(processed_SNVs_F123_3NP_chr[[chr]],truth_f123_chr[[chr]][,c("CHROM","POS","F123_GT")],by=c("CHROM","POS"))
    merged = merge(merged,F123_1NP_GAMIBHEAR_PHASE[[chr]][,c("CHROM","POS","GT")],by=c("CHROM","POS"))
    merged = merge(merged,F123_3NP_GAMIBHEAR_PHASE[[chr]][,c("CHROM","POS","GT")],by=c("CHROM","POS"))

    colnames(merged) = c("CHROM","POS","REF","ALT","REF_COUNT","ALT_COUNT","SAMPLE","F123_GT","GT_GAM_1NP","GT_GAM_3NP")
    merged$CAST_COUNT_F123=0
    merged$J129_COUNT_F123=0
    merged$CAST_COUNT_F123[merged$F123_GT=='0|1']=merged$REF_COUNT[merged$F123_GT=='0|1']
    merged$CAST_COUNT_F123[merged$F123_GT=='1|0']=merged$ALT_COUNT[merged$F123_GT=='1|0']
    merged$J129_COUNT_F123[merged$F123_GT=='0|1']=merged$ALT_COUNT[merged$F123_GT=='0|1']
    merged$J129_COUNT_F123[merged$F123_GT=='1|0']=merged$REF_COUNT[merged$F123_GT=='1|0']

    merged$CAST_COUNT_1NP=0
    merged$J129_COUNT_1NP=0
    merged$CAST_COUNT_1NP[merged$GT_GAM_1NP=='0|1']=merged$REF_COUNT[merged$GT_GAM_1NP=='0|1']
    merged$CAST_COUNT_1NP[merged$GT_GAM_1NP=='1|0']=merged$ALT_COUNT[merged$GT_GAM_1NP=='1|0']
    merged$J129_COUNT_1NP[merged$GT_GAM_1NP=='0|1']=merged$ALT_COUNT[merged$GT_GAM_1NP=='0|1']
    merged$J129_COUNT_1NP[merged$GT_GAM_1NP=='1|0']=merged$REF_COUNT[merged$GT_GAM_1NP=='1|0']

    merged$CAST_COUNT_3NP=0
    merged$J129_COUNT_3NP=0
    merged$CAST_COUNT_3NP[merged$GT_GAM_3NP=='0|1']=merged$REF_COUNT[merged$GT_GAM_3NP=='0|1']
    merged$CAST_COUNT_3NP[merged$GT_GAM_3NP=='1|0']=merged$ALT_COUNT[merged$GT_GAM_3NP=='1|0']
    merged$J129_COUNT_3NP[merged$GT_GAM_3NP=='0|1']=merged$ALT_COUNT[merged$GT_GAM_3NP=='0|1']
    merged$J129_COUNT_3NP[merged$GT_GAM_3NP=='1|0']=merged$REF_COUNT[merged$GT_GAM_3NP=='1|0']
    
    merged$TOTAL_COUNT=merged$REF_COUNT + merged$ALT_COUNT
    merged = merged[order(merged$POS),]
    
    return(merged)
})

#  1NP also
processed_SNVs_F123_1NP_chr =  split(processed_SNVs_F123_1NP,processed_SNVs_F123_1NP$CHROM)
processed_SNVs_F123_1NP_chr = processed_SNVs_F123_1NP_chr[c(paste('chr',1:19,sep=''))]

processed_SNVs_F123_1NP_chr_GAMIBHEAR = lapply(c(1:19),function(chr){
    merged = merge(processed_SNVs_F123_1NP_chr[[chr]],truth_f123_chr[[chr]][,c("CHROM","POS","F123_GT")],by=c("CHROM","POS"))
    merged = merge(merged,F123_1NP_GAMIBHEAR_PHASE[[chr]][,c("CHROM","POS","GT")],by=c("CHROM","POS"))
    merged = merge(merged,F123_3NP_GAMIBHEAR_PHASE[[chr]][,c("CHROM","POS","GT")],by=c("CHROM","POS"))
    
    colnames(merged) = c("CHROM","POS","REF","ALT","REF_COUNT","ALT_COUNT","SAMPLE","F123_GT","GT_GAM_1NP","GT_GAM_3NP")
    merged$CAST_COUNT_F123=0
    merged$J129_COUNT_F123=0
    merged$CAST_COUNT_F123[merged$F123_GT=='0|1']=merged$REF_COUNT[merged$F123_GT=='0|1']
    merged$CAST_COUNT_F123[merged$F123_GT=='1|0']=merged$ALT_COUNT[merged$F123_GT=='1|0']
    merged$J129_COUNT_F123[merged$F123_GT=='0|1']=merged$ALT_COUNT[merged$F123_GT=='0|1']
    merged$J129_COUNT_F123[merged$F123_GT=='1|0']=merged$REF_COUNT[merged$F123_GT=='1|0']
    
    merged$CAST_COUNT_1NP=0
    merged$J129_COUNT_1NP=0
    merged$CAST_COUNT_1NP[merged$GT_GAM_1NP=='0|1']=merged$REF_COUNT[merged$GT_GAM_1NP=='0|1']
    merged$CAST_COUNT_1NP[merged$GT_GAM_1NP=='1|0']=merged$ALT_COUNT[merged$GT_GAM_1NP=='1|0']
    merged$J129_COUNT_1NP[merged$GT_GAM_1NP=='0|1']=merged$ALT_COUNT[merged$GT_GAM_1NP=='0|1']
    merged$J129_COUNT_1NP[merged$GT_GAM_1NP=='1|0']=merged$REF_COUNT[merged$GT_GAM_1NP=='1|0']
    
    merged$CAST_COUNT_3NP=0
    merged$J129_COUNT_3NP=0
    merged$CAST_COUNT_3NP[merged$GT_GAM_3NP=='0|1']=merged$REF_COUNT[merged$GT_GAM_3NP=='0|1']
    merged$CAST_COUNT_3NP[merged$GT_GAM_3NP=='1|0']=merged$ALT_COUNT[merged$GT_GAM_3NP=='1|0']
    merged$J129_COUNT_3NP[merged$GT_GAM_3NP=='0|1']=merged$ALT_COUNT[merged$GT_GAM_3NP=='0|1']
    merged$J129_COUNT_3NP[merged$GT_GAM_3NP=='1|0']=merged$REF_COUNT[merged$GT_GAM_3NP=='1|0']
    
    merged$TOTAL_COUNT=merged$REF_COUNT + merged$ALT_COUNT
    merged = merged[order(merged$POS),]
    
    return(merged)
})

length(unique(processed_SNVs_F123_1NP$SAMPLE))
length(unique(processed_SNVs_F123_3NP$SAMPLE))
length(unique(processed_SNVs_H1_3NP_basic$SAMPLE))

saveRDS(processed_SNVs_F123_1NP_chr_GAMIBHEAR,"/fast/groups/ag_schwarz/Projects/project-gam/results/F123_1123_1NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds")
saveRDS(processed_SNVs_F123_3NP_chr_GAMIBHEAR,"/fast/groups/ag_schwarz/Projects/project-gam/results/F123_1019_3NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds")
saveRDS(processed_SNVs_H1_3NP_basic_chr_GAMIBHEAR,"/fast/groups/ag_schwarz/Projects/project-gam/results/H1_381_3NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds")

# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/results/F123_1123_1NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds /Users/jmarkow/Desktop/GAM/results/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/results/F123_1019_3NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds /Users/jmarkow/Desktop/GAM/results/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/results/H1_381_3NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds /Users/jmarkow/Desktop/GAM/results/

###############################################################################
################################################################################

processed_SNVs_F123_1NP_chr_GAMIBHEAR = readRDS("/Users/jmarkow/Desktop/GAM/results/F123_1123_1NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds")
processed_SNVs_F123_3NP_chr_GAMIBHEAR = readRDS("/Users/jmarkow/Desktop/GAM/results/F123_1019_3NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds")
processed_SNVs_H1_3NP_basic_chr_GAMIBHEAR = readRDS("/Users/jmarkow/Desktop/GAM/results/H1_381_3NP_processed_SNPs_merged_known_and_predicted_GT_110621.rds")

###
chromsize_mm10=read.table('/Users/jmarkow/Desktop/RNAseq/data/chromsize.tsv',header = T)
chromsize_hg38=read.table('/Users/jmarkow/Desktop/GAM/data/hg38_sorted_autosomes.bed',header = F)
chromsize_hg38 = chromsize_hg38[,c(1,3)]
colnames(chromsize_hg38) = c("chr","size")
###    

###
#  H1 3NPs
inch=0.393701
#pdf("/Users/jmarkow/Desktop/GAM/figures/H1_P1_P2_readcounts.pdf",width=30*inch,,height=15*inch)

#cex_value=1.75
#cex=1.75
par(mar=c(3,3,0.1,0.1)) #bltr # margins
par(mgp=c(1.8,0.5,0)) #titel,label,line # axis placement
cex_plot=1
cex_rel=1.2#0.8
par(cex=cex_plot)
par(las=1)
barwidth=0.1

colours_allele = c(rgb(251,176,59,maxColorValue = 255),rgb(158,0,93,maxColorValue = 255)) # P1, P2

#par(mfrow=c(1,2))
par(mfrow=c(2,1))

i=1
for (chr in processed_SNVs_H1_3NP_basic_chr_GAMIBHEAR){
    samples_raw_chr = unique(chr$SAMPLE)
    samples_sort_chr = samples_raw_chr[order(as.numeric(substring(samples_raw_chr, 7)))] 
    #samples_sort_chr = samples_sort_chr[1:20]
    for (sample in samples_sort_chr){
        # sample = "sample4"
        # chr = processed_SNVs_H1_3NP_basic_chr_GAMIBHEAR[[1]]
        plot_sample=chr[chr$SAMPLE==sample,]
        max_y = max(abs(plot_sample$TOTAL_COUNT))
        
        print(paste(i,sample,sep=" : "))
        # plot known SNPs frist
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_hg38[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,chromsize_hg38[i,2]+1e7,1e7),labels=seq(0,round((chromsize_hg38[i,2]+1e7)/1e6),10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$p1_dix),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$p2_dix, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_hg38[i,2]/2),y=-(max_y),labels = 'P1 Dixon',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_hg38[i,2]/2),y=max_y,labels = paste('P2 Dixon,chr',i,sample,sep=' '),cex=1.1,col = colours_allele[2])
        
        # plot predicted haplotype later
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_hg38[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,chromsize_hg38[i,2]+1e7,1e7),labels=seq(0,round((chromsize_hg38[i,2]+1e7)/1e6),10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$p1_gam),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$p2_gam, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_hg38[i,2]/2),y=-(max_y),labels = 'P1 GAM',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_hg38[i,2]/2),y=max_y,labels = 'P2 GAM',cex=1.1,col = colours_allele[2])
        
    }
    i=i+1
}
#dev.off()

################################################################################

###
#  F123 3NPs
inch=0.393701
#pdf("/Users/jmarkow/Desktop/GAM/figures/F123_3NP_CAST_J129_readcounts.pdf",width=30*inch,height=22.5*inch)

#cex_value=1.75
#cex=1.75
par(mar=c(3,3,0.1,0.1)) #bltr # margins
par(mgp=c(1.8,0.5,0)) #titel,label,line # axis placement
cex_plot=1
cex_rel=1.2#0.8
par(cex=cex_plot)
par(las=1)
barwidth=0.1

colours_allele = c(rgb(251,176,59,maxColorValue = 255),rgb(158,0,93,maxColorValue = 255)) # P1, P2

#par(mfrow=c(1,2))
par(mfrow=c(3,1))

i=1
for (chr in processed_SNVs_F123_3NP_chr_GAMIBHEAR){
    samples_raw_chr = unique(chr$SAMPLE)
    samples_sort_chr = samples_raw_chr[order(as.numeric(substring(samples_raw_chr, 7)))] 
    for (sample in samples_sort_chr){
        plot_sample=chr[chr$SAMPLE==sample,]
        max_y = max(abs(plot_sample$TOTAL_COUNT))
        
        print(paste(i,sample,sep=" : "))
        # plot known SNPs frist
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_mm10[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,round(chromsize_mm10[i,2]+1e7,-7),1e7),labels=seq(0,round(chromsize_mm10[i,2]+1e7,-7)/1e6,10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$CAST_COUNT_F123),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$J129_COUNT_F123, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_mm10[i,2]/2),y=-(max_y),labels = 'CAST F123',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_mm10[i,2]/2),y=max_y,labels = paste('J129 F123, chr',i,sample,sep=' '),cex=1.1,col = colours_allele[2])
        
        # plot predicted haplotype later
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_mm10[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,round(chromsize_mm10[i,2]+1e7,-7),1e7),labels=seq(0,round(chromsize_mm10[i,2]+1e7,-7)/1e6,10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$CAST_COUNT_1NP),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$J129_COUNT_1NP, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_mm10[i,2]/2),y=-(max_y),labels = 'CAST 1NP',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_mm10[i,2]/2),y=max_y,labels = 'J129 1NP',cex=1.1,col = colours_allele[2])
        
        # plot predicted haplotype later
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_mm10[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,round(chromsize_mm10[i,2]+1e7,-7),1e7),labels=seq(0,round(chromsize_mm10[i,2]+1e7,-7)/1e6,10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$CAST_COUNT_3NP),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$J129_COUNT_3NP, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_mm10[i,2]/2),y=-(max_y),labels = 'CAST 3NP',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_mm10[i,2]/2),y=max_y,labels = 'J129 3NP',cex=1.1,col = colours_allele[2])
        
    }
    i=i+1
}
#dev.off()




################################################################################


###
#  F123 1NPs
inch=0.393701
#pdf("/Users/jmarkow/Desktop/GAM/figures/F123_1NP_CAST_J129_readcounts.pdf",width=30*inch,height=22.5*inch)

#cex_value=1.75
#cex=1.75
par(mar=c(3,3,0.1,0.1)) #bltr # margins
par(mgp=c(1.8,0.5,0)) #titel,label,line # axis placement
cex_plot=1
cex_rel=1.2#0.8
par(cex=cex_plot)
par(las=1)
barwidth=0.1

colours_allele = c(rgb(251,176,59,maxColorValue = 255),rgb(158,0,93,maxColorValue = 255)) # P1, P2

#par(mfrow=c(1,2))
par(mfrow=c(3,1))

i=1
for (chr in processed_SNVs_F123_1NP_chr_GAMIBHEAR){
    samples_raw_chr = unique(chr$SAMPLE)
    samples_sort_chr = samples_raw_chr[order(as.numeric(substring(samples_raw_chr, 7)))] 
    for (sample in samples_sort_chr){
        plot_sample=chr[chr$SAMPLE==sample,]
        max_y = max(abs(plot_sample$TOTAL_COUNT))
        
        print(paste(i,sample,sep=" : "))
        # plot known SNPs frist
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_mm10[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,round(chromsize_mm10[i,2]+1e7,-7),1e7),labels=seq(0,round(chromsize_mm10[i,2]+1e7,-7)/1e6,10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$CAST_COUNT_F123),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$J129_COUNT_F123, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_mm10[i,2]/2),y=-(max_y),labels = 'CAST F123',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_mm10[i,2]/2),y=max_y,labels = paste('J129 F123, chr',i,sample,sep=' '),cex=1.1,col = colours_allele[2])
        
        # plot predicted haplotype later
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_mm10[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,round(chromsize_mm10[i,2]+1e7,-7),1e7),labels=seq(0,round(chromsize_mm10[i,2]+1e7,-7)/1e6,10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$CAST_COUNT_1NP),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$J129_COUNT_1NP, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_mm10[i,2]/2),y=-(max_y),labels = 'CAST 1NP',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_mm10[i,2]/2),y=max_y,labels = 'J129 1NP',cex=1.1,col = colours_allele[2])
        
        # plot predicted haplotype later
        plot(1, type="n", axes=T, 
             xlim=c(0,chromsize_mm10[i,2]),ylim=c(-(max_y),max_y), #,ylim=c(-(max_y),max_y*2),
             xlab="Genomic position (Mb)", 
             ylab="Read counts",
             main='',#paste('read count per haplotype at genomic position,\n',sample,chr,sep=' '),
             xaxt='n',yaxt='n',
             cex=cex_plot,
             cex.lab=cex_rel, cex.axis=cex_rel, #cex.main=cex_rel,
             las=1)
        axis(1,at=seq(0,round(chromsize_mm10[i,2]+1e7,-7),1e7),labels=seq(0,round(chromsize_mm10[i,2]+1e7,-7)/1e6,10),cex.axis=cex_rel)
        axis(2,at=seq(-150,150,50),labels=seq(-150,150,50),cex.axis=cex_rel)
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = -(plot_sample$CAST_COUNT_3NP),
             col=colours_allele[1],border=colours_allele[1])
        rect(xleft = plot_sample$POS - barwidth, 
             xright = plot_sample$POS + barwidth,
             ybottom = 0,
             ytop = plot_sample$J129_COUNT_3NP, 
             col=colours_allele[2],border=colours_allele[2])
        abline(h=0,lwd=0.5,col='grey70')
        text(x=(chromsize_mm10[i,2]/2),y=-(max_y),labels = 'CAST 3NP',cex=1.1,col = colours_allele[1]) 
        text(x=(chromsize_mm10[i,2]/2),y=max_y,labels = 'J129 3NP',cex=1.1,col = colours_allele[2])
        
    }
    i=i+1
}
#dev.off()




################################################################################
################################################################################

################################################################################
################################################################################
