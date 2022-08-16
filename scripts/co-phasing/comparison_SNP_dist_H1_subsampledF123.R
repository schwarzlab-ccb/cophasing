# was bisher geschah:
# cophasing H1, version 1, basic assignment without considering distance between reads and snps
# cophasing H1, version 2, dist cutoff around SNPs, 1kb, 50kb, 100kb, 200kb
# cophasing downsampled F123:
#   downsampled F123 SNPs to human SNP density, which is 1kb (10times), cophased those 10 iterations with version 2, dist cutoff: 50kb, 100kb, 200kb 
#
# now the subsampled F123 cophasing results look astonishingly good, which could be due to the SNP distribution. If we have a SNP every 1kb and the dist cutoff is min 50kb, we would expect good results
# compare SNP distance distribution between H1 and the 10 subsampled F123, maybe the distance between subsampled F123 is more steady, the distribution in H1 could be more brought/wide with lots of SNPs with an actual higher distance

# ------------------------------------------------------------------------------
### load known SNPs (independent of observation)
# H1 (clean for mappability and phasing)
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/H1_dixon_haplotypes_liftover_hg38.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/H1/Genotyping/intermediate/mappability/mapinter_v3 /Users/jmarkow/Desktop/GAM/Co-Phasing/data/hg38_mapp

H1_truth_in=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/H1_dixon_haplotypes_liftover_hg38.rds") #1630597       9
H1_truth=H1_truth[which(substr(H1_truth$GT,2,2)=='|'),] #1563859       9

library(data.table)
mapp=fread('/Users/jmarkow/Desktop/GAM/Co-Phasing/data/hg38_mapp')
mappdb=as.data.frame(mapp)
mappdb=mappdb[,c(1,2,7)]
colnames(mappdb)=c('CHROM','POS','MAPPABILITY')

H1_truth_mapp=merge(H1_truth,mappdb,by=c("CHROM","POS")) # 1554248      10
H1_truth_mapp=H1_truth_mapp[which(H1_truth_mapp$MAPPABILITY==1),] # 1453239      10

# > nrow(H1_truth_in)
# [1] 1630597
# > nrow(H1_truth)
# [1] 1563859
# > nrow(H1_truth_mapp)
# [1] 1453239

# F123 subsampled truth
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_human_SNPdensity_10times.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/
F123_subsampled_10=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/F123_mm10_GT_subsampled_to_human_SNPdensity_10times.rds")

# genome sizes
H1_size=read.table("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/hg38_sorted_autosomes.chrom.sizes",header = T)
F123_size=read.table("/fast/groups/ag_schwarz/Projects/project-gam/data/mm10_sorted_autosomes.chrom.sizes")
sum(H1_size$SIZE)
#[1] 2875001522
sum(F123_size$V2)
#[1] 2462745373

#sum(H1_size$SIZE)/nrow(H1_truth_mapp)
#2875001522/1453239 # theoretical SNP distance: genome size divided by number of SNPs
# [1] 1978.34

#nrow(F123_subsampled_10[[1]]) # all 10: 2462745 subsampled SNPs
#sum(F123_size$V2)/nrow(F123_subsampled_10[[1]])
#2462745373/2462745 # theoretical SNP distance: genome size divided by number of SNPs, I subsampled the F123 SNPs to achieve a SNP distance of 1kb
# 1000

### dist hists  ----------------------------------------------------------------
# H1
chroms=paste("chr",1:22,sep="")
H1_truth_mapp = H1_truth_mapp[order(H1_truth_mapp$POS),]
H1_truth_mapp = H1_truth_mapp[order(as.numeric(substring(H1_truth_mapp$CHROM, 4))),]
H1_truth_mapp$CHROM = factor(H1_truth_mapp$CHROM, levels = chroms)
H1_truth_mapp_chr=split(H1_truth_mapp,H1_truth_mapp$CHROM)
dist_list_H1=unlist(lapply(H1_truth_mapp_chr,function(chr){
    return(diff(chr$POS))
}))
# min(dist_list_H1)
# max(dist_list_H1)
# hist(dist_list_H1)
max(log10(dist_list_H1))
H1_hist=hist(log10(dist_list_H1),breaks=100,xlim = c(0,max(log10(dist_list_H1))))

# F123
chroms_F=paste("chr",1:19,sep="")
F123_sub_chr = lapply(F123_subsampled_10,function(iter){
    iter = iter[order(iter$POS),]
    iter = iter[order(as.numeric(substring(iter$CHROM, 4))),]
    iter$CHROM = factor(iter$CHROM, levels = chroms_F)
    iter_chr = split(iter,iter$CHROM)
    return(iter_chr)
})
dist_list_F123 = lapply(F123_sub_chr,function(iter){
    dist = unlist(lapply(iter,function(chr){
        return(diff(chr$POS))
    }))
    return(dist)
})

par(mfrow=c(4,3))
plot(H1_hist)
for (i in 1:length(dist_list_F123)){
    hist(log10(dist_list_F123[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i))
}

par(mfrow=c(1,1))
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between known SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,2e5))
abline(v=mean(log10(dist_list_H1)), col="gold") # 2.686611
for (i in 1:length(dist_list_F123)){
    hist(log10(dist_list_F123[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(dist_list_F123[[i]])), col=i) # 2.595...
}
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")

mean(dist_list_H1) # 1915.676
mean(dist_list_F123[[10]]) # 975
mean(log10(dist_list_H1)) # 2.686611
mean(log10(dist_list_F123[[10]])) # 2.595...


# boxplots
# par(mfrow=c(4,3))
# boxplot(log10(dist_list_H1),horizontal = T)
# for (i in 1:length(dist_list_F123)){
#     boxplot(log10(dist_list_F123[[i]]),
#          main = paste("SNP dist - F123 iter ",i),horizontal = T)
# }

### dist hist of observed SNPs  ------------------------------------------------
# H1
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/H1/AS/prelim/H1_mpileup_truth_merged_ref_alt_phased_checked.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/
H1_processed_SNVs_merged_truth_red=readRDS('/Users/jmarkow/Desktop/GAM/Co-Phasing/data/H1_mpileup_truth_merged_ref_alt_phased_checked.rds')
H1_processed_SNVs_merged_truth_red = H1_processed_SNVs_merged_truth_red[order(H1_processed_SNVs_merged_truth_red$POS),]
H1_processed_SNVs_merged_truth_red = H1_processed_SNVs_merged_truth_red[order(as.numeric(substring(H1_processed_SNVs_merged_truth_red$CHROM, 4))),]
H1_processed_SNVs_merged_truth_red$CHROM = factor(H1_processed_SNVs_merged_truth_red$CHROM, levels = chroms)
H1_processed_SNVs_merged_truth_red_chr=split(H1_processed_SNVs_merged_truth_red,H1_processed_SNVs_merged_truth_red$CHROM)
# split the observations by sample, see how far the distance of observed SNPs is
H1_processed_SNVs_merged_truth_red_chr_sample = lapply(H1_processed_SNVs_merged_truth_red_chr,function(chr){
    return(split(chr,chr$SAMPLE))
})

dist_list_obs_H1_sample = unlist(lapply(H1_processed_SNVs_merged_truth_red_chr_sample,function(chr){
    sample_dist = unlist(lapply(chr,function(sample){
        return(diff(sample$POS))
    }))
    return(sample_dist)
}))
length(dist_list_obs_H1_sample) # 1824600
max(dist_list_obs_H1_sample) # 203440810
mean(dist_list_obs_H1_sample) # 526705.2
median(dist_list_obs_H1_sample) # 20660
H1_hist_obs_sample = hist(log10(dist_list_obs_H1_sample),breaks=50,xlim = c(0,max(log10(dist_list_obs_H1_sample))))
abline(v=log10(mean(dist_list_obs_H1_sample)),col="red")
# abline(v=log10(median(dist_list_obs_H1_sample)),col="red")

# F123 subsampled
# on cluster
subsampled_F123_list = readRDS('/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_human_SNPdensity_10times.rds')
# libraries
library(GAMIBHEAR)
library(stringr)
library(data.table)

chrom_order=paste("chr",c(1:19),sep="")

# load the F123 mpileup VCFs again, to be sure the correct samples are assigned to the correct filenames
vcf_dir = '/fast/groups/ag_schwarz/Projects/project-gam/data/F123_bam_files_1123_1NP_mm10_Sept_2019_last_base_trimmed_mpileup'
files = list.files(path = vcf_dir, pattern = '*.mpileup.vcf$')
VCF_list = load_sample_VCF(files,vcf_dir)
processed_SNVs = process_sample_VCF(VCF_list)
#saveRDS(processed_SNVs,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/processed_F123_SNPs.rds")
# merge the SNPs with the 10 subsampled F123 genomes
processed_SNVs_merged_subsampled = lapply(subsampled_F123_list,function(i){
    obs_truth_merged = merge(processed_SNVs,i,by=c("CHROM","POS"))
    obs_truth_merged = obs_truth_merged[order(as.numeric(obs_truth_merged$POS)),]
    obs_truth_merged = obs_truth_merged[order(match(obs_truth_merged$CHROM,chrom_order)),]
    return(obs_truth_merged)
})

processed_SNVs_merged_subsampled_chr_sample_pos_dist = lapply(processed_SNVs_merged_subsampled,function(iter){
    per_chr = split(iter,iter$CHROM)
    per_chr_sample = unlist(lapply(per_chr,function(chr){
        per_sample = split(chr,chr$SAMPLE)
        sample_dist = unlist(lapply(per_sample,function(sample){
            return(diff(sample$POS))
        }))
        return(sample_dist)
    }))
    return(per_chr_sample)
})

saveRDS(processed_SNVs_merged_subsampled_chr_sample_pos_dist,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_obs_SNP_distance.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_obs_SNP_distance.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_chr_sample_pos_dist=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_obs_SNP_distance.rds")

F123_hist_obs_sample = hist(log10(processed_SNVs_merged_subsampled_chr_sample_pos_dist[[1]]),breaks=50)#,xlim = c(0,max(log10(dist_list_obs_H1_sample))))

par(mfrow=c(4,3))
plot(H1_hist_obs_sample)
for (i in 1:length(processed_SNVs_merged_subsampled_chr_sample_pos_dist)){
    hist(log10(processed_SNVs_merged_subsampled_chr_sample_pos_dist[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_obs_H1_sample))),
         main = paste("SNP dist - F123 iter ",i))
}

# distributions look the same again, but again, y-axis is different, check for number of SNPs observed total
unlist(lapply(processed_SNVs_merged_subsampled_chr_sample_pos_dist,length))
# 4721294 4716737 4722438 4721555 4717483 4722554 4722454 4717153 4719268 4718818
length(dist_list_obs_H1_sample)
# 1824600

par(mfrow=c(1,1))
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between oserved SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,3.5e5))
abline(v=mean(log10(dist_list_obs_H1_sample)), col="gold") # 
for (i in 1:length(processed_SNVs_merged_subsampled_chr_sample_pos_dist)){
    hist(log10(processed_SNVs_merged_subsampled_chr_sample_pos_dist[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_obs_H1_sample))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(processed_SNVs_merged_subsampled_chr_sample_pos_dist[[i]])), col=i) # 
}
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")


# how can the distribution be the same by shape




# check number of observed SNPs per sample in H1 and subsampled F123
# H1
H1_processed_SNVs_merged_truth_red_sample = split(H1_processed_SNVs_merged_truth_red,H1_processed_SNVs_merged_truth_red$SAMPLE)
H1_processed_SNVs_merged_truth_red_SNPs_per_sample = unlist(lapply(H1_processed_SNVs_merged_truth_red_sample,function(sample){
    return(nrow(sample))
}))

# F123
processed_SNVs_merged_subsampled_SNPs_per_sample = lapply(processed_SNVs_merged_subsampled,function(iter){
    per_sample = split(iter,iter$SAMPLE)
    sample_SNPs = unlist(lapply(per_sample,function(sample){
        return(nrow(sample))
        }))
    return(sample_SNPs)
})
saveRDS(processed_SNVs_merged_subsampled_SNPs_per_sample,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_obs_SNP_per_sample.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_obs_SNP_per_sample.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_SNPs_per_sample=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_obs_SNP_per_sample.rds")

par(mfrow=c(4,3))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=30)
for (i in 1:length(processed_SNVs_merged_subsampled_SNPs_per_sample)){
    hist(processed_SNVs_merged_subsampled_SNPs_per_sample[[i]],
         xlim=c(0,30000),
         breaks=30)
}

par(mfrow=c(1,1))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=60,
     ylim=c(0,150),
     col=adjustcolor("gold",alpha.f = 0.3),
     main="number of SNPs obs per sample",
     xlab="number of SNPs")
for (i in 1:length(processed_SNVs_merged_subsampled_SNPs_per_sample)){
    hist(processed_SNVs_merged_subsampled_SNPs_per_sample[[i]],
         xlim=c(0,30000),
         breaks=60,
         col=adjustcolor(i,alpha.f = 0.1),add=T)
}
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     breaks=60,
     col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")


################################################################################
################################################################################
################################################################################
# on the cluster, subsample F123 again to achieve the same SNP distance distribution as in H1

# load full F123
full_SNPs_mm10 = readRDS('/fast/groups/ag_schwarz/Projects/project-gam/prelim/clean_mm10_gt_merged_mappability_minmap1.rds')
nrow(full_SNPs_mm10)
# [1] 18150228
# subsample to 2kb distance
chromsize_mm10=read.table('/fast/groups/ag_schwarz/Projects/project-gam/RNAseq/data/Rieke_current/genome/chromsize.tsv',header = T)
sum(chromsize_mm10$size)/nrow(full_SNPs_mm10)
# 135.6867

# now we want to lower the SNP density to ~ 1 SNP every 2000bp (2kb), means we need x number of SNPs for the mouse genome length
nr_subsampled_SNPs = sum(chromsize_mm10$size)/2000 # 1231373; before with 1kb: 2462745
nr_subsampled_SNPs/nrow(full_SNPs_mm10)*100 # 6.784337; before with 1kb: 13.56867 % 

# I will sample indizes of the rows, so I can easily subset the SNP df
subsampled_F123_list_2kb = vector(mode = "list", length = 10)
for(seed in 1:10){
    set.seed(seed)
    SNPs_mm10_subsampled_indizes = sample(x = nrow(full_SNPs_mm10), size = nr_subsampled_SNPs, replace = FALSE)
    SNPs_mm10_subsampled = full_SNPs_mm10[sort(SNPs_mm10_subsampled_indizes),]
    subsampled_F123_list_2kb[[seed]] = SNPs_mm10_subsampled # 6.784334
    print(nrow(SNPs_mm10_subsampled)/nrow(full_SNPs_mm10)*100) #   0.129 0.871
    print(table(SNPs_mm10_subsampled$F123_GT)/nrow(SNPs_mm10_subsampled))
    SNPs_mm10_subsampled_chr = split(SNPs_mm10_subsampled, SNPs_mm10_subsampled$CHROM)
    mm10_SNP_subsampled_dist=unlist(lapply(SNPs_mm10_subsampled_chr,function(chr){
        dist=diff(chr$POS)
    }))
    print(mean(mm10_SNP_subsampled_dist)) # ca 1950
    print(median(mm10_SNP_subsampled_dist)) # 991
}

saveRDS(subsampled_F123_list_2kb,'/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_2kb_dist_10times.rds')
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_2kb_dist_10times.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/

# merge with processed SNPs 
processed_SNVs = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/processed_F123_SNPs.rds")

processed_SNVs_merged_subsampled_2kb = lapply(subsampled_F123_list_2kb,function(i){
    obs_truth_merged = merge(processed_SNVs,i,by=c("CHROM","POS"))
    obs_truth_merged = obs_truth_merged[order(as.numeric(obs_truth_merged$POS)),]
    obs_truth_merged = obs_truth_merged[order(match(obs_truth_merged$CHROM,chrom_order)),]
    return(obs_truth_merged)
})

processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist = lapply(processed_SNVs_merged_subsampled_2kb,function(iter){
    per_chr = split(iter,iter$CHROM)
    per_chr_sample = unlist(lapply(per_chr,function(chr){
        per_sample = split(chr,chr$SAMPLE)
        sample_dist = unlist(lapply(per_sample,function(sample){
            return(diff(sample$POS))
        }))
        return(sample_dist)
    }))
    return(per_chr_sample)
})

saveRDS(processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_2kb_obs_SNP_distance.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_2kb_obs_SNP_distance.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_2kb_obs_SNP_distance.rds")

# plot distance known SNPs (subsample to 2kb distance)
subsampled_F123_list_2kb=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/F123_mm10_GT_subsampled_to_2kb_dist_10times.rds")

chroms_F=paste("chr",1:19,sep="")
F123_sub_2kb_chr = lapply(subsampled_F123_list_2kb,function(iter){
    iter = iter[order(iter$POS),]
    iter = iter[order(as.numeric(substring(iter$CHROM, 4))),]
    iter$CHROM = factor(iter$CHROM, levels = chroms_F)
    iter_chr = split(iter,iter$CHROM)
    return(iter_chr)
})
dist_list_F123_2kb = lapply(F123_sub_2kb_chr,function(iter){
    dist = unlist(lapply(iter,function(chr){
        return(diff(chr$POS))
    }))
    return(dist)
})

par(mfrow=c(4,3))
plot(H1_hist)
for (i in 1:length(dist_list_F123_2kb)){
    hist(log10(dist_list_F123_2kb[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i))
}

par(mfrow=c(1,1))
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between known SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,2e5))
abline(v=mean(log10(dist_list_H1)), col="gold") # 2.686611
for (i in 1:length(dist_list_F123_2kb)){
    hist(log10(dist_list_F123_2kb[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(dist_list_F123_2kb[[i]])), col=i)
    #print(mean(log10(dist_list_F123_2kb[[i]])))  # 2.913
}
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")

mean(dist_list_H1) # 1915.676
mean(dist_list_F123_2kb[[10]]) # 1951.045
mean(log10(dist_list_H1)) # 2.686611
mean(log10(dist_list_F123_2kb[[10]])) # 2.913996

# plot distance observed SNPs
par(mfrow=c(1,1))
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between oserved SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,3.5e5))
abline(v=mean(log10(dist_list_obs_H1_sample)), col="gold") # 4.05915
for (i in 1:length(processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist)){
    hist(log10(processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_obs_H1_sample))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist[[i]])), col=i) # 4.31
}
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")
for (i in 1:length(processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist)){
    # print(mean(log10(processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist[[i]]))) # 4.31
    print(mean(processed_SNVs_merged_subsampled_2kb_chr_sample_pos_dist[[i]])) # 589602.9 - 599757.8
}
mean(dist_list_obs_H1_sample) # 526705.2

# number of SNPs
# F123 
processed_SNVs_merged_subsampled_2kb_SNPs_per_sample = lapply(processed_SNVs_merged_subsampled_2kb,function(iter){
    per_sample = split(iter,iter$SAMPLE)
    sample_SNPs = unlist(lapply(per_sample,function(sample){
        return(nrow(sample))
    }))
    return(sample_SNPs)
})
saveRDS(processed_SNVs_merged_subsampled_2kb_SNPs_per_sample,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_2kb_obs_SNP_per_sample.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_2kb_obs_SNP_per_sample.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_2kb_SNPs_per_sample=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_2kb_obs_SNP_per_sample.rds")

par(mfrow=c(4,3))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=30)
for (i in 1:length(processed_SNVs_merged_subsampled_2kb_SNPs_per_sample)){
     hist(processed_SNVs_merged_subsampled_2kb_SNPs_per_sample[[i]],
          xlim=c(0,30000),
          breaks=30)
}

par(mfrow=c(1,1))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=60,
     ylim=c(0,300),
     col=adjustcolor("gold",alpha.f = 0.3),
     main="number of SNPs obs per sample",
     xlab="number of SNPs")
for (i in 1:length(processed_SNVs_merged_subsampled_2kb_SNPs_per_sample)){
    hist(processed_SNVs_merged_subsampled_2kb_SNPs_per_sample[[i]],
         xlim=c(0,30000),
         breaks=30,
         col=adjustcolor(i,alpha.f = 0.1),add=T)
}
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     breaks=60,
     col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")

length(unique(H1_processed_SNVs_merged_truth_red$SAMPLE)) # 381

################################################################################
################################################################################
################################################################################
# on the cluster, subsample F123 again to achieve the same SNP distance distribution as in H1
# second try, go for 1978 exactly?
#sum(H1_size$SIZE)/nrow(H1_truth_mapp)
#2875001522/1453239 # theoretical SNP distance: genome size divided by number of SNPs
# [1] 1978.34
#mean(dist_list_H1) # 1915.676


sum(chromsize_mm10$size)/nrow(full_SNPs_mm10)
# 135.6867

# now we want to lower the SNP density to ~ 1 SNP every 2000bp (2kb), means we need x number of SNPs for the mouse genome length
nr_subsampled_SNPs = sum(chromsize_mm10$size)/1978 # 1245068; 1231373; before with 1kb: 2462745
nr_subsampled_SNPs/nrow(full_SNPs_mm10)*100 # 6.859795; 6.784337; before with 1kb: 13.56867 % 

# I will sample indizes of the rows, so I can easily subset the SNP df
subsampled_F123_list_1978 = vector(mode = "list", length = 10)
for(seed in 1:10){
    set.seed(seed)
    SNPs_mm10_subsampled_indizes = sample(x = nrow(full_SNPs_mm10), size = nr_subsampled_SNPs, replace = FALSE)
    SNPs_mm10_subsampled = full_SNPs_mm10[sort(SNPs_mm10_subsampled_indizes),]
    subsampled_F123_list_1978[[seed]] = SNPs_mm10_subsampled # 
    print(nrow(SNPs_mm10_subsampled)/nrow(full_SNPs_mm10)*100) #   
    print(table(SNPs_mm10_subsampled$F123_GT)/nrow(SNPs_mm10_subsampled))
    SNPs_mm10_subsampled_chr = split(SNPs_mm10_subsampled, SNPs_mm10_subsampled$CHROM)
    mm10_SNP_subsampled_dist=unlist(lapply(SNPs_mm10_subsampled_chr,function(chr){
        dist=diff(chr$POS)
    }))
    print(mean(mm10_SNP_subsampled_dist)) # 
    print(median(mm10_SNP_subsampled_dist)) # 
}
# [1] 6.859793
# 
# 0|1       1|0
# 0.1288604 0.8711396
# [1] 1928.493
# [1] 980
saveRDS(subsampled_F123_list_1978,'/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_1978_dist_10times.rds')
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_1978_dist_10times.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/

# merge with processed SNPs 
#processed_SNVs = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/processed_F123_SNPs.rds")

processed_SNVs_merged_subsampled_1978 = lapply(subsampled_F123_list_1978,function(i){
    obs_truth_merged = merge(processed_SNVs,i,by=c("CHROM","POS"))
    obs_truth_merged = obs_truth_merged[order(as.numeric(obs_truth_merged$POS)),]
    obs_truth_merged = obs_truth_merged[order(match(obs_truth_merged$CHROM,chrom_order)),]
    return(obs_truth_merged)
})

processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist = lapply(processed_SNVs_merged_subsampled_1978,function(iter){
    per_chr = split(iter,iter$CHROM)
    per_chr_sample = unlist(lapply(per_chr,function(chr){
        per_sample = split(chr,chr$SAMPLE)
        sample_dist = unlist(lapply(per_sample,function(sample){
            return(diff(sample$POS))
        }))
        return(sample_dist)
    }))
    return(per_chr_sample)
})

saveRDS(processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_1978_obs_SNP_distance.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_1978_obs_SNP_distance.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_1978_obs_SNP_distance.rds")

# plot distance known SNPs (subsample to 1978bp distance)
subsampled_F123_list_1978=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/F123_mm10_GT_subsampled_to_1978_dist_10times.rds")

chroms_F=paste("chr",1:19,sep="")
F123_sub_1978_chr = lapply(subsampled_F123_list_1978,function(iter){
    iter = iter[order(iter$POS),]
    iter = iter[order(as.numeric(substring(iter$CHROM, 4))),]
    iter$CHROM = factor(iter$CHROM, levels = chroms_F)
    iter_chr = split(iter,iter$CHROM)
    return(iter_chr)
})
dist_list_F123_1978 = lapply(F123_sub_1978_chr,function(iter){
    dist = unlist(lapply(iter,function(chr){
        return(diff(chr$POS))
    }))
    return(dist)
})

par(mfrow=c(4,3))
plot(H1_hist)
for (i in 1:length(dist_list_F123_1978)){
    hist(log10(dist_list_F123_1978[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i))
}

par(mfrow=c(1,1))
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between known SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,2e5))
abline(v=mean(log10(dist_list_H1)), col="gold") # 2.686611
for (i in 1:length(dist_list_F123_1978)){
    hist(log10(dist_list_F123_1978[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(dist_list_F123_1978[[i]])), col=i)
    print(mean(log10(dist_list_F123_1978[[i]])))  # 2.90859
}
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")

mean(dist_list_H1) # 1915.676
mean(dist_list_F123_1978[[10]]) # 1929.61
mean(log10(dist_list_H1)) # 2.686611
mean(log10(dist_list_F123_1978[[10]])) # 2.908903

# plot distance observed SNPs
par(mfrow=c(1,1))
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between oserved SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,3.5e5))
abline(v=mean(log10(dist_list_obs_H1_sample)), col="gold") # 4.05915
for (i in 1:length(processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist)){
    hist(log10(processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_obs_H1_sample))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist[[i]])), col=i) # 4.298
}
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")
for (i in 1:length(processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist)){
    print(mean(log10(processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist[[i]]))) # 4.29
    print(mean(processed_SNVs_merged_subsampled_1978_chr_sample_pos_dist[[i]])) # between 583970.6 and 594354.6
}
mean(dist_list_obs_H1_sample) # 526705.2

# number of SNPs
# F123 
processed_SNVs_merged_subsampled_1978_SNPs_per_sample = lapply(processed_SNVs_merged_subsampled_1978,function(iter){
    per_sample = split(iter,iter$SAMPLE)
    sample_SNPs = unlist(lapply(per_sample,function(sample){
        return(nrow(sample))
    }))
    return(sample_SNPs)
})
saveRDS(processed_SNVs_merged_subsampled_1978_SNPs_per_sample,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_2kb_obs_SNP_per_sample.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_1978_obs_SNP_per_sample.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_1978_SNPs_per_sample=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_1978_obs_SNP_per_sample.rds")

par(mfrow=c(4,3))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=30)
for (i in 1:length(processed_SNVs_merged_subsampled_1978_SNPs_per_sample)){
    hist(processed_SNVs_merged_subsampled_1978_SNPs_per_sample[[i]],
         xlim=c(0,30000),
         breaks=30)
}

par(mfrow=c(1,1))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=60,
     ylim=c(0,300),
     col=adjustcolor("gold",alpha.f = 0.3),
     main="number of SNPs obs per sample",
     xlab="number of SNPs")
for (i in 1:length(processed_SNVs_merged_subsampled_1978_SNPs_per_sample)){
    hist(processed_SNVs_merged_subsampled_1978_SNPs_per_sample[[i]],
         xlim=c(0,30000),
         breaks=30,
         col=adjustcolor(i,alpha.f = 0.1),add=T)
}
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     breaks=60,
     col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")

length(unique(H1_processed_SNVs_merged_truth_red$SAMPLE)) # 381




################################################################################
################################################################################
################################################################################
# number of observed SNPs in the total dataset? has the number of samples have an impact on that? must have, no?
lapply(processed_SNVs_merged_subsampled,function(iter){

    print(nrow(iter[!(duplicated(iter[c("CHROM","POS")])),]))
    #nrow(full_SNPs_mm10) # [1] 18150228 # from all F123 SNPs
    print(nrow(iter[!(duplicated(iter[c("CHROM","POS")])),])/18150228*100)
    # from subsampled F123 SNPs 2462745
    print(nrow(iter[!(duplicated(iter[c("CHROM","POS")])),])/2462745*100)
    
})
# representative:
# [1] 1592473 SNP observed in the dataset
# [1] 8.773846 % from all F123 SNPs
# [1] 64.66252 % from the subsampled F123 SNPs

lapply(processed_SNVs_merged_subsampled_2kb,function(iter){
    
    print(nrow(iter[!(duplicated(iter[c("CHROM","POS")])),]))
    #nrow(full_SNPs_mm10) # [1] 18150228
    print(nrow(iter[!(duplicated(iter[c("CHROM","POS")])),])/18150228*100)
    # from subsampled F123 SNPs 1231373
    print(nrow(iter[!(duplicated(iter[c("CHROM","POS")])),])/1231373*100)
    
})
# representative:
# [1] 796609
# [1] 4.388975
# [1] 64.69275

print(nrow(H1_processed_SNVs_merged_truth_red[!(duplicated(H1_processed_SNVs_merged_truth_red[c("CHROM","POS")])),]))
# 780308
print(nrow(H1_processed_SNVs_merged_truth_red[!(duplicated(H1_processed_SNVs_merged_truth_red[c("CHROM","POS")])),])/1453239*100)
# 53.6944

################################################################################
################################################################################
################################################################################
# why is there such a difference in phasing efficiency? Sasha shared hos result table
phasing_eff_table_subsampled_F123=fread("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/F123.subsampling.phasing.efficiencies.table")
dim(phasing_eff_table_subsampled_F123) #[1] 90  9
head(phasing_eff_table_subsampled_F123)

phasing_eff_table_subsampled_F123$res = as.numeric(substring(phasing_eff_table_subsampled_F123$Resolution, 5))
phasing_eff_table_subsampled_F123$res_col = factor(phasing_eff_table_subsampled_F123$res,levels=unique(phasing_eff_table_subsampled_F123$res),labels = c("black","red","gold"))
phasing_eff_table_subsampled_F123$dist_cutoff = as.numeric(substring(phasing_eff_table_subsampled_F123$Distance_cutoff, 13))
phasing_eff_table_subsampled_F123$dist_cutoff_col = factor(phasing_eff_table_subsampled_F123$dist_cutoff,levels=unique(phasing_eff_table_subsampled_F123$dist_cutoff),labels = c("lightblue1","hotpink4","yellow1"))
par(mfrow=c(2,3))
plot(x=jitter(phasing_eff_table_subsampled_F123$dist_cutoff),
     y=phasing_eff_table_subsampled_F123$Dual_efficiency,
     pch=19,
     col=adjustcolor(phasing_eff_table_subsampled_F123$res_col,alpha.f = 0.3))
plot(x=jitter(phasing_eff_table_subsampled_F123$res),
     y=phasing_eff_table_subsampled_F123$Dual_efficiency,
     pch=19,
     col=adjustcolor(phasing_eff_table_subsampled_F123$dist_cutoff_col,alpha.f = 0.3))
# higher dual with higher resolution (larger bins)
plot(x=jitter(phasing_eff_table_subsampled_F123$res),
     y=phasing_eff_table_subsampled_F123$CAST_efficiency,
     pch=19,
     col=adjustcolor(phasing_eff_table_subsampled_F123$dist_cutoff_col,alpha.f = 0.3))
points(x=jitter(phasing_eff_table_subsampled_F123$res),
     y=phasing_eff_table_subsampled_F123$S129_efficiency,
     pch=19,
     col=adjustcolor(phasing_eff_table_subsampled_F123$dist_cutoff_col,alpha.f = 0.3))
# no difference in the parents alleles,
# phasing efficiency goes down with larger bins
phasing_eff_table_subsampled_F123$total_eff = phasing_eff_table_subsampled_F123$S129_efficiency+phasing_eff_table_subsampled_F123$CAST_efficiency+phasing_eff_table_subsampled_F123$Dual_efficiency
phasing_eff_table_subsampled_F123$iter_col =  factor(substring(phasing_eff_table_subsampled_F123$Iteration, 10),levels=unique(substring(phasing_eff_table_subsampled_F123$Iteration, 10)),
                                                     labels = c("turquoise4","slategray1","violet","slategray2","violetred","slategray3","violetred1","slategray4","violetred3","violetred4"))
plot(x=jitter(phasing_eff_table_subsampled_F123$res),
       y=phasing_eff_table_subsampled_F123$total_eff,
       pch=19,
       col=adjustcolor(phasing_eff_table_subsampled_F123$iter_col,alpha.f = 0.8))
# iteratoins are less similar than distance cutoffs, 
# larger bins, less phased bins in total


# in general, the phasing efficiency is very high, in normal F123 traditional phasing we see between 75- 80% phasing efficiency.
# this might be a really good thing, as this was exactly what we wanted to achieve. 
# on the one hand, bins overlapping a SNP which maybe had not enough reads covering it to make the cut now includes many more reads in the surrounding and can be phased
# on the other hand, we might have removed contrary information due to subsampling:
# cophase the whole F123 dataset (1NP)?
################################################################################
################################################################################
################################################################################
# subsample again for the exact number of SNPs in our H1 truth
# > nrow(H1_truth_mapp)
# [1] 1453239

nr_subsampled_SNPs = 1453239 # 2kb:1231373; before with 1kb: 2462745
nr_subsampled_SNPs/nrow(full_SNPs_mm10)*100 # 8.006726 ###2kb:6.784337; before with 1kb: 13.56867 % 

# I will sample indizes of the rows, so I can easily subset the SNP df
subsampled_F123_list_H1 = vector(mode = "list", length = 10)
for(seed in 1:10){
    set.seed(seed)
    SNPs_mm10_subsampled_indizes = sample(x = nrow(full_SNPs_mm10), size = nr_subsampled_SNPs, replace = FALSE)
    SNPs_mm10_subsampled = full_SNPs_mm10[sort(SNPs_mm10_subsampled_indizes),]
    subsampled_F123_list_H1[[seed]] = SNPs_mm10_subsampled # 8.006726
    print(nrow(SNPs_mm10_subsampled)/nrow(full_SNPs_mm10)*100) #   0.129 0.871
    print(table(SNPs_mm10_subsampled$F123_GT)/nrow(SNPs_mm10_subsampled))
    SNPs_mm10_subsampled_chr = split(SNPs_mm10_subsampled, SNPs_mm10_subsampled$CHROM)
    mm10_SNP_subsampled_dist=unlist(lapply(SNPs_mm10_subsampled_chr,function(chr){
        dist=diff(chr$POS)
    }))
    print(mean(mm10_SNP_subsampled_dist)) # ca 1653
    print(median(mm10_SNP_subsampled_dist)) # 831
}

# [1] 8.006726
# 
# 0|1       1|0
# 0.1291859 0.8708141
# [1] 1653.02
# [1] 831
saveRDS(subsampled_F123_list_H1,'/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_H1_SNPamount_10times.rds')
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/F123_mm10_GT_subsampled_to_H1_SNPamount_10times.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/

# merge with processed SNPs 
#processed_SNVs = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/processed_F123_SNPs.rds")

processed_SNVs_merged_subsampled_H1 = lapply(subsampled_F123_list_H1,function(i){
    obs_truth_merged = merge(processed_SNVs,i,by=c("CHROM","POS"))
    obs_truth_merged = obs_truth_merged[order(as.numeric(obs_truth_merged$POS)),]
    obs_truth_merged = obs_truth_merged[order(match(obs_truth_merged$CHROM,chrom_order)),]
    return(obs_truth_merged)
})

processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist = lapply(processed_SNVs_merged_subsampled_H1,function(iter){
    per_chr = split(iter,iter$CHROM)
    per_chr_sample = unlist(lapply(per_chr,function(chr){
        per_sample = split(chr,chr$SAMPLE)
        sample_dist = unlist(lapply(per_sample,function(sample){
            return(diff(sample$POS))
        }))
        return(sample_dist)
    }))
    return(per_chr_sample)
})

saveRDS(processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_H1_SNPamount_obs_SNP_distance.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_H1_SNPamount_obs_SNP_distance.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_H1_SNPamount_obs_SNP_distance.rds")


# plot distance known SNPs (subsample to 1978bp distance)
subsampled_F123_list_H1=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/F123_mm10_GT_subsampled_to_H1_SNPamount_10times.rds")

chroms_F=paste("chr",1:19,sep="")
F123_sub_H1_chr = lapply(subsampled_F123_list_H1,function(iter){
    iter = iter[order(iter$POS),]
    iter = iter[order(as.numeric(substring(iter$CHROM, 4))),]
    iter$CHROM = factor(iter$CHROM, levels = chroms_F)
    iter_chr = split(iter,iter$CHROM)
    return(iter_chr)
})
dist_list_F123_H1 = lapply(F123_sub_H1_chr,function(iter){
    dist = unlist(lapply(iter,function(chr){
        return(diff(chr$POS))
    }))
    return(dist)
})

par(mfrow=c(4,3))
plot(H1_hist)
for (i in 1:length(dist_list_F123_H1)){
    hist(log10(dist_list_F123_H1[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i))
}

par(mfrow=c(1,1))
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between known SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,2e5))
abline(v=mean(log10(dist_list_H1)), col="gold") # 2.686611
for (i in 1:length(dist_list_F123_H1)){
    hist(log10(dist_list_F123_H1[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_H1))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(dist_list_F123_H1[[i]])), col=i)
    print(mean(log10(dist_list_F123_H1[[i]])))  # 2.83
}
plot(H1_hist,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")

mean(dist_list_H1) # 1915.676
mean(dist_list_F123_H1[[10]]) # 1653.207
mean(log10(dist_list_H1)) # 2.686611
mean(log10(dist_list_F123_H1[[10]])) # 2.83799
log10(mean(dist_list_H1)) # 3.282322
log10(mean(dist_list_F123_H1[[10]])) # 3.218327


# plot distance observed SNPs
par(mfrow=c(1,1))
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),
     main = "genomic distance between oserved SNPs\nin H1 and subsampled F123",
     xlab="log10(SNP distance)",
     ylim=c(0,3.5e5))
abline(v=mean(log10(dist_list_obs_H1_sample)), col="gold") # 4.05915
for (i in 1:length(processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist)){
    hist(log10(processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist[[i]]),breaks=50,
         xlim = c(0,max(log10(dist_list_obs_H1_sample))),
         main = paste("SNP dist - F123 iter ",i),
         col=adjustcolor(i,alpha.f = 0.1),
         add=T)
    abline(v=mean(log10(processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist[[i]])), col=i) # 4.298
}
plot(H1_hist_obs_sample,col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")
for (i in 1:length(processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist)){
    print(mean(log10(processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist[[i]]))) # 4.18
    print(mean(processed_SNVs_merged_subsampled_H1_chr_sample_pos_dist[[i]])) # between 512103.3 and 519929.4
}
mean(dist_list_obs_H1_sample) # 526705.2

# number of SNPs
# F123 
processed_SNVs_merged_subsampled_H1_SNPs_per_sample = lapply(processed_SNVs_merged_subsampled_H1,function(iter){
    per_sample = split(iter,iter$SAMPLE)
    sample_SNPs = unlist(lapply(per_sample,function(sample){
        return(nrow(sample))
    }))
    return(sample_SNPs)
})
saveRDS(processed_SNVs_merged_subsampled_H1_SNPs_per_sample,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_H1_obs_SNP_per_sample.rds")
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/prelim/subsampled_F123/subsampled_H1_obs_SNP_per_sample.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/
processed_SNVs_merged_subsampled_H1_SNPs_per_sample=readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/prelim/subsampled_H1_obs_SNP_per_sample.rds")

par(mfrow=c(4,3))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=30)
for (i in 1:length(processed_SNVs_merged_subsampled_H1_SNPs_per_sample)){
    hist(processed_SNVs_merged_subsampled_H1_SNPs_per_sample[[i]],
         xlim=c(0,30000),
         breaks=30)
}

par(mfrow=c(1,1))
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     xlim=c(0,30000),
     breaks=60,
     ylim=c(0,300),
     col=adjustcolor("gold",alpha.f = 0.3),
     main="number of SNPs obs per sample",
     xlab="number of SNPs")
for (i in 1:length(processed_SNVs_merged_subsampled_H1_SNPs_per_sample)){
    hist(processed_SNVs_merged_subsampled_H1_SNPs_per_sample[[i]],
         xlim=c(0,30000),
         breaks=30,
         col=adjustcolor(i,alpha.f = 0.1),add=T)
}
hist(H1_processed_SNVs_merged_truth_red_SNPs_per_sample,
     breaks=60,
     col=adjustcolor("gold",alpha.f = 0.3),add=T)
legend("topright",fill = c("gold",1:10),legend = c("H1",paste("F123 subsampled",1:10)),bty="n")

length(unique(H1_processed_SNVs_merged_truth_red$SAMPLE)) # 381

# does not really help either, as expectd


# for F123 we have 1123 sample 