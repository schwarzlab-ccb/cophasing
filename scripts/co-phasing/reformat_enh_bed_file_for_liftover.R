### reformat enhancer atlas file showing enhancer gene interactions for liftover


H1_enh_gene = read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction.bed",sep = "\t",header = F) # 227298      2
# this needs some parsing / cleaning / extracting:
# the data format in the file listed as "A:B-C_D$E$F$G$H I": 
# A - Chromosome of enhancer. 
# B - The starting position of enhancer.
# C - The ending position of enhancer. 
# D - Gene ensembl ID. 
# E - Gene ensembl Name. 
# F - chromosome of gene. 
# G - the position of the gene transcription start site. 
# H - The strand of DNA the gene located. 
# I - The predition score of the enhancer-gene interaction.
tmp = strsplit(as.character(H1_enh_gene$V1),"$",fixed=TRUE)
H1_enh_gene$gene_name = sapply(tmp, `[`, 2)
H1_enh_gene$gene_chr = sapply(tmp, `[`, 3)
H1_enh_gene$gene_TSS = sapply(tmp, `[`, 4)
H1_enh_gene$gene_strand = sapply(tmp, `[`, 5)
tmp_2 = sapply(tmp, `[`, 1)
tmp_3 = strsplit(as.character(tmp_2),"_",fixed=TRUE)
H1_enh_gene$gene_id = sapply(tmp_3, `[`, 2)
tmp_4 = sapply(tmp_3, `[`, 1)
tmp_5 = strsplit(as.character(tmp_4),"[:-]",fixed=FALSE)
H1_enh_gene$chr_enh = sapply(tmp_5, `[`, 1)
H1_enh_gene$start_enh = sapply(tmp_5, `[`, 2)
H1_enh_gene$end_enh = sapply(tmp_5, `[`, 3)
H1_enh_gene$merge_help = c(1:nrow(H1_enh_gene))

H1_enh_gene_enhcoord = H1_enh_gene[,c("chr_enh","start_enh","end_enh","gene_id","gene_name","gene_chr","gene_TSS","gene_strand","merge_help")] # 227298      9

H1_enh_gene$gene_end_fake = as.numeric(H1_enh_gene$gene_TSS)+1
H1_enh_gene_genecoord = H1_enh_gene[,c("gene_chr","gene_TSS","gene_end_fake","gene_strand","gene_id","gene_name","chr_enh","start_enh","end_enh","merge_help")] # 227298      10

# save as rds and write out as bed
saveRDS(H1_enh_gene_enhcoord,"/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_enhancercoordinates.rds")
saveRDS(H1_enh_gene_genecoord,"/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_genecoordinates.rds")

write.table(H1_enh_gene_enhcoord,"/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_enhancercoordinates.bed",sep = "\t",quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(H1_enh_gene_genecoord,"/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_genecoordinates.bed",sep = "\t",quote = FALSE, col.names = FALSE, row.names = FALSE)


################################################################################
################################################################################
############################## liftover ########################################
################################################################################
################################################################################
# merge the lifted positions again

H1_enh_gene_enhcoord_lift = read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_enhancercoordinates_liftover_hg19_to_hg38.bed",sep = "\t",header = F)
H1_enh_gene_genecoord_lift = read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_genecoordinates_liftover_hg19_to_hg38.bed",sep = "\t",header = F)

colnames(H1_enh_gene_enhcoord_lift) = c("chr_enh_enh","start_enh_enh","end_enh_enh","gene_id_enh","gene_name_enh","gene_chr_enh","gene_TSS_enh","gene_strand_enh","merge_help")
colnames(H1_enh_gene_genecoord_lift) = c("gene_chr_gene","gene_TSS_gene","gene_end_fake_gene","gene_strand_gene","gene_id_gene","gene_name_gene","chr_enh_gene","start_enh_gene","end_enh_gene","merge_help")

H1_enh_gene_lift = merge(H1_enh_gene_enhcoord_lift,
                         H1_enh_gene_genecoord_lift, 
                         by = "merge_help")

# dim(H1_enh_gene_enhcoord_lift) # 226905      9
# dim(H1_enh_gene_genecoord_lift) # 227216     10
# dim(H1_enh_gene_lift) # 226829     18

H1_enh_gene_lift = H1_enh_gene_lift[,c("chr_enh_enh", "start_enh_enh", "end_enh_enh", "gene_id_gene", "gene_name_gene", "gene_chr_gene", "gene_TSS_gene", "gene_strand_gene")]
colnames(H1_enh_gene_lift) = c("chr_enh", "start_enh", "end_enh", "gene_id", "gene_name", "gene_chr", "gene_TSS", "gene_strand")
saveRDS(H1_enh_gene_lift,"/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_liftover_hg38.rds")
# scp bihtext:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_liftover_hg38.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/

# merge with anno to compare the liftover gene position with our anno 
# load anno
annoV29_df = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.annotation.gtf.rds")
annoV29_df_g = annoV29_df[which(annoV29_df$type=="gene"),] # 58721
annoV29_df_g$gene_id_red = sapply(strsplit(annoV29_df_g$gene_id,".",fixed=TRUE), `[`, 1)

H1_enh_gene_lift_merge_annoV29_df_g = merge(H1_enh_gene_lift,
                                            annoV29_df_g[,c("seqnames", "start", "end", "width", "strand", "gene_id", "gene_type", "gene_name", "gene_id_red")],
                                            by.x = "gene_id",by.y="gene_id_red")

# dim(H1_enh_gene_lift) # 226829      8
# dim(annoV29_df_g) # 58721    26
# dim(H1_enh_gene_lift_merge_annoV29_df_g) # 217765     16
all(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh == H1_enh_gene_lift_merge_annoV29_df_g$seqnames) # FALSE
table(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh == H1_enh_gene_lift_merge_annoV29_df_g$seqnames)
# FALSE   TRUE
#    66 217699

#tmp = H1_enh_gene_lift_merge_annoV29_df_g[which(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh != H1_enh_gene_lift_merge_annoV29_df_g$seqnames),]
# most of those were moved from X to Y, so not interesting, or included from scaffolds into the chr ( chr22_KI270879v1_alt to chr22 )
# others are misc_RNAs with name and chr changes
#tmp[which(tmp$gene_name.x != tmp$gene_name.y),]

H1_enh_gene_lift_merge_annoV29_df_g$chr_enh_red = sapply(strsplit(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh,"_",fixed=TRUE), `[`, 1)
table(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh_red == H1_enh_gene_lift_merge_annoV29_df_g$seqnames)
# FALSE   TRUE
#    51 217714
H1_enh_gene_lift_merge_annoV29_df_g[which(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh_red != H1_enh_gene_lift_merge_annoV29_df_g$seqnames),]
# exclusively X to Y conversions and misc_RNA, can be filtered out

H1_enh_gene_lift_merge_annoV29_df_g = H1_enh_gene_lift_merge_annoV29_df_g[which(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh_red == H1_enh_gene_lift_merge_annoV29_df_g$seqnames),]
H1_enh_gene_lift_merge_annoV29_df_g = H1_enh_gene_lift_merge_annoV29_df_g[which(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh %in% paste("chr",1:22,sep="")),] # 215099     
H1_enh_gene_lift_merge_annoV29_df_g = H1_enh_gene_lift_merge_annoV29_df_g[which(H1_enh_gene_lift_merge_annoV29_df_g$gene_chr %in% paste("chr",1:22,sep="")),] # 215093

H1_enh_gene_lift_merge_annoV29_df_g$pos_diff = abs(H1_enh_gene_lift_merge_annoV29_df_g$gene_TSS-H1_enh_gene_lift_merge_annoV29_df_g$start)

min(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff) # 0
mean(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff) # 12080.03
median(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff) # 310
max(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff) # 124498771
head(sort(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff),n=20) # all 0
table(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff == 0)/nrow(H1_enh_gene_lift_merge_annoV29_df_g)*100
#     FALSE     TRUE
# 62.17869 37.82131
# TSS might not be by definition of start of the gene
table(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff < 1000)/nrow(H1_enh_gene_lift_merge_annoV29_df_g)*100
#    FALSE     TRUE
# 44.79039 55.20961
table(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff < 50000)/nrow(H1_enh_gene_lift_merge_annoV29_df_g)*100
#    FALSE      TRUE
# 5.129409 94.870591
# so almost 95% are within the 50kb bin that is the resolution of our maps
head(sort(H1_enh_gene_lift_merge_annoV29_df_g$pos_diff,decreasing = T),n=20)
# 124498771  30636129  22806045  22775513  18105424  18105424   3676844   3676844   3676844   3242636   3242636   3241485   3241485   3241485   2573154   2172910   2059619   1783651   1783651   1783651

table(H1_enh_gene_lift_merge_annoV29_df_g$gene_type == "protein_coding")/nrow(H1_enh_gene_lift_merge_annoV29_df_g)*100
#    FALSE     TRUE
# 25.26024 74.73976
H1_enh_gene_lift_merge_annoV29_df_g_protcod = H1_enh_gene_lift_merge_annoV29_df_g[which(H1_enh_gene_lift_merge_annoV29_df_g$gene_type == "protein_coding"),] #160760
min(H1_enh_gene_lift_merge_annoV29_df_g_protcod$pos_diff) # 0
mean(H1_enh_gene_lift_merge_annoV29_df_g_protcod$pos_diff) # 13504.38
median(H1_enh_gene_lift_merge_annoV29_df_g_protcod$pos_diff) # 1369
max(H1_enh_gene_lift_merge_annoV29_df_g_protcod$pos_diff) # 22806045 #
head(sort(H1_enh_gene_lift_merge_annoV29_df_g_protcod$pos_diff,decreasing = T),n=20)
# 22806045  2172910  2059619  1783651  1783651  1783651  1502149  1380351  1380351  1380351  1380351  1380351  1258197  1196948  1170575  1170575  1117542  1117218  1048732  1048732
# that is 456 - 21 bins

H1_enh_gene_lift_merge_annoV29_df_g$enh_mid = (H1_enh_gene_lift_merge_annoV29_df_g$start_enh + H1_enh_gene_lift_merge_annoV29_df_g$end_enh)/2
H1_enh_gene_lift_merge_annoV29_df_g$gene_mid = (H1_enh_gene_lift_merge_annoV29_df_g$start + H1_enh_gene_lift_merge_annoV29_df_g$end)/2

# sort
H1_enh_gene_lift_merge_annoV29_df_g = H1_enh_gene_lift_merge_annoV29_df_g[order(H1_enh_gene_lift_merge_annoV29_df_g$start_enh),]
H1_enh_gene_lift_merge_annoV29_df_g = H1_enh_gene_lift_merge_annoV29_df_g[order(as.numeric(substring(H1_enh_gene_lift_merge_annoV29_df_g$chr_enh, 4))),]

saveRDS(H1_enh_gene_lift_merge_annoV29_df_g,"/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_liftover_hg38_annoV29.rds")
#scp bihtext:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_liftover_hg38_annoV29.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/