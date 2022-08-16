# I recounted the mapped reads on gene level, now the results need to be merged with the anno info again
# also the ASE counts per SNPs can be reused, but need to be reaggregated on gene level

################################################################################
### load data

# load anno
annoV29_df = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.annotation.gtf.rds")

# load SNPs
H1_SNPs = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/H1_haplotypes_dixon/H1_dixon_haplotypes_liftover_hg38_rs_id_refalt_GT_chrompossort_only_phased.rds")

# load total counts
ENCFF675NTU_total = read.table("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/total_ENCFF675NTU/htseq_counts_downloaded_ENCFF675NTU_genelevel_introns_V29.table",nrows=58721)
ENCFF379NOY_total = read.table("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/total_ENCFF379NOY/htseq_counts_downloaded_ENCFF379NOY_genelevel_introns_V29.table",nrows=58721) 

# load ASE counts
ENCFF675NTU_ase = read.table("/fast/groups/ag_schwarz/Projects/project-gam/RNAseq/src/project-ponder-ase-new/postprocessing_ASEReadcounter_ENCFF675NTU_mindepth1/ase_gene_table/allele_counts_downloaded_ENCFF675NTU.sort.ase_genes.tsv",skip=1)
ENCFF379NOY_ase = read.table("/fast/groups/ag_schwarz/Projects/project-gam/RNAseq/src/project-ponder-ase-new/postprocessing_ASEReadcounter_ENCFF379NOY_mindepth1/ase_gene_table/allele_counts_downloaded_ENCFF379NOY.sort.ase_genes.tsv",skip=1)

################################################################################
### merge gene level total expression counts with annotation
ENCFF675NTU_total = ENCFF675NTU_total[,2:5]
ENCFF379NOY_total = ENCFF379NOY_total[,2:5]

colnames(ENCFF675NTU_total) = c("gene_id","gene_name","gene_type","total_counts_ENCFF675NTU")
colnames(ENCFF379NOY_total) = c("gene_id","gene_name","gene_type","total_counts_ENCFF379NOY")

length(unique(ENCFF675NTU_total$gene_id)) #[1] 58721

annoV29_df_g = annoV29_df[which(annoV29_df$type=="gene"),]
dim(annoV29_df_g) #[1] 58721    25

ENCFF675NTU_total_merge_anno = merge(ENCFF675NTU_total,annoV29_df_g,by="gene_id") # 58721    28
total_RNA = merge(ENCFF675NTU_total_merge_anno,ENCFF379NOY_total[,c(1,4)],by="gene_id") # 58721    29

total_RNA_red = total_RNA[,c("seqnames","start","end","width","strand","gene_id","gene_name.x","gene_type.x","level","total_counts_ENCFF675NTU","total_counts_ENCFF379NOY")]
colnames(total_RNA_red) = c("chr","start","end","width","strand","gene_id","gene_name","gene_type","level","total_counts_ENCFF675NTU","total_counts_ENCFF379NOY")

total_RNA_red = total_RNA_red[order(total_RNA_red$start),]
total_RNA_red = total_RNA_red[order(as.numeric(substring(total_RNA_red$chr, 4))),]

### calc TPM
tpm = function(counts, lengths) {
    rpk=counts/(lengths/1000)
    sum_rpk=sum(rpk)
    scaling=sum_rpk/1000000
    tpm=rpk/scaling
}
total_RNA_red$TPM_ENCFF675NTU = tpm(total_RNA_red$total_counts_ENCFF675NTU,total_RNA_red$width)
total_RNA_red$TPM_ENCFF379NOY = tpm(total_RNA_red$total_counts_ENCFF379NOY,total_RNA_red$width)

table(total_RNA_red$TPM_ENCFF675NTU >= 1,total_RNA_red$TPM_ENCFF379NOY >= 1)/nrow(total_RNA_red)*100
#            FALSE       TRUE
# FALSE 61.356244  1.783008
# TRUE   3.230531 33.630217

total_RNA_red_auto = total_RNA_red[which(total_RNA_red$chr%in%paste("chr",1:22,sep="")),]
nrow(total_RNA_red_auto) # 55746

total_RNA_red_auto_protcod = total_RNA_red_auto[which(total_RNA_red_auto$gene_type == "protein_coding"),]
nrow(total_RNA_red_auto_protcod) # 19019

table(total_RNA_red_auto_protcod$TPM_transcript_ENCFF675NTU >= 1,total_RNA_red_auto_protcod$TPM_transcript_ENCFF379NOY >= 1)/nrow(total_RNA_red_auto_protcod)*100
#           FALSE      TRUE
# FALSE 51.392971  0.317827
# TRUE   1.717633 46.571569

saveRDS(total_RNA_red_auto, "/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/H1_total_expression_ENCFF675NTU_ENCFF379NOY_gene_level_autosomes.rds")

cor(total_RNA_red_auto_protcod$TPM_ENCFF675NTU,total_RNA_red_auto_protcod$TPM_ENCFF379NOY) # 0.9994274

################################################################################
### aggregate ASE counts per gene level
### ENCFF675NTU 
################################################################################

colnames(ENCFF675NTU_ase) = c("chr","pos","gene_id","feature","strand","gene_name","start","end","variantID","refAllele","altAllele","refCount","altCount","totalCount","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs")

### add GT info to allele-specific readcounts 

ENCFF675NTU_ase_gt_merged=merge(x = ENCFF675NTU_ase,
                                y = H1_SNPs[,c("ID","H1")],
                                by.x = c("variantID"),
                                by.y = c("ID"), 
                                all.y = F)

# > dim(ENCFF379NOY_ase)
# [1] 172929     19
# > dim(H1_SNPs)
# [1] 1563859      10
# > dim(ENCFF675NTU_ase_gt_merged)
# [1] 196359     20


# sort
chrom_order=paste("chr",c(1:22),sep="")
ENCFF675NTU_ase_gt_merged = ENCFF675NTU_ase_gt_merged[order(ENCFF675NTU_ase_gt_merged$pos),]
ENCFF675NTU_ase_gt_merged = ENCFF675NTU_ase_gt_merged[order(match(ENCFF675NTU_ase_gt_merged$chr,chrom_order)),]
ENCFF675NTU_ase_gt_merged$chr = factor(ENCFF675NTU_ase_gt_merged$chr, levels = unique(ENCFF675NTU_ase_gt_merged$chr))

# phase the allele-specific read counts
ENCFF675NTU_ase_gt_merged$P1_count=0
ENCFF675NTU_ase_gt_merged$P2_count=0
ENCFF675NTU_ase_gt_merged$P1_count[which(ENCFF675NTU_ase_gt_merged$H1=='1|0')]=ENCFF675NTU_ase_gt_merged$altCount[which(ENCFF675NTU_ase_gt_merged$H1=='1|0')]
ENCFF675NTU_ase_gt_merged$P1_count[which(ENCFF675NTU_ase_gt_merged$H1=='0|1')]=ENCFF675NTU_ase_gt_merged$refCount[which(ENCFF675NTU_ase_gt_merged$H1=='0|1')]
ENCFF675NTU_ase_gt_merged$P2_count[which(ENCFF675NTU_ase_gt_merged$H1=='1|0')]=ENCFF675NTU_ase_gt_merged$refCount[which(ENCFF675NTU_ase_gt_merged$H1=='1|0')]
ENCFF675NTU_ase_gt_merged$P2_count[which(ENCFF675NTU_ase_gt_merged$H1=='0|1')]=ENCFF675NTU_ase_gt_merged$altCount[which(ENCFF675NTU_ase_gt_merged$H1=='0|1')]


### sort SNPs into transcript ranges ############################################

# split counts and anno by chromosome
ENCFF675NTU_ase_gt_merged_chr = split(ENCFF675NTU_ase_gt_merged,ENCFF675NTU_ase_gt_merged$chr) # no minimum read depth per SNP position
annoV29_df_g_chr = split(annoV29_df_g[,c("seqnames","start","end","width","gene_id")],annoV29_df_g$seqnames)

annoV29_df_g_chr = annoV29_df_g_chr[1:22]
# sort the SNPs into the transcript ranges
idx_chr=c(1:length(ENCFF675NTU_ase_gt_merged_chr)) # both 22, both are in the same order
ENCFF675NTU_ase_gt_anno_merged_chr=vector('list',length(ENCFF675NTU_ase_gt_merged_chr))

library('sqldf')
for (chr in idx_chr){
    print(chr)
    print(Sys.time())
    ase_chr=ENCFF675NTU_ase_gt_merged_chr[[chr]]
    anno_chr=annoV29_df_g_chr[[chr]]
    ENCFF675NTU_ase_gt_anno_merged_chr[[chr]]=sqldf("select * from ase_chr f1 left join anno_chr f2 on (f1.pos >= f2.start and f1.pos<= f2.end) ")
    #ASE_gt_anno_merged_chr[[chr]]=ASE_gt_anno_merged_chr[[chr]][which(!(is.na(ASE_gt_anno_merged_chr[[chr]]$feature))),]
    print('chr done')
    print(Sys.time())
}
ENCFF675NTU_ase_gt_anno_merged = do.call(rbind.data.frame, ENCFF675NTU_ase_gt_anno_merged_chr)

colnames(ENCFF675NTU_ase_gt_anno_merged) =  c("variantID","chr","SNP_pos","gene_id","feature","strand","gene_name","start_gene","end_gene","refAllele","altAllele","refCount","altCount","totalCount",
                                              "lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","H1","P1_count","P2_count","seqnames","start_g2","end_g2","width","gene_id_filter")
# filter out transcripts matched to the wrong gene id
dim(ENCFF675NTU_ase_gt_anno_merged) # [1] 255889     27
ENCFF675NTU_ase_gt_anno_merged = ENCFF675NTU_ase_gt_anno_merged[which(ENCFF675NTU_ase_gt_anno_merged$gene_id==ENCFF675NTU_ase_gt_anno_merged$gene_id_filter),]
dim(ENCFF675NTU_ase_gt_anno_merged) # [1] 196359     27

# all(ENCFF675NTU_ase_gt_anno_merged$start_gene == ENCFF675NTU_ase_gt_anno_merged$start_g2) # TRUE
# all(ENCFF675NTU_ase_gt_anno_merged$end_gene == ENCFF675NTU_ase_gt_anno_merged$end_g2) # TRUE
ENCFF675NTU_ase_gt_anno_merged = ENCFF675NTU_ase_gt_anno_merged[,c("variantID","chr","SNP_pos","gene_id","feature","strand","gene_name","start_gene","end_gene","width","refAllele","altAllele","refCount","altCount","totalCount","H1","P1_count","P2_count")]
annoV29_df_g_red = annoV29_df_g[,c("seqnames","start","end","width","strand","gene_id","gene_type","gene_name","level")]
ENCFF675NTU_ase_gt_anno_merged = merge(ENCFF675NTU_ase_gt_anno_merged,annoV29_df_g_red,by="gene_id",all.y=F)

ENCFF675NTU_ase_gt_anno_merged = ENCFF675NTU_ase_gt_anno_merged[,c("variantID","chr","SNP_pos","refAllele","altAllele","refCount","altCount","totalCount","H1","P1_count","P2_count",
                                                                   "start_gene","end_gene","width.x","gene_id","gene_name.x","gene_type","feature","strand.x")]
ENCFF675NTU_ase_gt_anno_merged = ENCFF675NTU_ase_gt_anno_merged[order(ENCFF675NTU_ase_gt_anno_merged$SNP_pos),]
ENCFF675NTU_ase_gt_anno_merged = ENCFF675NTU_ase_gt_anno_merged[order(match(ENCFF675NTU_ase_gt_anno_merged$chr,chrom_order)),]
colnames(ENCFF675NTU_ase_gt_anno_merged) = c("variantID","chr","SNP_pos","refAllele","altAllele","refCount","altCount","totalCount","H1","P1_count","P2_count","start_gene","end_gene","width","gene_id","gene_name","gene_type","feature","strand")

dim(ENCFF675NTU_ase_gt_anno_merged) #[1] 196359     19

saveRDS(ENCFF675NTU_ase_gt_anno_merged,"/fast/groups/ag_schwarz/Projects/project-gam/RNAseq/src/project-ponder-ase-new/postprocessing_ASEReadcounter_ENCFF675NTU_mindepth1/ASE_SNPs_allele_counts_merged_anno_V29_downloaded_ENCFF675NTU_gene_level.rds")

### aggregate SNPs over transcripts ############################################

### split by transcript id
ENCFF675NTU_ase_gt_anno_merged_splitg=split(ENCFF675NTU_ase_gt_anno_merged,ENCFF675NTU_ase_gt_anno_merged$gene_id,drop=TRUE)
length(ENCFF675NTU_ase_gt_anno_merged_splitg) # 15847

# bis hier gekommen

### aggregate
agg_g = vector('list',length(ENCFF675NTU_ase_gt_anno_merged_splitg))
for (idx in 1:length(ENCFF675NTU_ase_gt_anno_merged_splitg)){
    gene = ENCFF675NTU_ase_gt_anno_merged_splitg[[idx]]
    #print(idx)
    if (idx%%1000==0){
        percent = idx/length(ENCFF675NTU_ase_gt_anno_merged_splitg)*100
        print(percent)
    }
    
    number_snps = nrow(gene)
    
    if(number_snps==0){return(NULL)}
    agg = aggregate(cbind(gene$P1_count,gene$P2_count), 
                    by=list(gene$chr,
                            gene$start_gene,
                            gene$end_gene,
                            gene$width,
                            gene$strand,
                            gene$gene_id,
                            gene$gene_name,
                            gene$gene_type),FUN=sum)
    
    colnames(agg) = c('chrom','start_gene','end_gene','width_gene','strand','gene_id','gene_name','gene_type','P1_count_sum','P2_count_sum')
    agg$p_value_sum = binom.test(cbind(round(agg$P1_count_sum),round(agg$P2_count_sum)),p=0.5,alternative = c("two.sided"))$p.value
    agg$number_SNPs = number_snps
    agg$number_P1_alt=table(gene$H1)["1|0"]
    agg$number_P2_alt=table(gene$H1)["0|1"]
    ex_in_tab=table(gene$feature)
    if("intron" %in% names(ex_in_tab)){
        agg$number_intron_SNPs=ex_in_tab["intron"]
    }else{agg$number_intron_SNPs=0}
    if("exon" %in% names(ex_in_tab)){
        agg$number_exon_SNPs=ex_in_tab["exon"]
    }else{agg$number_exon_SNPs=0}
    agg_g[[idx]]=agg
    #return(agg)
}

agg_g_df=do.call(rbind.data.frame, agg_g)
dim(agg_g_df) # 15847    16

# calc additional stats
agg_g_df$p_adj_sum=p.adjust(agg_g_df$p_value_sum, method = "BH")
agg_g_df$sig_sum=ifelse(agg_g_df$p_adj_sum<=0.05,'TRUE','FALSE')
agg_g_df$ASE_ratio=agg_g_df$P1_count_sum/(agg_g_df$P2_count_sum+agg_g_df$P1_count_sum)
agg_g_df$log2foldchange=log2((agg_g_df$P1_count_sum+1)/(agg_g_df$P2_count_sum+1))
agg_g_df$number_P1_alt[which(is.na(agg_g_df$number_P1_alt))]=0
agg_g_df$number_P2_alt[which(is.na(agg_g_df$number_P2_alt))]=0
agg_g_df$SNP_ratio=agg_g_df$number_P1_alt/(agg_g_df$number_P1_alt+agg_g_df$number_P2_alt)

# sort
agg_g_df = agg_g_df[order(agg_g_df$start_gene),]
agg_g_df = agg_g_df[order(match(agg_g_df$chrom,chrom_order)),]
agg_g_df$chrom = factor(agg_g_df$chrom, levels = unique(agg_g_df$chrom))
saveRDS(agg_g_df,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/ASEReadcounter_ENCFF675NTU/ASE_incl_intro_anno_V29_downloaded_ENCFF675NTU_gene_level.rds")
agg_g_df_ENCFF675NTU = agg_g_df

# first glance
table(agg_g_df$sig_sum)
# FALSE  TRUE
# 15397   450
table(agg_g_df$sig_sum)/nrow(agg_g_df)*100
#     FALSE      TRUE
# 97.160346  2.839654

table(agg_g_df$log2foldchange>=1)/nrow(agg_g_df)*100
#    FALSE     TRUE
# 81.95242 18.04758

table(agg_g_df$log2foldchange>=1,agg_g_df$sig_sum)/nrow(agg_g_df)*100 # that is not much at all!
#            FALSE       TRUE
#  FALSE 80.153972  1.798448
#  TRUE  17.006373  1.041207

table(agg_g_df$log2foldchange>=1,agg_g_df$sig_sum)
#       FALSE  TRUE
# FALSE 12702   285
# TRUE   2695   165

agg_g_df_LFC_ASE = agg_g_df[which(agg_g_df$log2foldchange >= 1 & agg_g_df$sig_sum==TRUE),]
length(unique(agg_g_df_LFC_ASE$gene_id)) # 165
length(unique(agg_g_df_LFC_ASE$gene_name)) # 165
# that does still include genes that are not protein coding

length(unique(agg_g_df_LFC_ASE$gene_id[which(agg_g_df_LFC_ASE$gene_type=="protein_coding")])) # 124
################################################################################
################################################################################
### aggregate ASE counts per gene level
### ENCFF379NOY 
################################################################################

colnames(ENCFF379NOY_ase) = c("chr","pos","gene_id","feature","strand","gene_name","start","end","variantID","refAllele","altAllele","refCount","altCount","totalCount","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs")

### add GT info to allele-specific readcounts 

ENCFF379NOY_ase_gt_merged=merge(x = ENCFF379NOY_ase,
                                y = H1_SNPs[,c("ID","H1")],
                                by.x = c("variantID"),
                                by.y = c("ID"), 
                                all.y = F)

# > dim(ENCFF379NOY_ase)
# [1] 172929     19
# > dim(H1_SNPs)
# [1] 1563859      10
# > dim(ENCFF379NOY_ase_gt_merged)
# [1] 172929     20


# sort
chrom_order=paste("chr",c(1:22),sep="")
ENCFF379NOY_ase_gt_merged = ENCFF379NOY_ase_gt_merged[order(ENCFF379NOY_ase_gt_merged$pos),]
ENCFF379NOY_ase_gt_merged = ENCFF379NOY_ase_gt_merged[order(match(ENCFF379NOY_ase_gt_merged$chr,chrom_order)),]
ENCFF379NOY_ase_gt_merged$chr = factor(ENCFF379NOY_ase_gt_merged$chr, levels = unique(ENCFF379NOY_ase_gt_merged$chr))

# phase the allele-specific read counts
ENCFF379NOY_ase_gt_merged$P1_count=0
ENCFF379NOY_ase_gt_merged$P2_count=0
ENCFF379NOY_ase_gt_merged$P1_count[which(ENCFF379NOY_ase_gt_merged$H1=='1|0')]=ENCFF379NOY_ase_gt_merged$altCount[which(ENCFF379NOY_ase_gt_merged$H1=='1|0')]
ENCFF379NOY_ase_gt_merged$P1_count[which(ENCFF379NOY_ase_gt_merged$H1=='0|1')]=ENCFF379NOY_ase_gt_merged$refCount[which(ENCFF379NOY_ase_gt_merged$H1=='0|1')]
ENCFF379NOY_ase_gt_merged$P2_count[which(ENCFF379NOY_ase_gt_merged$H1=='1|0')]=ENCFF379NOY_ase_gt_merged$refCount[which(ENCFF379NOY_ase_gt_merged$H1=='1|0')]
ENCFF379NOY_ase_gt_merged$P2_count[which(ENCFF379NOY_ase_gt_merged$H1=='0|1')]=ENCFF379NOY_ase_gt_merged$altCount[which(ENCFF379NOY_ase_gt_merged$H1=='0|1')]


### sort SNPs into transcript ranges ############################################

# split counts and anno by chromosome
ENCFF379NOY_ase_gt_merged_chr = split(ENCFF379NOY_ase_gt_merged,ENCFF379NOY_ase_gt_merged$chr) # no minimum read depth per SNP position
annoV29_df_g_chr = split(annoV29_df_g[,c("seqnames","start","end","width","gene_id")],annoV29_df_g$seqnames)

annoV29_df_g_chr = annoV29_df_g_chr[1:22]
# sort the SNPs into the transcript ranges
idx_chr=c(1:length(ENCFF379NOY_ase_gt_merged_chr)) # both 22, both are in the same order
ENCFF379NOY_ase_gt_anno_merged_chr=vector('list',length(ENCFF379NOY_ase_gt_merged_chr))

library('sqldf')
for (chr in idx_chr){
    print(chr)
    print(Sys.time())
    ase_chr=ENCFF379NOY_ase_gt_merged_chr[[chr]]
    anno_chr=annoV29_df_g_chr[[chr]]
    ENCFF379NOY_ase_gt_anno_merged_chr[[chr]]=sqldf("select * from ase_chr f1 left join anno_chr f2 on (f1.pos >= f2.start and f1.pos<= f2.end) ")
    #ASE_gt_anno_merged_chr[[chr]]=ASE_gt_anno_merged_chr[[chr]][which(!(is.na(ASE_gt_anno_merged_chr[[chr]]$feature))),]
    print('chr done')
    print(Sys.time())
}
ENCFF379NOY_ase_gt_anno_merged = do.call(rbind.data.frame, ENCFF379NOY_ase_gt_anno_merged_chr)

colnames(ENCFF379NOY_ase_gt_anno_merged) =  c("variantID","chr","SNP_pos","gene_id","feature","strand","gene_name","start_gene","end_gene","refAllele","altAllele","refCount","altCount","totalCount",
                                              "lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs","H1","P1_count","P2_count","seqnames","start_g2","end_g2","width","gene_id_filter")
# filter out transcripts matched to the wrong gene id
dim(ENCFF379NOY_ase_gt_anno_merged) # [1] 226611     27
ENCFF379NOY_ase_gt_anno_merged = ENCFF379NOY_ase_gt_anno_merged[which(ENCFF379NOY_ase_gt_anno_merged$gene_id==ENCFF379NOY_ase_gt_anno_merged$gene_id_filter),]
dim(ENCFF379NOY_ase_gt_anno_merged) # [1] 172929     27

# all(ENCFF379NOY_ase_gt_anno_merged$start_gene == ENCFF379NOY_ase_gt_anno_merged$start_g2) # TRUE
# all(ENCFF379NOY_ase_gt_anno_merged$end_gene == ENCFF379NOY_ase_gt_anno_merged$end_g2) # TRUE
ENCFF379NOY_ase_gt_anno_merged = ENCFF379NOY_ase_gt_anno_merged[,c("variantID","chr","SNP_pos","gene_id","feature","strand","gene_name","start_gene","end_gene","width","refAllele","altAllele","refCount","altCount","totalCount","H1","P1_count","P2_count")]
annoV29_df_g_red = annoV29_df_g[,c("seqnames","start","end","width","strand","gene_id","gene_type","gene_name","level")]
ENCFF379NOY_ase_gt_anno_merged = merge(ENCFF379NOY_ase_gt_anno_merged,annoV29_df_g_red,by="gene_id",all.y=F)

ENCFF379NOY_ase_gt_anno_merged = ENCFF379NOY_ase_gt_anno_merged[,c("variantID","chr","SNP_pos","refAllele","altAllele","refCount","altCount","totalCount","H1","P1_count","P2_count",
                                                                   "start_gene","end_gene","width.x","gene_id","gene_name.x","gene_type","feature","strand.x")]
ENCFF379NOY_ase_gt_anno_merged = ENCFF379NOY_ase_gt_anno_merged[order(ENCFF379NOY_ase_gt_anno_merged$SNP_pos),]
ENCFF379NOY_ase_gt_anno_merged = ENCFF379NOY_ase_gt_anno_merged[order(match(ENCFF379NOY_ase_gt_anno_merged$chr,chrom_order)),]
colnames(ENCFF379NOY_ase_gt_anno_merged) = c("variantID","chr","SNP_pos","refAllele","altAllele","refCount","altCount","totalCount","H1","P1_count","P2_count","start_gene","end_gene","width","gene_id","gene_name","gene_type","feature","strand")

dim(ENCFF379NOY_ase_gt_anno_merged) #[1] 172929     19

saveRDS(ENCFF379NOY_ase_gt_anno_merged,"/fast/groups/ag_schwarz/Projects/project-gam/RNAseq/src/project-ponder-ase-new/postprocessing_ASEReadcounter_ENCFF379NOY_mindepth1/ASE_SNPs_allele_counts_merged_anno_V29_downloaded_ENCFF379NOY_gene_level.rds")

### aggregate SNPs over transcripts ############################################

### split by transcript id
ENCFF379NOY_ase_gt_anno_merged_splitg=split(ENCFF379NOY_ase_gt_anno_merged,ENCFF379NOY_ase_gt_anno_merged$gene_id,drop=TRUE)
length(ENCFF379NOY_ase_gt_anno_merged_splitg) # 15326

# bis hier gekommen

### aggregate
agg_g = vector('list',length(ENCFF379NOY_ase_gt_anno_merged_splitg))
for (idx in 1:length(ENCFF379NOY_ase_gt_anno_merged_splitg)){
    gene = ENCFF379NOY_ase_gt_anno_merged_splitg[[idx]]
    #print(idx)
    if (idx%%1000==0){
        percent = idx/length(ENCFF379NOY_ase_gt_anno_merged_splitg)*100
        print(percent)
    }
    
    number_snps = nrow(gene)
    
    if(number_snps==0){return(NULL)}
    agg = aggregate(cbind(gene$P1_count,gene$P2_count), 
                    by=list(gene$chr,
                            gene$start_gene,
                            gene$end_gene,
                            gene$width,
                            gene$strand,
                            gene$gene_id,
                            gene$gene_name,
                            gene$gene_type),FUN=sum)
    
    colnames(agg) = c('chrom','start_gene','end_gene','width_gene','strand','gene_id','gene_name','gene_type','P1_count_sum','P2_count_sum')
    agg$p_value_sum = binom.test(cbind(round(agg$P1_count_sum),round(agg$P2_count_sum)),p=0.5,alternative = c("two.sided"))$p.value
    agg$number_SNPs = number_snps
    agg$number_P1_alt=table(gene$H1)["1|0"]
    agg$number_P2_alt=table(gene$H1)["0|1"]
    ex_in_tab=table(gene$feature)
    if("intron" %in% names(ex_in_tab)){
        agg$number_intron_SNPs=ex_in_tab["intron"]
    }else{agg$number_intron_SNPs=0}
    if("exon" %in% names(ex_in_tab)){
        agg$number_exon_SNPs=ex_in_tab["exon"]
    }else{agg$number_exon_SNPs=0}
    agg_g[[idx]]=agg
    #return(agg)
}

agg_g_df=do.call(rbind.data.frame, agg_g)
dim(agg_g_df) # 15847    16

# calc additional stats
agg_g_df$p_adj_sum=p.adjust(agg_g_df$p_value_sum, method = "BH")
agg_g_df$sig_sum=ifelse(agg_g_df$p_adj_sum<=0.05,'TRUE','FALSE')
agg_g_df$ASE_ratio=agg_g_df$P1_count_sum/(agg_g_df$P2_count_sum+agg_g_df$P1_count_sum)
agg_g_df$log2foldchange=log2((agg_g_df$P1_count_sum+1)/(agg_g_df$P2_count_sum+1))
agg_g_df$number_P1_alt[which(is.na(agg_g_df$number_P1_alt))]=0
agg_g_df$number_P2_alt[which(is.na(agg_g_df$number_P2_alt))]=0
agg_g_df$SNP_ratio=agg_g_df$number_P1_alt/(agg_g_df$number_P1_alt+agg_g_df$number_P2_alt)

# sort
agg_g_df = agg_g_df[order(agg_g_df$start_gene),]
agg_g_df = agg_g_df[order(match(agg_g_df$chrom,chrom_order)),]
agg_g_df$chrom = factor(agg_g_df$chrom, levels = unique(agg_g_df$chrom))
saveRDS(agg_g_df,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/ASEReadcounter_ENCFF379NOY/ASE_incl_intro_anno_V29_downloaded_ENCFF379NOY_gene_level.rds")
agg_g_df_ENCFF379NOY = agg_g_df
# first glance
table(agg_g_df$sig_sum)
# FALSE  TRUE
# 14910   416
table(agg_g_df$sig_sum)/nrow(agg_g_df)*100
#     FALSE      TRUE
# 97.285658  2.714342

table(agg_g_df$log2foldchange>=1)/nrow(agg_g_df)*100
#    FALSE     TRUE
# 81.33238 18.66762

table(agg_g_df$log2foldchange>=1,agg_g_df$sig_sum)/nrow(agg_g_df)*100 # that is not much at all!
#           FALSE      TRUE
# FALSE 79.675062  1.657314
# TRUE  17.610596  1.057027

table(agg_g_df$log2foldchange>=1,agg_g_df$sig_sum)
#       FALSE  TRUE
# FALSE 12211   254
# TRUE   2699   162

agg_g_df_LFC_ASE = agg_g_df[which(agg_g_df$log2foldchange >= 1 & agg_g_df$sig_sum==TRUE),]
length(unique(agg_g_df_LFC_ASE$gene_id)) # 162
length(unique(agg_g_df_LFC_ASE$gene_name)) # 162
# that does still include genes that are not protein coding

length(unique(agg_g_df_LFC_ASE$gene_id[which(agg_g_df_LFC_ASE$gene_type=="protein_coding")])) # 123
################################################################################
################################################################################ merge total and ase
################################################################################
# merge with expressed
total_RNA_red_auto$expressed_ENCFF675NTU = FALSE
total_RNA_red_auto$expressed_ENCFF675NTU[which(total_RNA_red_auto$TPM_ENCFF675NTU >= 1)] = TRUE
total_RNA_red_auto$expressed_ENCFF379NOY = FALSE
total_RNA_red_auto$expressed_ENCFF379NOY[which(total_RNA_red_auto$TPM_ENCFF379NOY >= 1)] = TRUE
table(total_RNA_red_auto$expressed_ENCFF675NTU)/nrow(total_RNA_red_auto)*100
#    FALSE     TRUE
# 62.43677 37.56323
table(total_RNA_red_auto$expressed_ENCFF379NOY)/nrow(total_RNA_red_auto)*100
#    FALSE     TRUE
# 63.90772 36.09228


# ASE_ENCFF675NTU = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/ASEReadcounter_ENCFF675NTU/ASE_incl_intro_anno_V29_downloaded_ENCFF675NTU.rds")
# ASE_ENCFF379NOY = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/ASEReadcounter_ENCFF379NOY/ASE_incl_intro_anno_V29_downloaded_ENCFF379NOY.rds")
ASE_ENCFF675NTU = agg_g_df_ENCFF675NTU
ASE_ENCFF379NOY = agg_g_df_ENCFF379NOY


#### merge total gene expression back to anno, so we have all genes in one df
annoV29_total_exp_merge = merge(annoV29_df_g_red,
                                total_RNA_red_auto[,c("gene_id","total_counts_ENCFF675NTU","total_counts_ENCFF379NOY","TPM_ENCFF675NTU","TPM_ENCFF379NOY","expressed_ENCFF675NTU","expressed_ENCFF379NOY")],
                                by="gene_id",all.x=T) # 58721    15
annoV29_total_exp_merge$total_counts_ENCFF675NTU[which(is.na(annoV29_total_exp_merge$total_counts_ENCFF675NTU))] = 0 
annoV29_total_exp_merge$total_counts_ENCFF379NOY[which(is.na(annoV29_total_exp_merge$total_counts_ENCFF379NOY))] = 0 
annoV29_total_exp_merge$TPM_ENCFF675NTU[which(is.na(annoV29_total_exp_merge$TPM_ENCFF675NTU))] = 0 
annoV29_total_exp_merge$TPM_ENCFF379NOY[which(is.na(annoV29_total_exp_merge$TPM_ENCFF379NOY))] = 0 
annoV29_total_exp_merge$expressed_ENCFF675NTU[which(is.na(annoV29_total_exp_merge$expressed_ENCFF675NTU))] = FALSE
annoV29_total_exp_merge$expressed_ENCFF379NOY[which(is.na(annoV29_total_exp_merge$expressed_ENCFF379NOY))] = FALSE 





anno_total_ase = merge(annoV29_total_exp_merge,
                       ASE_ENCFF675NTU[,c("gene_id","P1_count_sum","P2_count_sum","p_value_sum","number_SNPs","number_P1_alt","number_P2_alt","number_intron_SNPs","number_exon_SNPs","p_adj_sum","sig_sum","ASE_ratio","log2foldchange","SNP_ratio")],
                       by="gene_id",all.x=T)

colnames(anno_total_ase) = c("gene_id","chrom","start_gene","end_gene","width_gene","strand","gene_type","gene_name","level","total_counts_ENCFF675NTU","total_counts_ENCFF379NOY","TPM_ENCFF675NTU","TPM_ENCFF379NOY","expressed_ENCFF675NTU","expressed_ENCFF379NOY","P1_count_sum_ENCFF675NTU","P2_count_sum_ENCFF675NTU","p_value_sum_ENCFF675NTU","number_SNPs_ENCFF675NTU","number_P1_alt_ENCFF675NTU","number_P2_alt_ENCFF675NTU","number_intron_SNPs_ENCFF675NTU","number_exon_SNPs_ENCFF675NTU","p_adj_sum_ENCFF675NTU","sig_sum_ENCFF675NTU","ASE_ratio_ENCFF675NTU","log2foldchange_ENCFF675NTU","SNP_ratio_ENCFF675NTU")

anno_total_ase = merge(anno_total_ase,
                       ASE_ENCFF379NOY[,c("gene_id","P1_count_sum","P2_count_sum","p_value_sum","number_SNPs","number_P1_alt","number_P2_alt","number_intron_SNPs","number_exon_SNPs","p_adj_sum","sig_sum","ASE_ratio","log2foldchange","SNP_ratio")],
                       by="gene_id",all.x=T)
colnames(anno_total_ase) = c("gene_id","chrom","start_gene","end_gene","width_gene","strand","gene_type","gene_name","level","total_counts_ENCFF675NTU","total_counts_ENCFF379NOY","TPM_ENCFF675NTU","TPM_ENCFF379NOY","expressed_ENCFF675NTU","expressed_ENCFF379NOY",
                             "P1_count_sum_ENCFF675NTU","P2_count_sum_ENCFF675NTU","p_value_sum_ENCFF675NTU","number_SNPs_ENCFF675NTU","number_P1_alt_ENCFF675NTU","number_P2_alt_ENCFF675NTU","number_intron_SNPs_ENCFF675NTU","number_exon_SNPs_ENCFF675NTU","p_adj_sum_ENCFF675NTU","sig_sum_ENCFF675NTU","ASE_ratio_ENCFF675NTU","log2foldchange_ENCFF675NTU","SNP_ratio_ENCFF675NTU",
                             "P1_count_sum_ENCFF379NOY","P2_count_sum_ENCFF379NOY","p_value_sum_ENCFF379NOY","number_SNPs_ENCFF379NOY","number_P1_alt_ENCFF379NOY","number_P2_alt_ENCFF379NOY","number_intron_SNPs_ENCFF379NOY","number_exon_SNPs_ENCFF379NOY","p_adj_sum_ENCFF379NOY","sig_sum_ENCFF379NOY","ASE_ratio_ENCFF379NOY","log2foldchange_ENCFF379NOY","SNP_ratio_ENCFF379NOY")
# keep the NA here for now, this makes it easier to understand which transcripts do not overlap SNPs

dim(anno_total_ase) # 58721    41
# sort
anno_total_ase = anno_total_ase[order(anno_total_ase$start_gene),]
anno_total_ase = anno_total_ase[order(match(anno_total_ase$chrom,chrom_order)),]

# add a quick column showing if transcript is expressed (TPM >=1) with sig ASE and LFC >=1
anno_total_ase$ENCFF675NTU_exp_ase_lfc = FALSE
anno_total_ase$ENCFF675NTU_exp_ase_lfc[which(anno_total_ase$expressed_ENCFF675NTU == TRUE & anno_total_ase$sig_sum_ENCFF675NTU == TRUE & abs(anno_total_ase$log2foldchange_ENCFF675NTU) >= 1)] = TRUE 
anno_total_ase$ENCFF379NOY_exp_ase_lfc = FALSE
anno_total_ase$ENCFF379NOY_exp_ase_lfc[which(anno_total_ase$expressed_ENCFF379NOY == TRUE & anno_total_ase$sig_sum_ENCFF379NOY == TRUE & abs(anno_total_ase$log2foldchange_ENCFF379NOY) >= 1)] = TRUE 

saveRDS(anno_total_ase,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/total_exp_and_ASE_incl_intro_anno_V29_ENCFF675NTU_ENCFF379NOY_gene_level.rds")


table(anno_total_ase$ENCFF675NTU_exp_ase_lfc)
# FALSE   TRUE
# 58425   296
table(anno_total_ase$ENCFF675NTU_exp_ase_lfc)/nrow(anno_total_ase)*100 # almost identical to the transcript values
# FALSE      TRUE
# 99.4959214  0.5040786

table(anno_total_ase$ENCFF379NOY_exp_ase_lfc)
# FALSE   TRUE
# 58432   289
table(anno_total_ase$ENCFF379NOY_exp_ase_lfc)/nrow(anno_total_ase)*100
# FALSE       TRUE
# 99.5078422  0.4921578

table(anno_total_ase$ENCFF675NTU_exp_ase_lfc,anno_total_ase$ENCFF379NOY_exp_ase_lfc)
#        FALSE   TRUE
# FALSE  58359     66
# TRUE      73    223

#exp_ase_genes = intersect(anno_total_ase$gene_name[which(anno_total_ase$ENCFF675NTU_exp_ase_lfc==TRUE)],anno_total_ase$gene_name[which(anno_total_ase$ENCFF379NOY_exp_ase_lfc==TRUE)])
################################################################################
# merge with imprinted
imprinted_genes_unroll = readRDS("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/imprinted_genes_h1_project_unroll_aliases.rds")

length(intersect(anno_total_ase$gene_name,imprinted_genes_unroll$Gene)) # 102

imp_exp_ase_genes = intersect(exp_ase_genes,imprinted_genes_unroll$Gene)
# "ZDBF2" "PEG10" "NDN"   "SNRPN" "SNURF" "PWAR6" 

anno_total_ase_imp = merge(anno_total_ase,
                           imprinted_genes_unroll[,c("Gene","Expressed.Allele","chr","start","end","start_plot","end_plot","col")],
                           by.y="Gene",by.x="gene_name",
                           all.x = T,all.y = F)

anno_total_ase_imp = anno_total_ase_imp[order(anno_total_ase_imp$start_gene),]
anno_total_ase_imp = anno_total_ase_imp[order(as.numeric(substring(anno_total_ase_imp$chrom, 4))),]
saveRDS(anno_total_ase_imp,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY_gene_level.rds")
anno_total_ase_only_imp = anno_total_ase_imp[which(!(is.na(anno_total_ase_imp$Expressed.Allele))),] # 103  50
table(anno_total_ase_only_imp$ENCFF675NTU_exp_ase_lfc) # 6 expressed allele-specificly
# FALSE  TRUE
#    97     6
table(anno_total_ase_only_imp$expressed_ENCFF675NTU) # roughly half expressed
# FALSE  TRUE
#    57    46

table(is.na(anno_total_ase_only_imp$ASE_ratio_ENCFF675NTU[which(anno_total_ase_only_imp$expressed_ENCFF675NTU==TRUE)]),exclude=NULL) # of those 46 that are expressed, 23 don't have SNPs
#  FALSE  TRUE
#     23    23
# those TRUE don't have SNPs in the transcript body

anno_total_ase_only_imp[anno_total_ase_only_imp$ENCFF675NTU_exp_ase_lfc==TRUE,]
unique(anno_total_ase_only_imp$gene_name[anno_total_ase_only_imp$ENCFF379NOY_exp_ase_lfc==TRUE,]) #"ZDBF2" "PEG10" "NDN"   "SNRPN" "SNURF" "PWAR6" "DNMT1" 
# all paternal, 
# 5 protcod, 1 lincRNA
# on chr 2,7,15
# almost all intron SNPs
# take care, the coordinates from the imprinted source does not exactly match the anno coordinates, actually they differ quite mutch, 17 or 44 50k bins:
# > 206274663-207129522
# [1] -854859
# > 94656325-94275636
# [1] 380689
# > 23685400-21471646
# [1] 2213754
# > 24978723-25058793
# [1] -80070
# > 24954986-22741227
# [1] 2213759
# > 25031873-25021872
# [1] 10001
# plot with the anno coordinates the allele-specific readcounts were aggregated over

### number of genes with (observed) SNP info
length(anno_total_ase_imp$gene_id[which(!(is.na(anno_total_ase_imp$number_SNPs_ENCFF675NTU)))]) # 15847
length(anno_total_ase_imp$gene_id[which(!(is.na(anno_total_ase_imp$number_SNPs_ENCFF379NOY)))]) # 15326

################################################################################
### comparison ENCODE data 

# merge with ENCODE gene quantification of these replicates 

ENCFF174OMR_g = read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/expression/RNA_Seq_downloaded_gene_quant/ENCFF174OMR/ENCFF174OMR.tsv",sep="\t") # 59526    17
ENCFF910OBU_g = read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/expression/RNA_Seq_downloaded_gene_quant/ENCFF910OBU/ENCFF910OBU.tsv",sep="\t") # 59526    17

# filter for entrys of genes marked with ENSG 
ENCFF174OMR_g_subset <- ENCFF174OMR_g[grep("ENSG", ENCFF174OMR_g$gene_id), ] 
# dim(ENCFF174OMR_g) # [1] 59526    17
# dim(ENCFF174OMR_g_subset) # [1] 58780    17
ENCFF910OBU_g_subset <- ENCFF910OBU_g[grep("ENSG", ENCFF910OBU_g$gene_id), ] 
# dim(ENCFF910OBU_g) # [1] 59526    17
# dim(ENCFF910OBU_g_subset) # [1] 58780    17


colnames(ENCFF174OMR_g_subset) = paste(colnames(ENCFF174OMR_g_subset),"_ENCFF174OMR",sep="")
anno_total_ase_imp_merged_ENCODE_g = merge(anno_total_ase_imp,ENCFF174OMR_g_subset[,c("gene_id_ENCFF174OMR","length_ENCFF174OMR","effective_length_ENCFF174OMR","expected_count_ENCFF174OMR","TPM_ENCFF174OMR","posterior_mean_count_ENCFF174OMR","pme_TPM_ENCFF174OMR")],
                                           by.x = "gene_id",by.y = "gene_id_ENCFF174OMR",all.x = TRUE)
colnames(ENCFF910OBU_g_subset) = paste(colnames(ENCFF910OBU_g_subset),"_ENCFF910OBU",sep="")
anno_total_ase_imp_merged_ENCODE_g = merge(anno_total_ase_imp_merged_ENCODE_g,ENCFF910OBU_g_subset[,c("gene_id_ENCFF910OBU","length_ENCFF910OBU","effective_length_ENCFF910OBU","expected_count_ENCFF910OBU","TPM_ENCFF910OBU","posterior_mean_count_ENCFF910OBU","pme_TPM_ENCFF910OBU")],
                                           by.x = "gene_id",by.y = "gene_id_ENCFF910OBU",all.x = TRUE)
dim(anno_total_ase_imp_merged_ENCODE_g) #  58721    62
# # why does this not have 58780 columns? bc in the V29 anno, there are only so many genes, I do not know why the downloaded gene quantifications have more genes than the anno
# # maybe it is just the value after the dot?
# diff_genes = setdiff(ENCFF910OBU_g_subset$gene_id_ENCFF910OBU,annoV29_df_g_red$gene_id) # 59
# gene_id_ENCFF910OBU = sapply(strsplit(as.character(ENCFF910OBU_g_subset$gene_id_ENCFF910OBU),".",fixed=TRUE), `[`, 1)
# gene_id_anno = sapply(strsplit(as.character(annoV29_df_g_red$gene_id),".",fixed=TRUE), `[`, 1)
# diff_genes_full_id = setdiff(gene_id_ENCFF910OBU,gene_id_anno) # 59
# gene_id_anno_2 = sapply(strsplit(as.character(unique(annoV29_df$gene_id)),".",fixed=TRUE), `[`, 1) # same
# diff_genes_full_id = setdiff(gene_id_ENCFF910OBU,gene_id_anno_2) # 59, same
# # nope, still the same number of genes missing
# # googling showed, those are genes with "No location information available", which might be located on scaffolds eg here: Un_KI270713.1

anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU = FALSE
anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR = FALSE
anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU[which(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF910OBU >= 1)] = TRUE
anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR[which(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF174OMR >= 1)] = TRUE

anno_total_ase_imp_merged_ENCODE_g$ENCFF910OBU_exp_ase_lfc = FALSE
anno_total_ase_imp_merged_ENCODE_g$ENCFF910OBU_exp_ase_lfc[which(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU == TRUE & anno_total_ase_imp_merged_ENCODE_g$sig_sum_ENCFF675NTU == TRUE & abs(anno_total_ase_imp_merged_ENCODE_g$log2foldchange_ENCFF675NTU) >= 1)] = TRUE
anno_total_ase_imp_merged_ENCODE_g$ENCFF174OMR_exp_ase_lfc = FALSE
anno_total_ase_imp_merged_ENCODE_g$ENCFF174OMR_exp_ase_lfc[which(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR == TRUE & anno_total_ase_imp_merged_ENCODE_g$sig_sum_ENCFF379NOY == TRUE & abs(anno_total_ase_imp_merged_ENCODE_g$log2foldchange_ENCFF379NOY) >= 1)] = TRUE

saveRDS(anno_total_ase_imp_merged_ENCODE_g,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY_ENCFF174OMR_ENCFF910OBU_gene_level.rds")
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY_ENCFF174OMR_ENCFF910OBU_gene_level.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/data/expression_data/

anno_total_ase_imp_merged_ENCODE_g = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/expression_data/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY_ENCFF174OMR_ENCFF910OBU_gene_level.rds")
# expression = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/expression_data/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY_ENCFF174OMR_ENCFF910OBU_gene_level.rds")
# expression = expression[order(expression$start_gene),]
# expression = expression[order(as.numeric(substring(expression$chrom, 4))),]
# write.table(expression,"/Users/jmarkow/Desktop/GAM/Co-Phasing/data/expression_data/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY_ENCFF174OMR_ENCFF910OBU_gene_level.tsv",
#           quote = FALSE, sep = "\t")
# NTU(ali) - IXN(tr_quant) - OBU(g_quant) - NWP(tr_quant) 
# NOY(ali) - VMY(tr_quant) - OMR(g_quant) - AFW(tr_quant) 

# TPMs genes
pairs(~TPM_ENCFF675NTU+TPM_ENCFF379NOY+TPM_ENCFF910OBU+TPM_ENCFF174OMR,data=anno_total_ase_imp_merged_ENCODE_g, 
      main="TPM comparison")

pairs(~log10(TPM_ENCFF675NTU+1)+log10(TPM_ENCFF379NOY+1)+log10(TPM_ENCFF910OBU+1)+log10(TPM_ENCFF174OMR+1),data=anno_total_ase_imp_merged_ENCODE_g, 
      main="TPM comparison")

plot(log10(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF675NTU+1),
           log10(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF910OBU+1),
     pch = 19, col = adjustcolor("black",alpha.f = 0.1),
     xlab = "log10(TPM_ENCFF675NTU+1)",
     ylab = "log10(TPM_ENCFF910OBU+1)")
abline(coef = c(0,1), col="red")

# raw counts
pairs(~total_counts_ENCFF675NTU+total_counts_ENCFF379NOY+expected_count_ENCFF910OBU+expected_count_ENCFF174OMR,data=anno_total_ase_imp_merged_ENCODE_g, 
      main="count comparison")

pairs(~log10(total_counts_ENCFF675NTU+1)+log10(total_counts_ENCFF379NOY+1)+log10(expected_count_ENCFF910OBU+1)+log10(expected_count_ENCFF174OMR+1),data=anno_total_ase_imp_merged_ENCODE_g, 
      main="count comparison")

plot(log10(anno_total_ase_imp_merged_ENCODE_g$total_counts_ENCFF675NTU+1),
     log10(anno_total_ase_imp_merged_ENCODE_g$expected_count_ENCFF910OBU+1),
     pch = 19, col = adjustcolor("black",alpha.f = 0.1),
     xlab = "log10(total_counts_ENCFF675NTU+1)",
     ylab = "log10(expected_count_ENCFF910OBU+1)")
abline(coef = c(0,1), col="red")

# length / width to consider
plot(log10(anno_total_ase_imp_merged_ENCODE_g$width_gene+1),
     log10(anno_total_ase_imp_merged_ENCODE_g$length_ENCFF174OMR+1),
     pch = 19, col = adjustcolor("black",alpha.f = 0.1),
     xlab = "log10(length_gene_ENCFF675NTU+1)",
     ylab = "log10(length_gene_ENCFF174OMR+1)")
abline(coef = c(0,1), col="red")

plot(anno_total_ase_imp_merged_ENCODE_g$width_gene,
     anno_total_ase_imp_merged_ENCODE_g$length_ENCFF174OMR,
     pch = 19, col = adjustcolor("black",alpha.f = 0.1),
     xlab = "length_gene_ENCFF675NTU",
     ylab = "length_gene_ENCFF174OMR")
abline(coef = c(0,1), col="red")


table(anno_total_ase_imp_merged_ENCODE_g$ENCFF174OMR_exp_ase_lfc,anno_total_ase_imp_merged_ENCODE_g$ENCFF379NOY_exp_ase_lfc)/nrow(anno_total_ase_imp_merged_ENCODE_g)*100
#             FALSE        TRUE
# FALSE 99.46356499  0.07152467
# TRUE   0.04427718  0.42063316

table(anno_total_ase_imp_merged_ENCODE_g$ENCFF910OBU_exp_ase_lfc,anno_total_ase_imp_merged_ENCODE_g$ENCFF675NTU_exp_ase_lfc)/nrow(anno_total_ase_imp_merged_ENCODE_g)*100
#             FALSE        TRUE
# FALSE 99.45164422  0.06811873
# TRUE   0.04427718  0.43595988

table(anno_total_ase_imp_merged_ENCODE_g$ENCFF910OBU_exp_ase_lfc,anno_total_ase_imp_merged_ENCODE_g$ENCFF174OMR_exp_ase_lfc)/nrow(anno_total_ase_imp_merged_ENCODE_g)*100
#            FALSE       TRUE
# FALSE 99.4073670  0.1123959
# TRUE   0.1277226  0.3525144

table(anno_total_ase_imp_merged_ENCODE_g$ENCFF675NTU_exp_ase_lfc,anno_total_ase_imp_merged_ENCODE_g$ENCFF379NOY_exp_ase_lfc)/nrow(anno_total_ase_imp_merged_ENCODE_g)*100
#            FALSE       TRUE
# FALSE 99.3835255  0.1123959
# TRUE   0.1243167  0.3797619

# There is more difference between the replicates than between the methods!


library("ggVennDiagram")

expr_tr_list <- list(
    ENCFF675NTU = anno_total_ase_imp_merged_ENCODE_g$gene_id[which(anno_total_ase_imp_merged_ENCODE_g$ENCFF675NTU_exp_ase_lfc)], # 296
    ENCFF910OBU = anno_total_ase_imp_merged_ENCODE_g$gene_id[which(anno_total_ase_imp_merged_ENCODE_g$ENCFF910OBU_exp_ase_lfc)], # 282
    ENCFF379NOY = anno_total_ase_imp_merged_ENCODE_g$gene_id[which(anno_total_ase_imp_merged_ENCODE_g$ENCFF379NOY_exp_ase_lfc)], # 289
    ENCFF174OMR = anno_total_ase_imp_merged_ENCODE_g$gene_id[which(anno_total_ase_imp_merged_ENCODE_g$ENCFF174OMR_exp_ase_lfc)]) # 273


# lapply(expr_tr_list,length)
ggVennDiagram(expr_tr_list, label_alpha = 0, label_color="white") 

################################################################################
### volcano plot - gene level

#################### ENCFF675NTU
data=anno_total_ase_imp_merged_ENCODE_g[anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE,] # 20940    66
agg='p_adj_sum_ENCFF675NTU' 
title='Allele-specific Gene Expression hESC H1 - ENCFF675NTU'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001

par(mfrow=c(1,1))

plot(x=data$log2foldchange_ENCFF675NTU,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=700,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF675NTU>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF675NTU<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=400,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=400,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Paternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Maternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
# data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # -
# data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "PWAR6" "SNURF"

#################### ENCFF910OBU
data=anno_total_ase_imp_merged_ENCODE_g[anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU==TRUE,] # 15278    66
agg='p_adj_sum_ENCFF675NTU'# bc this is the data it came from
title='Allele-specific Gene Expression hESC H1 - ENCFF910OBU'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001

par(mfrow=c(1,1))

plot(x=data$log2foldchange_ENCFF675NTU,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=700,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF675NTU>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF675NTU<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=400,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=400,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "KCNQ1"(chr11) "ATP10A" (chr15), both protcod
data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "PWAR6" "SNURF"

data[which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff)),] # "KCNQ1"(chr11) "ATP10A" (chr15), both protcod
data[which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff)),] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "PWAR6" "SNURF"

#################### ENCFF379NOY
data=anno_total_ase_imp_merged_ENCODE_g[anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE,] # 20120    66
agg='p_adj_sum_ENCFF379NOY'
title='Allele-specific Gene Expression hESC H1 - ENCFF379NOY'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001


plot(x=data$log2foldchange_ENCFF379NOY,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=700,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF379NOY>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF379NOY<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=400,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=400,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Paternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Maternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # 
data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "PWAR6" "SNURF"


#################### ENCFF174OMR
data=anno_total_ase_imp_merged_ENCODE_g[anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR==TRUE,] # 15001    66
agg='p_adj_sum_ENCFF379NOY'# bc this is the data it came from
title='Allele-specific Gene Expression hESC H1 - ENCFF174OMR'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001

par(mfrow=c(1,1))

plot(x=data$log2foldchange_ENCFF379NOY,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=700,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF379NOY>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF379NOY<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=400,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=400,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Paternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Maternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "KCNQ1"
data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "PWAR6" "SNURF"

######################################## 
######################################## 
######################################## 
#### only protein coding genes
#################### ENCFF379NOY
data=anno_total_ase_imp_merged_ENCODE_g[which(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE & anno_total_ase_imp_merged_ENCODE_g$gene_type== "protein_coding"),] # 11305    66
agg='p_adj_sum_ENCFF379NOY'
title='Allele-specific protein coding Gene Expression hESC H1 - ENCFF379NOY'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001


plot(x=data$log2foldchange_ENCFF379NOY,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=700,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF379NOY>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF379NOY<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=400,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=400,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Paternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Maternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # 
data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF379NOY_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "SNURF"

#################### ENCFF675NTU
data=anno_total_ase_imp_merged_ENCODE_g[which(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE & anno_total_ase_imp_merged_ENCODE_g$gene_type== "protein_coding"),] # 11569    66
agg='p_adj_sum_ENCFF675NTU' 
title='Allele-specific protein coding Gene Expression hESC H1 - ENCFF675NTU'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001

par(mfrow=c(1,1))

plot(x=data$log2foldchange_ENCFF675NTU,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=400,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF675NTU>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF675NTU<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=200,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=200,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Paternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Maternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # -
data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF675NTU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "SNURF"

#################### ENCFF910OBU
data=anno_total_ase_imp_merged_ENCODE_g[which(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU==TRUE & anno_total_ase_imp_merged_ENCODE_g$gene_type== "protein_coding"),] # 12909    66
agg='p_adj_sum_ENCFF675NTU'# bc this is the data it came from
title='Allele-specific protein coding Gene Expression hESC H1 - ENCFF910OBU'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001

par(mfrow=c(1,1))

plot(x=data$log2foldchange_ENCFF675NTU,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF675NTU[which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF675NTU<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=400,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF675NTU>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF675NTU<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=200,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=200,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF675NTU[which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "KCNQ1"(chr11) "ATP10A" (chr15), both protcod
data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "SNURF"

data[which((data$Expressed.Allele=='Maternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff)),] # "KCNQ1"(chr11) "ATP10A" (chr15), both protcod
data[which((data$Expressed.Allele=='Paternal')&(data$ENCFF910OBU_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff)),] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "SNURF"

#################### ENCFF174OMR
data=anno_total_ase_imp_merged_ENCODE_g[which(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR==TRUE & anno_total_ase_imp_merged_ENCODE_g$gene_type== "protein_coding"),] # 12787    66
agg='p_adj_sum_ENCFF379NOY'# bc this is the data it came from
title='Allele-specific protein coding Gene Expression hESC H1 - ENCFF174OMR'
col_mom='indianred2'
col_dad='deepskyblue'
p_val_cutoff=0.05 #0.01 #0.001

par(mfrow=c(1,1))

plot(x=data$log2foldchange_ENCFF379NOY,y=-log(data[,agg]),xlab='Log2foldchange',ylab='-log(adjusted p-value)', pch=19,main=title,col='grey75')
# adj p value cut off
abline(h=-log(p_val_cutoff))
text(x=-6,y=-10,label=paste('p-adj >= ', p_val_cutoff))

# L2FC 
abline(v=c(-1,1))
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY>=1)&(data[,agg]<=p_val_cutoff))]),col=col_dad,pch=19)
points(x=data$log2foldchange_ENCFF379NOY[which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$log2foldchange_ENCFF379NOY<=(-1))&(data[,agg]<=p_val_cutoff))]),col=col_mom,pch=19)

text(x=5,y=700,label='|Log2foldchange| >= 1')
nr_ase_genes_dad=table((data$log2foldchange_ENCFF379NOY>=1),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
nr_ase_genes_mom=table((data$log2foldchange_ENCFF379NOY<=(-1)),(data[,agg]<=p_val_cutoff))['TRUE','TRUE']
text(x=-5,y=400,label=nr_ase_genes_mom,col = col_mom)
text(x=5,y=400,label=nr_ase_genes_dad,col = col_dad)

# imprinted

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Paternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Paternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='blue')

points(x=data$log2foldchange_ENCFF379NOY[which((data$Expressed.Allele=='Maternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))],
       y=-log(data[,agg][which((data$Expressed.Allele=='Maternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE))]),
       pch=21,col='black',bg='red4')
legend('topleft',pch=c(19,19),col=c(col_dad,col_mom,'blue','red4'),legend=c('Paternal Allele expressed','Maternal Allele expressed','Imprinted Paternal Allele','Imprinted Maternal Allele'),bty='n')
data$gene_name[which((data$Expressed.Allele=='Maternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "KCNQ1"
data$gene_name[which((data$Expressed.Allele=='Paternal')&(data$ENCFF174OMR_exp_ase_lfc==TRUE)&(data[,agg]<=p_val_cutoff))] # "SNRPN" "NDN"   "ZDBF2" "PEG10" "SNURF"


