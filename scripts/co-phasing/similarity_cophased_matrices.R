# structural similarity index 
# compare downsampled cophased F123 matrices with traditionally phased F123 matrices
    # complete matrix, all vs all (CAST vs J129, cophased vs traditionally phased)
        # copmare the results with the results from python
    # use a sliding window along the diagonal, exclude most of the NAs, take the mean of resulting values
    # with the sliding window, maybe also compare different chromosomes of the same type (haplotype, phased)

# in H1 conda env:
# conda install -c conda-forge r-spatialpack 

library("SpatialPack")

#### load all matrices
dir_matrices="/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices"

### CAST - traditional phasing
CAST_orig = vector(mode="list",length=19)

### CAST - subsampled cophasing
CAST_coph = vector(mode="list",length=19)

### J129 - traditional phasing
J129_orig = vector(mode="list",length=19)

### J129 - subsampled cophasing
J129_coph = vector(mode="list",length=19)


for (chr in 1:19){
    # chr = 6
    print(chr)
    CAST_orig[[chr]] = data.matrix(read.table(gzfile(paste(dir_matrices,"/CAST.original.res_50K/segtable.res_50000.CAST.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    CAST_coph[[chr]] = data.matrix(read.table(gzfile(paste(dir_matrices,"/CAST.dist_cutoff_100K.res_50K/segtable.dist_cutoff_100000.res_50000.CAST.iteration_1.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    J129_orig[[chr]] = data.matrix(read.table(gzfile(paste(dir_matrices,"/S129.original.res_50K/segtable.res_50000.S129.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    J129_coph[[chr]] = data.matrix(read.table(gzfile(paste(dir_matrices,"/S129.dist_cutoff_100K.res_50K/segtable.dist_cutoff_100000.res_50000.S129.iteration_1.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
}
#chr = 19
#CAST_orig[[chr]] = data.matrix(read.table(gzfile(paste(dir_matrices,"/CAST.original.res_50K/segtable.res_50000.CAST.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))

# replace Na with 0
CAST_orig = lapply(CAST_orig,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
CAST_coph = lapply(CAST_coph,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
J129_orig = lapply(J129_orig,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
J129_coph = lapply(J129_coph,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})


SSIM(CAST_orig[[chr]], CAST_coph[[chr]], alpha = 1, beta = 1, gamma = 1, eps = c(0.01, 0.03), L = 1)
# *** caught segfault ***
#     address 0x55c959537000, cause 'memory not mapped'
# 
# Traceback:
#     1: SSIM(CAST_orig[[chr]], CAST_coph[[chr]], alpha = 1, beta = 1,     gamma = 1, eps = c(0.01, 0.03), L = 1)
# 
# Possible actions:
#     1: abort (with core dump, if enabled)
# 2: normal R exit
# 3: exit R without saving workspace
# 4: exit R saving workspace
# Selection:
    

chr = 19
CAST_orig = data.matrix(read.table(gzfile(paste(dir_matrices,"/CAST.original.res_50K/segtable.res_50000.CAST.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
CAST_coph = data.matrix(read.table(gzfile(paste(dir_matrices,"/CAST.dist_cutoff_100K.res_50K/segtable.dist_cutoff_100000.res_50000.CAST.iteration_1.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
J129_orig = data.matrix(read.table(gzfile(paste(dir_matrices,"/S129.original.res_50K/segtable.res_50000.S129.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
J129_coph = data.matrix(read.table(gzfile(paste(dir_matrices,"/S129.dist_cutoff_100K.res_50K/segtable.dist_cutoff_100000.res_50000.S129.iteration_1.NPMI.chr",chr,".txt.gz",sep="")),header=T,row.names=1,sep = "\t"))

# replace Na with 0
CAST_orig_0 = CAST_orig
CAST_coph_0 = CAST_coph
J129_orig_0 = J129_orig
J129_coph_0 = J129_coph

CAST_orig_0[is.na(CAST_orig_0)] <- 0
CAST_coph_0[is.na(CAST_coph_0)] <- 0
J129_orig_0[is.na(J129_orig_0)] <- 0
J129_coph_0[is.na(J129_coph_0)] <- 0


SSIM(CAST_orig_0, CAST_orig_0, L = 1) # 1
SSIM(CAST_orig_0, CAST_coph_0, L = 1) 
# $SSIM
# [1] 0.7722388
# $comps
# luminance  contrast structure
# 0.9945795 0.9999937 0.7764525
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.19960403 0.02310282 0.22157276 0.02326943 0.01790221

# SSIM(CAST_orig_0, CAST_coph_0, L = 1)  == SSIM(CAST_coph_0, CAST_orig_0, L = 1) 

SSIM(J129_orig_0, J129_coph_0, L = 1)
# $SSIM
# [1] 0.7407316
# $comps
# luminance  contrast structure
# 0.9958924 0.9999986 0.7437878
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.18703702 0.02198259 0.20482551 0.02190926 0.01620779

SSIM(CAST_orig_0, J129_orig_0, L = 1)
# $SSIM
# [1] 0.5351005
# $comps
# luminance  contrast structure
# 0.9978921 0.9996973 0.5363932
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.19960403 0.02310282 0.18703702 0.02198259 0.01187940


SSIM(CAST_coph_0, J129_coph_0, L = 1)
# $SSIM
# [1] 0.5393301
# $comps
# luminance  contrast structure
# 0.9969229 0.9995556 0.5412353
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.22157276 0.02326943 0.20482551 0.02190926 0.01201417

### results are very different from the results calculated in python, where similar matrices had values around 57 and opposing matrices around 20
# here we have more details in which part of ssim has most influence, and as expected, the structure part shows the most differences



# different NA replacements

# replace Na with - 1 # does not work, bc NA get introduce in results (luminance) due to drastic negative influence
# replace Na with 3, stick mith L = 1
CAST_orig_3 = CAST_orig
CAST_coph_3 = CAST_coph
J129_orig_3 = J129_orig
J129_coph_3 = J129_coph

CAST_orig_3[is.na(CAST_orig_3)] <- -3
CAST_coph_3[is.na(CAST_coph_3)] <- -3
J129_orig_3[is.na(J129_orig_3)] <- -3
J129_coph_3[is.na(J129_coph_3)] <- -3


SSIM(CAST_orig_3, CAST_coph_3, L = 1) 
# 0.6350533 
# luminance  contrast structure
# 0.8674943 0.9942115 0.7363169

SSIM(J129_orig_3, J129_coph_3, L = 1)
# 0.6545031
# luminance  contrast structure
# 0.9256665 0.9963665 0.7096399

SSIM(CAST_orig_3, J129_orig_3, L = 1)
# 0.718988
# luminance  contrast structure
# 0.9909434 0.9997231 0.7257600

SSIM(CAST_coph_3, J129_coph_3, L = 1)
# 0.6688805
# luminance  contrast structure
# 0.9612822 0.9989422 0.6965580

SSIM(CAST_coph_3, J129_coph_3, L = 3)
# 0.669718
# luminance  contrast structure
# 0.9614673 0.9989448 0.6972941

# almost the same with L = 1 and L = 3 

# let's see if and how differences between different chromosomes is possible

# CAST orig vs cophased
SSIM_CAST = vector(mode="numeric",length=19)
luminance_CAST = vector(mode="numeric",length=19)
contrast_CAST = vector(mode="numeric",length=19)
structure_CAST = vector(mode="numeric",length=19)
# J129 orig vs cophased
SSIM_J129 = vector(mode="numeric",length=19)
luminance_J129 = vector(mode="numeric",length=19)
contrast_J129 = vector(mode="numeric",length=19)
structure_J129 = vector(mode="numeric",length=19)
# CAST vs J129 orig
SSIM_orig = vector(mode="numeric",length=19)
luminance_orig = vector(mode="numeric",length=19)
contrast_orig = vector(mode="numeric",length=19)
structure_orig = vector(mode="numeric",length=19)
# CAST vs J129 cophased
SSIM_coph = vector(mode="numeric",length=19)
luminance_coph = vector(mode="numeric",length=19)
contrast_coph = vector(mode="numeric",length=19)
structure_coph = vector(mode="numeric",length=19)



for (chr in 1:19){
    # chr = 10
    print(chr)
    ssim_all_cast = SSIM(CAST_orig[[chr]], CAST_coph[[chr]], L = 1) 
    SSIM_CAST[[chr]] = ssim_all_cast$SSIM
    luminance_CAST[[chr]] = ssim_all_cast$comps[1]
    contrast_CAST[[chr]] = ssim_all_cast$comps[2]
    structure_CAST[[chr]] = ssim_all_cast$comps[3]
    
    ssim_all_J129 = SSIM(J129_orig[[chr]], J129_coph[[chr]], L = 1) 
    SSIM_J129[[chr]] = ssim_all_J129$SSIM
    luminance_J129[[chr]] = ssim_all_J129$comps[1]
    contrast_J129[[chr]] = ssim_all_J129$comps[2]
    structure_J129[[chr]] = ssim_all_J129$comps[3]
    
    ssim_all_orig = SSIM(CAST_orig[[chr]], J129_orig[[chr]], L = 1) 
    SSIM_orig[[chr]] = ssim_all_orig$SSIM
    luminance_orig[[chr]] = ssim_all_orig$comps[1]
    contrast_orig[[chr]] = ssim_all_orig$comps[2]
    structure_orig[[chr]] = ssim_all_orig$comps[3]
    
    ssim_all_coph = SSIM(CAST_coph[[chr]], J129_coph[[chr]], L = 1) 
    SSIM_coph[[chr]] = ssim_all_coph$SSIM
    luminance_coph[[chr]] = ssim_all_coph$comps[1]
    contrast_coph[[chr]] = ssim_all_coph$comps[2]
    structure_coph[[chr]] = ssim_all_coph$comps[3]
}

saveRDS(list(SSIM_CAST,
             luminance_CAST,
             contrast_CAST,
             structure_CAST,
             SSIM_J129,
             luminance_J129,
             contrast_J129,
             structure_J129,
             SSIM_orig,
             luminance_orig,
             contrast_orig,
             structure_orig,
             SSIM_coph,
             luminance_coph,
             contrast_coph,
             structure_coph,
             SSIM_orig - SSIM_coph,
             SSIM_CAST - SSIM_J129),
        paste(dir_matrices,"SSIM_chr.rds",sep="/"))


# > SSIM_CAST
# [1] 0.7576193 0.7325226 0.7292646 0.7756147 0.7817201 0.7877008 0.8054583 0.7965067 0.7605551 0.7821187 0.8115980 0.7646442 0.7724997 0.7658144 0.7737549 0.7717518 0.7374430 0.8468244 0.7722388
# > SSIM_J129
# [1] 0.7409315 0.7240106 0.7343470 0.7570627 0.7733497 0.7621517 0.7737374 0.8131075 0.7668890 0.8055505 0.8272674 0.7280773 0.7603534 0.7607608 0.7595115 0.7819587 0.7413909 0.8226142 0.7407316
# > SSIM_orig
# [1] 0.3538783 0.4226191 0.3882641 0.4373973 0.4045846 0.4135557 0.4654670 0.4412108 0.4219609 0.3967572 0.4443722 0.3973478 0.4046012 0.4605636 0.4123515 0.3992182 0.4517728 0.4602127 0.5351005
# > SSIM_coph
# [1] 0.3553838 0.4230228 0.3860251 0.4526337 0.4232937 0.4101892 0.4744933 0.4524217 0.4237220 0.4157512 0.4639843 0.4183253 0.4108065 0.4724577 0.4094314 0.4115624 0.4582353 0.4815624 0.5393301
# >
#     > luminance_CAST
# [1] 0.9968194 0.9933010 0.9980556 0.9954432 0.9975994 0.9965446 0.9984183 0.9965657 0.9950355 0.9968661 0.9976567 0.9983570 0.9979616 0.9959338 0.9968626 0.9981226 0.9961923 0.9983189 0.9945795
# > luminance_J129
# [1] 0.9968230 0.9964263 0.9961757 0.9956712 0.9974524 0.9961185 0.9977205 0.9972842 0.9960735 0.9976458 0.9974867 0.9982142 0.9979802 0.9956724 0.9972577 0.9971287 0.9956455 0.9987789 0.9958924
# > luminance_orig
# [1] 0.9970119 0.9997048 0.9997472 0.9995966 0.9987498 0.9996135 0.9997787 0.9997893 0.9995011 0.9987292 0.9992522 0.9982320 0.9887599 0.9997120 0.9999752 0.9999430 0.9989810 0.9993599 0.9978921
# > luminance_coph
# [1] 0.9970155 0.9984504 0.9999964 0.9995246 0.9986430 0.9994620 0.9994771 0.9999360 0.9997898 0.9992064 0.9991547 0.9980839 0.9888032 0.9997768 0.9999983 0.9999925 0.9986881 0.9990165 0.9969229
# >
#     > contrast_CAST
# [1] 0.9996997 0.9994273 0.9997063 0.9997958 0.9998257 0.9997724 0.9994183 0.9997655 0.9998077 0.9997222 0.9998283 0.9996919 0.9994602 0.9992150 0.9997463 0.9998267 0.9997507 0.9997150 0.9999937
# > contrast_J129
# [1] 0.9998474 0.9996582 0.9998186 0.9996467 0.9998198 0.9999527 0.9995088 0.9998264 0.9999506 0.9998755 0.9997257 0.9997497 0.9997260 0.9993786 0.9997834 0.9998021 0.9998005 0.9996807 0.9999986
# > contrast_orig
# [1] 0.9996552 0.9999924 0.9993930 0.9995935 0.9993373 0.9999993 0.9999488 0.9999856 0.9997624 0.9999680 0.9998884 0.9999893 0.9986912 0.9999791 0.9997252 0.9999870 0.9999702 0.9999397 0.9996973
# > contrast_coph
# [1] 0.9998153 0.9999928 0.9991989 0.9997549 0.9993259 0.9999177 0.9999169 0.9999649 0.9995048 1.0000000 0.9998034 0.9999750 0.9991297 0.9999413 0.9996837 0.9999797 0.9999492 0.9999234 0.9995556
# >
#     > structure_CAST
# [1] 0.7602650 0.7378855 0.7309000 0.7793243 0.7837379 0.7906120 0.8072039 0.7994390 0.7644968 0.7847955 0.8136441 0.7661386 0.7744957 0.7695452 0.7763871 0.7733374 0.7404462 0.8484922 0.7764525
# > structure_J129
# [1] 0.7434064 0.7268558 0.7372999 0.7606228 0.7754647 0.7651577 0.7758862 0.8154634 0.7699501 0.8075520 0.8295794 0.7295625 0.7621011 0.7645424 0.7617651 0.7843656 0.7447820 0.8238831 0.7437878
# > structure_orig
# [1] 0.3550613 0.4227471 0.3885982 0.4377518 0.4053597 0.4137160 0.4655939 0.4413101 0.4222718 0.3972748 0.4447544 0.3980558 0.4097369 0.4607059 0.4124750 0.3992462 0.4522471 0.4605353 0.5363932
# > structure_coph
# [1] 0.3565135 0.4236824 0.3863360 0.4529601 0.4241548 0.4104438 0.4747810 0.4524666 0.4240211 0.4160815 0.4644682 0.4191389 0.4158202 0.4725910 0.4095616 0.4115739 0.4588606 0.4820734 0.5412353

# difference between (CAST vs J129 orig) and (CAST vs J129 cophased) --> these should be super similar

# abs(SSIM_CAST - SSIM_J129)
# [1] 0.016687744 0.008511969 0.005082419 0.018552033 0.008370377 0.025549079 0.031720969 0.016600860 0.006333874 0.023431857 0.015669400 0.036566894 0.012146323 0.005053675 0.014243359 0.010206935 0.003947950 0.024210143 0.031507123
# abs(SSIM_orig - SSIM_coph)
# [1] 0.0015054437 0.0004036754 0.0022390497 0.0152364224 0.0187090232 0.0033665274 0.0090263209 0.0112109182 0.0017611384 0.0189940197 0.0196120944 0.0209774804 0.0062053479 0.0118940579 0.0029201017 0.0123441815 0.0064625783 0.0213496460 0.0042295518


# make the plots now, see about the comparison later 
#scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/

ssim_first_direct = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr.rds")
SSIM_CAST = ssim_first_direct[[1]]
luminance_CAST = ssim_first_direct[[2]]
contrast_CAST = ssim_first_direct[[3]]
structure_CAST = ssim_first_direct[[4]]
SSIM_J129 = ssim_first_direct[[5]]
luminance_J129 = ssim_first_direct[[6]]
contrast_J129 = ssim_first_direct[[7]]
structure_J129 = ssim_first_direct[[8]]
SSIM_orig = ssim_first_direct[[9]]
luminance_orig = ssim_first_direct[[10]]
contrast_orig = ssim_first_direct[[11]]
structure_orig = ssim_first_direct[[12]]
SSIM_coph = ssim_first_direct[[13]]
luminance_coph = ssim_first_direct[[14]]
contrast_coph = ssim_first_direct[[15]]
structure_coph = ssim_first_direct[[16]]
SSIM_diff_orig_coph = ssim_first_direct[[17]]
SSIM_diff_CAST_J129 = ssim_first_direct[[18]]

boxplot(list(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph),
        ylab = "Structural similarity index SSIM",
        main = "",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        boxcol="grey75",
        medcol="grey75",
        whiskcol="grey75",
        staplecol="grey75",
        outcol="grey75",
        outbg="grey75",
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)

points(x=jitter(rep(1:4,each=length(SSIM_CAST))),
     y=cbind(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph),
     pch=19)




# 4 plots for SSIM and the 3 contributers luminance, contrast and stucture
par(mfrow = c(2,2))
### SSIM
boxplot(list(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph),#,SSIM_diff_orig_coph,SSIM_diff_CAST_J129),
        ylab = "Structural similarity index SSIM",
        main = "",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        boxcol="grey75",
        medcol="grey75",
        whiskcol="grey75",
        staplecol="grey75",
        outcol="grey75",
        outbg="grey75",
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)

points(x=jitter(rep(1:4,each=length(SSIM_CAST))),
       y=cbind(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph),
       pch=19)

### luminance 
boxplot(list(luminance_CAST,luminance_J129,luminance_orig,luminance_coph),
        ylab = "Structural similarity index luminance",
        main = "",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        boxcol="grey75",
        medcol="grey75",
        whiskcol="grey75",
        staplecol="grey75",
        outcol="grey75",
        outbg="grey75",
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)

points(x=jitter(rep(1:4,each=length(luminance_CAST))),
       y=cbind(luminance_CAST,luminance_J129,luminance_orig,luminance_coph),
       pch=19)

### contrast
boxplot(list(contrast_CAST,contrast_J129,contrast_orig,contrast_coph),
        ylab = "Structural similarity index contrast",
        main = "",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        boxcol="grey75",
        medcol="grey75",
        whiskcol="grey75",
        staplecol="grey75",
        outcol="grey75",
        outbg="grey75",
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)

points(x=jitter(rep(1:4,each=length(contrast_CAST))),
       y=cbind(contrast_CAST,contrast_J129,contrast_orig,contrast_coph),
       pch=19)

### structure
boxplot(list(structure_CAST,structure_J129,structure_orig,structure_coph),
        ylab = "Structural similarity index structure",
        main = "",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        boxcol="grey75",
        medcol="grey75",
        whiskcol="grey75",
        staplecol="grey75",
        outcol="grey75",
        outbg="grey75",
        cex.lab=1.3,
        cex.axis=1.3)

points(x=jitter(rep(1:4,each=length(structure_CAST))),
       y=cbind(structure_CAST,structure_J129,structure_orig,structure_coph),
       pch=19)

# difference between original and cophased SSIM and CAST and J129 SSIM
# SSIM_orig - SSIM_coph
# SSIM_CAST - SSIM_J129

par(mfrow=c(1,1))
boxplot(list(SSIM_diff_orig_coph,SSIM_diff_CAST_J129),
        ylab = "Difference in SSIM",
        main = "",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        #ylim = c(0,1),
        boxcol="grey75",
        medcol="grey75",
        whiskcol="grey75",
        staplecol="grey75",
        outcol="grey75",
        outbg="grey75",
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:2,labels = c("SSIM\norig - cophased","SSIM\nCAST - S129"), tick = FALSE, cex.axis=1)

points(x=jitter(rep(1:2,each=length(SSIM_diff_orig_coph))),
       y=cbind(SSIM_diff_orig_coph,SSIM_diff_CAST_J129),
       pch=19)

################################################################################
################################################################################
################################################################################
# calc a baseline aka how similar are contact matrices from different chromosomes
# this seems to go into memory issues quite quickly, so let's just compare CAST chr 19 vs all others for now

SSIM(CAST_coph[[1]], J129_coph[[2]], L = 1) 
# works, all chr vs all chr is an additional 171 comparisons, times 8?
# this always crashes at 40 something bc of memory
# also, with different sizes, it keeps changimg. So comparing same size matrices, the results are reproduceable, btu with different sizes, they are not, they keep changing
# # CAST orig vs cophased
# SSIM_CAST_all_chr = vector(mode="numeric",length=171)
# # J129 orig vs cophased
# SSIM_J129_all_chr = vector(mode="numeric",length=171)
# # CAST vs J129 orig
# SSIM_orig_all_chr = vector(mode="numeric",length=171)
# # CAST vs J129 cophased
# SSIM_coph_all_chr = vector(mode="numeric",length=171)
# 
# CAST orig
SSIM_CAST_orig_19_vs_all_chr = vector(mode="numeric",length=18)
# J129 orig
SSIM_J129_orig_19_vs_all_chr = vector(mode="numeric",length=18)
# CAST cophased
SSIM_CAST_coph_19_vs_all_chr = vector(mode="numeric",length=18)
# J129 cophased
SSIM_J129_coph_19_vs_all_chr = vector(mode="numeric",length=18)


for (i in 1:18){
    print(i)
    # SSIM_CAST_all_chr[i] = SSIM(CAST_orig[[j]], CAST_coph[[k]], L = 1)
    # 
    # SSIM_J129_all_chr[i] = SSIM(J129_orig[[j]], J129_coph[[k]], L = 1)
    # 
    # SSIM_orig_all_chr[i] = SSIM(CAST_orig[[j]], J129_orig[[k]], L = 1)
    # 
    # SSIM_coph_all_chr[i] = SSIM(CAST_coph[[j]], J129_coph[[k]], L = 1)
    
    SSIM_CAST_orig_19_vs_all_chr[i] = SSIM(CAST_orig[[i]], CAST_orig[[19]], L = 1)
    
    SSIM_J129_orig_19_vs_all_chr[i] = SSIM(J129_orig[[i]], J129_orig[[19]], L = 1)
    
    SSIM_CAST_coph_19_vs_all_chr[i] = SSIM(CAST_coph[[i]], CAST_coph[[19]], L = 1)
    
    SSIM_J129_coph_19_vs_all_chr[i] = SSIM(J129_coph[[i]], J129_coph[[19]], L = 1)
}
# failed as well

# SSIM(CAST_orig[[18]], CAST_orig[[19]], L = 1)
# $SSIM
# [1] 0.1312057
# $comps
# luminance  contrast structure
# 0.9965020 0.9968348 0.1320844
# 
# SSIM(CAST_orig[[17]], CAST_orig[[18]], L = 1)
# $SSIM
# [1] 0.07462757
# $comps
# luminance   contrast  structure
# 0.97887150 0.99984542 0.07625016

# SSIM(CAST_orig[[16]], CAST_orig[[17]], L = 1)
# $SSIM
# [1] 0.07717305
# $comps
# luminance   contrast  structure
# 0.96710304 0.99994651 0.07980244

# SSIM(CAST_orig[[16]], CAST_orig[[15]], L = 1)
# $SSIM
# [1] 0.0798612
# $comps
# luminance   contrast  structure
# 0.99328824 0.99997883 0.08040254

# makes no sense so far, results are constantly low but not reproduceable with different size matrices

################################################################################
################################################################################
################################################################################
# sliding window along diagonal, take mean
# compare with top10 most common and top5 most differential contacts

# take windows of 10 MB? this is 50kb resolution, so 200 bins make 10 MB
# dim(CAST_orig[[19]])
# [1] 1229 1229

windowsize = 200
windowamount = nrow(CAST_orig[[19]]) - windowsize

ssim_sliding = vector(mode="numeric",length=windowamount)
for (i in  1:windowamount){
    if(i%%100 == 0){print(i)}
        ssim_sliding[i] =  SSIM(CAST_orig[[19]][i:(i+windowsize),i:(i+windowsize)], CAST_coph[[19]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
}

# > mean(ssim_sliding)
# [1] 0.8094989
# > SSIM(CAST_orig[[19]], CAST_coph[[19]], L = 1)$SSIM
# [1] 0.7722388


ssim_sliding_CAST = vector(mode="list",length=19)
ssim_sliding_J129 = vector(mode="list",length=19)
ssim_sliding_orig = vector(mode="list",length=19)
ssim_sliding_coph = vector(mode="list",length=19)

windowsize = 200

for (chr in 1:19){
    windowamount = nrow(CAST_orig[[chr]]) - windowsize
    
    ssim_sliding_CAST[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_J129[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_orig[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_coph[[chr]] = vector(mode="numeric",length=windowamount)
    
    for (i in  1:windowamount){
        if(i%%100 == 0){print(paste(chr,i))}
        ssim_sliding_CAST[[chr]][i] =  SSIM(CAST_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], CAST_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_J129[[chr]][i] =  SSIM(J129_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_orig[[chr]][i] =  SSIM(CAST_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_coph[[chr]][i] =  SSIM(CAST_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
    }
}

saveRDS(list(ssim_sliding_CAST,
             ssim_sliding_J129,
             ssim_sliding_orig,
             ssim_sliding_coph),
        paste(dir_matrices,"SSIM_chr_sliding_200bins_10mb.rds",sep="/"))

ssim_sliding_CAST_10mb = ssim_sliding_CAST
ssim_sliding_J129_10mb = ssim_sliding_J129
ssim_sliding_orig_10mb = ssim_sliding_orig
ssim_sliding_coph_10mb = ssim_sliding_coph

# again with window size 1Mb, aka 20 windows (maybe too small)

ssim_sliding_CAST = vector(mode="list",length=19)
ssim_sliding_J129 = vector(mode="list",length=19)
ssim_sliding_orig = vector(mode="list",length=19)
ssim_sliding_coph = vector(mode="list",length=19)

windowsize = 20

for (chr in 1:19){
    windowamount = nrow(CAST_orig[[chr]]) - windowsize
    
    ssim_sliding_CAST[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_J129[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_orig[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_coph[[chr]] = vector(mode="numeric",length=windowamount)
    
    for (i in  1:windowamount){
        if(i%%100 == 0){print(paste(chr,i))}
        ssim_sliding_CAST[[chr]][i] =  SSIM(CAST_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], CAST_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_J129[[chr]][i] =  SSIM(J129_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_orig[[chr]][i] =  SSIM(CAST_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_coph[[chr]][i] =  SSIM(CAST_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
    }
}

saveRDS(list(ssim_sliding_CAST,
             ssim_sliding_J129,
             ssim_sliding_orig,
             ssim_sliding_coph),
        paste(dir_matrices,"SSIM_chr_sliding_20bins_1mb.rds",sep="/"))

ssim_sliding_CAST_1mb = ssim_sliding_CAST
ssim_sliding_J129_1mb = ssim_sliding_J129
ssim_sliding_orig_1mb = ssim_sliding_orig
ssim_sliding_coph_1mb = ssim_sliding_coph

# again with the border to NA, so exactly the area around the diagonal that is actually calculated
ssim_sliding_CAST = vector(mode="list",length=19)
ssim_sliding_J129 = vector(mode="list",length=19)
ssim_sliding_orig = vector(mode="list",length=19)
ssim_sliding_coph = vector(mode="list",length=19)

windowsize = 60

for (chr in 1:19){
    windowamount = nrow(CAST_orig[[chr]]) - windowsize
    
    ssim_sliding_CAST[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_J129[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_orig[[chr]] = vector(mode="numeric",length=windowamount)
    ssim_sliding_coph[[chr]] = vector(mode="numeric",length=windowamount)
    
    for (i in  1:windowamount){
        if(i%%100 == 0){print(paste(chr,i))}
        ssim_sliding_CAST[[chr]][i] =  SSIM(CAST_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], CAST_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_J129[[chr]][i] =  SSIM(J129_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_orig[[chr]][i] =  SSIM(CAST_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_orig[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        ssim_sliding_coph[[chr]][i] =  SSIM(CAST_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], J129_coph[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
    }
}

saveRDS(list(ssim_sliding_CAST,
             ssim_sliding_J129,
             ssim_sliding_orig,
             ssim_sliding_coph),
        paste(dir_matrices,"SSIM_chr_sliding_60bins_3mb.rds",sep="/"))

ssim_sliding_CAST_3mb = ssim_sliding_CAST
ssim_sliding_J129_3mb = ssim_sliding_J129
ssim_sliding_orig_3mb = ssim_sliding_orig
ssim_sliding_coph_3mb = ssim_sliding_coph


# get the mean for all of those
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_200bins_10mb.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_20bins_1mb.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_60bins_3mb.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/


ssim_sliding_10mb = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_200bins_10mb.rds")
ssim_sliding_CAST_10mb = ssim_sliding_10mb[[1]]
ssim_sliding_J129_10mb = ssim_sliding_10mb[[2]]
ssim_sliding_orig_10mb = ssim_sliding_10mb[[3]]
ssim_sliding_coph_10mb = ssim_sliding_10mb[[4]]

ssim_sliding_1mb = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_20bins_1mb.rds")
ssim_sliding_CAST_1mb = ssim_sliding_1mb[[1]]
ssim_sliding_J129_1mb = ssim_sliding_1mb[[2]]
ssim_sliding_orig_1mb = ssim_sliding_1mb[[3]]
ssim_sliding_coph_1mb = ssim_sliding_1mb[[4]]

ssim_sliding_3mb = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_60bins_3mb.rds")
ssim_sliding_CAST_3mb = ssim_sliding_3mb[[1]]
ssim_sliding_J129_3mb = ssim_sliding_3mb[[2]]
ssim_sliding_orig_3mb = ssim_sliding_3mb[[3]]
ssim_sliding_coph_3mb = ssim_sliding_3mb[[4]]

boxplot(list(ssim_sliding_CAST_10mb[[19]],
             ssim_sliding_CAST_3mb[[19]],
             ssim_sliding_CAST_1mb[[19]]))

par(mfrow=c(3,1))
boxplot(ssim_sliding_CAST_10mb)
abline(h=mean(SSIM_CAST),col="red") # 0.7750342
abline(h=mean(unlist(ssim_sliding_CAST_10mb)),col="cyan") # 0.8263076
boxplot(ssim_sliding_CAST_3mb)
abline(h=mean(SSIM_CAST),col="red")
abline(h=mean(unlist(ssim_sliding_CAST_3mb)),col="cyan") # 0.835688
boxplot(ssim_sliding_CAST_1mb)
abline(h=mean(SSIM_CAST),col="red")
abline(h=mean(unlist(ssim_sliding_CAST_1mb)),col="cyan") # 0.8406515

par(mfrow=c(2,2))
boxplot(list(SSIM_CAST,unlist(ssim_sliding_CAST_10mb),unlist(ssim_sliding_CAST_3mb),unlist(ssim_sliding_CAST_1mb)),
        ylab = "Structural similarity index SSIM",
        main = "CAST: original vs cophased",
        xaxt = "n",
        xlab = "Window size",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("Chromosome","10 Mb","3 Mb","1 Mb"), tick = FALSE, cex.axis=1)
boxplot(list(SSIM_J129,unlist(ssim_sliding_J129_10mb),unlist(ssim_sliding_J129_3mb),unlist(ssim_sliding_J129_1mb)),
        ylab = "Structural similarity index SSIM",
        main = "S129: original vs cophased",
        xaxt = "n",
        xlab = "Window size",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("Chromosome","10 Mb","3 Mb","1 Mb"), tick = FALSE, cex.axis=1)
boxplot(list(SSIM_orig,unlist(ssim_sliding_orig_10mb),unlist(ssim_sliding_orig_3mb),unlist(ssim_sliding_orig_1mb)),
        ylab = "Structural similarity index SSIM",
        main = "original: CAST vs S129",
        xaxt = "n",
        xlab = "Window size",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("Chromosome","10 Mb","3 Mb","1 Mb"), tick = FALSE, cex.axis=1)
boxplot(list(SSIM_coph,unlist(ssim_sliding_coph_10mb),unlist(ssim_sliding_coph_3mb),unlist(ssim_sliding_coph_1mb)),
        ylab = "Structural similarity index SSIM",
        main = "cophased: CAST vs S129",
        xaxt = "n",
        xlab = "Window size",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3)
axis(1,at = 1:4,labels = c("Chromosome","10 Mb","3 Mb","1 Mb"), tick = FALSE, cex.axis=1)

### window sizes next to each other
par(mfrow=c(2,2))
boxplot(list(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph),
        ylab = "Structural similarity index SSIM",
        main = "Chromosome",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3,
        notch=T)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)
ar(mfrow=c(1,1))
boxplot(list(unlist(ssim_sliding_CAST_10mb),unlist(ssim_sliding_J129_10mb),unlist(ssim_sliding_orig_10mb),unlist(ssim_sliding_coph_10mb)),
        ylab = "Structural similarity index SSIM",
        main = "10 Mb windows",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3,
        notch=T)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)
boxplot(list(unlist(ssim_sliding_CAST_3mb),unlist(ssim_sliding_J129_3mb),unlist(ssim_sliding_orig_3mb),unlist(ssim_sliding_coph_3mb)),
        ylab = "Structural similarity index SSIM",
        main = "3 Mb windows",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3,
        notch=T)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)
boxplot(list(unlist(ssim_sliding_CAST_1mb),unlist(ssim_sliding_J129_1mb),unlist(ssim_sliding_orig_1mb),unlist(ssim_sliding_coph_1mb)),
        ylab = "Structural similarity index SSIM",
        main = "1 Mb windows",
        xaxt = "n",
        xlab = "Comparison of originally phased and cophased contact matrices",
        ylim = c(0,1),
        cex.lab=1.3,
        cex.axis=1.3,
        notch=T)
axis(1,at = 1:4,labels = c("CAST\norig vs cophased","S129\norig vs cophased","original\nCAST vs S129","cophased\nCAST vs S129"), tick = FALSE, cex.axis=1)

# where is it that they don't match up
par(mfrow=c(4,1))
plot(y = ssim_sliding_CAST_10mb[[19]],
     x = c(1:length(ssim_sliding_CAST_10mb[[19]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (10 Mb)",
     ylab = "SSIM")
plot(y = ssim_sliding_CAST_1mb[[19]],
     x = c(1:length(ssim_sliding_CAST_1mb[[19]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "SSIM")
plot(y = ssim_sliding_J129_10mb[[19]],
     x = c(1:length(ssim_sliding_CAST_10mb[[19]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (10 Mb)",
     ylab = "SSIM")
plot(y = ssim_sliding_J129_1mb[[19]],
     x = c(1:length(ssim_sliding_CAST_1mb[[19]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "SSIM")
# looks promising
# try to plot a heatmap with it
# mount cluster
# sshfs -o allow_other,defer_permissions,follow_symlinks jmarkow_m@bihtext:/ mnt/bih_cluster
dir_mount="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices"

# CAST_orig_19 = data.matrix(read.table(gzfile(paste(dir_mount,"/CAST.original.res_50K/segtable.res_50000.CAST.NPMI.chr19.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
# CAST_coph_19 = data.matrix(read.table(gzfile(paste(dir_mount,"/CAST.dist_cutoff_100K.res_50K/segtable.dist_cutoff_100000.res_50000.CAST.iteration_1.NPMI.chr19.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
## J129_orig_19 = data.matrix(read.table(gzfile(paste(dir_mount,"/S129.original.res_50K/segtable.res_50000.S129.NPMI.chr19.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
## J129_coph_19 = data.matrix(read.table(gzfile(paste(dir_mount,"/S129.dist_cutoff_100K.res_50K/segtable.dist_cutoff_100000.res_50000.S129.iteration_1.NPMI.chr19.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))


CAST_orig_19_red = CAST_orig_19[1:100,1:100]
CAST_orig_19_red[is.na(CAST_orig_19_red)] <- 0
CAST_orig_19_NA = CAST_orig_19 
CAST_orig_19_NA[is.na(CAST_orig_19_NA)] <- 0

CAST_coph_19_red = CAST_coph_19[1:100,1:100]
CAST_coph_19_red[is.na(CAST_coph_19_red)] <- 0
CAST_coph_19_NA = CAST_coph_19 
CAST_coph_19_NA[is.na(CAST_coph_19_NA)] <- 0

library(RColorBrewer)
pal <- colorRampPalette(c("white", "black"))
npmi_col_base = brewer.pal(n = 11, name = "RdYlBu")
npmi_col = colorRampPalette(c(rev(npmi_col_base)))(100) 

heatmap(CAST_orig_19_red,col = pal(20))
heatmap(CAST_orig_19_NA[1:300,1:300],Colv = NA, Rowv = NA, col = pal(20)) 
heatmap(CAST_orig_19_NA[1:500,1:500],Colv = NA, Rowv = NA, col = npmi_col) 
heatmap(CAST_orig_19_NA[500:900,500:900],Colv = NA, Rowv = NA, col = npmi_col) 
heatmap(CAST_coph_19_NA[500:900,500:900],Colv = NA, Rowv = NA, col = npmi_col) 

heatmap(CAST_orig_19_NA[400:600,400:600],Colv = NA, Rowv = NA, col = npmi_col) 
heatmap(CAST_coph_19_NA[400:600,400:600],Colv = NA, Rowv = NA, col = npmi_col) 

heatmap(CAST_orig_19_NA[900:1000,900:1000],Colv = NA, Rowv = NA, col = npmi_col) 
heatmap(CAST_coph_19_NA[900:1000,900:1000],Colv = NA, Rowv = NA, col = npmi_col) 

#heatmap(CAST_orig_19_NA,Colv = NA, Rowv = NA, col = pal(20)) 
################################################################################
# are those siginificant?
# first of all, are they normally distributed? (p.value > 0.05 means normal distribution)
shapiro.test(SSIM_CAST) # yes
qqnorm(SSIM_CAST) # does not look too good
hist(SSIM_CAST)
shapiro.test(SSIM_J129) # yes
shapiro.test(SSIM_orig) # yes
shapiro.test(SSIM_coph) # yes

shapiro.test(ssim_sliding_CAST_10mb[[19]]) # p-value < 2.2e-16, however, All normality tests are too sensitive to sample size. 
qqnorm(ssim_sliding_CAST_10mb[[19]])  # does not look too good
hist(ssim_sliding_CAST_10mb[[19]])
shapiro.test(ssim_sliding_J129_10mb[[19]]) # p-value < 2.2e-16
shapiro.test(ssim_sliding_orig_10mb[[19]]) # p-value < 2.2e-16
shapiro.test(ssim_sliding_coph_10mb[[19]]) # p-value < 2.2e-16

shapiro.test(ssim_sliding_CAST_1mb[[19]]) # p-value < 2.2e-16
qqnorm(ssim_sliding_CAST_1mb[[19]])  # does not look too good
shapiro.test(ssim_sliding_J129_1mb[[19]]) # p-value < 2.2e-16
shapiro.test(ssim_sliding_orig_1mb[[19]]) # p-value < 2.2e-16
shapiro.test(ssim_sliding_coph_1mb[[19]]) # p-value < 2.2e-16


pairwise.t.test(c(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph), rep(1:4,each=length(SSIM_CAST)))
#   1      2      3   
# 2 0.86   -      -   
# 3 <2e-16 <2e-16 -   
# 4 <2e-16 <2e-16 0.86
# cool, 1 vs 3 and 4 is sig, 2 vs 3 and 4 as well

# on the cluster:
pairwise.t.test(c(unlist(ssim_sliding_CAST_10mb),unlist(ssim_sliding_J129_10mb),unlist(ssim_sliding_orig_10mb),unlist(ssim_sliding_coph_10mb)), rep(1:4,each=length(unlist(ssim_sliding_CAST_10mb)))) # p.adjust.method do not change anything
#   1      2      3
# 2 <2e-16 -      -
# 3 <2e-16 <2e-16 -
# 4 <2e-16 <2e-16 <2e-16

pairwise.t.test(c(unlist(ssim_sliding_CAST_3mb),unlist(ssim_sliding_J129_3mb),unlist(ssim_sliding_orig_3mb),unlist(ssim_sliding_coph_3mb)), rep(1:4,each=length(unlist(ssim_sliding_CAST_3mb))))
#   1      2      3
# 2 1e-13  -      -
# 3 <2e-16 <2e-16 -
# 4 <2e-16 <2e-16 <2e-16

pairwise.t.test(c(unlist(ssim_sliding_CAST_1mb),unlist(ssim_sliding_J129_1mb),unlist(ssim_sliding_orig_1mb),unlist(ssim_sliding_coph_1mb)), rep(1:4,each=length(unlist(ssim_sliding_CAST_1mb))),paired = T)
#   1      2      3     
# 2 <2e-16 -      -     
# 3 <2e-16 <2e-16 -     
# 4 <2e-16 <2e-16 <2e-16
pairwise.t.test(c(unlist(ssim_sliding_CAST_1mb),unlist(ssim_sliding_J129_1mb),unlist(ssim_sliding_orig_1mb),unlist(ssim_sliding_coph_1mb)), rep(1:4,each=length(unlist(ssim_sliding_CAST_1mb))),paired = F)
#   1       2       3
# 2 1.5e-07 -       -
# 3 < 2e-16 < 2e-16 -
# 4 < 2e-16 < 2e-16 < 2e-16
# 
# P value adjustment method: holm

t.test(c(SSIM_CAST,SSIM_J129),rep(1:2,each=length(unlist(ssim_sliding_CAST_10mb))))
# Welch Two Sample t-test
# 
# data:  c(SSIM_CAST, SSIM_J129) and rep(1:2, each = length(unlist(ssim_sliding_CAST_10mb)))
# t = -144.03, df = 46.432, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#     -0.7391469 -0.7187766
# sample estimates:
#     mean of x mean of y 
# 0.7710383 1.5000000 

# all of this does not make much sense




one.way <- aov(c(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph) ~ rep(1:4,each=length(SSIM_CAST)))

summary(one.way)
# does not make much sense either 

ssim_chrom = data.frame(cbind(c(SSIM_CAST,SSIM_J129,SSIM_orig,SSIM_coph), rep(c("C","J","O","P"),each=length(SSIM_CAST))),stringsAsFactors = T)
colnames(ssim_chrom) = c("SSIM","comparison")
class(ssim_chrom$SSIM) <- "numeric"

ssim_10mb = data.frame(cbind(c(unlist(ssim_sliding_CAST_10mb),unlist(ssim_sliding_J129_10mb),unlist(ssim_sliding_orig_10mb),unlist(ssim_sliding_coph_10mb)), rep(c("C","J","O","P"),each=length(unlist(ssim_sliding_CAST_10mb)))),stringsAsFactors = T)
colnames(ssim_10mb) = c("SSIM","comparison")
class(ssim_10mb$SSIM) <- "numeric"

ssim_3mb = data.frame(cbind(c(unlist(ssim_sliding_CAST_3mb),unlist(ssim_sliding_J129_3mb),unlist(ssim_sliding_orig_3mb),unlist(ssim_sliding_coph_3mb)), rep(c("C","J","O","P"),each=length(unlist(ssim_sliding_CAST_3mb)))),stringsAsFactors = T)
colnames(ssim_3mb) = c("SSIM","comparison")
class(ssim_3mb$SSIM) <- "numeric"

ssim_1mb = data.frame(cbind(c(unlist(ssim_sliding_CAST_1mb),unlist(ssim_sliding_J129_1mb),unlist(ssim_sliding_orig_1mb),unlist(ssim_sliding_coph_1mb)), rep(c("C","J","O","P"),each=length(unlist(ssim_sliding_CAST_1mb)))),stringsAsFactors = T)
colnames(ssim_1mb) = c("SSIM","comparison")
class(ssim_1mb$SSIM) <- "numeric"

# Kruskal-Wallis rank sum test
kruskal.test(SSIM ~ comparison, data = ssim_chrom)
# Kruskal-Wallis rank sum test
# data:  SSIM by comparison
# Kruskal-Wallis chi-squared = 56.835, df = 3, p-value = 2.787e-12

kruskal.test(SSIM ~ comparison, data = ssim_10mb)
# Kruskal-Wallis chi-squared = 105500, df = 3, p-value < 2.2e-16

kruskal.test(SSIM ~ comparison, data = ssim_3mb)
# Kruskal-Wallis chi-squared = 68729, df = 3, p-value < 2.2e-16

kruskal.test(SSIM ~ comparison, data = ssim_1mb)
# Kruskal-Wallis chi-squared = 54582, df = 3, p-value < 2.2e-16


pairwise.wilcox.test(ssim_chrom$SSIM, ssim_chrom$comparison,p.adjust.method = "BH")
# Pairwise comparisons using Wilcoxon rank sum test 
# data:  ssim_chrom$SSIM and ssim_chrom$comparison 
#   C       J       O   
# J 0.31    -       -   
# O 8.5e-11 8.5e-11 -   
# P 8.5e-11 8.5e-11 0.31

pairwise.wilcox.test(ssim_10mb$SSIM, ssim_10mb$comparison,p.adjust.method = "BH")
#   C      J      O     
# J <2e-16 -      -     
# O <2e-16 <2e-16 -     
# P <2e-16 <2e-16 <2e-16

pairwise.wilcox.test(ssim_3mb$SSIM, ssim_3mb$comparison,p.adjust.method = "BH")
#   C       J       O      
# J 3.3e-11 -       -      
# O < 2e-16 < 2e-16 -      
# P < 2e-16 < 2e-16 < 2e-16

pairwise.wilcox.test(ssim_1mb$SSIM, ssim_1mb$comparison,p.adjust.method = "BH")
#   C       J       O      
# J 6.6e-05 -       -      
# O < 2e-16 < 2e-16 -      
# P < 2e-16 < 2e-16 < 2e-16

# checked also with notches, talked to Maja and Stella, both do not know what to do (and if there needs to be something done) when handling so many points, as everything seems to be significant with that pile of data :P

################################################################################
################################################################################
################################################################################
# calculate SSIM for haplotype-specific ‚H1 maps as well, correlate with most common and most differential 

# load all the H1 NMPI matrices
# mount cluster
# sshfs -o allow_other,defer_permissions,follow_symlinks jmarkow_m@bihtext:/ mnt/bih_cluster
# dir_mount_H1="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/H1_phased/Curated/at50Kb"
dir_H1="/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/H1_phased/Curated/at50Kb"

NPMI_1 = vector(mode="list",length=22)
NPMI_2 = vector(mode="list",length=22)
zscores_common = vector(mode="list",length=22)
zscores_diff =vector(mode="list",length=22)
zscores_diff_1 = vector(mode="list",length=22)
zscores_diff_2 = vector(mode="list",length=22)

for (chr in 1:22){
    # chr = 6
    print(chr)
    #chr=1
    NPMI_1[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/NPMI_Genome1/H1.Genome1_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    NPMI_2[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/NPMI_Genome2/H1.Genome2_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    
    zscores_common[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_common/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_diff/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_1[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_2[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
}

# replace NA with 0
NPMI_1 = lapply(NPMI_1,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
NPMI_2 = lapply(NPMI_2,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
zscores_common = lapply(zscores_common,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
zscores_diff = lapply(zscores_diff,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
zscores_diff_1 = lapply(zscores_diff_1,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})
zscores_diff_2 = lapply(zscores_diff_2,function(chrom){
    chrom[is.na(chrom)] <- 0
    return(chrom)})

# so I want sliding windows:
    # SSIM NPMI_1 and NPMI_2
    # sum abs(most common)
    # sum abs(most differential)

ssim_sliding_NPMI = vector(mode="list",length=22)
sum_sliding_zmostdiff1 = vector(mode="list",length=22)
sum_sliding_zmostdiff2 = vector(mode="list",length=22)
sum_sliding_zdiff = vector(mode="list",length=22)
sum_sliding_zcom = vector(mode="list",length=22)

windowsize = 200

for (chr in 1:22){
    windowamount = nrow(NPMI_1[[chr]]) - windowsize
    
    ssim_sliding_NPMI[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff1[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff2[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zdiff[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zcom[[chr]] = vector(mode="numeric",length=windowamount)
    
    for (i in  1:windowamount){
        if(i%%100 == 0){print(paste(chr,i))}
        ssim_sliding_NPMI[[chr]][i] =  SSIM(NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)], NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        sum_sliding_zmostdiff1[[chr]][i] =  sum(abs(zscores_diff_1[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zmostdiff2[[chr]][i] =  sum(abs(zscores_diff_2[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zdiff[chr]][i] =  sum(abs(zscores_diff[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zcom[[chr]][i] =  sum(abs(zscores_common[[chr]][i:(i+windowsize),i:(i+windowsize)]))
    }
}

saveRDS(list(ssim_sliding_NPMI,
             sum_sliding_zmostdiff1,
             sum_sliding_zmostdiff2,
             sum_sliding_zdiff
             sum_sliding_zcom),
        paste(dir_matrices,"SSIM_chr_sliding_200bins_10mb_H1_zscore_diff_common_sum.rds",sep="/"))

ssim_sliding_NPMI_10mb = ssim_sliding_NPMI
sum_sliding_zmostdiff1_10mb = sum_sliding_zmostdiff1
sum_sliding_zmostdiff2_10mb = sum_sliding_zmostdiff2
sum_sliding_zdiff10mb = sum_sliding_zmostdiff
sum_sliding_zcom_10mb = sum_sliding_zcom

####

ssim_sliding_NPMI = vector(mode="list",length=22)
sum_sliding_zmostdiff1 = vector(mode="list",length=22)
sum_sliding_zmostdiff2 = vector(mode="list",length=22)
sum_sliding_zmostdiff = vector(mode="list",length=22)
sum_sliding_zcom = vector(mode="list",length=22)

windowsize = 20

for (chr in 1:22){
    windowamount = nrow(NPMI_1[[chr]]) - windowsize
    
    ssim_sliding_NPMI[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff1[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff2[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zdiff[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zcom[[chr]] = vector(mode="numeric",length=windowamount)
    
    for (i in  1:windowamount){
        if(i%%100 == 0){print(paste(chr,i))}
        ssim_sliding_NPMI[[chr]][i] =  SSIM(NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)], NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        sum_sliding_zmostdiff1[[chr]][i] =  sum(abs(zscores_diff_1[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zmostdiff2[[chr]][i] =  sum(abs(zscores_diff_2[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zdiff[chr]][i] =  sum(abs(zscores_diff[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zcom[[chr]][i] =  sum(abs(zscores_common[[chr]][i:(i+windowsize),i:(i+windowsize)]))
    }
}

saveRDS(list(ssim_sliding_NPMI,
             sum_sliding_zmostdiff1,
             sum_sliding_zmostdiff2,
             sum_sliding_zdiff
             sum_sliding_zcom),
        paste(dir_matrices,"SSIM_chr_sliding_20bins_1mb_H1_zscore_diff_common_sum.rds",sep="/"))

ssim_sliding_NPMI_1mb = ssim_sliding_NPMI
sum_sliding_zmostdiff1_1mb = sum_sliding_zmostdiff1
sum_sliding_zmostdiff2_1mb = sum_sliding_zmostdiff2
sum_sliding_zdiff1mb = sum_sliding_zmostdiff
sum_sliding_zcom_1mb = sum_sliding_zcom

####

ssim_sliding_NPMI = vector(mode="list",length=22)
sum_sliding_zmostdiff1 = vector(mode="list",length=22)
sum_sliding_zmostdiff2 = vector(mode="list",length=22)
sum_sliding_zmostdiff = vector(mode="list",length=22)
sum_sliding_zcom = vector(mode="list",length=22)


windowsize = 60

for (chr in 1:22){
    windowamount = nrow(NPMI_1[[chr]]) - windowsize
    
    ssim_sliding_NPMI[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff1[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff2[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zdiff[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zcom[[chr]] = vector(mode="numeric",length=windowamount)
    
    for (i in  1:windowamount){
        if(i%%100 == 0){print(paste(chr,i))}
        ssim_sliding_NPMI[[chr]][i] =  SSIM(NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)], NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)], L = 1)$SSIM
        sum_sliding_zmostdiff1[[chr]][i] =  sum(abs(zscores_diff_1[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zmostdiff2[[chr]][i] =  sum(abs(zscores_diff_2[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zdiff[chr]][i] =  sum(abs(zscores_diff[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        sum_sliding_zcom[[chr]][i] =  sum(abs(zscores_common[[chr]][i:(i+windowsize),i:(i+windowsize)]))
    }
}

saveRDS(list(ssim_sliding_NPMI,
             sum_sliding_zmostdiff1,
             sum_sliding_zmostdiff2,
             sum_sliding_zdiff
             sum_sliding_zcom),
        paste(dir_matrices,"SSIM_chr_sliding_60bins_3mb_H1_zscore_diff_common_sum.rds",sep="/"))

ssim_sliding_NPMI_3mb = ssim_sliding_NPMI
sum_sliding_zmostdiff1_3mb = sum_sliding_zmostdiff1
sum_sliding_zmostdiff2_3mb = sum_sliding_zmostdiff2
sum_sliding_zdiff3mb = sum_sliding_zmostdiff
sum_sliding_zcom_3mb = sum_sliding_zcom

# prepare plots
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_200bins_10mb_H1_zscore_diff_common_sum.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_20bins_1mb_H1_zscore_diff_common_sum.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_60bins_3mb_H1_zscore_diff_common_sum.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/

ssim_sliding_10mb = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_200bins_10mb_H1_zscore_diff_common_sum.rds")
ssim_sliding_NPMI_10mb = ssim_sliding_10mb[[1]]
sum_sliding_zmostdiff1_10mb = ssim_sliding_10mb[[2]]
sum_sliding_zmostdiff2_10mb = ssim_sliding_10mb[[3]]
sum_sliding_zdiff10mb = ssim_sliding_10mb[[4]]
sum_sliding_zcom_10mb = ssim_sliding_10mb[[5]]

ssim_sliding_1mb = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_20bins_1mb_H1_zscore_diff_common_sum.rds")
ssim_sliding_NPMI_1mb = ssim_sliding_1mb[[1]]
sum_sliding_zmostdiff1_1mb = ssim_sliding_1mb[[2]]
sum_sliding_zmostdiff2_1mb = ssim_sliding_1mb[[3]]
sum_sliding_zdiff1mb = ssim_sliding_1mb[[4]]
sum_sliding_zcom_1mb = ssim_sliding_1mb[[5]]

ssim_sliding_3mb = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_60bins_3mb_H1_zscore_diff_common_sum.rds")
ssim_sliding_NPMI_3mb = ssim_sliding_3mb[[1]]
sum_sliding_zmostdiff1_3mb = ssim_sliding_3mb[[2]]
sum_sliding_zmostdiff2_3mb = ssim_sliding_3mb[[3]]
sum_sliding_zdiff3mb = ssim_sliding_3mb[[4]]
sum_sliding_zcom_3mb = ssim_sliding_3mb[[5]]

par(mfrow=c(5,1))
plot(y = ssim_sliding_NPMI_10mb[[22]],
     x = c(1:length(ssim_sliding_NPMI_10mb[[22]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (10 Mb)",
     ylab = "SSIM")
plot(y = sum_sliding_zmostdiff1_10mb[[22]],
     x = c(1:length(sum_sliding_zmostdiff1_10mb[[22]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "zscore diff 1")
plot(y = sum_sliding_zmostdiff2_10mb[[22]],
     x = c(1:length(sum_sliding_zmostdiff2_10mb[[22]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "zscore diff 2")
plot(y = sum_sliding_zdiff10mb[[22]],
     x = c(1:length(sum_sliding_zdiff10mb[[22]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "zscore most diff")
plot(y = sum_sliding_zcom_10mb[[22]],
     x = c(1:length(sum_sliding_zcom_10mb[[22]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (10 Mb)",
     ylab = "zscore most common")


par(mfrow=c(5,1))
plot(y = ssim_sliding_NPMI_10mb[[1]],
     x = c(1:length(ssim_sliding_NPMI_10mb[[1]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (10 Mb)",
     ylab = "SSIM")
plot(y = sum_sliding_zmostdiff1_10mb[[1]],
     x = c(1:length(sum_sliding_zmostdiff1_10mb[[1]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "zscore diff 1")
plot(y = sum_sliding_zmostdiff2_10mb[[1]],
     x = c(1:length(sum_sliding_zmostdiff2_10mb[[1]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "zscore diff 2")
plot(y = sum_sliding_zdiff10mb[[1]],
     x = c(1:length(sum_sliding_zdiff10mb[[1]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (1 Mb)",
     ylab = "zscore most diff")
plot(y = sum_sliding_zcom_10mb[[1]],
     x = c(1:length(sum_sliding_zcom_10mb[[1]])),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "sliding window (10 Mb)",
     ylab = "zscore most common")

for (chr in 1:22){
    par(mfrow=c(5,1))
    plot(y = ssim_sliding_NPMI_10mb[[chr]],
         x = c(1:length(ssim_sliding_NPMI_10mb[[chr]])),
         pch = 19,
         main = paste(chr),
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (10 Mb)",
         ylab = "SSIM")
    plot(y = sum_sliding_zdiff10mb[[chr]],
         x = c(1:length(sum_sliding_zdiff10mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (10 Mb)",
         ylab = "zscore diff")
    plot(y = sum_sliding_zcom_10mb[[chr]],
         x = c(1:length(sum_sliding_zcom_10mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (10 Mb)",
         ylab = "zscore most common")
    plot(y = sum_sliding_zmostdiff1_10mb[[chr]],
         x = c(1:length(sum_sliding_zmostdiff1_10mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (10 Mb)",
         ylab = "zscore most diff 1")
    plot(y = sum_sliding_zmostdiff2_10mb[[chr]],
         x = c(1:length(sum_sliding_zmostdiff2_10mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (10 Mb)",
         ylab = "zscore most diff 2")
}


for (chr in 1:22){
    par(mfrow=c(5,1))
    plot(y = ssim_sliding_NPMI_3mb[[chr]],
         x = c(1:length(ssim_sliding_NPMI_3mb[[chr]])),
         pch = 19,
         main = paste(chr),
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (3 Mb)",
         ylab = "SSIM")
    plot(y = sum_sliding_zdiff3mb[[chr]],
         x = c(1:length(sum_sliding_zdiff3mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (3 Mb)",
         ylab = "zscore diff")
    plot(y = sum_sliding_zcom_3mb[[chr]],
         x = c(1:length(sum_sliding_zcom_3mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (3 Mb)",
         ylab = "zscore most common")
    plot(y = sum_sliding_zmostdiff1_3mb[[chr]],
         x = c(1:length(sum_sliding_zmostdiff1_3mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (3 Mb)",
         ylab = "zscore most diff 1")
    plot(y = sum_sliding_zmostdiff2_3mb[[chr]],
         x = c(1:length(sum_sliding_zmostdiff2_3mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (3 Mb)",
         ylab = "zscore most diff 2")
}

for (chr in 1:22){
    par(mfrow=c(5,1))
    plot(y = ssim_sliding_NPMI_1mb[[chr]],
         x = c(1:length(ssim_sliding_NPMI_1mb[[chr]])),
         pch = 19,
         main = paste(chr),
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (1 Mb)",
         ylab = "SSIM")
    plot(y = sum_sliding_zdiff1mb[[chr]],
         x = c(1:length(sum_sliding_zdiff1mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (1 Mb)",
         ylab = "zscore diff")
    plot(y = sum_sliding_zcom_1mb[[chr]],
         x = c(1:length(sum_sliding_zcom_1mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (1 Mb)",
         ylab = "zscore most common")
    plot(y = sum_sliding_zmostdiff1_1mb[[chr]],
         x = c(1:length(sum_sliding_zmostdiff1_1mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (1 Mb)",
         ylab = "zscore most diff 1")
    plot(y = sum_sliding_zmostdiff2_1mb[[chr]],
         x = c(1:length(sum_sliding_zmostdiff2_1mb[[chr]])),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "sliding window (1 Mb)",
         ylab = "zscore most diff 2")
}


par(mfrow=c(5,5))
for (chr in 1:22){
    plot(y = sum_sliding_zmostdiff1_10mb[[chr]],
         x = sum_sliding_zmostdiff2_10mb[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "zscore diff 2",
         ylab = "zscore diff 1")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_10mb[[chr]],
         y = sum_sliding_zcom_10mb[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "SSIM NPMI matrices",
         ylab = "zscore most common")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_3mb[[chr]],
         y = sum_sliding_zcom_3mb[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "SSIM NPMI matrices",
         ylab = "zscore most common")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_1mb[[chr]],
         y = sum_sliding_zcom_1mb[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "SSIM NPMI matrices",
         ylab = "zscore most common")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = sum_sliding_zdiff10mb[[chr]],
         y = sum_sliding_zcom_10mb[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlab = "zscore",
         ylab = "zscore most common")
}


#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
# SSIM can't handle NA, and also, it doesn't really care about the actual structure, as long as both matrices have the same dimensions
# so we change all Matrices into linear representations, 
# we get all NA positions from M1, from M2, merge them, remove all positions that are NA in either M1 and/or M2, also remove those in the zscore, diff, common eg matrices, 
# compute again

### first, proof of principle

a = matrix(data=c(1:9),nrow=3)/10
b = a
b[2,2] = 0.1
SSIM(a,b,L=1)
# $SSIM
# [1] 0.8908136
# $comps
# luminance  contrast structure
# 0.9956836 0.9944041 0.8997100
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.50000000 0.06666667 0.45555556 0.08246914 0.06666667

SSIM(as.vector(a),as.vector(b),L=1)
# $SSIM
# [1] NaN
# $comps
# luminance  contrast structure
# 1       NaN       NaN

al = matrix(data=a,nrow=1)
bl = matrix(data=b,nrow=1)
SSIM(al,bl,L=1)
# $SSIM
# [1] 0.8908136
# $comps
# luminance  contrast structure
# 0.9956836 0.9944041 0.8997100
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.50000000 0.06666667 0.45555556 0.08246914 0.06666667

a[3,3] = NA
b[1,2] = NA

SSIM(a,b,L=1)
# Error in SSIM(a, b, L = 1) : NA/NaN/Inf in foreign function call (arg 1)
al = matrix(data=a,nrow=1)
bl = matrix(data=b,nrow=1)
aNAidx = which(is.na(a))
bNAidx = which(is.na(b))
allNAidx = union(aNAidx,bNAidx)
al_noNA=matrix(al[-allNAidx],nrow=1)
bl_noNA=matrix(bl[-allNAidx],nrow=1)
SSIM(al_noNA,bl_noNA,L=1)
# $SSIM
# [1] 0.8470747
# $comps
# luminance  contrast structure
# 0.9911528 0.9939989 0.8597955
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.45714286 0.05959184 0.40000000 0.07428571 0.05714286

al_zeroNA = al
al_zeroNA[allNAidx] = 0
bl_zeroNA = bl
bl_zeroNA[allNAidx] = 0
SSIM(al_zeroNA,bl_zeroNA,L=1)
# $SSIM
# [1] 0.8983667
# $comps
# luminance  contrast structure
# 0.9911544 0.9998451 0.9065246
# $stats
# x.bar      x.var      y.bar      y.var        cov
# 0.35555556 0.08246914 0.31111111 0.08543210 0.07604938

# load maps again, linearize them and remove all entries where genome 1 and genome 2 have NAs in the NPMI map
# srun --pty --partition long --mem=300G --time 7-00 bash -i
# conda activate H1


library("SpatialPack")

dir_H1="/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/H1_phased/Curated/at50Kb"
dir_matrices="/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices"

NPMI_1 = vector(mode="list",length=22)
NPMI_2 = vector(mode="list",length=22)
zscores_common = vector(mode="list",length=22)
zscores_diff =vector(mode="list",length=22)
zscores_diff_1 = vector(mode="list",length=22)
zscores_diff_2 = vector(mode="list",length=22)

for (chr in 1:22){
    # chr = 6
    print(chr)
    #chr=1
    NPMI_1[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/NPMI_Genome1/H1.Genome1_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    NPMI_2[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/NPMI_Genome2/H1.Genome2_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    
    zscores_common[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_common/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_diff/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_1[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_2[[chr]] = data.matrix(read.table(gzfile(paste(dir_H1,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
}



# NPMI_1_0 = lapply(NPMI_1,function(chrom){
#     chrom[is.na(chrom)] <- 0
#     return(chrom)})
# NPMI_2_0 = lapply(NPMI_2,function(chrom){
#     chrom[is.na(chrom)] <- 0
#     return(chrom)})

NPMI_SSIM = vector(mode="numeric",length=22)
NPMI_SSIM_noNA = vector(mode="numeric",length=22)

for (chr in 1:22){

    print(chr)

    NPMI_SSIM[[chr]] =  SSIM(NPMI_1_0[[chr]], NPMI_2_0[[chr]], L = 1)$SSIM
    
    NA_idx_1 = which(is.na(NPMI_1[[chr]]))
    NA_idx_2 = which(is.na(NPMI_2[[chr]]))
    allNAidx = union(NA_idx_1,NA_idx_2)
    NPMI_1_lin =  matrix(data=NPMI_1[[chr]][-allNAidx],nrow=1)
    NPMI_2_lin =  matrix(data=NPMI_2[[chr]][-allNAidx],nrow=1)
    NPMI_SSIM_noNA[[chr]] =  SSIM(NPMI_1_lin, NPMI_2_lin, L = 1)$SSIM
}

saveRDS(list(NPMI_SSIM,
             NPMI_SSIM_noNA),
        paste(dir_matrices,"SSIM_chr_H1_NA_noNA.rds",sep="/"))


### sliding window first, then linearize, to avoid  hassle with the coordinates. 

ssim_sliding_NPMI = vector(mode="list",length=22)
sum_sliding_zmostdiff1 = vector(mode="list",length=22)
sum_sliding_zmostdiff2 = vector(mode="list",length=22)
sum_sliding_zmostdiff = vector(mode="list",length=22)
sum_sliding_zcom = vector(mode="list",length=22)
nonNAvalues = vector(mode="list",length=22)

# windowsize = 80 # the range is 4Mb, so 80 bins, not 3Mb = 60bins
# windowsize = 20 
 windowsize = 200

for (chr in 1:22){
    windowamount = nrow(NPMI_1[[chr]]) - windowsize

    ssim_sliding_NPMI[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff1[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff2[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zmostdiff[[chr]] = vector(mode="numeric",length=windowamount)
    sum_sliding_zcom[[chr]] = vector(mode="numeric",length=windowamount)
    nonNAvalues[[chr]] = vector(mode="numeric",length=windowamount)

    for (i in  1:windowamount){
        if(i%%100 == 0){print(paste(chr,i))}
        
        # get NA in slide 
        NA_idx_1 = which(is.na(NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        NA_idx_2 = which(is.na(NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)]))
        allNAidx = union(NA_idx_1,NA_idx_2)
        nonNAvalues[[chr]][i] = (length(NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)]) - length(allNAidx))
        
        if(length(allNAidx)==length(NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)])){
            ssim_sliding_NPMI[[chr]][i] =  NA
            sum_sliding_zmostdiff1[[chr]][i] =  NA
            sum_sliding_zmostdiff2[[chr]][i] =  NA
            sum_sliding_zmostdiff[[chr]][i] =  NA
            sum_sliding_zcom[[chr]][i] =  NA
        } else {
            if(length(allNAidx)==0){
                # linearize the window, then calc SSIM
                NPMI_1_linslide =  matrix(data=NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)],nrow=1)
                NPMI_2_linslide =  matrix(data=NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)],nrow=1)
                zscores_diff_1_linslide =  matrix(data=zscores_diff_1[[chr]][i:(i+windowsize),i:(i+windowsize)],nrow=1)
                zscores_diff_2_linslide =  matrix(data=zscores_diff_2[[chr]][i:(i+windowsize),i:(i+windowsize)],nrow=1)
                zscores_diff_linslide =  matrix(data=zscores_diff[[chr]][i:(i+windowsize),i:(i+windowsize)],nrow=1)
                zscores_common_linslide =  matrix(data=zscores_common[[chr]][i:(i+windowsize),i:(i+windowsize)],nrow=1)
                
                ssim_sliding_NPMI[[chr]][i] =  SSIM(NPMI_1_linslide, NPMI_2_linslide, L = 1)$SSIM
                sum_sliding_zmostdiff1[[chr]][i] =  sum(abs(zscores_diff_1_linslide),na.rm=T)
                sum_sliding_zmostdiff2[[chr]][i] =  sum(abs(zscores_diff_2_linslide),na.rm=T)
                sum_sliding_zmostdiff[[chr]][i] =  sum(abs(zscores_diff_linslide),na.rm=T)
                sum_sliding_zcom[[chr]][i] =  sum(abs(zscores_common_linslide),na.rm=T)
            } else {
                # linearize the window, then calc SSIM
                NPMI_1_linslide =  matrix(data=NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
                NPMI_2_linslide =  matrix(data=NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
                zscores_diff_1_linslide =  matrix(data=zscores_diff_1[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
                zscores_diff_2_linslide =  matrix(data=zscores_diff_2[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
                zscores_diff_linslide =  matrix(data=zscores_diff[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
                zscores_common_linslide =  matrix(data=zscores_common[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
                
                ssim_sliding_NPMI[[chr]][i] =  SSIM(NPMI_1_linslide, NPMI_2_linslide, L = 1)$SSIM
                sum_sliding_zmostdiff1[[chr]][i] =  sum(abs(zscores_diff_1_linslide),na.rm=T)
                sum_sliding_zmostdiff2[[chr]][i] =  sum(abs(zscores_diff_2_linslide),na.rm=T)
                sum_sliding_zmostdiff[[chr]][i] =  sum(abs(zscores_diff_linslide),na.rm=T)
                sum_sliding_zcom[[chr]][i] =  sum(abs(zscores_common_linslide),na.rm=T)
            }
        }
    }
}

saveRDS(list(ssim_sliding_NPMI,
             sum_sliding_zmostdiff1,
             sum_sliding_zmostdiff2,
             sum_sliding_zmostdiff,
             sum_sliding_zcom,
             nonNAvalues),
        #paste(dir_matrices,"SSIM_chr_sliding_80bins_4mb_H1_zscore_diff_common_sum_noNA.rds",sep="/"))
        #paste(dir_matrices,"SSIM_chr_sliding_20bins_1mb_H1_zscore_diff_common_sum_noNA.rds",sep="/"))
        paste(dir_matrices,"SSIM_chr_sliding_200bins_10mb_H1_zscore_diff_common_sum_noNA.rds",sep="/"))


# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_80bins_4mb_H1_zscore_diff_common_sum_noNA.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_20bins_1mb_H1_zscore_diff_common_sum_noNA.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_sliding_200bins_10mb_H1_zscore_diff_common_sum_noNA.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/F123/subsampled_cophased_F123_matrices/SSIM_chr_H1_NA_noNA.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/

H1_SSIM = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_H1_NA_noNA.rds")

ssim_sliding_10mb_noNA = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_200bins_10mb_H1_zscore_diff_common_sum_noNA.rds")
ssim_sliding_NPMI_10mb_noNA = ssim_sliding_10mb_noNA[[1]]
sum_sliding_zmostdiff1_10mb_noNA = ssim_sliding_10mb_noNA[[2]]
sum_sliding_zmostdiff2_10mb_noNA = ssim_sliding_10mb_noNA[[3]]
sum_sliding_zdiff10mb_noNA = ssim_sliding_10mb_noNA[[4]]
sum_sliding_zcom_10mb_noNA = ssim_sliding_10mb_noNA[[5]]
nonNAvalues_10mb = ssim_sliding_10mb_noNA[[6]]

ssim_sliding_1mb_noNA = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_20bins_1mb_H1_zscore_diff_common_sum_noNA.rds")
ssim_sliding_NPMI_1mb_noNA = ssim_sliding_1mb_noNA[[1]]
sum_sliding_zmostdiff1_1mb_noNA = ssim_sliding_1mb_noNA[[2]]
sum_sliding_zmostdiff2_1mb_noNA = ssim_sliding_1mb_noNA[[3]]
sum_sliding_zdiff1mb_noNA = ssim_sliding_1mb_noNA[[4]]
sum_sliding_zcom_1mb_noNA = ssim_sliding_1mb_noNA[[5]]
nonNAvalues_1mb = ssim_sliding_1mb_noNA[[6]]

ssim_sliding_4mb_noNA = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_80bins_4mb_H1_zscore_diff_common_sum_noNA.rds")
ssim_sliding_NPMI_4mb_noNA = ssim_sliding_4mb_noNA[[1]]
sum_sliding_zmostdiff1_4mb_noNA = ssim_sliding_4mb_noNA[[2]]
sum_sliding_zmostdiff2_4mb_noNA = ssim_sliding_4mb_noNA[[3]]
sum_sliding_zdiff4mb_noNA = ssim_sliding_4mb_noNA[[4]]
sum_sliding_zcom_4mb_noNA = ssim_sliding_4mb_noNA[[5]]
nonNAvalues_4mb = ssim_sliding_4mb_noNA[[6]]

par(mar=c(4.1,4.1,1.1,0.1)) #bltr # margins


sum_sliding_zcom_1mb_noNA

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(y = sum_sliding_zmostdiff1_10mb_noNA[[chr]],
         x = sum_sliding_zmostdiff2_10mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         main = paste("chr",chr, ": 10Mb",sep=""),
         xlab = "zscore diff 2",
         ylab = "zscore diff 1")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(y = sum_sliding_zmostdiff1_4mb_noNA[[chr]],
         x = sum_sliding_zmostdiff2_4mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         main = paste("chr",chr, ": 4Mb",sep=""),
         xlab = "zscore diff 2",
         ylab = "zscore diff 1")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(y = sum_sliding_zmostdiff1_1mb_noNA[[chr]],
         x = sum_sliding_zmostdiff2_1mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         main = paste("chr",chr, ": 1Mb",sep=""),
         xlab = "zscore diff 2",
         ylab = "zscore diff 1")
}

par(mfrow=c(2,2))
plot(y = unlist(sum_sliding_zmostdiff1_1mb_noNA),
     x = unlist(sum_sliding_zmostdiff2_1mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     main = "1Mb",
     xlab = "zscore diff 2",
     ylab = "zscore diff 1")
plot(y = unlist(sum_sliding_zmostdiff1_4mb_noNA),
     x = unlist(sum_sliding_zmostdiff2_4mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     main = "4Mb",
     xlab = "zscore diff 2",
     ylab = "zscore diff 1")
plot(y = unlist(sum_sliding_zmostdiff1_10mb_noNA),
     x = unlist(sum_sliding_zmostdiff2_10mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     main = "10Mb",
     xlab = "zscore diff 2",
     ylab = "zscore diff 1")


### NPMI vs szcore most common

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_10mb_noNA[[chr]],
         y = sum_sliding_zcom_10mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlim= c(0.4,1),
         main = paste("chr",chr, ": 10Mb",sep=""),
         xlab = "SSIM NPMI matrices",
         ylab = "zscore most common")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_4mb_noNA[[chr]],
         y = sum_sliding_zcom_4mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlim= c(0.4,1),
         main = paste("chr",chr, ": 4Mb",sep=""),
         xlab = "SSIM NPMI matrices",
         ylab = "zscore most common")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_1mb_noNA[[chr]],
         y = sum_sliding_zcom_1mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlim= c(0.4,1),
         main = paste("chr",chr, ": 1Mb",sep=""),
         xlab = "SSIM NPMI matrices",
         ylab = "zscore most common")
}

par(mfrow=c(2,2))
plot(x = unlist(ssim_sliding_NPMI_1mb_noNA),
     y = unlist(sum_sliding_zcom_1mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlim= c(0.4,1),
     main = "1Mb",
     xlab = "SSIM NPMI matrices",
     ylab = "zscore most common")
plot(x = unlist(ssim_sliding_NPMI_4mb_noNA),
     y = unlist(sum_sliding_zcom_4mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlim= c(0.4,1),
     main = "4Mb",
     xlab = "SSIM NPMI matrices",
     ylab = "zscore most common")
plot(x = unlist(ssim_sliding_NPMI_10mb_noNA),
     y = unlist(sum_sliding_zcom_10mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlim= c(0.4,1),
     main = "10Mb",
     xlab = "SSIM NPMI matrices",
     ylab = "zscore most common")


### NPMI vs szcore diff

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_10mb_noNA[[chr]],
         y = sum_sliding_zdiff10mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlim= c(0.4,1),
         main = paste("chr",chr, ": 10Mb",sep=""),
         xlab = "SSIM NPMI matrices",
         ylab = "zscore difference")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_4mb_noNA[[chr]],
         y = sum_sliding_zdiff4mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlim= c(0.4,1),
         main = paste("chr",chr, ": 4Mb",sep=""),
         xlab = "SSIM NPMI matrices",
         ylab = "zscore difference")
}

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = ssim_sliding_NPMI_1mb_noNA[[chr]],
         y = sum_sliding_zdiff1mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         xlim= c(0.4,1),
         main = paste("chr",chr, ": 1Mb",sep=""),
         xlab = "SSIM NPMI matrices",
         ylab = "zscore difference")
}

par(mfrow=c(2,2))
plot(x = unlist(ssim_sliding_NPMI_1mb_noNA),
     y = unlist(sum_sliding_zdiff1mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlim= c(0.4,1),
     main = "1Mb",
     xlab = "SSIM NPMI matrices",
     ylab = "zscore difference")
plot(x = unlist(ssim_sliding_NPMI_4mb_noNA),
     y = unlist(sum_sliding_zdiff4mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlim= c(0.4,1),
     main = "4Mb",
     xlab = "SSIM NPMI matrices",
     ylab = "zscore difference")
plot(x = unlist(ssim_sliding_NPMI_10mb_noNA),
     y = unlist(sum_sliding_zdiff10mb_noNA),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlim= c(0.4,1),
     main = "10Mb",
     xlab = "SSIM NPMI matrices",
     ylab = "zscore difference")

par(mfrow=c(1,1))
pairs(~unlist(ssim_sliding_NPMI_1mb_noNA) + 
          unlist(sum_sliding_zdiff1mb_noNA) + 
          unlist(sum_sliding_zcom_1mb_noNA) + 
          unlist(sum_sliding_zmostdiff1_1mb_noNA) +  
          unlist(sum_sliding_zmostdiff2_1mb_noNA) +
          unlist(nonNAvalues_1mb),
      labels = c("SSIM","zscore differences","10% most common","5% most different G1", "5% most different G2", "amount non NA in sliding window" ),
      pch = 19,
      col = adjustcolor("black", alpha.f = 0.5),
      cex = 0.1)

pairs(~unlist(ssim_sliding_NPMI_4mb_noNA) + 
          unlist(sum_sliding_zdiff4mb_noNA) + 
          unlist(sum_sliding_zcom_4mb_noNA) + 
          unlist(sum_sliding_zmostdiff1_4mb_noNA) +  
          unlist(sum_sliding_zmostdiff2_4mb_noNA) +
          unlist(nonNAvalues_4mb),
      labels = c("SSIM","zscore differences","10% most common","5% most different G1", "5% most different G2", "amount non NA in sliding window" ),
      pch = 19,
      col = adjustcolor("black", alpha.f = 0.5),
      cex = 0.1)

pairs(~unlist(ssim_sliding_NPMI_10mb_noNA) + 
          unlist(sum_sliding_zdiff10mb_noNA) + 
          unlist(sum_sliding_zcom_10mb_noNA) + 
          unlist(sum_sliding_zmostdiff1_10mb_noNA) +  
          unlist(sum_sliding_zmostdiff2_10mb_noNA) +
          unlist(nonNAvalues_10mb),
      labels = c("SSIM","zscore differences","10% most common","5% most different G1", "5% most different G2", "amount non NA in sliding window" ),
      pch = 19,
      col = adjustcolor("black", alpha.f = 0.5),
      cex = 0.1)

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = sum_sliding_zdiff10mb_noNA[[chr]],
         y = sum_sliding_zcom_10mb_noNA[[chr]],
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         main = paste("chr",chr, ": 10Mb",sep=""),
         xlab = "zscore",
         ylab = "zscore most common")
}
par(mfrow=c(5,1))
plot(y = ssim_sliding_NPMI_1mb_noNA[[chr]],
     x = c(1:length(ssim_sliding_NPMI_1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,1))
# points(y = ssim_sliding_NPMI_4mb_noNA[[chr]],
#        x = c(40:(length(ssim_sliding_NPMI_4mb_noNA[[chr]])+39)), pch = 19, col = "red")
# points(y = ssim_sliding_NPMI_10mb_noNA[[chr]],
#        x = c(200:(length(ssim_sliding_NPMI_10mb_noNA[[chr]])+199)), pch = 19, col = "blue")
# not really, go and figure out the actual indices

plot(y = ssim_sliding_NPMI_1mb_noNA[[chr]],
     x = c(1:length(ssim_sliding_NPMI_1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,1))
points(y = ssim_sliding_NPMI_4mb_noNA[[chr]],
       x = c(1:length(ssim_sliding_NPMI_4mb_noNA[[chr]])), pch = 19, col = "red")
points(y = ssim_sliding_NPMI_10mb_noNA[[chr]],
       x = c(1:length(ssim_sliding_NPMI_10mb_noNA[[chr]])), pch = 19, col = "blue")

# plot those together in one plot, plot the values in the middle of the sliding window
# we had windowsizes of 20, 80 and 200
# meaning with 50kb resolution the windows were 1Mb, 4Mb and 10Mb wide
# so for the 1Mb, the first value is at 500k, the next 550k, next 600k etc
# for the 4Mb sliding windows, the first value is at 2Mb, the next at 2050k, the next at 2100k
# for the 10mb, first value at 5mb, next at 5050k, next at 5100k
# chr = 19
# length(ssim_sliding_NPMI_1mb_noNA[[chr]]) # 1153
# length(ssim_sliding_NPMI_4mb_noNA[[chr]]) # 1093
# length(ssim_sliding_NPMI_10mb_noNA[[chr]]) # 973
names(ssim_sliding_NPMI_1mb_noNA[[chr]]) = seq(from=500000,by=50000,length.out = length(ssim_sliding_NPMI_1mb_noNA[[chr]]) )
names(ssim_sliding_NPMI_4mb_noNA[[chr]]) = seq(from=2000000,by=50000,length.out = length(ssim_sliding_NPMI_4mb_noNA[[chr]]) )
names(ssim_sliding_NPMI_10mb_noNA[[chr]]) = seq(from=5000000,by=50000,length.out = length(ssim_sliding_NPMI_10mb_noNA[[chr]]) )

chr = 20
plot(y = ssim_sliding_NPMI_1mb_noNA[[chr]],
     x = as.numeric(names(ssim_sliding_NPMI_1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,1))
points(y = ssim_sliding_NPMI_4mb_noNA[[chr]],
       x = as.numeric(names(ssim_sliding_NPMI_4mb_noNA[[chr]])), pch = 19, col = "red")
points(y = ssim_sliding_NPMI_10mb_noNA[[chr]],
       x = as.numeric(names(ssim_sliding_NPMI_10mb_noNA[[chr]])), pch = 19, col = "blue")
### add those to the large plots in script visualize_diff_contacts_zscores_npmi_maps_phased_H1.R


#############################################################
#############################################################
#############################################################
# # 05.10.21: when plotting the sliding SSIM together with the NPMI maps, 
# # the 1Mb SSIM values are very sparse and there is so far no good explanation for this, 
# # so let's check what is going on here
# 
# # also the points do not start where they should, eg the 10Mb should start at 5Mb, but this could be explained by mising data
# # also save the number of non NA values in the square, just to check if this biases the computed SSIM values a lot
# 
# # values for the 1Mb sliding window starts at 500kb (middle of 1Mb)
# head(ssim_sliding_NPMI_1mb_noNA[[22]])
# # NPMI map has first values at maybe 17Mb, with 50kb resolution this is column 
# NPMI_1[[22]][1:5,1:5]
# NPMI_1[[22]][340:350,340:350]
# # earlier even
# which(!(is.na(diag(NPMI_1[[22]]))))[1] # 334
# NPMI_1[[22]][330:340,330:340] #chr22.16650000.16700000
# which(names(ssim_sliding_NPMI_1mb_noNA[[22]])==16650000) # 324. makes sense, they start 10 bins in
# ssim_sliding_NPMI_1mb_noNA[[22]][320:350]
# # already here many are NaN: 17200000  17250000    17350000  17400000  17450000  17500000  17550000  17600000 17650000, why?
# windowsize = 20 
# chr = 22
# windowamount = nrow(NPMI_1[[chr]]) - windowsize
# i = 335
# # get NA in slide 
#         NA_idx_1 = which(is.na(NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)]))
#         NA_idx_2 = which(is.na(NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)]))
#         allNAidx = union(NA_idx_1,NA_idx_2) # if there are no NAs in both regions, this is empty: allNAidx: integer(0)
#     # which triggers an error in the next step, bc [-allNAidx] with an empty allNAidx only gives: numeric(0)
#         NPMI_1_linslide =  matrix(data=NPMI_1[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
#         NPMI_2_linslide =  matrix(data=NPMI_2[[chr]][i:(i+windowsize),i:(i+windowsize)][-allNAidx],nrow=1)
# ### I corrected it in the code above, is fine now
# ssim_sliding_NPMI_1mb_noNA[[22]][320:350] # does not contain any NaNs anymore 


#############################################################
#############################################################
# 05.10.21
# influence of amount of missing data (NA) in sliding window
par(mfrow=c(5,5))
for (chr in 1:22){
    plot(y = ssim_sliding_NPMI_1mb_noNA[[chr]],
         x = nonNAvalues_1mb[[chr]],
#         x = (max(nonNAvalues_1mb[[chr]])-nonNAvalues_1mb[[chr]]),
         pch = 19,
         col = adjustcolor("black", alpha.f = 0.5),
         cex = 0.1,
         ylab = "SSIM",
         xlab = "amount non NA values in sliding window")
    #lines(lowess(nonNAvalues_1mb[[chr]], ssim_sliding_NPMI_1mb_noNA[[chr]]), col = "blue")
    abline(lm(ssim_sliding_NPMI_1mb_noNA[[chr]] ~ nonNAvalues_1mb[[chr]]), col = "red")
}
par(mfrow=c(1,1))
par(mar=c(4.1,4.1,0.1,0.1)) #bltr # margins
plot(y = unlist(ssim_sliding_NPMI_1mb_noNA),
     x = unlist(nonNAvalues_1mb),
     pch = 19,
     col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     ylab = "SSIM",
     xlab = "amount non NA values in sliding window")
abline(lm(unlist(ssim_sliding_NPMI_1mb_noNA) ~ unlist(nonNAvalues_1mb)), col = "red")

# clear negative correlation! more NAs make area more similar