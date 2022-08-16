# plot phased H1 NPMI matrices and Zscores
# 16.06.2021


# It s easiest to mount the cluster and loed the gz files directly
# locally, change to "GAM", then:
# sshfs -o allow_other,defer_permissions,follow_symlinks jmarkow_m@bihtext:/ mnt/bih_cluster
dir="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/H1_phased/Curated/at50Kb"

# plot NPMI, chr 22 for now
NPMI_22_1 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome1/H1.Genome1_chr22_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
NPMI_22_2 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome2/H1.Genome2_chr22_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))

dim(NPMI_22_1) # 1017 1017
table(is.na(NPMI_22_1))
# FALSE   TRUE 
# 84381 949908


library(gplots)
library(RColorBrewer)
heatmap.2(NPMI_22_1, dendrogram = "none", Rowv = F, Colv = F, trace = "none")
heatmap.2(NPMI_22_2, dendrogram = "none", Rowv = F, Colv = F, trace = "none")

# plot Zscores, chr 22 for now
zscores_common_22 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_common/H1.Genome1_H1.Genome2_chr22_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
zscores_diff_22 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_diff/H1.Genome1_H1.Genome2_chr22_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
zscores_diff_1_22 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr22_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
zscores_diff_2_22 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr22_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
dim(zscores_common_22) # 1017 1017
dim(zscores_diff_1_22) # 1017 1017
dim(zscores_diff_2_22) # 1017 1017


heatmap.2(zscores_common_22, dendrogram = "none", Rowv = F, Colv = F, trace = "none")
heatmap.2(zscores_diff_1_22, dendrogram = "none", Rowv = F, Colv = F, trace = "none")
heatmap.2(zscores_diff_2_22x, dendrogram = "none", Rowv = F, Colv = F, trace = "none")


### testing how to plot the matrix itself as colour and then rotate by 45 degrees

test_npmi_2 = NPMI_22_2[400:499,400:499]
test_npmi_1 = NPMI_22_1[400:499,400:499]
test_zdif_2 = zscores_diff_2_22[400:499,400:499]
test_zdif_1 = zscores_diff_1_22[400:499,400:499]
test_zcom = zscores_common_22[400:499,400:499]
test_zdiff = zscores_diff_22[400:499,400:499]

heatmap.2(test_npmi_2, dendrogram = "none", Rowv = F, Colv = F, trace = "none")

tickpos = vector(mode  ="numeric", length = 25)
for (i in 1:25){
    mb = i*10000000
    tickpos[i] = sqrt(mb^2 + mb^2)/5e4
}

plotTriMatrix <- function(x) {
    x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ii <- cut(x, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
    colors <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ii]
    xxx=matrix(colors,nrow=nrow(x),ncol=ncol(x))
    
    ## clear lower triangle
    xxx[lower.tri(xxx)] <- NA
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(0, 60), xlab="", ylab="", asp=1)
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}

plotTriMatrix(test_npmi_2)

plotTriMatrix_w <- function(x) {
    x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ii <- cut(x, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
    colors <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ii]
    xxx=matrix(colors,nrow=nrow(x),ncol=ncol(x))
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}

plotTriMatrix_w(test_npmi_2)

plotTriMatrix_w_both <- function(x,y) {
    x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ix <- cut(x, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
    colorsx <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ix]
    xxx=matrix(colorsx,nrow=nrow(x),ncol=ncol(x))
    
    iy <- cut(y, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
    colorsy <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[iy]
    yyy=matrix(colorsy,nrow=nrow(y),ncol=ncol(y))
    
    ## clear lower triangle
    xxx[lower.tri(xxx)] <- yyy[lower.tri(yyy)]
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1,  xaxt="none")
#    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    abline(h=0)
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}

plotTriMatrix_w(test_npmi_2)
plotTriMatrix_w(test_npmi_1)
plotTriMatrix_w_both(test_npmi_2,test_npmi_1)
plotTriMatrix_w(test_zcom)
# plotTriMatrix_w(test_zdiff)
# plotTriMatrix_w_both(test_zdif_2,test_zdif_1)
# those do not work yet bc the have negative values 

# also figure out the axis and plot SNP positions!
colnames(NPMI_22_2)[500] #chr22.24950000.25000000
colnames(NPMI_22_2)[400] #chr22.19950000.20000000
# this is 100 * 50kb resolution = 5000000
#sqrt(5e6^2 + 5e6^2) #7071068 , times 2 I guess bc of 50kb resolution

# add the axis names 
chromsize_hg38=read.table('/Users/jmarkow/Desktop/GAM/data/hg38_sorted_autosomes.bed',header = F)
chromsize_hg38 = chromsize_hg38[,c(1,3)]
colnames(chromsize_hg38) = c("chr","size")
###    

plotTriMatrix_w_both(NPMI_22_1,NPMI_22_2)
chromsize_hg38[22,] # 50818468
colnames(NPMI_22_1)[1] #chr22.0.50000
colnames(NPMI_22_1)[ncol(NPMI_22_1)] #chr22.50800000.50818468
sqrt(chromsize_hg38[22,2]^2 + chromsize_hg38[22,2]^2)/5e4 # 1437.363 # perfect
sqrt(5e7^2 + 5e7^2)/5e4 # 1414.214 corresponds to 50Mb

# largest chr: 1
chromsize_hg38[1,2] # 248956422
# 248956422
# 249000000
#  50818468
# count in 10mb steps, so we have 5 marks on chr 22, 24 on chr 1
tickpos = vector(mode  ="numeric", length = 25)
for (i in 1:25){
    mb = i*10000000
    tickpos[i] = sqrt(mb^2 + mb^2)/5e4
}

# add the SNPs

SNPs_H1 = readRDS("/Users/jmarkow/Desktop/GAM/data/H1_truths_mapp_filtered_H1_GT.rds")
nrow(SNPs_H1) #1506858
SNPs_H1_possort = SNPs_H1[order(SNPs_H1$POS),]
SNPs_H1 = SNPs_H1_possort[order(as.numeric(substring(SNPs_H1_possort$CHROM, 4))),]

plotTriMatrix_w_both(NPMI_22_1,NPMI_22_2)
snpchr = SNPs_H1[SNPs_H1$CHROM=="chr22",2]
SNP_pos = sqrt(snpchr^2 + snpchr^2)/5e4

points(y=rep(0,times=length(snpchr)),
       x=SNP_pos,
       cex=0.1,pch=19)

# with 1 SNP every 1 kb, but 50kb resolution, we have quite a lot of bins with SNPs inside, maybe I can make a heatmap for that too?
inch=0.393701
#pdf("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_chr22.pdf",width=200*inch,height=30*inch)
plotTriMatrix_w_both(NPMI_22_1,NPMI_22_2)
points(y=rep(0,times=length(snpchr)),
       x=SNP_pos,
       cex=0.05,pch=19)
#dev.off()
# seems like y and x coordinates of SNPs are off 

# make those comparative npmi matrix plots for all chromosomes
dir="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/H1_phased/Curated/at50Kb"

# add the axis names 
chromsize_hg38=read.table('/Users/jmarkow/Desktop/GAM/data/hg38_sorted_autosomes.bed',header = F)
chromsize_hg38 = chromsize_hg38[,c(1,3)]
colnames(chromsize_hg38) = c("chr","size")

tickpos = vector(mode  ="numeric", length = 25)
for (i in 1:25){
    mb = i*10000000
    tickpos[i] = sqrt(mb^2 + mb^2)/5e4
}

SNPs_H1 = readRDS("/Users/jmarkow/Desktop/GAM/data/H1_truths_mapp_filtered_H1_GT.rds")
SNPs_H1_possort = SNPs_H1[order(SNPs_H1$POS),]
SNPs_H1 = SNPs_H1_possort[order(as.numeric(substring(SNPs_H1_possort$CHROM, 4))),]

inch=0.393701

for (chr in 1:22){
    
    #chr=1
    NPMI_1 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome1/H1.Genome1_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    NPMI_2 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome2/H1.Genome2_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    
    zscores_common = data.matrix(read.table(gzfile(paste(dir,"/ZScore_common/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff = data.matrix(read.table(gzfile(paste(dir,"/ZScore_diff/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_1 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_2 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    
    
    pdf(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_chr",chr,".pdf",sep=""),width=200*inch,height=30*inch) #90
    par(mfrow=c(4,1))
    plotTriMatrix_w_both(NPMI_1,NPMI_2)
    
    snpchr = SNPs_H1[substring(SNPs_H1$CHROM, 4)==chr,2]
    SNP_pos = sqrt(snpchr^2 + snpchr^2)/5e4
    points(y=rep(0,times=length(snpchr)),
           x=SNP_pos,
           cex=0.05,pch=19)
    
    plotTriMatrix_zscores_w(zscores_diff)
    plotTriMatrix_zscores_w(zscores_common)
    plotTriMatrix_zcores_w_both(zscores_diff_1,zscores_diff_2)
    
    dev.off()
    
}

# the matrices are not centered at 0, also the dimension of the plot, the height/width ratio needs to be adapted per chromosome

# for now, try to plot the zscores properly, inlcuding the negative values

# test_zdif_2 = zscores_diff_2_22[400:499,400:499]
# test_zdif_1 = zscores_diff_1_22[400:499,400:499]
# test_zcom = zscores_common_22[400:499,400:499]
# test_zdiff = zscores_diff_22[400:499,400:499]

min(test_zdiff,na.rm = T) # -4.28689
max(test_zdiff,na.rm = T) #  4.510395
x = test_zdiff

x[1,1:3]
ii[1:3]
plotTriMatrix_zscores_w <- function(x) {
    #x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ii <- cut(x, breaks = seq(-5, 5, len = 100), 
              include.lowest = TRUE)
    colors <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ii]
    xxx=matrix(colors,nrow=nrow(x),ncol=ncol(x))
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}
# 
# color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
#     scale = (length(lut)-1)/(max-min)
#     
#     dev.new(width=1.75, height=5)
#     plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
#     axis(2, ticks, las=1)
#     for (i in 1:(length(lut)-1)) {
#         y = (i-1)/scale + min
#         rect(0,y,10,y+1/scale, col=lut[i], border=NA)
#     }
# }
# color.bar(colorRampPalette(c("light green", "yellow", "orange", "red"))(100), -1)

plotTriMatrix_zscores_w(test_zdiff)
plotTriMatrix_zscores_w(test_zcom)

plotTriMatrix_zcores_w_both <- function(x,y) {
    #x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ix <- cut(x, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    colorsx <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ix]
    xxx=matrix(colorsx,nrow=nrow(x),ncol=ncol(x))
    
    iy <- cut(y, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    colorsy <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[iy]
    yyy=matrix(colorsy,nrow=nrow(y),ncol=ncol(y))
    
    ## clear lower triangle
    xxx[lower.tri(xxx)] <- yyy[lower.tri(yyy)]
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1,  xaxt="none")
    #    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    abline(h=0)
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}

plotTriMatrix_zcores_w_both(test_zdif_1,test_zdif_2)
plotTriMatrix_zcores_w_both(zscores_diff_1_22,zscores_diff_2_22)

### make a colour code bar for the number of SNPs per bin
bins_all =  seq(0,250000000, 50000)
bins <- cut(snpchr, breaks = bins_all)
length(bins)
length(snpchr)
head(bins)
head(snpchr)
snp_dens=table(bins)
length(snp_dens)
tail(snp_dens)
barplot(snp_dens) # nice 
snp_dens[1000:1010]
#barplot(log(snp_dens+1))
colors_snps <- colorRampPalette(c("white","orchid4"))(max(snp_dens))[snp_dens]
s_dens_pos = sqrt(bins_all^2 + bins_all^2)/5e4
#s_dens_pos = s_dens_pos[1:length(colors_snps)]
points(x = s_dens_pos[1:length(colors_snps)],
       y = rep(100,times=length(colors_snps)),
       pch = 19,
       col = colors_snps) 


# diverge0 <- function(p, ramp) {
#     # p: a trellis object resulting from rasterVis::levelplot
#     # ramp: the name of an RColorBrewer palette (as character), a character 
#     #       vector of colour names to interpolate, or a colorRampPalette.
#     require(RColorBrewer)
#     require(rasterVis)
#     if(length(ramp)==1 && is.character(ramp) && ramp %in% 
#        row.names(brewer.pal.info)) {
#         ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
#     } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
#         ramp <- colorRampPalette(ramp)
#     } else if(!is.function(ramp)) 
#         stop('ramp should be either the name of a RColorBrewer palette, ', 
#              'a vector of colours to be interpolated, or a colorRampPalette.')
#     rng <- range(p$legend[[1]]$args$key$at)
#     s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
#     i <- findInterval(rng[which.min(abs(rng))], s)
#     zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
#     p$legend[[1]]$args$key$at <- s[zlim]
#     p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
#     p
# }
# try = diverge0(snp_dens, ramp=colorRampPalette(c("white","orchid4"))(max(snp_dens)))
# install.packages("rasterVis")


# col5 <- colorRampPalette(c("white","orchid4"))
# color_levels=max(snp_dens) #the number of colors to use
# max_absolute_value=max(snp_dens) #what is the maximum absolute value of raster?
# color_sequence=seq(-max_absolute_value,max_absolute_value,length.out=color_levels+1)
# plot(img, col=col5(n=color_levels), breaks=color_sequence, axes=FALSE)
# points(x = s_dens_pos[1:length(colors_snps)],
#        y = rep(200,times=length(colors_snps)),
#        pch = 19,
#        col = col5(n=color_levels), breaks=color_sequence)
#        
# plot(mask_data, add=T)
# 
# heatmap.2(snp_dens)
plotTriMatrix_zcores_w_both(zscores_diff_1_22,zscores_diff_2_22)
chr=22
snpchr = SNPs_H1[substring(SNPs_H1$CHROM, 4)==chr,2]
bins_all =  seq(0,250000000, 50000)
bins <- cut(snpchr, breaks = bins_all)
snp_dens=table(bins)
s_dens_pos = sqrt(bins_all^2 + bins_all^2)/5e4

my.colors<-colorRampPalette(c("white","orchid4")) #creates a function my.colors which interpolates n colors between blue, white and red
my_colors2 = colorRampPalette(c("white","orchid4"))(max(snp_dens)+1) # +1 for the 0 values
#colorRampPalette(c("white","orchid4"))(max(snp_dens))[snp_dens]
#color.df<-data.frame(val=seq(0,max(snp_dens),1), colorname=my.colors(max(snp_dens)+1)) 
color_df = data.frame(val=0:max(snp_dens),colorname=my_colors2,stringsAsFactors = F)
all(color_df$colorname==my_colors2) #TRUE
snp_dens_df = data.frame(val=snp_dens,order=c(1:length(snp_dens)))
data_col <- merge(snp_dens_df, color_df,by.x="val.Freq",by.y="val")
data_col_sort = data_col[order(data_col$order),]

points(x = s_dens_pos[-length(s_dens_pos)],
       y = rep(60,times=nrow(snp_dens_df)),
       pch = 19,
       col = data_col_sort$colorname,cex=0.1)

points(x=s_dens_pos[401],y=100,col="red",pch=19)
reg1.with.color[401:411,]
barplot(snp_dens)
which(snp_dens>0)[1:10]
max(snp_dens)
which(snp_dens==max(snp_dens))
snp_dens_df[688,]
s_dens_pos[688]
xpos=sqrt((3.435e+07)^2 + (3.435e+07)^2)/5e4
points(x=xpos,y=159,pch=19,col="#8B4789")#reg1.with.color_sort[688,"colorname"])

points(x = s_dens_pos[-length(s_dens_pos)],
       y = rep(300,times=nrow(snp_dens_df)),
       pch = 19,
       col = reg1.with.color_sort$colorname)

plot(x=color_df$val,y=color_df$val,col=color_df$colorname,pch=19)
plot(x=color_df$val,y=color_df$val,col=my_colors2,pch=19)
all(color_df$colorname==my_colors2) #TRUE


4000000/50000 # 80 bins
sqrt(4000000^2 + 4000000^2)/5e4
# ok so we have 80 bins, 113 measures on the diagonal, as we are centering them at 0, we can cut off at 60 both sides for sure, now we have to deal with the width and height of the plot
# I would assume a ratio of chr length / 4000000
chromsize_hg38$size/4000000
# [1] 62.23911 60.54838 49.57389 47.55364 45.38456 42.70149 39.83649 36.28466 34.59868 33.44936 33.77166 33.31883 28.59108 26.76093 25.49780 22.58459 20.81436 20.09332 14.65440
# [20] 16.11104 11.67750 12.70462
chromsize_hg38$size/(4000000*4)

200/90


plotTriMatrix_zcores_w_both <- function(x,y) {
    #x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ix <- cut(x, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    colorsx <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ix]
    xxx=matrix(colorsx,nrow=nrow(x),ncol=ncol(x))
    
    iy <- cut(y, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    colorsy <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[iy]
    yyy=matrix(colorsy,nrow=nrow(y),ncol=ncol(y))
    
    ## clear lower triangle
    xxx[lower.tri(xxx)] <- yyy[lower.tri(yyy)]
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1,  xaxt="none")
    #    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    #abline(h=0)
    
    snpchr = SNPs_H1[substring(SNPs_H1$CHROM, 4)==chr,2]
    SNP_pos = sqrt(snpchr^2 + snpchr^2)/5e4
    points(y=rep(0,times=length(snpchr)),
           x=SNP_pos,
           cex=0.05,pch=19)
    
    snpchr = SNPs_H1[substring(SNPs_H1$CHROM, 4)==chr,2]
    bins_all =  seq(0,250000000, 50000)
    bins <- cut(snpchr, breaks = bins_all)
    snp_dens=table(bins)
    my_colors2 = colorRampPalette(c("white","orchid4"))(max(snp_dens)+1) # do this max over all so it is comparable
    color_df = data.frame(val=0:max(snp_dens),colorname=my_colors2,stringsAsFactors = F)
    all(color_df$colorname==my_colors2) #TRUE
    snp_dens_df = data.frame(val=snp_dens,order=c(1:length(snp_dens)))
    data_col <- merge(snp_dens_df, color_df,by.x="val.Freq",by.y="val")
    data_col_sort = data_col[order(data_col$order),]
    
    s_dens_pos = sqrt(bins_all^2 + bins_all^2)/5e4
    
    points(x = s_dens_pos[-length(s_dens_pos)],
           y = rep(60,times=nrow(snp_dens_df)),
           pch = 19,
           col = data_col_sort$colorname,cex=0.1)
    
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}

# try to get white to the middle
# no, do not use white, bc we need to see the mids (scale is for the zscore differences, from - over 0 to +)
colors_test <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(100)[1:100]

plot(x = 1:100,
     y = rep(10,times=100),
     ylim = c(0,10),
     pch=19,
     col=colors_test)
colors_test2 <- colorRampPalette(c("steelblue4","lightblue","white","darkgoldenrod1", "orangered"))(100)[1:100]
points(x = 1:100,
       y = rep(9,times=100),
       pch=19,
       col=colors_test2)
colors_test3 <- colorRampPalette(c("steelblue4","lightblue","white","darkgoldenrod1", "red4"))(100)[1:100]
points(x = 1:100,
       y = rep(8,times=100),
       pch=19,
       col=colors_test3)
colors_test4 <- colorRampPalette(c("steelblue4","deepskyblue","lightblue","white","darkgoldenrod1", "orangered", "red4"))(100)[1:100]
points(x = 1:100,
       y = rep(7,times=100),
       pch=19,
       col=colors_test4)
colors_test5 <- colorRampPalette(c("darkblue","steelblue4","lightblue","white","darkgoldenrod1", "orangered", "red4"))(100)[1:100]
points(x = 1:100,
       y = rep(6,times=100),
       pch=19,
       col=colors_test5)
colors_test6 <- colorRampPalette(c("midnightblue","steelblue4","lightblue","white","darkgoldenrod1", "orangered", "red4"))(100)[1:100]
points(x = 1:100,
       y = rep(5,times=100),
       pch=19,
       col=colors_test6)
colors_test7 <- colorRampPalette(c("steelblue4","steelblue2","wheat1","darkgoldenrod1", "orangered", "red4"))(100)[1:100]
points(x = 1:100,
       y = rep(5,times=100),
       pch=19,
       col=colors_test7)
colors_test8 <- colorRampPalette(c("steelblue4","steelblue2","wheat1", "orangered", "red4"))(100)[1:100]
points(x = 1:100,
       y = rep(4,times=100),
       pch=19,
       col=colors_test8)


# find max number of SNPs per bin over all chromosomes, so the colour code is comparable over chromosomes

SNPs_H1_chr = split(SNPs_H1,SNPs_H1$CHROM)
SNPs_H1_chr = SNPs_H1_chr[order(as.numeric(substring(names(SNPs_H1_chr), 4)))]

bins_all =  seq(0,250000000, 50000)
bins_all_mean = seq(25000,250000000, 50000)
bin_pos = sqrt(bins_all_mean^2 + bins_all_mean^2)/5e4

SNPs_H1_chr_binned = lapply(SNPs_H1_chr,function(chr){
    bins = cut(chr$POS, breaks = bins_all)
    return(table(bins))
})
max_snps_per_bin = max(unlist(lapply(SNPs_H1_chr_binned,max)))
# check how the SNP dens is distributed
hist(unlist(SNPs_H1_chr_binned),breaks=c(0:max_snps_per_bin),xlim=c(1,max_snps_per_bin),ylim=c(0,3000))
#hist(unlist(SNPs_H1_chr_binned),breaks=c(0:max_snps_per_bin),xlim=c(1,max_snps_per_bin),ylim=c(0,20))
# majority of bins has artificial 0 (bc every chr has binnes up until the length of chr 1, so bins exceed chr size)
# majority of (real) bins has up to 50 maybe 100 snps, so maybe a notched col scale would be good
# less than 5 bins have more than 15 SNPs (each)

# colour code for the SNPs density
snp_dens_col = colorRampPalette(c("white","orchid4"))(max_snps_per_bin+1) # do this max over all so it is comparable
snp_dens_col_log = colorRampPalette(c("white","orchid4"))(ceiling(log(max_snps_per_bin))+1)# I put the +1 outsiede the log, bc later, I want to assign all values a col, and then assign 0 the white

snp_dens_col_base = brewer.pal(n = 11, name = "BrBG")
snp_dens_col_base = snp_dens_col_base[6:11]
snp_dens_col_base[1] = "white"
snp_dens_col = colorRampPalette(c(snp_dens_col_base))(max_snps_per_bin+1)
snp_dens_col_log = colorRampPalette(c(snp_dens_col_base))(ceiling(log(max_snps_per_bin))+1)# I put the +1 outsiede the log, bc later, I want to assign all values a col, and then assign 0 the white

# create dataframe with x position and colour (snp dens) for each chromosome
color_spn_dense_df = data.frame(val=0:max(max_snps_per_bin),colorname = snp_dens_col,stringsAsFactors = F)
color_spn_dense_log_df = data.frame(val=0:ceiling(log(max(max_snps_per_bin))),colorname = snp_dens_col_log,stringsAsFactors = F)

# since the majority of bins has very few snps, but the log does not really cut it, I decided to make a gradient up until 100 SNPs and set everything above to the max colour 
snp_dens_col_100lim = colorRampPalette(c(snp_dens_col_base))(101)
color_spn_dense_100lim_df = data.frame(val=0:max(max_snps_per_bin),colorname = c(snp_dens_col_100lim,rep(snp_dens_col_100lim[length(snp_dens_col_100lim)],times=(max(max_snps_per_bin)-100))),stringsAsFactors = F)


# #test plot
# test=color_spn_dense_df
# test$val_log = ceiling(log(test$val+1))
# test = merge(test,color_spn_dense_log_df, by.x="val_log",by.y="val")
# test = test[order(test$val),]
# plot(x=test$val,
#      y=rep(0,times=nrow(test)),
#      col=test$colorname.x,pch=19)
# points(x=test$val,
#        y=rep(0.5,times=nrow(test)),
#        col=test$colorname.y,pch=19)
# points(x=color_spn_dense_100lim_df$val,
#        y=rep(-0.5,times=nrow(color_spn_dense_100lim_df)),
#        col=color_spn_dense_100lim_df$colorname,pch=19)
# 
# hist(unlist(SNPs_H1_chr_binned),breaks=c(0:max_snps_per_bin),xlim=c(1,max_snps_per_bin),ylim=c(0,3000))
# points(x=color_spn_dense_100lim_df$val,
#        y=rep(1000,times=nrow(color_spn_dense_100lim_df)),
#        col=color_spn_dense_100lim_df$colorname,pch=19)
# points(x=test$val,
#      y=rep(1500,times=nrow(test)),
#      col=test$colorname.x,pch=19)

SNPs_H1_chr_binned_col = lapply(SNPs_H1_chr_binned,function(chr){
    # chr = SNPs_H1_chr_binned[[1]]
    chr_df = data.frame(val=chr,order=c(1:length(chr)))
    chr_df$val_log = ceiling(log(chr_df$val.Freq+1))
    # hist(chr_df$val_log,breaks=seq((-0.5),6.5,1))
    chr_df$val_log[which(chr_df$val.Freq==0)] = 0
    # hist(chr_df$val_log,breaks=seq((-0.5),6.5,1))
    
    df = merge(chr_df, color_spn_dense_df,by.x="val.Freq",by.y="val")
    df = merge(df, color_spn_dense_log_df,by.x="val_log",by.y="val")
    df = merge(df, color_spn_dense_100lim_df,by.x="val.Freq",by.y="val")
    
    df_sort = df[order(df$order),]
    df_sort$plot_pos = bin_pos
    df_sort$bin_mids = bins_all_mean
    
    df_sort = df_sort[,c(4,8,9,3,1,5,7,2,6)]
    colnames(df_sort) = c("order","plot_pos","bin_mids","bin","value","colour","colour_100lim","value_log","colour_log")
    return(df_sort)
})

# plot(NA, type="n", xlim=c(0, 7042), ylim=c(0,22), xlab="", ylab="")
# # sqrt(chromsize_hg38[1,2]^2 + chromsize_hg38[1,2]^2)/5e4 #7041.551
# for (chr in 1:22){
#     points(x = SNPs_H1_chr_binned_col[[chr]]$plot_pos,
#            y = rep(chr,times=nrow(SNPs_H1_chr_binned_col[[chr]])),
#            pch = 19,
#            col = SNPs_H1_chr_binned_col[[chr]]$colour,cex=0.1)
#     points(x = SNPs_H1_chr_binned_col[[chr]]$plot_pos,
#            y = rep(chr-0.5,times=nrow(SNPs_H1_chr_binned_col[[chr]])),
#            pch = 19,
#            col = SNPs_H1_chr_binned_col[[chr]]$colour_log,cex=0.1)
# }


points(x = SNPs_H1_chr_binned_col[[chr]]$plot_pos,
       y = rep(60,times=nrow(SNPs_H1_chr_binned_col[[chr]])),
       pch = 19,
       col = SNPs_H1_chr_binned_col[[chr]]$colour,cex=0.1)


#cmap='YlGnBu_r' and we normally set vmin=0 and vmax=np.percentile (99)

#test_col = colorRampPalette(c("white","orchid4"))(max_snps_per_bin+1) 
#display.brewer.pal(11, "RdYlBu")
npmi_col_base = brewer.pal(n = 11, name = "RdYlBu")
npmi_col = colorRampPalette(c(rev(npmi_col_base)))(100) 

plot_NPMI_w_both <- function(x,y,range=c(0,d)){
    x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ix <- cut(x, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
#    colorsx <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ix]
    colorsx = npmi_col[ix]
    xxx=matrix(colorsx,nrow=nrow(x),ncol=ncol(x))
    
    iy <- cut(y, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
    colorsy = npmi_col[iy]
    yyy=matrix(colorsy,nrow=nrow(y),ncol=ncol(y))
    
    ## clear lower triangle
    xxx[lower.tri(xxx)] <- yyy[lower.tri(yyy)]
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    #    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    plot(NA, type="n", xlim=range, ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
   # abline(h=0)
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}
testi_1 = NPMI_1[400:500,400:500]
testi_2 = NPMI_2[400:500,400:500]
plot_NPMI_w_both(testi_1,testi_2)
points(x = SNPs_H1_chr_binned_col[[1]]$plot_pos,
       y = rep(60,times=nrow(SNPs_H1_chr_binned_col[[1]])),
       pch = 19,
       col = SNPs_H1_chr_binned_col[[1]]$colour,cex=1)
points(x = SNPs_H1_chr_binned_col[[1]]$plot_pos,
                  y = rep(-60,times=nrow(SNPs_H1_chr_binned_col[[1]])),
                  pch = 19,
                  col = SNPs_H1_chr_binned_col[[1]]$colour_log,cex=1)

# looking good, no on to the Zscores
# this might be trickier, bc we have negative and 0 values in here
# let's see if it just works like that
# 
 zscore_col_base = brewer.pal(n = 11, name = "PuOr")
# zscore_col = colorRampPalette(c(rev(zscore_col_base)))(100) 
# 
# plot(x=1:100,
#      y=rep(1,times=100),
#      ylim = c(1,4),
#      col=zscore_col[1:100],
#      pch=19)
# 
 zscore_col_base[6] = "#FFFFBF"
 zscore_col = colorRampPalette(c(rev(zscore_col_base)))(100) 
# 
# points(x=1:100,
#      y=rep(2,times=100),
#      col=zscore_col[1:100],
#      pch=19)

# zscore_col_base = brewer.pal(n = 11, name = "BrBG")
# zscore_col = colorRampPalette(c(rev(zscore_col_base)))(100) 
# 
# points(x=1:100,
#        y=rep(3,times=100),
#        col=zscore_col[1:100],
#        pch=19)
# 
# zscore_col_base[6] = "#FFFFBF"
# zscore_col = colorRampPalette(c(rev(zscore_col_base)))(100) # looks like mold

# points(x=1:100,
#        y=rep(4,times=100),
#        col=zscore_col[1:100],
#        pch=19)


plot_zscores_w <- function(x,title = "") {
    #x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ii <- cut(x, breaks = seq(-5, 5, len = 100), 
              include.lowest = TRUE)
   # colors <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ii]
    colors = zscore_col[ii]
    xxx=matrix(colors,nrow=nrow(x),ncol=ncol(x))
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab=title, asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}
# testi_dif = zscores_diff[400:500,400:500]
# testi_com = zscores_common[400:500,400:500]
# testi_z1 = zscores_diff_1[400:500,400:500]
# testi_z2 = zscores_diff_2[400:500,400:500]
# plotTriMatrix_zscores_w(zscores_diff)
# plotTriMatrix_zscores_w(zscores_common)
# plotTriMatrix_zcores_w_both(zscores_diff_1,zscores_diff_2)

plot_zscores_w(testi_dif)
plot_zscores_w(testi_com)

plot_zcores_w_both <- function(x,y,range=c(0,d),title = "") {
    #x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ix <- cut(x, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    #colorsx <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ix]
    colorsx = zscore_col[ix]
    xxx=matrix(colorsx,nrow=nrow(x),ncol=ncol(x))
    
    iy <- cut(y, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    #colorsy <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[iy]
    colorsy = zscore_col[iy]
    yyy=matrix(colorsy,nrow=nrow(y),ncol=ncol(y))
    
    ## clear lower triangle
    xxx[lower.tri(xxx)] <- yyy[lower.tri(yyy)]
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    #plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1,  xaxt="none")
    plot(NA, type="n", xlim=range, ylim=c(-60, 60), xlab="", ylab=title, asp=1, yaxt="none", xaxt="none")
    #    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    #abline(h=0)

    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}

plot_zcores_w_both(testi_z1,testi_z2)
points(x = SNPs_H1_chr_binned_col[[1]]$plot_pos,
       y = rep(0,times=nrow(SNPs_H1_chr_binned_col[[1]])),
       pch = 19,
       col = SNPs_H1_chr_binned_col[[1]]$colour,cex=1)

# try to incorporate SNPs at diagonal
plot_zcores_w_both_snp_diag <- function(x,y,snps) {
    #x[which(x<0)]=0
    #x[which(is.na(x))]=0
    x=zscores_diff_1
    y=zscores_diff_2
    #snps=SNPs_H1_chr_binned_col[[chr]]$colour_100lim[1:nrow(zscores_diff_1)]
    
    ix <- cut(x, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    #colorsx <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ix]
    colorsx = zscore_col[ix]
    xxx=matrix(colorsx,nrow=nrow(x),ncol=ncol(x))
    
    iy <- cut(y, breaks = seq(-5,5, len = 100), 
              include.lowest = TRUE)
    #colorsy <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[iy]
    colorsy = zscore_col[iy]
    yyy=matrix(colorsy,nrow=nrow(y),ncol=ncol(y))
    
    ## replace lower triangle
    xxx[lower.tri(xxx)] <- yyy[lower.tri(yyy)]
    ## fill with SNP heatmap at the diagonal
    diag(xxx) = snps
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    #abline(h=0)
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}

plot_zcores_w_both_snp_diag(zscores_diff_1,zscores_diff_2,SNPs_H1_chr_binned_col[[1]]$colour[1:nrow(zscores_diff_1)])

### plot all together
for (chr in 1:22){
    # chr = 6
    print(chr)
    #chr=1
    NPMI_1 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome1/H1.Genome1_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    NPMI_2 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome2/H1.Genome2_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    
    zscores_common = data.matrix(read.table(gzfile(paste(dir,"/ZScore_common/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff = data.matrix(read.table(gzfile(paste(dir,"/ZScore_diff/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_1 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_2 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    
    height_plot = (200/(chromsize_hg38[chr,2]/8000000))*4.3
    pdf(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_chr",chr,".pdf",sep=""),width=200*inch,height=height_plot*inch) #90
    #pdf(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_chr_test.pdf",sep=""),width=200*inch,height=height_plot*inch) 
    par(mar=c(3,2.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(5,1))
    plot_NPMI_w_both(NPMI_1,NPMI_2)
    # points(x = SNPs_H1_chr_binned_col[[chr]]$plot_pos,
    #        y = rep(75,times=nrow(SNPs_H1_chr_binned_col[[chr]])),
    #        pch = 19,
    #        col = SNPs_H1_chr_binned_col[[chr]]$colour,cex=0.1)
    # points(x = SNPs_H1_chr_binned_col[[chr]]$plot_pos,
    #        y = rep(-45,times=nrow(SNPs_H1_chr_binned_col[[chr]])),
    #        pch = 19,
    #        col = SNPs_H1_chr_binned_col[[chr]]$colour_log,cex=0.1)

    plot_zscores_w(zscores_diff)
    plot_zscores_w(zscores_common)
    plot_zcores_w_both_snp_diag(zscores_diff_1,zscores_diff_2,SNPs_H1_chr_binned_col[[chr]]$colour_100lim[1:nrow(zscores_diff_1)])
    # plot_zcores_w_both(zscores_diff_1,zscores_diff_2)
    # points(x = SNPs_H1_chr_binned_col[[chr]]$plot_pos,
    #        y = rep(15,times=nrow(SNPs_H1_chr_binned_col[[chr]])),
    #        pch = 19,
    #        col = SNPs_H1_chr_binned_col[[chr]]$colour,cex=0.1)
    # plot(x=SNPs_H1_chr_binned_col[[chr]]$plot_pos,
    #      y=SNPs_H1_chr_binned_col[[chr]]$value,
    #      col=SNPs_H1_chr_binned_col[[chr]]$colour,
    #      pch=19)
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="number of SNPs per window")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    # it seems here the coordinates are off, recalc!
    dev.off()
    
}

#################################################################################
#################################################################################
#################################################################################
# add imprinted genes : iprinted in human ESC!

#################################################################################
#################################################################################
#################################################################################
# check the bins with higesht SNP density
#nr_bins_chr = unlist(lapply(SNPs_H1_chr_binned,length))
SNPs_H1_binned = unlist(SNPs_H1_chr_binned)
head(rev(sort(SNPs_H1_binned)),30)
# chr6.(3.26e+07,3.265e+07]   chr13.(6.305e+07,6.31e+07]   chr11.(4.83e+07,4.835e+07]   chr16.(1.255e+07,1.26e+07] 
# 356                          344                          340                          337 
# chr9.(1.294e+08,1.2945e+08]      chr8.(4.3e+06,4.35e+06]    chr2.(3.775e+07,3.78e+07]   chr11.(4.835e+07,4.84e+07] 
# 282                          245                          245                          241 
# chr13.(8.575e+07,8.58e+07]    chr3.(9.815e+07,9.82e+07]    chr7.(5.28e+07,5.285e+07]     chr3.(9.9e+07,9.905e+07] 
# 237                          229                          223                          222 
# chr8.(3.35e+06,3.4e+06]        chr5.(8e+06,8.05e+06]    chr17.(6.6e+07,6.605e+07]   chr10.(2.605e+07,2.61e+07] 
# 214                          209                          208                          208 
# chr8.(5.75e+06,5.8e+06]      chr8.(4.9e+06,4.95e+06]   chr13.(5.885e+07,5.89e+07]   chr10.(8.345e+07,8.35e+07] 
# 205                          201                          199                          199 
# chr10.(5.69e+07,5.695e+07]    chr8.(8.84e+07,8.845e+07]    chr5.(8.51e+07,8.515e+07]    chr6.(3.265e+07,3.27e+07] 
# 197                          196                          196                          193 
# chr11.(1.034e+08,1.0345e+08]    chr2.(3.77e+07,3.775e+07]   chr16.(7.64e+07,7.645e+07]     chr11.(5.35e+06,5.4e+06] 
# 192                          191                          190                          189 
# chr6.(6.06e+07,6.065e+07]    chr8.(1.42e+07,1.425e+07] 
# 189                          186 

chr_high_snp = 3
pos_high_snp = 9.9e+07 + 25000
plot_pos_high_snp = SNPs_H1_chr_binned_col[[chr_high_snp]][which(SNPs_H1_chr_binned_col[[chr_high_snp]]$bin_mids==pos_high_snp),2]

plot_NPMI_w_both(NPMI_1[[chr_high_snp]],NPMI_2[[chr_high_snp]],range = c(plot_pos_high_snp-60,plot_pos_high_snp+60))

plot_zcores_w_both(zscores_diff_1[[chr_high_snp]],zscores_diff_2[[chr_high_snp]],range = c(plot_pos_high_snp-60,plot_pos_high_snp+60))

points(x=SNPs_H1_chr_binned_col[[chr_high_snp]]$plot_pos,
       y=rep(0,times=nrow(SNPs_H1_chr_binned_col[[chr_high_snp]])),#SNPs_H1_chr_binned_col$chr6$value,
       col=SNPs_H1_chr_binned_col[[chr_high_snp]]$colour_100lim,
       pch=19,cex=0.5)

# chr6 : 3.26e+07 nothing
# chr13 : 6.305e+07 maybe a TAD change downstream, but nothing in this window
# chr11 : 4.83e+07 only phased data for one copy, can't say
# chr16 : 1.255e+07 looking ok
# chr9 : 1.294e+08 yes there is a change right there, but again more like next to it
# chr8 : 4.3e+06 not too muhc, but wuite many high snp dense regions
# chr2 : 3.775e+0 above it, but notas anchor point
# chr11 : 4.835e+07 just as the bin next to it, see above
# chr13 : 8.575e+07 zscore diff maybe, but real map not so much
# chr3 : 9.815e+07 nope
# chr7 : 5.28e+07 there is difference, but maybe not bc of the SNPs
# chr3 : 9.9e+07  maybe, but not really
#################################################################################
#################################################################################
#################################################################################
# plots for lab meeting
# chr 6, 75Mb - 84Mb
75000000/50000 # 1540
84000000/50000 # 1640
SNPs_H1_binned[1500]
SNPs_H1_binned[1660]

# chr = 6
print(chr)
#chr=1
NPMI_1 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome1/H1.Genome1_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
NPMI_2 = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome2/H1.Genome2_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))

zscores_common = data.matrix(read.table(gzfile(paste(dir,"/ZScore_common/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
zscores_diff = data.matrix(read.table(gzfile(paste(dir,"/ZScore_diff/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
zscores_diff_1 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
zscores_diff_2 = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))


NPMI_1_p=NPMI_1[1500:1660,1500:1660]
NPMI_2_p=NPMI_2[1500:1660,1500:1660]
zscores_common_p = zscores_common[1500:1660,1500:1660]
zscores_diff_p = zscores_diff[1500:1660,1500:1660]
zscores_diff_1_p = zscores_diff_1[1500:1660,1500:1660]
zscores_diff_2_p = zscores_diff_2[1500:1660,1500:1660]

plot_NPMI_w_both(NPMI_1_p,NPMI_2_p)

plot_zscores_w(zscores_diff_p)
plot_zscores_w(zscores_common_p)
plot_zcores_w_both(zscores_diff_1_p,zscores_diff_2_p)
#plot_zcores_w_both_snp_diag(zscores_diff_1_p,zscores_diff_2_p,SNPs_H1_chr_binned_col[[chr]]$colour[1:nrow(zscores_diff_1)])
snps = SNPs_H1_chr_binned_col[[6]][1500:1660,]
points(x=snps$plot_pos-snps$plot_pos[1],
     y=rep(-60,times=nrow(snps)),
     col=snps$colour,
     pch=19,cex=0.5)
plot(x=snps$plot_pos,
     y=snps$value,
     col=snps$colour,
     pch=19)

SNPs_H1_chr_binned_col = lapply(SNPs_H1_chr_binned_col,function(chr){
    chr$pos = bins_all_mean
    return(chr)
})
snps = SNPs_H1_chr_binned_col[[6]][1500:1660,]
plot(x = snps$pos,
     y = snps$value,
     col=snps$colour_log,
     pch=19,cex=0.5)
bar=barplot(snps$value,col=1,ylab = "number of SNPs",ylim=c(-15,160))
points(x=bar,
       y=rep(-10,times=nrow(snps)),
       col=snps$colour_log,
       pch=19,cex=0.5)
points(x=bar,
       y=rep(-5,times=nrow(snps)),
       col=snps$colour,
       pch=19,cex=0.5)

#################################################################################
#################################################################################
#################################################################################
# correlation snp density, diff contacts
# rowsum maps
dir="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/H1_phased/Curated/at50Kb"

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
    NPMI_1[[chr]] = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome1/H1.Genome1_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    NPMI_2[[chr]] = data.matrix(read.table(gzfile(paste(dir,"/NPMI_Genome2/H1.Genome2_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))

    zscores_common[[chr]] = data.matrix(read.table(gzfile(paste(dir,"/ZScore_common/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff[[chr]] = data.matrix(read.table(gzfile(paste(dir,"/ZScore_diff/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_1[[chr]] = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
    zscores_diff_2[[chr]] = data.matrix(read.table(gzfile(paste(dir,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
}

diff_sum = lapply(zscores_diff,function(chr){
    return(colSums(abs(chr),na.rm=T))}) # this gets all contacts, upstream and downstream of the bin

diff_top_sum = lapply(1:22,function(chr){
    x1 = zscores_diff_1[[chr]]
    x2 = zscores_diff_2[[chr]]
    cs1 = colSums(abs(x1),na.rm=T)
    cs2 = colSums(abs(x2),na.rm=T)
    return(cs1+cs2)})

par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = as.numeric(SNPs_H1_chr_binned[[chr]][1:length(diff_sum[[chr]])]),
         y = diff_sum[[chr]],
         pch = 19,
         xlab = "SNP density",
         ylab = "differential contact score sum")
}

max_diff_top_sum = ceiling(max(unlist(lapply(diff_top_sum,max))))
par(mfrow=c(5,5))
for (chr in 1:22){
    plot(x = as.numeric(SNPs_H1_chr_binned[[chr]][1:length(diff_top_sum[[chr]])]),
         y = diff_top_sum[[chr]],
         pch = 19,
         xlab = "SNP density",
         ylab = "differential contact score sum",
         xlim = c(0,max_snps_per_bin),
         ylim = c(0,max_diff_top_sum)) 
}

SNPs_H1_binned_all_chr = unlist(lapply(1:22,function(chr){
    tmp = as.numeric(SNPs_H1_chr_binned[[chr]][1:length(diff_top_sum[[chr]])])
    return(tmp)
})) # 57509
diff_top_sum_all_chr = unlist(diff_top_sum) # 57509
diff_sum_all_chr = unlist(diff_sum) # 57509

par(mfrow=c(1,1))
plot(x = SNPs_H1_binned_all_chr,
     y = diff_top_sum_all_chr,
     pch = 19,
     xlab = "SNP density",
     ylab = "differential contact score top 5% sum",
     xlim = c(0,max_snps_per_bin),
     ylim = c(0,max_diff_top_sum),
     col = adjustcolor("black",alpha.f = 0.1)) 
lines(lowess(SNPs_H1_binned_all_chr, diff_top_sum_all_chr), col = "steelblue4")
lines(lowess(SNPs_H1_binned_all_chr, diff_top_sum_all_chr, f = 1/4), col = "steelblue2")
abline(lm(diff_top_sum_all_chr ~ SNPs_H1_binned_all_chr, data = mtcars), col = "lightblue")

max_diff_sum = max(diff_sum_all_chr)
plot(x = SNPs_H1_binned_all_chr,
     y = diff_sum_all_chr,
     pch = 19,
     xlab = "SNP density",
     ylab = "differential contact score sum",
     xlim = c(0,max_snps_per_bin),
     ylim = c(0,max_diff_sum),
     col = adjustcolor("black",alpha.f = 0.1)) 
lines(lowess(SNPs_H1_binned_all_chr, diff_sum_all_chr), col = "steelblue4")
lines(lowess(SNPs_H1_binned_all_chr, diff_sum_all_chr, f = 1/4), col = "steelblue2")
abline(lm(diff_sum_all_chr ~ SNPs_H1_binned_all_chr, data = mtcars), col = "lightblue")

# there is nothing whatsoever


### calc correlation
cor(SNPs_H1_binned_all_chr,diff_top_sum_all_chr) # 0.1741943
cov(SNPs_H1_binned_all_chr,diff_top_sum_all_chr) # 107.7419

cor(SNPs_H1_binned_all_chr,diff_sum_all_chr) # 0.3732369
cov(SNPs_H1_binned_all_chr,diff_sum_all_chr) # 614.1111


cor(SNPs_H1_binned_all_chr, diff_top_sum_all_chr, method = c("pearson")) # 0.1741943
cor(SNPs_H1_binned_all_chr, diff_top_sum_all_chr, method = c("kendall")) # 0.2180486
cor(SNPs_H1_binned_all_chr, diff_top_sum_all_chr, method = c("spearman")) # 0.3274759
cor(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("pearson")) # 0.3732369
cor(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("kendall")) # 0.2559206
cor(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("spearman")) # 0.3769934

cor.test(SNPs_H1_binned_all_chr, diff_top_sum_all_chr, method = c("pearson")) # 
# Pearson's product-moment correlation
# 
# data:  SNPs_H1_binned_all_chr and diff_top_sum_all_chr
# t = 42.421, df = 57507, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.166258 0.182108
# sample estimates:
# cor 
# 0.1741943

cor.test(SNPs_H1_binned_all_chr, diff_top_sum_all_chr, method = c("kendall")) 
# Kendall's rank correlation tau
# 
# data:  SNPs_H1_binned_all_chr and diff_top_sum_all_chr
# z = 76.166, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# 0.2180486 

cor.test(SNPs_H1_binned_all_chr, diff_top_sum_all_chr, method = c("spearman")) 
# Spearman's rank correlation rho
# 
# data:  SNPs_H1_binned_all_chr and diff_top_sum_all_chr
# S = 2.1319e+13, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.3274759 

cor.test(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("pearson")) 
# Pearson's product-moment correlation
# 
# data:  SNPs_H1_binned_all_chr and diff_sum_all_chr
# t = 96.476, df = 57507, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.3661809 0.3802500
# sample estimates:
# cor 
# 0.3732369

cor.test(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("kendall")) 
# Kendall's rank correlation tau
# 
# data:  SNPs_H1_binned_all_chr and diff_sum_all_chr
# z = 89.428, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau
# 0.2559206

cor.test(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("spearman")) 
# 	Kendall's rank correlation tau
# 
# data:  SNPs_H1_binned_all_chr and diff_sum_all_chr
# z = 89.428, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#     tau 
# 0.2559206 

# library("car")
# scatterplot(diff_sum_all_chr ~ SNPs_H1_binned_all_chr, 
#             smooth = TRUE, grid = FALSE, frame = FALSE,
#             pch = 19,
#             xlab = "SNP density",
#             ylab = "differential contact score sum",
#             xlim = c(0,max_snps_per_bin),
#             ylim = c(0,max_diff_sum),
#             col = adjustcolor("black",alpha.f = 0.1))
# scatterplot(diff_top_sum_all_chr ~ SNPs_H1_binned_all_chr, 
#             smooth = TRUE, grid = FALSE, frame = FALSE,
#             pch = 19,
#             xlab = "SNP density",
#             ylab = "differential contact score sum",
#             xlim = c(0,max_snps_per_bin),
#             ylim = c(0,max_diff_sum),
#             col = adjustcolor("black",alpha.f = 0.1))
# nice but does not really help


### fazit, checking bins with highest SNP denisity by eye, and corelation show no relation of SNP density and diff contact, damn it
#################################################################################
#################################################################################
#################################################################################
# genes iprinted in human ES cells
# link data
# https://www.geneimprint.com/site/genes-by-species
# link paper
# load data
imprinted_genes = read.csv("/Users/jmarkow/Desktop/GAM/imprinted_genes_h1_project.csv",header = T) # all imprinted, 126 genes
# remove thise genes with no very clear plottable location, eg as 8q24.22
imprinted_genes$Location = as.character(imprinted_genes$Location)
keep_id = grepl(x = imprinted_genes$Location, pattern = ":", fixed = TRUE)
imprinted_genes = imprinted_genes[keep_id,]
dim(imprinted_genes) # 122 5
imprinted_genes_pos_tmp  = strsplit(imprinted_genes$Location,"\\D")
imprinted_genes_pos_chr  = unlist(lapply(imprinted_genes_pos_tmp, `[[`, 1))
imprinted_genes_pos_start  = unlist(lapply(imprinted_genes_pos_tmp, `[[`, 2))
imprinted_genes_pos_end  = unlist(lapply(imprinted_genes_pos_tmp, `[[`, 3))
imprinted_genes$chr = as.numeric(imprinted_genes_pos_chr)
imprinted_genes$start = as.numeric(imprinted_genes_pos_start)
imprinted_genes$end = as.numeric(imprinted_genes_pos_end)

imprinted_genes$start[94] = 8669
imprinted_genes$end[94] = 16689

imprinted_genes$start_plot = sqrt(imprinted_genes$start^2 + imprinted_genes$start^2)/5e4
imprinted_genes$end_plot = sqrt(imprinted_genes$end^2 + imprinted_genes$end^2)/5e4

unique(imprinted_genes$Expressed.Allele)
imprinted_genes = imprinted_genes[which(imprinted_genes$Expressed.Allele=="Paternal" | imprinted_genes$Expressed.Allele=="Maternal"),]
imprinted_genes$col = "maroon3" #orchid3
imprinted_genes$col[which(imprinted_genes$Expressed.Allele=="Paternal")] = "royalblue3"

saveRDS(imprinted_genes,"/Users/jmarkow/Desktop/GAM/imprinted_genes_h1_project.rds")
# imprinted_genes = readRDS("/Users/jmarkow/Desktop/GAM/imprinted_genes_h1_project.rds")

# unroll the aliases
imprinted_genes$Gene = as.character(imprinted_genes$Gene)
imprinted_genes$Aliases = as.character(imprinted_genes$Aliases)
imprinted_genes_split = split(imprinted_genes,imprinted_genes$Gene,drop=T)
imprinted_genes_split_unroll = lapply(imprinted_genes_split,function(imp){
    # imp = imprinted_genes_split[[1]]
    # imp = imprinted_genes_split[["LRRTM1"]]
    if(imp$Aliases==""){return(imp)}
    ali = gsub("\\s", "", imp$Aliases)
    ali = unlist(strsplit(ali,","))
    for (a in 1:length(ali)){
        imp[(a+1),] = imp
        imp[(a+1),1] = ali[a]
    }
    return(imp)
})
#imprinted_genes_split_unroll[["DLX5"]]
#imprinted_genes_split_unroll[["LRRTM1"]]
imprinted_genes_unroll =  do.call(rbind,imprinted_genes_split_unroll)
dim(imprinted_genes_unroll) # 472  11
saveRDS(imprinted_genes_unroll,"/Users/jmarkow/Desktop/GAM/imprinted_genes_h1_project_unroll_aliases.rds")

# scp /Users/jmarkow/Desktop/GAM/imprinted_genes_h1_project.csv bih:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1
# scp /Users/jmarkow/Desktop/GAM/imprinted_genes_h1_project.rds bih:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1
# scp /Users/jmarkow/Desktop/GAM/imprinted_genes_h1_project_unroll_aliases.rds bih:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1

###

pdf("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/imprinted_genes_cophased_H1_contact_map.pdf",width=25*inch,height=50*inch)
par(mfrow=c(2,1))
for (gene in 1:nrow(imprinted_genes)){
    print(gene)
    chr_high_snp = imprinted_genes$chr[gene]
    start_gene = imprinted_genes$start_plot[gene]
    end_gene = imprinted_genes$end_plot[gene]
    

    plot_NPMI_w_both(NPMI_1[[chr_high_snp]],NPMI_2[[chr_high_snp]],range = c(start_gene-60,start_gene+60))
    rect(start_gene,-1,end_gene,1,col=imprinted_genes$col[gene])
    points(x=SNPs_H1_chr_binned_col[[chr_high_snp]]$plot_pos,
           y=rep(0,times=nrow(SNPs_H1_chr_binned_col[[chr_high_snp]])),
           col=SNPs_H1_chr_binned_col[[chr_high_snp]]$colour_100lim,
           pch=19,cex=0.5)

    
    plot_zcores_w_both(zscores_diff_1[[chr_high_snp]],zscores_diff_2[[chr_high_snp]],range = c(start_gene-60,start_gene+60))
    rect(start_gene,-1,end_gene,1,col=imprinted_genes$col[gene])
    points(x=SNPs_H1_chr_binned_col[[chr_high_snp]]$plot_pos,
           y=rep(0,times=nrow(SNPs_H1_chr_binned_col[[chr_high_snp]])),
           col=SNPs_H1_chr_binned_col[[chr_high_snp]]$colour_100lim,
           pch=19,cex=0.5)


}
dev.off()
# the matrix is being plotted off centered and I have no clue why
# check around 47-56, all in the same area (that we don't have good data of -.-)

#png("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/imprinted_genes_cophased_H1_contact_map.png",width=2000,height=8000,res=350) #height=4000*nrow(imprinted_genes) #4448000

#png("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/imprinted_genes_cophased_H1_contact_map2.png",width=2000,height=4448000,res=350) #height=4000*nrow(imprinted_genes) #4448000
#par(mfrow=c(2*nrow(imprinted_genes),1))

for (gene in 1:nrow(imprinted_genes)){
    #gene = 2
    print(gene)
    gene_name = imprinted_genes$Gene[gene]
    chr_high_snp = imprinted_genes$chr[gene]
    start_gene = imprinted_genes$start_plot[gene]
    end_gene = imprinted_genes$end_plot[gene]
    
    png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/imprinted_genes/imprinted_genes_cophased_H1_contact_map_",gene_name,".png",sep=""),width=2000,height=4000,res=350) #height=4000*nrow(imprinted_genes) #4448000
    par(mfrow=c(2,1))
    par(mar=c(3,0.1,0.1,0.1)) #bltr # margins

    plot_NPMI_w_both(NPMI_1[[chr_high_snp]],NPMI_2[[chr_high_snp]],range = c(start_gene-60,start_gene+60))
    rect(start_gene,-1,end_gene,1,col=imprinted_genes$col[gene])
    points(x=SNPs_H1_chr_binned_col[[chr_high_snp]]$plot_pos,
           y=rep(0,times=nrow(SNPs_H1_chr_binned_col[[chr_high_snp]])),
           col=SNPs_H1_chr_binned_col[[chr_high_snp]]$colour_100lim,
           pch=19,cex=0.33)
    text(x=start_gene,y=62,labels = paste(gene_name,", chr",chr_high_snp," : ",imprinted_genes$start[gene], " - ",imprinted_genes$end[gene],sep=""))
    
    
    plot_zcores_w_both(zscores_diff_1[[chr_high_snp]],zscores_diff_2[[chr_high_snp]],range = c(start_gene-60,start_gene+60))
    rect(start_gene,-1,end_gene,1,col=imprinted_genes$col[gene])
    points(x=SNPs_H1_chr_binned_col[[chr_high_snp]]$plot_pos,
           y=rep(0,times=nrow(SNPs_H1_chr_binned_col[[chr_high_snp]])),
           col=SNPs_H1_chr_binned_col[[chr_high_snp]]$colour_100lim,
           pch=19,cex=0.33)
    
    dev.off()
}

################################################################################
################################################################################
################################################################################
# next to SNP density, I want to consider gene density as well
# I downloaded 3 different versions: latest is hg38V38, then I also downloaded V29 bc this is the version used for the available RNA-Seq data, and the V29 with ucsc names
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/ucsc_names/gencode.v29.primary_assembly.annotation.gtf.gz


### read in zipped annotation
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("rtracklayer")

gtf <- rtracklayer::import(gzfile("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz"))
# dim(gtf)
# NULL
# > class(gtf)
# [1] "GRanges"
# attr(,"package")
# [1] "GenomicRanges"

gtf_df=as.data.frame(gtf)
# dim(gtf_df)
# [1] 2742734      25
saveRDS(gtf_df,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.rds")

### clean up the gtf file:
# chr 1 - 22, keep X,Y just bc
length(unique(gtf_df$seqnames)) #[1] 47
chr_keep = paste("chr",c(1:22,"X","Y"),sep="")
gtf_df_chrfilter = gtf_df[which(gtf_df$seqnames%in%chr_keep),]
dim(gtf_df_chrfilter) #[1] 2741874      25
length(unique(gtf_df_chrfilter$seqnames)) #24

# clean up columns I don't need, that are uninformative
table(gtf_df_chrfilter$score, exclude = NULL)
     #    <NA>
     # 2741874
table(gtf_df_chrfilter$phase, exclude = NULL)
    #      0       1       2    <NA>
    # 551973  146113  213659 1830129

gtf_df_chrfilter=within(gtf_df_chrfilter,rm("score"))
saveRDS(gtf_df_chrfilter,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.cleaned.rds")
# this still has pseudogenes etc


#### extract genes

gtf_df_genes = gtf_df_chrfilter[which(gtf_df_chrfilter$type=="gene"),]

dim(gtf_df_genes)
#[1] 58684    24
table(gtf_df_chrfilter$gene_type)

# 3prime_overlapping_ncRNA          antisense      bidirectional_promoter_lncRNA                          IG_C_gene                    IG_C_pseudogene
#                      148              46636                               1468                                296                                 33
# IG_D_gene                         IG_J_gene                    IG_J_pseudogene                      IG_pseudogene                          IG_V_gene
#       152                                76                                  9                                  3                               1169
# IG_V_pseudogene                     lincRNA                       macro_lncRNA                              miRNA                           misc_RNA
#             664                       61986                                  3                               5643                               6639
# non_coding           polymorphic_pseudogene               processed_pseudogene               processed_transcript                     protein_coding
#          6                             1796                              32225                              19045                            2502462
# pseudogene                         ribozyme                               rRNA                    rRNA_pseudogene                             scaRNA
#         76                               24                                156                               1500                                147
# scRNA                        sense_intronic                  sense_overlapping                             snoRNA                              snRNA
#     3                                  3577                               1506                               2853                               5700
# sRNA                                    TEC                          TR_C_gene                          TR_D_gene                          TR_J_gene
#   15                                   3250                                 70                                 16                                316
# TR_J_pseudogene                   TR_V_gene                    TR_V_pseudogene   transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#              12                         826                                123                               3834                               4552
# transcribed_unprocessed_pseudogene             translated_processed_pseudogene                 unitary_pseudogene             unprocessed_pseudogene                           vaultRNA
#                              18968                                           8                                482                              13398                                  3

table(gtf_df_genes$gene_type)

# 3prime_overlapping_ncRNA          antisense      bidirectional_promoter_lncRNA                          IG_C_gene                    IG_C_pseudogene
#                       32               5587                                 73                                 14                                  9
# IG_D_gene                         IG_J_gene                    IG_J_pseudogene                      IG_pseudogene                          IG_V_gene
#       37                                 18                                  3                                  1                                144
# IG_V_pseudogene                     lincRNA                       macro_lncRNA                              miRNA                           misc_RNA
#             188                        7635                                  1                               1881                               2213
# non_coding           polymorphic_pseudogene               processed_pseudogene               processed_transcript                     protein_coding
#          2                               41                              10196                                601                              19927
# pseudogene                         ribozyme                               rRNA                    rRNA_pseudogene                             scaRNA
#         18                                8                                 52                                500                                 49
# scRNA                        sense_intronic                  sense_overlapping                             snoRNA                              snRNA
#     1                                   895                                180                                943                               1900
# sRNA                                    TEC                          TR_C_gene                          TR_D_gene                          TR_J_gene
#    5                                   1060                                  6                                  4                                 79
# TR_J_pseudogene                   TR_V_gene                    TR_V_pseudogene   transcribed_processed_pseudogene     transcribed_unitary_pseudogene
#               4                         106                                 33                                481                                124
# transcribed_unprocessed_pseudogene             translated_processed_pseudogene                 unitary_pseudogene             unprocessed_pseudogene                           vaultRNA
#                                886                                           2                                 95                               2649                                  1
# https://www.gencodegenes.org/pages/biotypes.html
# note: there are no lncRNA in here:
# the website says:
# lncRNA - Generic long non-coding RNA biotype that replaced the following biotypes: 
# 3prime_overlapping_ncRNA, antisense, bidirectional_promoter_lncRNA, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic and sense_overlapping.

# so I might need to edit this, orcheck the version of the annotation / ref genome that I am using here, see if it fits.

# TEC - To be Experimentally Confirmed. 


table(gtf_df_genes$gene_type=="protein_coding")#/nrow(gtf_df_genes)*100
# FALSE  TRUE
# 38757 19927
# FALSE     TRUE
# 66.04356 33.95644



saveRDS(gtf_df_chrfilter,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.cleaned.rds")

####

annoV38 <- rtracklayer::import(gzfile("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v38.annotation.gtf.gz"))
annoV38_df = as.data.frame(annoV38)
# dim(annoV38_df)
# [1] 3150424      26
saveRDS(annoV38_df,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v38.annotation.gtf.rds")

length(unique(annoV38_df$seqnames)) #25
annoV38_df = annoV38_df[which(annoV38_df$seqnames%in%chr_keep),]
dim(annoV38_df) #[1] 3150281      26

table(annoV38_df$gene_type)
# IG_C_gene                     IG_C_pseudogene                          IG_D_gene                          IG_J_gene                    IG_J_pseudogene
#       281                                  33                                152                                 76                                  9
# IG_pseudogene                       IG_V_gene                    IG_V_pseudogene                             lncRNA                              miRNA
#             3                            1134                                661                             241995                               5637
# misc_RNA               polymorphic_pseudogene               processed_pseudogene                     protein_coding                         pseudogene
#     6636                                 1818                              32114                            2800821                                 64
# ribozyme                                 rRNA                    rRNA_pseudogene                             scaRNA                              scRNA
#       24                                  141                               1491                                147                                  3
# snoRNA                                  snRNA                               sRNA                                TEC                          TR_C_gene
#   2829                                   5703                                 15                               3218                                 70
# TR_D_gene                           TR_J_gene                    TR_J_pseudogene                          TR_V_gene                    TR_V_pseudogene
#        16                                 316                                 12                                823                                123
# transcribed_processed_pseudogene                  transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene    translated_processed_pseudogene  translated_unprocessed_pseudogene
#                             4307                                            5376                              20542                                  8                                  9
# unitary_pseudogene                                        unprocessed_pseudogene                          vault_RNA
#                554                                                         13117                                  3

annoV38_df_genes = annoV38_df[which(annoV38_df$type=="gene"),]

dim(annoV38_df_genes)
# [1] 60612    26
table(annoV38_df_genes$gene_type)
# IG_C_gene                     IG_C_pseudogene                          IG_D_gene                          IG_J_gene                    IG_J_pseudogene
#        14                                   9                                 37                                 18                                  3
# IG_pseudogene                       IG_V_gene                    IG_V_pseudogene                             lncRNA                              miRNA
#             1                             145                                187                              16888                               1879
# misc_RNA               polymorphic_pseudogene               processed_pseudogene                     protein_coding                         pseudogene
#     2212                                   49                              10163                              19942                                 15
# ribozyme                                 rRNA                    rRNA_pseudogene                             scaRNA                              scRNA
#        8                                   47                                497                                 49                                  1
# snoRNA                                  snRNA                               sRNA                                TEC                          TR_C_gene
#    943                                   1901                                  5                               1056                                  6
# TR_D_gene                           TR_J_gene                    TR_J_pseudogene                          TR_V_gene                    TR_V_pseudogene
#         4                                  79                                  4                                106                                 33
# transcribed_processed_pseudogene                  transcribed_unitary_pseudogene transcribed_unprocessed_pseudogene    translated_processed_pseudogene  translated_unprocessed_pseudogene
#                              502                                             143                                950                                  2                                  1
# unitary_pseudogene                                        unprocessed_pseudogene                          vault_RNA
#                 98                                                          2614                                  1


table(annoV38_df_genes$gene_type=="protein_coding")/nrow(annoV38_df_genes)*100
# FALSE  TRUE
# 40670 19942 # 15 more compared to V29 UCSC
# FALSE     TRUE
# 67.09892 32.90108
table(annoV38_df_genes$gene_type=="protein_coding" | annoV38_df_genes$gene_type=="lncRNA")/nrow(annoV38_df_genes)*100
# FALSE     TRUE
# 39.23645 60.76355

table(annoV38_df_genes$tag, exclude = NULL)
#      fragmented_locus             ncRNA_host      overlapping_locus                    PAR         pseudo_consens reference_genome_error              retrogene         semi_processed                   <NA>
#                     2                   1844                  10583                     44                   7347                     14                    388                      8                  40382

saveRDS(annoV38_df_genes[,1:15],"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v38.annotation.gtf.allgenes.rds")

annoV38_df_genes_lncRNA_protcod = annoV38_df_genes[which(annoV38_df_genes$gene_type=="protein_coding" | annoV38_df_genes$gene_type=="lncRNA"),1:15] #[1] 36830    15
saveRDS(annoV38_df_genes_lncRNA_protcod,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v38.annotation.gtf.lncRNA_protcod_genes.rds")

###
annoV38_df_transcripts = annoV38_df[which(annoV38_df$type == "transcript"),] #236975     26
annoV38_df_transcripts_lncRNA_protcod = annoV38_df_transcripts[which((annoV38_df_transcripts$gene_type=="protein_coding" | annoV38_df_transcripts$gene_type=="lncRNA")&(annoV38_df_transcripts$transcript_type=="protein_coding" | annoV38_df_transcripts$transcript_type=="lncRNA")),] # 133780     26

saveRDS(annoV38_df_transcripts_lncRNA_protcod,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v38.annotation.gtf.lncRNA_protcod_transcripts.rds")


###### V29


####

annoV29 <- rtracklayer::import(gzfile("/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.annotation.gtf.gz"))
annoV29_df = as.data.frame(annoV29)
dim(annoV29_df) #2742017      25
saveRDS(annoV29_df,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.annotation.gtf.rds")

length(unique(annoV29_df$seqnames)) #
annoV29_df = annoV29_df[which(annoV29_df$seqnames%in%chr_keep),]
dim(annoV29_df) #[1] 2741874      25

annoV29_df_genes = annoV29_df[which(annoV29_df$type=="gene"),]
dim(annoV29_df_genes) #58684    25

table(annoV29_df_genes$gene_type=="protein_coding")/nrow(annoV29_df_genes)*100
#    FALSE     TRUE
# 66.04356 33.95644
#table(annoV29_df_genes$gene_type=="protein_coding" | annoV29_df_genes$gene_type=="lncRNA")/nrow(annoV29_df_genes)*100
# there is no single "lncRNA" in this dataset

table(annoV29_df_genes$tag, exclude = NULL)
#      fragmented_locus             ncRNA_host                 orphan      overlapping_locus                    PAR         pseudo_consens reference_genome_error              retrogene         semi_processed                   <NA>
#                     4                   1681                      1                  10391                     45                   7361                     17                    381                      6                  38797

saveRDS(annoV29_df_genes[,1:14],"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.annotation.gtf.allgenes.rds")

# I didn't the next block, bc there is no single lnc RNA type in thsi version
# annoV29_df_genes_lncRNA_protcod = annoV29_df_genes[which(annoV29_df_genes$gene_type=="protein_coding" | annoV29_df_genes$gene_type=="lncRNA"),1:14] #[1] 36830    15
# saveRDS(annoV29_df_genes_lncRNA_protcod,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.annotation.gtf.lncRNA_protcod_genes.rds")
# 
# ###
# annoV29_df_transcripts = annoV29_df[which(annoV29_df$type == "transcript"),] #236975     26
# annoV29_df_transcripts_lncRNA_protcod = annoV29_df_transcripts[which((annoV29_df_transcripts$gene_type=="protein_coding" | annoV29_df_transcripts$gene_type=="lncRNA")&(annoV29_df_transcripts$transcript_type=="protein_coding" | annoV29_df_transcripts$transcript_type=="lncRNA")),] #
# 
# saveRDS(annoV29_df_transcripts_lncRNA_protcod,"/fast/groups/ag_schwarz/Projects/project-gam/H1/data/anno/gencode.v29.annotation.gtf.lncRNA_protcod_transcripts.rds")


####### bin genes, all and lnc/protcod only
# mount the clustert, do it locally
dir="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/"

annoV38_df_genes_lncRNA_protcod = readRDS(paste(dir,"/H1/data/anno/gencode.v38.annotation.gtf.lncRNA_protcod_genes.rds",sep="")) # 36830    15
annoV38_df_genes = readRDS(paste(dir,"/H1/data/anno/gencode.v38.annotation.gtf.allgenes.rds",sep="")) # 60612    15
chr_keep = paste("chr",c(1:22,"X","Y"),sep="")

# take code from SNP density and use for gene density

annoV38_df_genes_chr = split(annoV38_df_genes,annoV38_df_genes$seqnames)
annoV38_df_genes_chr = annoV38_df_genes_chr[chr_keep]

bins_all =  seq(0,250000000, 50000)
bins_all_mean = seq(25000,250000000, 50000)
bin_pos = sqrt(bins_all_mean^2 + bins_all_mean^2)/5e4

annoV38_df_genes_chr_binned = lapply(annoV38_df_genes_chr,function(chr){
    bins = cut(chr$start, breaks = bins_all)
    return(table(bins))
})
max_genes_per_bin = max(unlist(lapply(annoV38_df_genes_chr_binned,max)))
# check how the SNP dens is distributed
hist(unlist(annoV38_df_genes_chr_binned),breaks=c(0:max_genes_per_bin),xlim=c(1,max_genes_per_bin),ylim=c(0,10000))

# colour code for the SNPs density
gene_dens_col = colorRampPalette(c("white","maroon"))(max_genes_per_bin+1) # do this max over all so it is comparable
gene_dens_col_log = colorRampPalette(c("white","maroon"))(ceiling(log(max_genes_per_bin))+1)# I put the +1 outsiede the log, bc later, I want to assign all values a col, and then assign 0 the white

# gene_dens_col_base = brewer.pal(n = 11, name = "BrBG")
# gene_dens_col_base = gene_dens_col_base[6:11]
# gene_dens_col_base[1] = "white"
# gene_dens_col = colorRampPalette(c(gene_dens_col_base))(max_genes_per_bin+1)
# gene_dens_col_log = colorRampPalette(c(gene_dens_col_base))(ceiling(log(max_genes_per_bin))+1)# I put the +1 outsiede the log, bc later, I want to assign all values a col, and then assign 0 the white

# create dataframe with x position and colour (gene dens) for each chromosome
color_gene_dense_df = data.frame(val=0:max(max_genes_per_bin),colorname = gene_dens_col,stringsAsFactors = F)
color_gene_dense_log_df = data.frame(val=0:ceiling(log(max(max_genes_per_bin))),colorname = gene_dens_col_log,stringsAsFactors = F)

# since the majority of bins has very few genes, but the log does not really cut it, I decided to make a gradient up until 100 genes and set everything above to the max colour 
gene_dens_col_10lim = colorRampPalette(c(gene_dens_col))(11)
color_gene_dense_10lim_df = data.frame(val=0:max(max_genes_per_bin),colorname = c(gene_dens_col_10lim,rep(gene_dens_col_10lim[length(gene_dens_col_10lim)],times=(max(max_genes_per_bin)-10))),stringsAsFactors = F)


annoV38_df_genes_chr_binned_col = lapply(annoV38_df_genes_chr_binned,function(chr){
    # chr = annoV38_df_genes_chr_binned[[1]]
    chr_df = data.frame(val=chr,order=c(1:length(chr)))
    chr_df$val_log = ceiling(log(chr_df$val.Freq+1))
    # hist(chr_df$val_log,breaks=seq((-0.5),6.5,1))
    chr_df$val_log[which(chr_df$val.Freq==0)] = 0
    # hist(chr_df$val_log,breaks=seq((-0.5),6.5,1))
    
    df = merge(chr_df, color_gene_dense_df,by.x="val.Freq",by.y="val")
    df = merge(df, color_gene_dense_log_df,by.x="val_log",by.y="val")
    df = merge(df, color_gene_dense_10lim_df,by.x="val.Freq",by.y="val")
    
    df_sort = df[order(df$order),]
    df_sort$plot_pos = bin_pos
    df_sort$bin_mids = bins_all_mean
    
    df_sort = df_sort[,c(4,8,9,3,1,5,7,2,6)]
    colnames(df_sort) = c("order","plot_pos","bin_mids","bin","value","colour","colour_10lim","value_log","colour_log")
    return(df_sort)
})

annoV38_df_genes_chr_binned_col = lapply(annoV38_df_genes_chr_binned_col,function(chr){
    chr$pos = bins_all_mean
    return(chr)
})

plot(annoV38_df_genes_chr_binned_col$chr1$bin_mids,annoV38_df_genes_chr_binned_col$chr1$value,col=annoV38_df_genes_chr_binned_col$chr1$colour_10lim,pch=19)
# plot(annoV38_df_genes_chr_binned_col$chr1$bin_mids,annoV38_df_genes_chr_binned_col$chr1$value,col=annoV38_df_genes_chr_binned_col$chr1$colour,pch=19)
# plot(annoV38_df_genes_chr_binned_col$chr1$bin_mids,annoV38_df_genes_chr_binned_col$chr1$value_log,col=annoV38_df_genes_chr_binned_col$chr1$colour_log,pch=19)

# bins_all =  seq(0,250000000, 50000)
# snpchr = SNPs_H1[substring(SNPs_H1$CHROM, 4)==chr,2]
# 
# bins <- cut(snpchr, breaks = bins_all)
# length(bins)
# length(snpchr)
# head(bins)
# head(snpchr)
# snp_dens=table(bins)


################################################################################
# plot all chr again, with SNP density, gene density, impirnted genes, as png so there is proper positioning

imprinted_genes_chr = split(imprinted_genes,imprinted_genes$chr)
names(imprinted_genes_chr) = paste("chr",names(imprinted_genes_chr),sep="")



for (chr in 1:22){
#   chr=22
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
   
   
    png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_chr",chr,".png",sep=""),width=width_plot,height=1000)#,res=350) #1000 #
        par(mar=c(3,2.1,0.1,0.1)) #bltr # margins
        par(mfrow=c(6,1))
        plot_NPMI_w_both(NPMI_1_chr,NPMI_2_chr)
    
        plot_zscores_w(zscores_diff_chr)
        plot_zscores_w(zscores_common_chr)
        plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr)
    
        plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
             y=SNPs_H1_chr_binned_col[[chr]]$value,
             col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
             pch=19,xlim=c(0,chromsize_hg38[chr,2]),
             xaxt="none",xlab="",
             ylab="number of SNPs per window")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
        
        plot(x=annoV38_df_genes_chr_binned_col[[chr]]$bin_mids,
             y=annoV38_df_genes_chr_binned_col[[chr]]$value,
             col=annoV38_df_genes_chr_binned_col[[chr]]$colour_10lim,
             pch=19,xlim=c(0,chromsize_hg38[chr,2]),
             xaxt="none",xlab="",
             ylab="number of genes per window")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
        
        chrname = paste("chr",chr,sep="")
        if( chrname%in% names(imprinted_genes_chr)){
            rect(imprinted_genes_chr[[chrname]]$start,-0.5,imprinted_genes_chr[[chrname]]$end,0.5,col=imprinted_genes_chr[[chrname]]$col)
        }
    dev.off()
}


plot(annoV38_df_genes_chr_binned_col$chr1$bin_mids,annoV38_df_genes_chr_binned_col$chr1$value,col=annoV38_df_genes_chr_binned_col$chr1$colour_10lim,pch=19)

### compare V29 and V38 gene annotations
dir="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/"

annoV29_df_genes = readRDS(paste(dir,"/H1/data/anno/gencode.v29.annotation.gtf.allgenes.rds",sep="")) #  58684    14
annoV38_df_genes = readRDS(paste(dir,"/H1/data/anno/gencode.v38.annotation.gtf.allgenes.rds",sep="")) #  60612    15
anno_ver_comp_genes = merge(annoV38_df_genes,annoV29_df_genes,by="gene_id") # 37281    28  # wow, that is not many genes left considering what we put in!
anno_ver_comp_genes_name = merge(annoV38_df_genes,annoV29_df_genes,by="gene_name") #  38048    28 #  a few more, but wow
head(anno_ver_comp_genes_name[which(anno_ver_comp_genes_name$gene_id.x!=anno_ver_comp_genes_name$gene_id.y),]) # ok, seems like gene_id includes "versions" as in V38 some are .01 higher than in V29
# let's merge by gene_id without .version
annoV29_df_genes$gene_id_full =  sapply(strsplit(annoV29_df_genes$gene_id,"\\."), `[`, 1)
annoV38_df_genes$gene_id_full =  sapply(strsplit(annoV38_df_genes$gene_id,"\\."), `[`, 1)
anno_ver_comp_genes_id_basic = merge(annoV38_df_genes,annoV29_df_genes,by="gene_id_full") # 58433    30 # that looks much better
nrow(anno_ver_comp_genes_id_basic[which(anno_ver_comp_genes_id_basic$gene_name.x!=anno_ver_comp_genes_id_basic$gene_name.y),]) # there are quite some gene name changes going on here! 20983

# ok, check how many start and end pos change 
table(anno_ver_comp_genes_id_basic$start.x==anno_ver_comp_genes_id_basic$start.y)
# FALSE  TRUE 
# 9995 48438 
table(anno_ver_comp_genes_id_basic$end.x==anno_ver_comp_genes_id_basic$end.y)
# FALSE  TRUE 
# 9950 48483 
table((anno_ver_comp_genes_id_basic$end.x==anno_ver_comp_genes_id_basic$end.y)&(anno_ver_comp_genes_id_basic$start.x==anno_ver_comp_genes_id_basic$start.y))
# FALSE  TRUE 
# 14779 43654 

plot(x = anno_ver_comp_genes_id_basic$start.x,
     y = anno_ver_comp_genes_id_basic$start.y,
     pch=19,
     col=adjustcolor("black",alpha.f = 0.1))

plot(x = anno_ver_comp_genes_id_basic$end.x,
     y = anno_ver_comp_genes_id_basic$end.y,
     pch=19,
     col=adjustcolor("black",alpha.f = 0.1))

anno_ver_comp_genes_id_basic$start_diff = anno_ver_comp_genes_id_basic$start.x - anno_ver_comp_genes_id_basic$start.y
anno_ver_comp_genes_id_basic$end_diff = anno_ver_comp_genes_id_basic$end.x - anno_ver_comp_genes_id_basic$end.y
max(anno_ver_comp_genes_id_basic$start_diff)
max(anno_ver_comp_genes_id_basic$end_diff)
anno_ver_comp_genes_id_basic[which.max(anno_ver_comp_genes_id_basic$end_diff),] # that is a chr jump
table(anno_ver_comp_genes_id_basic$seqnames.x==anno_ver_comp_genes_id_basic$seqnames.y)
# FALSE  TRUE 
# 89 58344
anno_ver_comp_genes_id_basic_same_chr = anno_ver_comp_genes_id_basic[which(anno_ver_comp_genes_id_basic$seqnames.x==anno_ver_comp_genes_id_basic$seqnames.y),]
plot(x = anno_ver_comp_genes_id_basic_same_chr$start.x,
     y = anno_ver_comp_genes_id_basic_same_chr$start.y,
     pch=19,
     col=adjustcolor("black",alpha.f = 0.1))

plot(x = anno_ver_comp_genes_id_basic_same_chr$end.x,
     y = anno_ver_comp_genes_id_basic_same_chr$end.y,
     pch=19,
     col=adjustcolor("black",alpha.f = 0.1))

#hist(anno_ver_comp_genes_id_basic_same_chr$start_diff[which(anno_ver_comp_genes_id_basic_same_chr$start_diff!=0)],breaks=seq(min(anno_ver_comp_genes_id_basic_same_chr$start_diff)-1000,max(anno_ver_comp_genes_id_basic_same_chr$start_diff)+1000,1000))
hist(log10(abs(anno_ver_comp_genes_id_basic_same_chr$start_diff))[which(anno_ver_comp_genes_id_basic_same_chr$start_diff!=0)]*sign(anno_ver_comp_genes_id_basic_same_chr$start_diff[which(anno_ver_comp_genes_id_basic_same_chr$start_diff!=0)]),breaks=seq(-6,8,0.1),
     xlab = "log10(position difference)",
     main = "start pos difference between V38 and V29 of hg38")

hist(log10(abs(anno_ver_comp_genes_id_basic_same_chr$end_diff))[which(anno_ver_comp_genes_id_basic_same_chr$end_diff!=0)]*sign(anno_ver_comp_genes_id_basic_same_chr$end_diff[which(anno_ver_comp_genes_id_basic_same_chr$end_diff!=0)]),breaks=seq(-6.5,8,0.1),
     xlab = "log10(position difference)",
     main = "end pos difference between V38 and V29 of hg38")
# seems like it is mirrored? what could that mean?
anno_ver_comp_genes_id_basic_same_chr$pos_diff = anno_ver_comp_genes_id_basic_same_chr$start_diff - anno_ver_comp_genes_id_basic_same_chr$end_diff

#hist(log10(abs(anno_ver_comp_genes_id_basic_same_chr$pos_diff))[which(anno_ver_comp_genes_id_basic_same_chr$pos_diff!=0)]*sign(anno_ver_comp_genes_id_basic_same_chr$pos_diff[which(anno_ver_comp_genes_id_basic_same_chr$pos_diff!=0)]),breaks=seq(-6.5,8,0.1),
hist(log10(abs(anno_ver_comp_genes_id_basic_same_chr$pos_diff))[which(anno_ver_comp_genes_id_basic_same_chr$pos_diff!=0)]*sign(anno_ver_comp_genes_id_basic_same_chr$pos_diff[which(anno_ver_comp_genes_id_basic_same_chr$pos_diff!=0)]),breaks=seq(-6.5,8,0.1),
     xlab = "log10(position difference)",
     main = "pos difference between V38 and V29 of hg38")

##############################################################################
### check correlation between gene dense regions and top differential contacts

annoV38_df_genes_chr_binned_all_chr = unlist(lapply(1:22,function(chr){
    tmp = as.numeric(annoV38_df_genes_chr_binned[[chr]][1:length(diff_top_sum[[chr]])])
    return(tmp)}))

    
par(mfrow=c(1,1))
plot(x = annoV38_df_genes_chr_binned_all_chr,
     y = diff_top_sum_all_chr,
     pch = 19,
     xlab = "gene density",
     ylab = "differential contact score top 5% sum",
     xlim = c(0,max_genes_per_bin),
     ylim = c(0,max_diff_top_sum),
     col = adjustcolor("black",alpha.f = 0.1)) 
lines(lowess(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr), col = "steelblue4")
lines(lowess(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr, f = 1/4), col = "steelblue2")
abline(lm(diff_top_sum_all_chr ~ annoV38_df_genes_chr_binned_all_chr), col = "lightblue")

plot(x = annoV38_df_genes_chr_binned_all_chr,
     y = diff_sum_all_chr,
     pch = 19,
     xlab = "gene density",
     ylab = "differential contact score sum",
     xlim = c(0,max_genes_per_bin),
     ylim = c(0,max_diff_sum),
     col = adjustcolor("black",alpha.f = 0.1)) 
lines(lowess(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr), col = "steelblue4")
lines(lowess(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr, f = 1/4), col = "steelblue2")
abline(lm(diff_sum_all_chr ~ annoV38_df_genes_chr_binned_all_chr), col = "lightblue")

# there is nothing whatsoever


### calc correlation
cor(annoV38_df_genes_chr_binned_all_chr,diff_top_sum_all_chr) # -0.01870123
cov(annoV38_df_genes_chr_binned_all_chr,diff_top_sum_all_chr) # -0.6382972

cor(annoV38_df_genes_chr_binned_all_chr,diff_sum_all_chr) #  -0.02510523
cov(annoV38_df_genes_chr_binned_all_chr,diff_sum_all_chr) # -2.279441


cor(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr, method = c("pearson")) # -0.01870123
cor(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr, method = c("kendall")) # 0.01606614
cor(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr, method = c("spearman")) # 0.02147341
cor(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr, method = c("pearson")) # -0.02510523
cor(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr, method = c("kendall")) # -0.0006712014
cor(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr, method = c("spearman")) # 0.0004004468

cor.test(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr, method = c("pearson"))  
# Pearson's product-moment correlation
# 
# data:  annoV38_df_genes_chr_binned_all_chr and diff_top_sum_all_chr
# t = -4.4855, df = 57507, p-value = 7.29e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.02687013 -0.01052984
# sample estimates:
# cor 
# -0.01870123 

cor.test(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr, method = c("kendall")) 
# Kendall's rank correlation tau
# 
# data:  annoV38_df_genes_chr_binned_all_chr and diff_top_sum_all_chr
# z = 5.0172, p-value = 5.244e-07
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# 0.01606614 

# cor.test(annoV38_df_genes_chr_binned_all_chr, diff_top_sum_all_chr, method = c("spearman")) 
# Spearman's rank correlation rho
# 
# data:  annoV38_df_genes_chr_binned_all_chr and diff_top_sum_all_chr
# S = 3.1019e+13, p-value = 2.604e-07
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.02147341 

cor.test(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr, method = c("pearson")) 
# Pearson's product-moment correlation
# 
# data:  annoV38_df_genes_chr_binned_all_chr and diff_sum_all_chr
# t = -6.0223, df = 57507, p-value = 1.73e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.03327141 -0.01693570
# sample estimates:
# cor 
# -0.02510523

cor.test(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr, method = c("kendall")) 
# Kendall's rank correlation tau
# 
# data:  annoV38_df_genes_chr_binned_all_chr and diff_sum_all_chr
# z = -0.20968, p-value = 0.8339
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# -0.0006712014

cor.test(annoV38_df_genes_chr_binned_all_chr, diff_sum_all_chr, method = c("spearman")) 
# Spearman's rank correlation rho
# 
# data:  annoV38_df_genes_chr_binned_all_chr and diff_sum_all_chr
# S = 3.1687e+13, p-value = 0.9235
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.0004004468 

################################################################################
### total expression
################################################################################
dir="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster"

total_RNA_red_auto = readRDS(paste(dir,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/H1_total_expression_ENCFF675NTU_ENCFF379NOY.rds",sep="")) #[1] 199269     18
total_RNA_red_auto_protcod = total_RNA_red_auto[which(total_RNA_red_auto$gene_type == "protein_coding" | total_RNA_red_auto$transcript_type == "protein_coding" ),]

total_RNA_red_auto$mids = total_RNA_red_auto$start + (total_RNA_red_auto$end - total_RNA_red_auto$start)/2
total_RNA_red_auto$log10_TPM_transcript_ENCFF675NTU = log10(total_RNA_red_auto$TPM_transcript_ENCFF675NTU +1)
total_RNA_red_auto$log10_TPM_transcript_ENCFF379NOY = log10(total_RNA_red_auto$TPM_transcript_ENCFF379NOY +1)

total_RNA_red_auto_chr = split(total_RNA_red_auto,total_RNA_red_auto$chr,drop=T)

for (chr in 1:22){
    #   chr=22
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    
    
    png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_total_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_total_expression_chr",chr,".png",sep=""),width=width_plot,height=1167)#,res=350) #1000 #
    par(mar=c(3,2.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(7,1))
    plot_NPMI_w_both(NPMI_1_chr,NPMI_2_chr)
    
    plot_zscores_w(zscores_diff_chr)
    plot_zscores_w(zscores_common_chr)
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr)
    
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="number of SNPs per window")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    plot(x=annoV38_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV38_df_genes_chr_binned_col[[chr]]$value,
         col=annoV38_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="number of genes per window")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    plot(x=total_RNA_red_auto_chr[[chr]]$mids,
         y=total_RNA_red_auto_chr[[chr]]$log10_TPM_transcript_ENCFF675NTU,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,-0.5,imprinted_genes_chr[[chrname]]$end,0.5,col=imprinted_genes_chr[[chrname]]$col)
    }
    dev.off()
}

# merge expression with imprinted genes 

imprinted_genes_expression = merge(imprinted_genes,total_RNA_red_auto,by.x="Gene",by.y="gene_name",all.x = T,all.y = F)
imprinted_genes_expressed = imprinted_genes_expression[which(imprinted_genes_expression$TPM_transcript_ENCFF379NOY>=1),] #179 transcripts
length(unique(imprinted_genes_expressed$Gene)) # 34 of 112
imprinted_genes_expressed_split = split(imprinted_genes_expressed,imprinted_genes_expressed$Gene,drop=T)
imprinted_genes_expressed_split_maxexp = lapply(imprinted_genes_expressed_split,function(g){
    return(g[which.max(g$TPM_transcript_ENCFF675NTU),])
})
imprinted_genes_expressed_maxexp = do.call(rbind,imprinted_genes_expressed_split_maxexp)

################################################################################
### ASE and total expression
################################################################################
anno_total_ase_imp = readRDS(paste(dir,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/ASE_H1/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY.rds",sep=""))
anno_total_ase_imp$mids = anno_total_ase_imp$start_transcript + (anno_total_ase_imp$end_transcript - anno_total_ase_imp$start_transcript)/2
anno_total_ase_imp$log10_TPM_transcript_ENCFF675NTU = log10(anno_total_ase_imp$TPM_transcript_ENCFF675NTU +1)
anno_total_ase_imp$log10_TPM_transcript_ENCFF379NOY = log10(anno_total_ase_imp$TPM_transcript_ENCFF379NOY +1)
unique(anno_total_ase_imp$chrom)
anno_total_ase_imp_chr = split(anno_total_ase_imp,anno_total_ase_imp$chrom,drop=T)
names(anno_total_ase_imp_chr)

for (chr in 1:22){
    #   chr=22
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    
    
    png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr",chr,".png",sep=""),width=width_plot,height=1500)#,res=350) #1000 #
    par(mar=c(3,2.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(9,1))
    plot_NPMI_w_both(NPMI_1_chr,NPMI_2_chr)
    
    plot_zscores_w(zscores_diff_chr)
    plot_zscores_w(zscores_common_chr)
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr)
    
    # SNP density
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="number of SNPs per window")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # gene density
    plot(x=annoV38_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV38_df_genes_chr_binned_col[[chr]]$value,
         col=annoV38_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="number of genes per window")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    #par(mfrow=c(3,1))
    # total expression
    plot(x=anno_total_ase_imp_chr[[chr]]$mids,
         y=anno_total_ase_imp_chr[[chr]]$log10_TPM_transcript_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_chr[[chr]]$ENCFF379NOY_exp_ase_lfc+1,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # ASE ratio
    plot(x=anno_total_ase_imp_chr[[chr]]$mids,
         y=anno_total_ase_imp_chr[[chr]]$ASE_ratio_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_chr[[chr]]$ENCFF379NOY_exp_ase_lfc+1,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="ASE")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # LFC
    plot(x=anno_total_ase_imp_chr[[chr]]$mids,
         y=anno_total_ase_imp_chr[[chr]]$log2foldchange_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_chr[[chr]]$ENCFF379NOY_exp_ase_lfc+1,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="LFC")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # imprinted genes
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,-0.5,imprinted_genes_chr[[chrname]]$end,0.5,col=imprinted_genes_chr[[chrname]]$col)
    }
    dev.off()
}

################################################################################
### plot the AB compartments as well 
################################################################################
comp=read.csv("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/compartments/210222.Compartments.AB.H1.phased.first.PCA.csv",header = T, sep="\t")

### preprocess ----------------------------------------------------------------
# check size
dim(comp) # 11749     7

# sort by chromosome and position
comp = comp[order(comp$start),]
comp = comp[order(as.numeric(substring(comp$chrom, 4))),]
comp$chrom = factor(comp$chrom, levels = unique(comp$chrom))

# add 1bp to every start position so that intervals are exclusive and not overlapping
comp$start = comp$start + 1 
comp$mid = comp$start + (comp$end - comp$start)/2

# check chromosome set
head(comp)
unique(comp$chrom)
# no need to exclude sex or mitochondria chromosomes

# combine called compartments of both chromosomes into one value
comp$combined_comp_raw = paste(comp$Genome1_AB1, comp$Genome2_AB1, sep = '')
comp$combined_comp_raw = factor(comp$combined_comp_raw, levels = c("AA", "AB", "BA", "BB"))

table(comp$combined_comp_raw) / nrow(comp) * 100
#       AA       AB       BA       BB
# 32.75172 18.59733 21.46566 27.18529

# reduce he compartment calls to AA / AB / BB
comp$combined_comp = comp$combined_comp_raw
comp$combined_comp[which(comp$combined_comp == "BA")] = "AB"
comp$combined_comp = factor(comp$combined_comp, levels = c("AA", "AB", "BB"))

# how many genomic windows are differentially called?
comp_dist = table(comp$combined_comp) / nrow(comp) * 100
# AA       AB       BB 
# 32.75172 40.06298 27.18529

comp_chr = split(comp,comp$chrom)
plot(x=comp_chr[[chr]]$mid,
     y=comp_chr[[chr]]$Genome1_score1,
     type = "l",
     xlim=c(0,chromsize_hg38[chr,2]),
     #ylim = c(-7,7),
     xaxt="none",xlab="",
     ylab="genome score")
lines(x=comp_chr[[chr]]$mid,
      y=comp_chr[[chr]]$Genome2_score1,col="grey80")
points(x=comp_chr[[chr]]$mid,
      y=comp_chr[[chr]]$Genome1_score1,
      pch=19, col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome1_AB1],alpha.f=1))
points(x=comp_chr[[chr]]$mid,
       y=comp_chr[[chr]]$Genome2_score1,
       pch=19,  col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome2_AB1],alpha.f=1))
points(x=comp_chr[[chr]]$mid,
       y=rep(0,times=nrow(comp_chr[[chr]])),
       pch=19, col=adjustcolor(c("limegreen","lightskyblue2","maroon")[comp_chr[[chr]]$combined_comp],alpha.f=1))
legend("bottomleft",col=c("limegreen","lightskyblue2","maroon"), pch=19,legend = c("A/AA","AB","B/BB"),bty="n")
axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))


for (chr in 1:22){
    #   chr=22
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    
    if(chr==22){
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr",chr,".png",sep=""),width=width_plot,height=1500)#,res=350) #1000 #
    }else{
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr",chr,".png",sep=""),width=width_plot,height=1666)#,res=350) #1000 #
    }
    par(mar=c(3,2.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(10,1))
    if(chr==22){par(mfrow=c(9,1))}
    plot_NPMI_w_both(NPMI_1_chr,NPMI_2_chr)
    
    plot_zscores_w(zscores_diff_chr)
    plot_zscores_w(zscores_common_chr)
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr)
    
    # SNP density
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="number of SNPs per window")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # gene density
    plot(x=annoV38_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV38_df_genes_chr_binned_col[[chr]]$value,
         col=annoV38_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="number of genes per window")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))

    # AB compartments    # somehow I only have this up until chr 21
    if(chr<22){
        plot(x=comp_chr[[chr]]$mid,
             y=comp_chr[[chr]]$Genome1_score1,
             type = "l",
             xlim=c(0,chromsize_hg38[chr,2]),
             #ylim = c(-7,7),
             xaxt="none",xlab="",
             ylab="genome score")
        lines(x=comp_chr[[chr]]$mid,
              y=comp_chr[[chr]]$Genome2_score1,col="grey80")
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome1_score1,
               pch=19, col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome1_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome2_score1,
               pch=19,  col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome2_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=rep(0,times=nrow(comp_chr[[chr]])),
               pch=19, col=adjustcolor(c("limegreen","lightskyblue2","maroon")[comp_chr[[chr]]$combined_comp],alpha.f=1))
        legend("bottomleft",col=c("limegreen","lightskyblue2","maroon"), pch=19,legend = c("A/AA","AB","B/BB"),bty="n")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    }

    # total expression
    plot(x=anno_total_ase_imp_chr[[chr]]$mids,
         y=anno_total_ase_imp_chr[[chr]]$log10_TPM_transcript_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_chr[[chr]]$ENCFF379NOY_exp_ase_lfc+1,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # ASE ratio
    plot(x=anno_total_ase_imp_chr[[chr]]$mids,
         y=anno_total_ase_imp_chr[[chr]]$ASE_ratio_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_chr[[chr]]$ENCFF379NOY_exp_ase_lfc+1,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="ASE")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # LFC
    plot(x=anno_total_ase_imp_chr[[chr]]$mids,
         y=anno_total_ase_imp_chr[[chr]]$log2foldchange_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_chr[[chr]]$ENCFF379NOY_exp_ase_lfc+1,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="LFC")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # imprinted genes
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,-0.5,imprinted_genes_chr[[chrname]]$end,0.5,col=imprinted_genes_chr[[chrname]]$col)
    }
    dev.off()
}

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# I tried correlating SNV density and gene density per bin with sum of zscore of all or only top 5% most differential contacts
# that did not bring the results I hoped for
# now I try to consider a larger window around the bin

### expression maybe?
# just plain expression of gene with diff contact around it
# ase of gene with diff contact around it

### gene dens with diff contact around it
### SNP density with diff contact around it

### calc the zscore value for this "area"

diff_sum = lapply(zscores_diff,function(chr){
    return(colSums(abs(chr),na.rm=T))}) # this gets all contacts, upstream and downstream of the bin

diff_top_sum = lapply(1:22,function(chr){
    x1 = zscores_diff_1[[chr]]
    x2 = zscores_diff_2[[chr]]
    cs1 = colSums(abs(x1),na.rm=T)
    cs2 = colSums(abs(x2),na.rm=T)
    return(cs1+cs2)})


#################################################################################################################################
#################################################################################################################################
# dir_c = "/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/H1_phased/Curated/at50Kb"
# chr = 2
# NPMI_1_2 = data.matrix(read.table(gzfile(paste(dir_c,"/NPMI_Genome1/H1.Genome1_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
# NPMI_2_2 = data.matrix(read.table(gzfile(paste(dir_c,"/NPMI_Genome2/H1.Genome2_chr",chr,"_NPMI_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
# 
# zscores_common_2 = data.matrix(read.table(gzfile(paste(dir_c,"/ZScore_common/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_common_top10%_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
# zscores_diff_2 = data.matrix(read.table(gzfile(paste(dir_c,"/ZScore_diff/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_4000000_cutoff.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
# zscores_diff_1_2 = data.matrix(read.table(gzfile(paste(dir_c,"/ZScore_Top5_Genome1/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome1.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))
# zscores_diff_2_2 = data.matrix(read.table(gzfile(paste(dir_c,"/ZScore_Top5_Genome2/H1.Genome1_H1.Genome2_chr",chr,"_ZScore_diff_top5per_4000000_cutoff_H1.Genome2.txt.gz",sep="")),header=T,row.names=1,sep = "\t"))

chr = 2
NPMI_1_2 = NPMI_1[[chr]]
NPMI_2_2 = NPMI_2[[chr]]

zscores_common_2 = zscores_common[[chr]]
zscores_diff_2 = zscores_diff[[chr]]
zscores_diff_1_2 = zscores_diff_1[[chr]]
zscores_diff_2_2 = zscores_diff_2[[chr]]


# we are working with 50kb bins, 
# dim(NPMI_1_2) is [1] 4844 4844
# how large is the area we are looking at?
# at 10mb, which coordinates in the cluster? 10000000/50000 = 200
# 239mb ? 239000000/50000 = 4780 ja das passt ungefähr
NPMI_1_2[200,200:300]
# non NA goes from 10000000 to 14000000, that is 4000000 (4e+06 aka 4Mb) 4000000/50000 = 80 50kb bins (from the diagonal, so -80:+80)
# that would be the max range

range =  80

sum_zscores_common_2_range_80 = vector(mode="numeric",length=nrow(zscores_common_2))
sum_zscores_diff_2_range_80 = vector(mode="numeric",length=nrow(zscores_diff_2))
sum_zscores_top5_diff_2_range_80 = vector(mode="numeric",length=nrow(zscores_diff_1_2))

for (bin in (range+1):(nrow(zscores_diff_2)-(range+1))){
    sum_zscores_common_2_range_80[bin] = sum(abs(zscores_common_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T)
    sum_zscores_diff_2_range_80[bin] = sum(abs(zscores_diff_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T)
    
    #sum_zscores_top5_diff_2_range_80[bin] = sum(abs(zscores_diff_1_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T) + sum(abs(zscores_diff_2_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T)
}


# 1 MB range 1000000/50000 = 20 bins
range =  20

sum_zscores_common_2_range_20_1mb = vector(mode="numeric",length=nrow(zscores_common_2))
sum_zscores_diff_2_range_20_1mb = vector(mode="numeric",length=nrow(zscores_diff_2))
sum_zscores_top5_diff_2_range_20_1mb = vector(mode="numeric",length=nrow(zscores_diff_1_2))

for (bin in (range+1):(nrow(zscores_diff_2)-(range+1))){
    sum_zscores_common_2_range_20_1mb[bin] = sum(abs(zscores_common_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T)
    sum_zscores_diff_2_range_20_1mb[bin] = sum(abs(zscores_diff_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T)
    
    #sum_zscores_top5_diff_2_range_20_1mb[bin] = sum(abs(zscores_diff_1_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T) + sum(abs(zscores_diff_2_2[(bin+range):(bin-range),(bin+range):(bin-range)]),na.rm=T)
}

# gene expression is not yet binned, in 50kb bins, try corr with gene and SNP density
SNPs_H1 = readRDS("/Users/jmarkow/Desktop/GAM/data/H1_truths_mapp_filtered_H1_GT.rds")
nrow(SNPs_H1) #1506858
SNPs_H1_possort = SNPs_H1[order(SNPs_H1$POS),]
SNPs_H1 = SNPs_H1_possort[order(as.numeric(substring(SNPs_H1_possort$CHROM, 4))),]

SNPs_H1_chr = split(SNPs_H1,SNPs_H1$CHROM)
SNPs_H1_chr = SNPs_H1_chr[order(as.numeric(substring(names(SNPs_H1_chr), 4)))]

bins_all =  seq(0,250000000, 50000)
bins_all_mean = seq(25000,250000000, 50000)
bin_pos = sqrt(bins_all_mean^2 + bins_all_mean^2)/5e4

SNPs_H1_chr_binned = lapply(SNPs_H1_chr,function(chr){
    bins = cut(chr$POS, breaks = bins_all)
    return(table(bins))
})
max_snps_per_bin = max(unlist(lapply(SNPs_H1_chr_binned,max)))

plot(sum_zscores_diff_2_range_20_1mb,
     SNPs_H1_chr_binned[[2]][1:length(sum_zscores_common_2_range_20_1mb)])

cor.test(SNPs_H1_chr_binned[[2]][1:length(sum_zscores_common_2_range_20_1mb)], sum_zscores_common_2_range_20_1mb, method = c("pearson")) 
# Pearson's product-moment correlation
# 
# data:  SNPs_H1_binned_all_chr and diff_sum_all_chr
# t = 96.476, df = 57507, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.3661809 0.3802500
# sample estimates:
# cor 
# 0.3732369

cor.test(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("kendall")) 
# Kendall's rank correlation tau
# 
# data:  SNPs_H1_binned_all_chr and diff_sum_all_chr
# z = 89.428, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau
# 0.2559206

cor.test(SNPs_H1_binned_all_chr, diff_sum_all_chr, method = c("spearman")) 

################################################################################
################################################################################
### gene level expression data
anno_total_ase_imp_merged_ENCODE_g = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/expression_data/total_exp_and_ASE_incl_intro_anno_V29_imprinted_ENCFF675NTU_ENCFF379NOY_ENCFF174OMR_ENCFF910OBU_gene_level.rds")
anno_total_ase_imp_merged_ENCODE_g$mids = anno_total_ase_imp_merged_ENCODE_g$start_gene + (anno_total_ase_imp_merged_ENCODE_g$end_gene - anno_total_ase_imp_merged_ENCODE_g$start_gene)/2
anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF675NTU = log10(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF675NTU +1)
anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF379NOY = log10(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF379NOY +1)
anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF910OBU = log10(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF910OBU +1)
anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF174OMR = log10(anno_total_ase_imp_merged_ENCODE_g$TPM_ENCFF174OMR +1)

anno_total_ase_imp_merged_ENCODE_g$col = "black"
anno_total_ase_imp_merged_ENCODE_g$col[which(anno_total_ase_imp_merged_ENCODE_g$ENCFF379NOY_exp_ase_lfc==TRUE)] = "red"
anno_total_ase_imp_merged_ENCODE_g$col[which(anno_total_ase_imp_merged_ENCODE_g$gene_type!="protein_coding")] = adjustcolor(anno_total_ase_imp_merged_ENCODE_g$col[which(anno_total_ase_imp_merged_ENCODE_g$gene_type!="protein_coding")], alpha.f = 0.1)
# table(anno_total_ase_imp_merged_ENCODE_g$col)
# #0000001A #FF00001A     black       red 
#     38715        66     19717       223  

anno_total_ase_imp_merged_ENCODE_g_chr = split(anno_total_ase_imp_merged_ENCODE_g,anno_total_ase_imp_merged_ENCODE_g$chrom)
anno_total_ase_imp_merged_ENCODE_g_chr = anno_total_ase_imp_merged_ENCODE_g_chr[order(as.numeric(substring(names(anno_total_ase_imp_merged_ENCODE_g_chr), 4)))]
anno_total_ase_imp_merged_ENCODE_g_chr = anno_total_ase_imp_merged_ENCODE_g_chr[1:22]

################################################################################
################################################################################
### enhancers (H1, from enhanceratlas.org)
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_enhanceratlas.bed /Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_SNPs_enhanceratlas.bed /Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction.bed /Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/

H1_enh = read.csv("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/H1_enhancers_enhanceratlas.bed",sep = "\t",header = F) # 58821     4
colnames(H1_enh) = c("chr","start_enh","end_enh","enh_signal")
H1_enh$enh_mid = (H1_enh$start_enh + H1_enh$end_enh)/2

H1_enh_chr = split(H1_enh,H1_enh$chr)
H1_enh_chr = H1_enh_chr[order(as.numeric(substring(names(H1_enh_chr), 4)))]
H1_enh_chr = H1_enh_chr[1:22]

# load enhancers with SNPs and merge them with the H1 dixon truth to colour them in regarding theri haplotype (colour alt allele)
H1_enh_SNPs =  read.csv("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/H1_enhancers_SNPs_enhanceratlas.bed",sep = "\t",header = F) # 4408    8 #
# the SNP region was defined as the risk site and its flanking two base pairs
colnames(H1_enh_SNPs) = c("chr","start_enh","end_enh","enh_signal","chr_snp","start_snp","end_snp","SNP_ID")
H1_enh_SNPs_hap = merge(H1_enh_SNPs, SNPs_H1,by.x = "SNP_ID",by.y="ID") # 974  14, that is not at all much, given we have haplotype info of 1506858 SNPs 

################################################################################
################################################################################
### enhancer gene relationships (H1, Enhancer atlas)
# 
# H1_enh_gene = read.csv("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers/H1_enhancers_gene_interaction.bed",sep = "\t",header = F) # 227298      2
# # this needs some parsing / cleaning / extracting:
# # the data format in the file listed as "A:B-C_D$E$F$G$H I": 
# # A - Chromosome of enhancer. 
# # B - The starting position of enhancer.
# # C - The ending position of enhancer. 
# # D - Gene ensembl ID. 
# # E - Gene ensembl Name. 
# # F - chromosome of gene. 
# # G - the position of the gene transcription start site. 
# # H - The strand of DNA the gene located. 
# # I - The predition score of the enhancer-gene interaction.
# tmp = strsplit(as.character(H1_enh_gene$V1),"$",fixed=TRUE)
# H1_enh_gene$gene_name = sapply(tmp, `[`, 2)
# H1_enh_gene$gene_chr = sapply(tmp, `[`, 3)
# H1_enh_gene$gene_TSS = sapply(tmp, `[`, 4)
# H1_enh_gene$gene_strand = sapply(tmp, `[`, 5)
# tmp_2 = sapply(tmp, `[`, 1)
# tmp_3 = strsplit(as.character(tmp_2),"_",fixed=TRUE)
# H1_enh_gene$gene_id = sapply(tmp_3, `[`, 2)
# tmp_4 = sapply(tmp_3, `[`, 1)
# tmp_5 = strsplit(as.character(tmp_4),"[:-]",fixed=FALSE)
# H1_enh_gene$chr_enh = sapply(tmp_5, `[`, 1)
# H1_enh_gene$start_enh = sapply(tmp_5, `[`, 2)
# H1_enh_gene$end_enh = sapply(tmp_5, `[`, 3)
# H1_enh_gene = H1_enh_gene[,c("chr_enh","start_enh","end_enh","gene_id","gene_name","gene_chr","gene_TSS","gene_strand")] # 227298      8
# H1_enh_gene$mid_y = abs(as.numeric(H1_enh_gene$start_enh) - as.numeric(H1_enh_gene$gene_TSS))/2
# H1_enh_gene$mid_x = pmin(as.numeric(H1_enh_gene$start_enh), as.numeric(H1_enh_gene$gene_TSS)) + H1_enh_gene$mid_y

#dirc="/Users/jmarkow/Desktop/GAM/mnt/bih_cluster/fast/groups/ag_schwarz/Projects/project-gam/"
#H1_enh_gene = readRDS(paste(dirc,"/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/anno_H1/H1_enhancers_gene_interaction_liftover_hg38.rds",sep=""))
H1_enh_gene = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/enhancers//H1_enhancers_gene_interaction_liftover_hg38_annoV29.rds") # 215093     20
colnames(H1_enh_gene) = c("gene_id_atlas", "chr_enh", "start_enh", "end_enh", "gene_name_atlas", "gene_chr", "gene_TSS", "gene_strand_atlas", "chr_anno", "gene_start_anno", "gene_end_anno", "width", "gene_strand_anno", "gene_id_anno", "gene_type", "gene_name_anno", "chr_enh_red", "pos_diff", "enh_mid", "gene_mid_anno")
H1_enh_gene = H1_enh_gene[,c("chr_enh", "start_enh", "end_enh", "enh_mid", "gene_name_atlas", "gene_id_atlas", "gene_chr", "gene_TSS", "gene_strand_atlas", "chr_anno", "gene_start_anno", "gene_end_anno", "width",  "gene_mid_anno", "gene_strand_anno", "gene_id_anno", "gene_type", "gene_name_anno", "pos_diff")]
H1_enh_gene$mid_y = abs(as.numeric(H1_enh_gene$enh_mid) - as.numeric(H1_enh_gene$gene_mid_anno))/2
H1_enh_gene$mid_x = pmin(as.numeric(H1_enh_gene$enh_mid), as.numeric(H1_enh_gene$gene_mid_anno)) + H1_enh_gene$mid_y
H1_enh_gene$col = "grey75"
H1_enh_gene$col[which(H1_enh_gene$gene_type=="protein_coding")] = "black"

# all(H1_enh_gene$chr_enh==H1_enh_gene$gene_chr) # TRUE

H1_enh_gene_chr = split(H1_enh_gene,H1_enh_gene$chr_enh)
H1_enh_gene_chr = H1_enh_gene_chr[order(as.numeric(substring(names(H1_enh_gene_chr), 4)))]
#H1_enh_gene_chr = H1_enh_gene_chr[1:22]


# plot(x =  c(H1_enh_gene$start_enh[1:5],H1_enh_gene$gene_TSS[1:5]),
#      y =  rep(0,times=10),
#      ylim = c(-903169.5,903169.5))
# apply(H1_enh_gene[1:5,],1,function(edge){
#     x=c(edge['start_enh'],edge['gene_TSS'],edge['mid_x'])
#     y=c(0,0,edge['mid_y'])
#     polygon(x,y)})
################################################################################
# gene density for V29



annoV29_df_genes = readRDS(paste(dirc,"/H1/data/anno/gencode.v29.annotation.gtf.rds",sep="")) # 2742017      25
annoV29_df_genes = annoV29_df_genes[which(annoV29_df_genes$type=="gene"),] # 58721    25
chr_keep = paste("chr",c(1:22),sep="")

annoV29_df_genes_chr = split(annoV29_df_genes,annoV29_df_genes$seqnames)
annoV29_df_genes_chr = annoV29_df_genes_chr[chr_keep]

bins_all =  seq(0,250000000, 50000)
bins_all_mean = seq(25000,250000000, 50000)
bin_pos = sqrt(bins_all_mean^2 + bins_all_mean^2)/5e4

annoV29_df_genes_chr_binned = lapply(annoV29_df_genes_chr,function(chr){
    bins = cut(chr$start, breaks = bins_all)
    return(table(bins))
})
max_genes_per_bin = max(unlist(lapply(annoV29_df_genes_chr_binned,max)))

# colour code for the SNPs density
gene_dens_col = colorRampPalette(c("white","maroon"))(max_genes_per_bin+1) # do this max over all so it is comparable
gene_dens_col_log = colorRampPalette(c("white","maroon"))(ceiling(log(max_genes_per_bin))+1)# I put the +1 outsiede the log, bc later, I want to assign all values a col, and then assign 0 the white

# create dataframe with x position and colour (gene dens) for each chromosome
color_gene_dense_df = data.frame(val=0:max(max_genes_per_bin),colorname = gene_dens_col,stringsAsFactors = F)
color_gene_dense_log_df = data.frame(val=0:ceiling(log(max(max_genes_per_bin))),colorname = gene_dens_col_log,stringsAsFactors = F)

# since the majority of bins has very few genes, but the log does not really cut it, I decided to make a gradient up until 100 genes and set everything above to the max colour 
gene_dens_col_10lim = colorRampPalette(c(gene_dens_col))(11)
color_gene_dense_10lim_df = data.frame(val=0:max(max_genes_per_bin),colorname = c(gene_dens_col_10lim,rep(gene_dens_col_10lim[length(gene_dens_col_10lim)],times=(max(max_genes_per_bin)-10))),stringsAsFactors = F)

annoV29_df_genes_chr_binned_col = lapply(annoV29_df_genes_chr_binned,function(chr){
    # chr = annoV29_df_genes_chr_binned[[1]]
    chr_df = data.frame(val=chr,order=c(1:length(chr)))
    chr_df$val_log = ceiling(log(chr_df$val.Freq+1))
    # hist(chr_df$val_log,breaks=seq((-0.5),6.5,1))
    chr_df$val_log[which(chr_df$val.Freq==0)] = 0
    # hist(chr_df$val_log,breaks=seq((-0.5),6.5,1))
    
    df = merge(chr_df, color_gene_dense_df,by.x="val.Freq",by.y="val")
    df = merge(df, color_gene_dense_log_df,by.x="val_log",by.y="val")
    df = merge(df, color_gene_dense_10lim_df,by.x="val.Freq",by.y="val")
    
    df_sort = df[order(df$order),]
    df_sort$plot_pos = bin_pos
    df_sort$bin_mids = bins_all_mean
    
    df_sort = df_sort[,c(4,8,9,3,1,5,7,2,6)]
    colnames(df_sort) = c("order","plot_pos","bin_mids","bin","value","colour","colour_10lim","value_log","colour_log")
    return(df_sort)
})

annoV29_df_genes_chr_binned_col = lapply(annoV29_df_genes_chr_binned_col,function(chr){
    chr$pos = bins_all_mean
    return(chr)
})


annoV29_df_genes_chr_binned_all_chr = unlist(lapply(1:22,function(chr){
    tmp = as.numeric(annoV29_df_genes_chr_binned[[chr]][1:length(diff_top_sum[[chr]])])
    return(tmp)}))

################################################################################
### plot big ass plot with ASE on gene level, enhancers and enhancer gene interactions from enhancer atlas

tickpos = vector(mode  ="numeric", length = 25)
for (i in 1:25){
    mb = i*10000000
    tickpos[i] = sqrt(mb^2 + mb^2)/5e4
}

plot_NPMI_w_both_t <- function(x,y,range = c(0,d), title = "NPMI both genomes"){
    x[which(x<0)]=0
    #x[which(is.na(x))]=0
    
    ix <- cut(x, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
    #    colorsx <- colorRampPalette(c("steelblue4","lightblue","white","wheat1","darkgoldenrod1", "orangered","red4"))(99)[ix]
    colorsx = npmi_col[ix]
    xxx=matrix(colorsx,nrow=nrow(x),ncol=ncol(x))
    
    iy <- cut(y, breaks = seq(0, 1, len = 100), 
              include.lowest = TRUE)
    colorsy = npmi_col[iy]
    yyy=matrix(colorsy,nrow=nrow(y),ncol=ncol(y))
    
    ## clear lower triangle
    xxx[lower.tri(xxx)] <- yyy[lower.tri(yyy)]
    
    ## calculate diag
    nr <- nrow(xxx)
    nc <- ncol(xxx)
    d <- sqrt(nr^2 + nc^2)
    d2 <- 0.5 * d
    
    ## empty plot area
    #    plot(NA, type="n", xlim=c(0, d), ylim=c(-60, 60), xlab="", ylab="", asp=1, yaxt="none", xaxt="none")
    plot(NA, type="n", xlim = range, ylim=c(-60, 60), xlab="", ylab=title, asp=1, yaxt="none", xaxt="none") # don't set xlim according to chromsize, bc the scaling is different, 
    axis(1,at=tickpos,labels=paste(seq(10,250,10),"Mb",sep=""))
    # abline(h=0)
    
    ## plot matrix and rotate 45
    rasterImage(xxx,
                xleft=d2, xright=d2+nc, ybottom=-d2, ytop=-d2+nr,
                interpolate=FALSE, angle=45)
}


for (chr in 1:22){
    #   chr=2
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    # I use 1666 height for 10 panels, so 166.6 per panel
    panel_height = 1000/6
    panels = 11
    if(chr==22){
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*(panels-1))#,res=350) #1000 # # with 9 panels 1500
    }else{
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667
    }
    par(mar=c(3,4.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(panels,1))
    if(chr==22){par(mfrow=c((panels-1),1))}
    plot_NPMI_w_both_t(NPMI_1_chr,NPMI_2_chr)
    abline(h=0,col="grey75")
    
    plot_zscores_w(zscores_diff_chr,title="contact difference")
    abline(h=0,col="grey75")
    
    plot_zscores_w(zscores_common_chr,title="most common contacts")
    abline(h=0,col="grey75")
    
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr,title="most differential contacts")
    abline(h=0,col="grey75")
    
    
    # enhancers and enhancer contact
    # H1_enh_chr
            # H1_enh$width=abs(H1_enh$start_enh-H1_enh$end_enh)
            # max(H1_enh$width)
            # [1] 36050
            # so there is actually no reason to plot them as rectabgles, they are so small
            # rect(H1_enh_chr[[chr]]$start_enh,-0.5,H1_enh_chr[[chr]]$end_enh,0.5)
    # enhancer contacts
    max_y=max(abs(H1_enh_gene_chr[[chr]]$mid_y))
    plot(x =  c(H1_enh_gene_chr[[chr]]$enh_mid,H1_enh_gene_chr[[chr]]$gene_mid_anno),
         y =  rep(0,times=2*nrow(H1_enh_gene_chr[[chr]])),
         pch = 19,
         col = c(rep("orange",times=nrow(H1_enh_gene_chr[[chr]])),rep("cyan2",times=nrow(H1_enh_gene_chr[[chr]]))),
         ylim = c((-0.2)*max_y,max_y),
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="H1 enhancers")
    apply(H1_enh_gene_chr[[chr]],1,function(edge){
        x=c(edge['enh_mid'],edge['gene_mid_anno'],edge['mid_x'])
        y=c(0,0,edge['mid_y'])
        polygon(x,y,border=edge['col'])})
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    points(x=H1_enh_chr[[chr]]$enh_mid,
         y=rep((-0.1)*max_y,times=nrow(H1_enh_chr[[chr]])),
         pch=19, 
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="H1 enhancers")
    legend("bottomleft",col=c("black","orange","cyan2"), pch=19,legend = c("all enhancers","interacting enhancers","interacting genes"),bty="n")
    legend("left",col=c("black","grey75"),lty=1,lwd=3,legend=c("protein coding gene", "other gene type"),bty="n")
    
    
    # imprinted genes
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,(-0.15)*max_y,imprinted_genes_chr[[chrname]]$end,(-0.2)*max_y,col=imprinted_genes_chr[[chrname]]$col)
    }

    # total expression gene level
    plot(x=anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$log10_TPM_ENCFF379NOY,
         pch=19, col = anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("gene expression in log10(TPM+1)","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # ASE ratio
    plot(x=anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$ASE_ratio_ENCFF379NOY,
         pch=19, col = anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="ASE ratio")
    abline(h=c(1/3,0.5,2/3),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("ASE ratio","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # LFC
    plot(x=anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$log2foldchange_ENCFF379NOY,
         pch=19, col = anno_total_ase_imp_merged_ENCODE_g_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="L2FC")
    abline(h=c(-1,0,1),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("Log2foldchange","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # AB compartments    # somehow I only have this up until chr 21
    if(chr<22){
        plot(x=comp_chr[[chr]]$mid,
             y=comp_chr[[chr]]$Genome1_score1,
             type = "l",
             xlim=c(0,chromsize_hg38[chr,2]),
             #ylim = c(-7,7),
             xaxt="none",xlab="",
             ylab="AB compartments")#ylab="genome score")
        lines(x=comp_chr[[chr]]$mid,
              y=comp_chr[[chr]]$Genome2_score1,col="grey80")
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome1_score1,
               pch=19, col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome1_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome2_score1,
               pch=19,  col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome2_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=rep(0,times=nrow(comp_chr[[chr]])),
               pch=19, col=adjustcolor(c("limegreen","lightskyblue2","maroon")[comp_chr[[chr]]$combined_comp],alpha.f=1))
        legend("bottomleft",col=c("limegreen","lightskyblue2","maroon"), pch=19,legend = c("A/AA","AB","B/BB"),bty="n")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    }
    
    # SNP density
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="SNPs per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("SNP density per 50kb bin"),bty="n")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # gene density # this is 38, let'S take 29 since all other results are based on this
    plot(x=annoV29_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV29_df_genes_chr_binned_col[[chr]]$value,
         col=annoV29_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="genes per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("Gene density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))

    dev.off()
}


################################################################################
### plot big ass plot 
### filter genes for TPM, ASE 
### filter enhancer contact for expressed genes

#H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g = merge(H1_enh_gene,anno_total_ase_imp_merged_ENCODE_g,by = "gene_name")
dim(H1_enh_gene) #215093 was curated earlier now, this was the value before: 227298
dim(anno_total_ase_imp_merged_ENCODE_g) #58721
#dim(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g) #170055

anno_total_ase_imp_merged_ENCODE_g$gene_id_red = sapply(strsplit(anno_total_ase_imp_merged_ENCODE_g$gene_id,".",fixed=TRUE), `[`, 1)
H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g = merge(H1_enh_gene,anno_total_ase_imp_merged_ENCODE_g,by.x = "gene_id_atlas",by.y="gene_id_red")
dim(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g) #bc of previous filters:215093 , before 218143 # not ideal but much better

H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp = H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g[which(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE & H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g$sig_sum_ENCFF379NOY==TRUE),] # 3920
H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr = split(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp,H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp$chrom) 
H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr = H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[order(as.numeric(substring(names(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr), 4)))]
H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr = H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[1:22]

H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp = H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g[H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE,] # 167774
H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr = split(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp,H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp$chrom) #all(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp$chrom==H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp$chr_enh) # TRUE
H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr = H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[order(as.numeric(substring(names(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr), 4)))]
H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr = H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[1:22]

anno_total_ase_imp_merged_ENCODE_g_exp = anno_total_ase_imp_merged_ENCODE_g[anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE,]# 20120  72
anno_total_ase_imp_merged_ENCODE_g_exp_chr = split(anno_total_ase_imp_merged_ENCODE_g_exp,anno_total_ase_imp_merged_ENCODE_g_exp$chrom)
anno_total_ase_imp_merged_ENCODE_g_exp_chr = anno_total_ase_imp_merged_ENCODE_g_exp_chr[order(as.numeric(substring(names(anno_total_ase_imp_merged_ENCODE_g_exp_chr), 4)))]
anno_total_ase_imp_merged_ENCODE_g_exp_chr = anno_total_ase_imp_merged_ENCODE_g_exp_chr[1:22]

anno_total_ase_imp_merged_ENCODE_g_ase_exp = anno_total_ase_imp_merged_ENCODE_g[which(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE & anno_total_ase_imp_merged_ENCODE_g$sig_sum_ENCFF379NOY==TRUE),] # 379  72
anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr = split(anno_total_ase_imp_merged_ENCODE_g_ase_exp,anno_total_ase_imp_merged_ENCODE_g_ase_exp$chrom)
anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr = anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[order(as.numeric(substring(names(anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr), 4)))]
anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr = anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[1:22]

# only expressed genes 
for (chr in 1:22){
    #   chr=2
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    # I use 1666 height for 10 panels, so 166.6 per panel
    panel_height = 1000/6
    panels = 11
    if(chr==22){
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*(panels-1))#,res=350) #1000 # # with 9 panels 1500
    }else{
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667
    }
    par(mar=c(3,4.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(panels,1))
    if(chr==22){par(mfrow=c((panels-1),1))}
    plot_NPMI_w_both_t(NPMI_1_chr,NPMI_2_chr)
    abline(h=0,col="grey75")
    
    plot_zscores_w(zscores_diff_chr,title="contact difference")
    abline(h=0,col="grey75")
    
    plot_zscores_w(zscores_common_chr,title="most common contacts")
    abline(h=0,col="grey75")
    
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr,title="most differential contacts")
    abline(h=0,col="grey75")
    
    
    # enhancers and enhancer contact
    # H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr
    # enhancer contacts
    max_y=max(abs(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mid_y))
    plot(x =  c(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$enh_mid,H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$gene_mid_anno),
         y =  rep(0,times=2*nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]])),
         pch = 19,
         col = c(rep("orange",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]])),rep("cyan2",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]))),
         ylim = c((-0.2)*max_y,max_y),
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="H1 enhancers")
    apply(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]],1,function(edge){
        x=c(edge['enh_mid'],edge['gene_mid_anno'],edge['mid_x'])
        y=c(0,0,edge['mid_y'])
        polygon(x,y,border=edge['col.x'])})
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    points(x=H1_enh_chr[[chr]]$enh_mid,
           y=rep((-0.1)*max_y,times=nrow(H1_enh_chr[[chr]])),
           pch=19, 
           xlim=c(0,chromsize_hg38[chr,2]),
           xaxt="none",xlab="",
           ylab="H1 enhancers")
    legend("bottomleft",col=c("black","orange","cyan2"), pch=19,legend = c("all enhancers","interacting enhancers","interacting genes"),bty="n")
    legend("left",col=c("black","grey75"),lty=1,lwd=3,legend=c("protein coding gene", "other gene type"),bty="n")

    # imprinted genes
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,(-0.15)*max_y,imprinted_genes_chr[[chrname]]$end,(-0.2)*max_y,col=imprinted_genes_chr[[chrname]]$col)
    }
    
    # total expression gene level
    plot(x=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$log10_TPM_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    abline(h=0.3,col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("gene expression in log10(TPM+1)","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # ASE ratio
    plot(x=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$ASE_ratio_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="ASE ratio")
    abline(h=c(1/3,0.5,2/3),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("ASE ratio","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    #table(is.na(anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$ASE_ratio_ENCFF379NOY),anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$col)
    
    # LFC
    plot(x=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$log2foldchange_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="L2FC")
    abline(h=c(-1,0,1),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("Log2foldchange","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    
    # AB compartments    # somehow I only have this up until chr 21
    if(chr<22){
        plot(x=comp_chr[[chr]]$mid,
             y=comp_chr[[chr]]$Genome1_score1,
             type = "l",
             xlim=c(0,chromsize_hg38[chr,2]),
             #ylim = c(-7,7),
             xaxt="none",xlab="",
             ylab="AB compartments")#ylab="genome score")
        lines(x=comp_chr[[chr]]$mid,
              y=comp_chr[[chr]]$Genome2_score1,col="grey80")
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome1_score1,
               pch=19, col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome1_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome2_score1,
               pch=19,  col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome2_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=rep(0,times=nrow(comp_chr[[chr]])),
               pch=19, col=adjustcolor(c("limegreen","lightskyblue2","maroon")[comp_chr[[chr]]$combined_comp],alpha.f=1))
        legend("bottomleft",col=c("limegreen","lightskyblue2","maroon"), pch=19,legend = c("A/AA","AB","B/BB"),bty="n")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    }
    
    # SNP density
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="SNPs per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("SNP density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # gene density # this is 38, let'S take 29 since all other results are based on this
    plot(x=annoV29_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV29_df_genes_chr_binned_col[[chr]]$value,
         col=annoV29_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19,xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="genes per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("Gene density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    dev.off()
}

# expressed and significant ase genes

for (chr in 1:22){
    #   chr=15
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    # I use 1666 height for 10 panels, so 166.6 per panel
    panel_height = 1000/6
    panels = 11
    if(chr==22){
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_ase_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*(panels-1))#,res=350) #1000 # # with 9 panels 1500
    }else{
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_ase_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667
    } #        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_ase_chr",chr,"larger.png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667

    par(mar=c(3,4.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(panels,1))
    if(chr==22){par(mfrow=c((panels-1),1))}
    plot_NPMI_w_both_t(NPMI_1_chr,NPMI_2_chr)
    abline(h=0,col="grey75")
    
    plot_zscores_w(zscores_diff_chr,title="contact difference")
    abline(h=0,col="grey75")
    
    plot_zscores_w(zscores_common_chr,title="most common contacts")
    abline(h=0,col="grey75")
    
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr,title="most differential contacts")
    abline(h=0,col="grey75")
    
    
    # enhancers and enhancer contact
    # H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr
    # enhancer contacts
    max_y=max(abs(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mid_y))
    plot(x =  c(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$enh_mid,H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$gene_mid_anno),
         y =  rep(0,times=2*nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]])),
         pch = 19,
         col = c(rep("orange",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]])),rep("cyan2",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]))),
         ylim = c((-0.2)*max_y,max_y),
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="H1 enhancers")
    apply(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]],1,function(edge){
        x=c(edge['enh_mid'],edge['gene_mid_anno'],edge['mid_x'])
        y=c(0,0,edge['mid_y'])
        polygon(x,y,border=edge['col.x'])})
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    points(x=H1_enh_chr[[chr]]$enh_mid,
           y=rep((-0.1)*max_y,times=nrow(H1_enh_chr[[chr]])),
           pch=19, 
           xlim=c(0,chromsize_hg38[chr,2]),
           xaxt="none",xlab="",
           ylab="H1 enhancers")
    legend("bottomleft",col=c("black","orange","cyan2"), pch=19,legend = c("all enhancers","interacting enhancers","interacting genes"),bty="n")
    legend("left",col=c("black","grey75"),lty=1,lwd=3,legend=c("protein coding gene", "other gene type"),bty="n")
    
    # imprinted genes
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,(-0.15)*max_y,imprinted_genes_chr[[chrname]]$end,(-0.2)*max_y,col=imprinted_genes_chr[[chrname]]$col)
    }
    
    # total expression gene level
    plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$log10_TPM_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    abline(h=0.3,col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("gene expression in log10(TPM+1)","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # ASE ratio
    plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$ASE_ratio_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="ASE ratio")
    abline(h=c(1/3,0.5,2/3),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("ASE ratio","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))

    # LFC
    plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$log2foldchange_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="L2FC")
    abline(h=c(-1,0,1),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("Log2foldchange","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))

    
    
    
    # AB compartments    # somehow I only have this up until chr 21
    if(chr<22){
        plot(x=comp_chr[[chr]]$mid,
             y=comp_chr[[chr]]$Genome1_score1,
             type = "l",
             xlim=c(0,chromsize_hg38[chr,2]),
             #ylim = c(-7,7),
             xaxt="none",xlab="",
             ylab="AB compartments")#ylab="genome score")
        lines(x=comp_chr[[chr]]$mid,
              y=comp_chr[[chr]]$Genome2_score1,col="grey80")
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome1_score1,
               pch=19, col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome1_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome2_score1,
               pch=19,  col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome2_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=rep(0,times=nrow(comp_chr[[chr]])),
               pch=19, col=adjustcolor(c("limegreen","lightskyblue2","maroon")[comp_chr[[chr]]$combined_comp],alpha.f=1))
        legend("bottomleft",col=c("limegreen","lightskyblue2","maroon"), pch=19,legend = c("A/AA","AB","B/BB"),bty="n")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    }
    
    # SNP density
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19, cex=2, xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="SNPs per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("SNP density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # gene density # this is 38, let'S take 29 since all other results are based on this
    plot(x=annoV29_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV29_df_genes_chr_binned_col[[chr]]$value,
         col=annoV29_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19, cex=2, xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="genes per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("Gene density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    dev.off()
}

# 22 has no protein coding genes that are expressed and have ase?
# tmp22 = anno_total_ase_imp_merged_ENCODE_g_chr[[22]]
# dim(tmp22) # 1353   73
# plot(tmp22$log10_TPM_ENCFF379NOY,tmp22$ASE_ratio_ENCFF379NOY,
#      pch=19,col=tmp22$col)
# abline(v=log10(2))
# table(tmp22$expressed_ENCFF379NOY,tmp22$sig_sum_ENCFF379NOY,tmp22$gene_type=="protein_coding")
# tmp22_t = tmp22[which(tmp22$expressed_ENCFF379NOY == TRUE & tmp22$sig_sum_ENCFF675NTU == TRUE),]
# 
# tmp15 = anno_total_ase_imp_merged_ENCODE_g_chr[[15]]
# tmp15 = tmp15[order(tmp15$start_gene),]
# tmp15 = tmp15[170:300,]
# table(tmp15$expressed_ENCFF379NOY == TRUE)
# # FALSE  TRUE 
# #    80    51 
# table(tmp15$gene_type == "protein_coding")
# # FALSE  TRUE 
# #  122     9 
# #tmp15 = tmp15[which(tmp15$expressed_ENCFF379NOY == TRUE),]
# table(tmp15$expressed_ENCFF379NOY == TRUE,tmp15$gene_type == "protein_coding")

# plot(x=tmp15$mids,
#      y=tmp15$log10_TPM_ENCFF379NO,
#      pch=19, col=tmp15$col,
#      xlim=c(20000000,30000000),
#      xaxt="none",xlab="",
#      ylab="total gene expression")
# abline(v=c(24000000,26000000))
# abline(h=log10(2))
# 
# plot(x=anno_total_ase_imp_merged_ENCODE_g_chr[[15]]$mids,
#      y=anno_total_ase_imp_merged_ENCODE_g_chr[[15]]$log10_TPM_ENCFF379NOY,
#      pch=19, col=anno_total_ase_imp_merged_ENCODE_g_chr[[15]]$col,
#      xlim=c(20000000,30000000),
#      xaxt="none",xlab="",
#      ylab="total gene expression")
# abline(v=c(24000000,26000000))
# abline(h=log10(2))


#############################################################
#############################################################
#############################################################
# add similarity measures and sliding sums into the plot
# take care to plot the sliding windows from the middle of the window on

H1_SSIM = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_H1_NA_noNA.rds")

ssim_sliding_10mb_noNA = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_200bins_10mb_H1_zscore_diff_common_sum_noNA.rds")
ssim_sliding_10mb_noNA = lapply(ssim_sliding_10mb_noNA,function(measure){
    measure_out = lapply(measure,function(chrom){
        names(chrom) = seq(from=5000000,by=50000,length.out = length(chrom))
        return(chrom)
    })
    return(measure_out)
})
ssim_sliding_NPMI_10mb_noNA = ssim_sliding_10mb_noNA[[1]]
sum_sliding_zmostdiff1_10mb_noNA = ssim_sliding_10mb_noNA[[2]]
sum_sliding_zmostdiff2_10mb_noNA = ssim_sliding_10mb_noNA[[3]]
sum_sliding_zdiff10mb_noNA = ssim_sliding_10mb_noNA[[4]]
sum_sliding_zcom_10mb_noNA = ssim_sliding_10mb_noNA[[5]]
nonNAvalues_10mb = ssim_sliding_10mb_noNA[[6]]

ssim_sliding_1mb_noNA = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_20bins_1mb_H1_zscore_diff_common_sum_noNA.rds")
ssim_sliding_1mb_noNA = lapply(ssim_sliding_1mb_noNA,function(measure){
    measure_out = lapply(measure,function(chrom){
        names(chrom) = seq(from=500000,by=50000,length.out = length(chrom))
        return(chrom)
    })
    return(measure_out)
})
#names(ssim_sliding_NPMI_1mb_noNA[[chr]]) = seq(from=500000,by=50000,length.out = length(ssim_sliding_NPMI_1mb_noNA[[chr]]) )
ssim_sliding_NPMI_1mb_noNA = ssim_sliding_1mb_noNA[[1]]
sum_sliding_zmostdiff1_1mb_noNA = ssim_sliding_1mb_noNA[[2]]
sum_sliding_zmostdiff2_1mb_noNA = ssim_sliding_1mb_noNA[[3]]
sum_sliding_zdiff1mb_noNA = ssim_sliding_1mb_noNA[[4]]
sum_sliding_zcom_1mb_noNA = ssim_sliding_1mb_noNA[[5]]
nonNAvalues_1mb = ssim_sliding_1mb_noNA[[6]]

ssim_sliding_4mb_noNA = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/ssim/SSIM_chr_sliding_80bins_4mb_H1_zscore_diff_common_sum_noNA.rds")
ssim_sliding_4mb_noNA = lapply(ssim_sliding_4mb_noNA,function(measure){
    measure_out = lapply(measure,function(chrom){
        names(chrom) = seq(from=2000000,by=50000,length.out = length(chrom))
        return(chrom)
    })
    return(measure_out)
})
#names(ssim_sliding_NPMI_4mb_noNA[[chr]]) = seq(from=2000000,by=50000,length.out = length(ssim_sliding_NPMI_4mb_noNA[[chr]]) )
ssim_sliding_NPMI_4mb_noNA = ssim_sliding_4mb_noNA[[1]]
sum_sliding_zmostdiff1_4mb_noNA = ssim_sliding_4mb_noNA[[2]]
sum_sliding_zmostdiff2_4mb_noNA = ssim_sliding_4mb_noNA[[3]]
sum_sliding_zdiff4mb_noNA = ssim_sliding_4mb_noNA[[4]]
sum_sliding_zcom_4mb_noNA = ssim_sliding_4mb_noNA[[5]]
nonNAvalues_4mb = ssim_sliding_4mb_noNA[[6]]

chr = 12
par(mfrow=c(5,1))
plot(y = ssim_sliding_NPMI_1mb_noNA[[chr]],
     x = as.numeric(names(ssim_sliding_NPMI_1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,1))
points(y = ssim_sliding_NPMI_4mb_noNA[[chr]],
       x = as.numeric(names(ssim_sliding_NPMI_4mb_noNA[[chr]])), pch = 19, col = "red")
points(y = ssim_sliding_NPMI_10mb_noNA[[chr]],
       x = as.numeric(names(ssim_sliding_NPMI_10mb_noNA[[chr]])), pch = 19, col = "blue")

plot(y = sum_sliding_zmostdiff1_1mb_noNA[[chr]],
     x = as.numeric(names(sum_sliding_zmostdiff1_1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,max(sum_sliding_zmostdiff1_10mb_noNA[[chr]])))
points(y = sum_sliding_zmostdiff1_4mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zmostdiff1_4mb_noNA[[chr]])), pch = 19, col = "red")
points(y = sum_sliding_zmostdiff1_10mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zmostdiff1_10mb_noNA[[chr]])), pch = 19, col = "blue")

plot(y = sum_sliding_zmostdiff2_1mb_noNA[[chr]],
     x = as.numeric(names(sum_sliding_zmostdiff2_1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,max(sum_sliding_zmostdiff2_10mb_noNA[[chr]])))
points(y = sum_sliding_zmostdiff2_4mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zmostdiff2_4mb_noNA[[chr]])), pch = 19, col = "red")
points(y = sum_sliding_zmostdiff2_10mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zmostdiff2_10mb_noNA[[chr]])), pch = 19, col = "blue")

plot(y = sum_sliding_zdiff1mb_noNA[[chr]],
     x = as.numeric(names(sum_sliding_zdiff1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,max(sum_sliding_zdiff10mb_noNA[[chr]])))
points(y = sum_sliding_zdiff4mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zdiff4mb_noNA[[chr]])), pch = 19, col = "red")
points(y = sum_sliding_zdiff10mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zdiff10mb_noNA[[chr]])), pch = 19, col = "blue")


plot(y = sum_sliding_zcom_1mb_noNA[[chr]],
     x = as.numeric(names(sum_sliding_zcom_1mb_noNA[[chr]])),
     pch = 19,
     ylim=c(0,max(sum_sliding_zcom_10mb_noNA[[chr]])))
points(y = sum_sliding_zcom_4mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zcom_4mb_noNA[[chr]])), pch = 19, col = "red")
points(y = sum_sliding_zcom_10mb_noNA[[chr]],
       x = as.numeric(names(sum_sliding_zcom_10mb_noNA[[chr]])), pch = 19, col = "blue")

### add those to the large plots


for (chr in 1:22){
    #   chr=15
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    # I use 1666 height for 10 panels, so 166.6 per panel
    panel_height = 1000/6
    panels = 16
    if(chr==22){
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_ase_SSIM_slidingsum_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*(panels-1))#,res=350) #1000 # # with 9 panels 1500
    }else{
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_ase_SSIM_slidingsum_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667
    } #        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_ase_chr",chr,"larger.png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667
    
    par(mar=c(3,4.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(panels,1))
    if(chr==22){par(mfrow=c((panels-1),1))}
    plot_NPMI_w_both_t(NPMI_1_chr,NPMI_2_chr)
    abline(h=0,col="grey75")
    
    plot(y = ssim_sliding_NPMI_1mb_noNA[[chr]],
         x = as.numeric(names(ssim_sliding_NPMI_1mb_noNA[[chr]])),
         pch = 19,
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,1))
    points(y = ssim_sliding_NPMI_4mb_noNA[[chr]],
           x = as.numeric(names(ssim_sliding_NPMI_4mb_noNA[[chr]])), pch = 19, col = "hotpink3")
    points(y = ssim_sliding_NPMI_10mb_noNA[[chr]],
           x = as.numeric(names(ssim_sliding_NPMI_10mb_noNA[[chr]])), pch = 19, col = "deepskyblue3")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("black","hotpink3","deepskyblue3"), pch=19,legend = c("1Mb SSIM","4Mb SSIM","10Mb SSIM"),bty="n")

    plot_zscores_w(zscores_diff_chr,title="contact difference")
    abline(h=0,col="grey75")
    
    plot(y = sum_sliding_zdiff1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zdiff1mb_noNA[[chr]])),
         pch = 19,
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zdiff10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zdiff4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zdiff4mb_noNA[[chr]])), pch = 19, col = "hotpink3")
    points(y = sum_sliding_zdiff10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zdiff10mb_noNA[[chr]])), pch = 19, col = "deepskyblue3")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("black","hotpink3","deepskyblue3"), pch=19,legend = c("1Mb sum differences","4Mb  sum differences","10Mb  sum differences"),bty="n")
    
    plot_zscores_w(zscores_common_chr,title="most common contacts")
    abline(h=0,col="grey75")
    
    plot(y = sum_sliding_zcom_1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zcom_1mb_noNA[[chr]])),
         pch = 19,
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zcom_10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zcom_4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zcom_4mb_noNA[[chr]])), pch = 19, col = "hotpink3")
    points(y = sum_sliding_zcom_10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zcom_10mb_noNA[[chr]])), pch = 19, col = "deepskyblue3")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("black","hotpink3","deepskyblue3"), pch=19,legend = c("1Mb sum 10% most common","4Mb  sum 10% most common","10Mb  sum 10% most common"),bty="n")
    
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr,title="most differential contacts")
    abline(h=0,col="grey75")
    
    plot(y = sum_sliding_zmostdiff1_1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zmostdiff1_1mb_noNA[[chr]])),
         pch = 19, col = "burlywood1",
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zmostdiff1_10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zmostdiff1_4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff1_4mb_noNA[[chr]])), pch = 19, col = "darkorange2")
    points(y = sum_sliding_zmostdiff1_10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff1_10mb_noNA[[chr]])), pch = 19, col = "darkorange4")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("burlywood1","darkorange2","darkorange4"), pch=19,legend = c("1Mb sum 5% most different Genome 1","4Mb  sum 5% most different Genome 1","10Mb  sum 5% most different Genome 1"),bty="n")
    
    plot(y = sum_sliding_zmostdiff2_1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zmostdiff2_1mb_noNA[[chr]])),
         pch = 19,col = "lavender",
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zmostdiff2_10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zmostdiff2_4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff2_4mb_noNA[[chr]])), pch = 19, col = "mediumpurple2")
    points(y = sum_sliding_zmostdiff2_10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff2_10mb_noNA[[chr]])), pch = 19, col = "purple4")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("lavender","mediumpurple2","purple4"), pch=19,legend = c("1Mb sum 5% most different Genome 2","4Mb  sum 5% most different Genome 2","10Mb  sum 5% most different Genome 2"),bty="n")
    
    # enhancers and enhancer contact
    # H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr
    # enhancer contacts
    max_y=max(abs(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mid_y))
    plot(x =  c(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$enh_mid,H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$gene_mid_anno),
         y =  rep(0,times=2*nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]])),
         pch = 19,
         col = c(rep("orange",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]])),rep("cyan2",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]))),
         ylim = c((-0.2)*max_y,max_y),
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="H1 enhancers")
    apply(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]],1,function(edge){
        x=c(edge['enh_mid'],edge['gene_mid_anno'],edge['mid_x'])
        y=c(0,0,edge['mid_y'])
        polygon(x,y,border=edge['col.x'])})
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    points(x=H1_enh_chr[[chr]]$enh_mid,
           y=rep((-0.1)*max_y,times=nrow(H1_enh_chr[[chr]])),
           pch=19, 
           xlim=c(0,chromsize_hg38[chr,2]),
           xaxt="none",xlab="",
           ylab="H1 enhancers")
    legend("bottomleft",col=c("black","orange","cyan2"), pch=19,legend = c("all enhancers","interacting enhancers","interacting genes"),bty="n")
    legend("left",col=c("black","grey75"),lty=1,lwd=3,legend=c("protein coding gene", "other gene type"),bty="n")
    
    # imprinted genes
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,(-0.15)*max_y,imprinted_genes_chr[[chrname]]$end,(-0.2)*max_y,col=imprinted_genes_chr[[chrname]]$col)
    }
    
    # total expression gene level
    plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$log10_TPM_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    abline(h=0.3,col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("gene expression in log10(TPM+1)","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # ASE ratio
    plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$ASE_ratio_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="ASE ratio")
    abline(h=c(1/3,0.5,2/3),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("ASE ratio","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # LFC
    plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$log2foldchange_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="L2FC")
    abline(h=c(-1,0,1),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("Log2foldchange","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    
    
    
    # AB compartments    # somehow I only have this up until chr 21
    if(chr<22){
        plot(x=comp_chr[[chr]]$mid,
             y=comp_chr[[chr]]$Genome1_score1,
             type = "l",
             xlim=c(0,chromsize_hg38[chr,2]),
             #ylim = c(-7,7),
             xaxt="none",xlab="",
             ylab="AB compartments")#ylab="genome score")
        lines(x=comp_chr[[chr]]$mid,
              y=comp_chr[[chr]]$Genome2_score1,col="grey80")
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome1_score1,
               pch=19, col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome1_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome2_score1,
               pch=19,  col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome2_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=rep(0,times=nrow(comp_chr[[chr]])),
               pch=19, col=adjustcolor(c("limegreen","lightskyblue2","maroon")[comp_chr[[chr]]$combined_comp],alpha.f=1))
        legend("bottomleft",col=c("limegreen","lightskyblue2","maroon"), pch=19,legend = c("A/AA","AB","B/BB"),bty="n")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    }
    
    # SNP density
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19, cex=2, xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="SNPs per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("SNP density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # gene density # this is 38, let'S take 29 since all other results are based on this
    plot(x=annoV29_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV29_df_genes_chr_binned_col[[chr]]$value,
         col=annoV29_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19, cex=2, xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="genes per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("Gene density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    dev.off()
}

#############################################################
#############################################################
# 05.10.21
# compare SSIM with AB compartments
# AB compartments are called in 200k bins
# integrate with SSIM in 1Mb sliding windows, still 50k bins


require(data.table)
comp_dt=as.data.table(comp)
#lapply(seq_along(x), function(i) paste(names(x)[[i]], x[[i]]))

ssim_sliding_NPMI_10mb_noNA_df = lapply(seq_along(ssim_sliding_NPMI_10mb_noNA), function(chrom_idx){
    chrom = cbind(ssim_sliding_NPMI_10mb_noNA[[chrom_idx]],
                  paste("chr",chrom_idx,sep=""),
                  as.numeric(names(ssim_sliding_NPMI_10mb_noNA[[chrom_idx]])),
                  as.numeric(names(ssim_sliding_NPMI_10mb_noNA[[chrom_idx]]))+50000)
    colnames(chrom) = c("SSIM_1Mb","chrom","start","end")
    rownames(chrom) <- NULL
    return(chrom)
})
ssim_sliding_NPMI_10mb_noNA_df = do.call(rbind,ssim_sliding_NPMI_10mb_noNA_df)
SSIM_1Mb_dt = as.data.table(ssim_sliding_NPMI_10mb_noNA_df)
SSIM_1Mb_dt$SSIM_1Mb = as.numeric(SSIM_1Mb_dt$SSIM_1Mb)
SSIM_1Mb_dt$start = as.numeric(SSIM_1Mb_dt$start)
SSIM_1Mb_dt$start = SSIM_1Mb_dt$start+1
SSIM_1Mb_dt$end = as.numeric(SSIM_1Mb_dt$end)

setkey(SSIM_1Mb_dt, chrom, start, end)
comp_SSIM_1Mb_dt = foverlaps(comp_dt, SSIM_1Mb_dt, type="any", mult="all") 

boxplot(comp_SSIM_1Mb_dt$SSIM_1Mb~comp_SSIM_1Mb_dt$combined_comp,
        outline=F,
        col=c("limegreen","lightskyblue2","maroon"),
        main="",
        ylab="SSIM",
        xlab="compartment",
        notch=T)

# YAY, there is a difference! 

comp_SSIM_1Mb_dt$col=c("limegreen","lightskyblue2","maroon")[comp_SSIM_1Mb_dt$combined_comp]
my_comparisons <- list( c("AA", "AB"), c("AA", "BB"), c("AB", "BB"))
library("ggsignif")
library("rstatix")
require("gridExtra")
library("ggplot2")
library("ggpubr")
ggplot(comp_SSIM_1Mb_dt, aes(x=combined_comp, y=SSIM_1Mb,fill=col)) +  ylab("SSIM") + xlab("Compartment Level") + 
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

# YAY, there is a significant difference! 

# group strength: 
comp_SSIM_1Mb_dt_noSSIMNA = comp_SSIM_1Mb_dt[which(!(is.na(comp_SSIM_1Mb_dt$SSIM_1Mb))),]
# > dim(comp_SSIM_1Mb_dt)
# [1] 44670    14
# > dim(comp_SSIM_1Mb_dt_noSSIMNA)
# [1] 43903    14
table(comp_SSIM_1Mb_dt_noSSIMNA$combined_comp)
# AA    AB    BB 
# 14293 17600 12010
table(comp_SSIM_1Mb_dt_noSSIMNA$combined_comp)/nrow(comp_SSIM_1Mb_dt_noSSIMNA)*100
# AA       AB       BB 
# 32.55586 40.08838 27.35576 

# gene expression and compartments
table(anno_total_ase_imp_merged_ENCODE_g_ase_exp$expressed_ENCFF675NTU,exclude = NULL) # TRUE 379 
table(anno_total_ase_imp_merged_ENCODE_g_ase_exp$expressed_ENCFF379NOY,exclude = NULL) # TRUE 379 
table(anno_total_ase_imp_merged_ENCODE_g_ase_exp$expressed_ENCFF174OMR,exclude = NULL) # FALSE 50 TRUE 329
table(anno_total_ase_imp_merged_ENCODE_g_ase_exp$expressed_ENCFF910OBU,exclude = NULL) # FALSE 51 TRUE 328
# so this table I am using right now contains only expressed genes as to my previous analysis

table(anno_total_ase_imp_merged_ENCODE_g_ase_exp$sig_sum_ENCFF675NTU,exclude = NULL) # FALSE 88 TRUE 291
table(anno_total_ase_imp_merged_ENCODE_g_ase_exp$sig_sum_ENCFF379NOY,exclude = NULL) # TRUE 379 
# and only those that actually already show ase, this is biased

# I should take all genes first, check gene expression in compartments, 
# then only take expressed genes and check ASE and L2FC of those
dim(anno_total_ase_imp_merged_ENCODE_g) # [1] 58721    72
table(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU,exclude = NULL) # FALSE 37781 TRUE 20940
table(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY,exclude = NULL) # FALSE 38601 TRUE 20120 
table(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR,exclude = NULL) # FALSE 43720 TRUE 15001
table(anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU,exclude = NULL) # FALSE 43443 TRUE 15278
# so this table I am using right now contains only expressed genes as to my previous analysis

table(anno_total_ase_imp_merged_ENCODE_g$sig_sum_ENCFF675NTU,exclude = NULL) # FALSE 15397 TRUE 450 NA 42874
table(anno_total_ase_imp_merged_ENCODE_g$sig_sum_ENCFF379NOY,exclude = NULL) # FALSE 14910 TRUE 416 NA 43395

# all genes
anno_total_ase_imp_merged_ENCODE_g_sort = anno_total_ase_imp_merged_ENCODE_g[order(anno_total_ase_imp_merged_ENCODE_g$start_gene),] # 58721
anno_total_ase_imp_merged_ENCODE_g_sort = anno_total_ase_imp_merged_ENCODE_g_sort[order(as.numeric(substring(anno_total_ase_imp_merged_ENCODE_g_sort$chrom, 4))),]

anno_total_ase_imp_merged_ENCODE_g_sort_dt = as.data.table(anno_total_ase_imp_merged_ENCODE_g_sort)
setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
comp_anno_total_ase_imp_merged_ENCODE_g = foverlaps(comp_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 53776, without X chromosome etc 


# expressed ASE genes
anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort = anno_total_ase_imp_merged_ENCODE_g_ase_exp[order(anno_total_ase_imp_merged_ENCODE_g_ase_exp$start_gene),] # 379
anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort = anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort[order(as.numeric(substring(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort$chrom, 4))),]

anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt = as.data.table(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort)
setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(comp_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) 


# all genes
comp_anno_total_ase_imp_merged_ENCODE_g$colAB=c("limegreen","lightskyblue2","maroon")[comp_anno_total_ase_imp_merged_ENCODE_g$combined_comp]
my_comparisons <- list( c("AA", "AB"), c("AA", "BB"), c("AB", "BB"))
ggplot(comp_anno_total_ase_imp_merged_ENCODE_g, aes(x=combined_comp, y=log10_TPM_ENCFF675NTU,fill=colAB)) +  ylab("log10_TPM_ENCFF675NTU") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g, aes(x=combined_comp, y=log10_TPM_ENCFF379NOY,fill=colAB)) +  ylab("log10_TPM_ENCFF379NOY") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g, aes(x=combined_comp, y=log10_TPM_ENCFF910OBU,fill=colAB)) +  ylab("log10_TPM_ENCFF910OBU") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g, aes(x=combined_comp, y=log10_TPM_ENCFF174OMR,fill=colAB)) +  ylab("log10_TPM_ENCFF174OMR") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

# nope none of them show great signal -.- 
# consider only expressed ones 

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),], aes(x=combined_comp, y=log10_TPM_ENCFF675NTU,fill=colAB)) +  ylab("log10_TPM_ENCFF675NTU") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),], aes(x=combined_comp, y=log10_TPM_ENCFF379NOY,fill=colAB)) +  ylab("log10_TPM_ENCFF379NOY") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU==TRUE),], aes(x=combined_comp, y=log10_TPM_ENCFF910OBU,fill=colAB)) +  ylab("log10_TPM_ENCFF910OBU") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR==TRUE),], aes(x=combined_comp, y=log10_TPM_ENCFF174OMR,fill=colAB)) +  ylab("log10_TPM_ENCFF174OMR") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

# better, but stil a no :(

# ase of those
# all genes
ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),], aes(x=combined_comp, y=ASE_ratio_ENCFF675NTU,fill=colAB)) +  ylab("ASE ENCFF675NTU") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),], aes(x=combined_comp, y=ASE_ratio_ENCFF379NOY,fill=colAB)) +  ylab("ASE ENCFF379NOY") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

 # ase genes
ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=ASE_ratio_ENCFF675NTU,fill=colAB)) +  ylab("ASE ENCFF675NTU") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=ASE_ratio_ENCFF379NOY,fill=colAB)) +  ylab("ASE ENCFF379NOY") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

# lfc of those

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),], aes(x=combined_comp, y=log2foldchange_ENCFF675NTU,fill=colAB)) +  ylab("LFC") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),], aes(x=combined_comp, y=log2foldchange_ENCFF379NOY,fill=colAB)) +  ylab("LFC") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)




# now think hard how you want to compare those, only genes, or number of genes in comp
# well, how is the expression in a,b,ab comp,
# how is the ase in a,b,ab comp
comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp$colAB=c("limegreen","lightskyblue2","maroon")[comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp$combined_comp]
my_comparisons <- list( c("AA", "AB"), c("AA", "BB"), c("AB", "BB"))
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=TPM_ENCFF675NTU,fill=col)) +  ylab("TPM") + xlab("Compartment Level") + 
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)
# 
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=TPM_ENCFF379NOY,fill=col)) +  ylab("TPM") + xlab("Compartment Level") + 
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)
# 
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=TPM_ENCFF174OMR,fill=col)) +  ylab("TPM") + xlab("Compartment Level") + 
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)
# 
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=TPM_ENCFF910OBU,fill=col)) +  ylab("TPM") + xlab("Compartment Level") + 
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)
# 
# # log
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=log10_TPM_ENCFF675NTU,fill=col)) +  ylab("log10_TPM") + xlab("Compartment Level") + 
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)
# 
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=log10_TPM_ENCFF379NOY,fill=col)) +  ylab("log10_TPM") + xlab("Compartment Level") + 
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=log10_TPM_ENCFF910OBU,fill=col)) +  ylab("log10_TPM") + xlab("Compartment Level") + 
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=log10_TPM_ENCFF174OMR,fill=col)) +  ylab("log10_TPM") + xlab("Compartment Level") + 
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() + 
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

# OBU and OMR are promising, but not at all significant


# ASE

# take care to get the whole df, not just the expressed ones? just to have a look you know
ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=ASE_ratio_ENCFF675NTU,fill=col)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp, y=ASE_ratio_ENCFF379NOY,fill=col)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons)

# ASE but difference between AB and BA compartments
my_comparisons2 <- list( c("AA", "AB"), c("AA", "BA"),c("AA", "BB"),c("AB", "BA"), c("AB", "BB"),c("BA", "BB"))
ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp_raw, y=ASE_ratio_ENCFF675NTU,fill=col)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp_raw, y=ASE_ratio_ENCFF379NOY,fill=col)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)


my_comparisons2 <- list( c("AA", "AB"), c("AA", "BA"),c("AA", "BB"),c("AB", "BA"), c("AB", "BB"),c("BA", "BB"))
ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),], aes(x=combined_comp_raw, y=ASE_ratio_ENCFF675NTU,fill=colAB)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressedENCFF379NOY==TRUE),], aes(x=combined_comp_raw, y=ASE_ratio_ENCFF379NOY,fill=colAB)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)

my_comparisons2 <- list( c("AA", "AB"), c("AA", "BA"),c("AA", "BB"),c("AB", "BA"), c("AB", "BB"),c("BA", "BB"))
ggdata = comp_anno_total_ase_imp_merged_ENCODE_g[which((comp_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE) & (comp_anno_total_ase_imp_merged_ENCODE_g$gene_type=="protein_coding")),]
ggplot(ggdata, aes(x=combined_comp_raw, y=ASE_ratio_ENCFF675NTU,fill=colAB)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)

ggplot(comp_anno_total_ase_imp_merged_ENCODE_g[which(comp_anno_total_ase_imp_merged_ENCODE_g$expressedENCFF379NOY==TRUE),], aes(x=combined_comp_raw, y=ASE_ratio_ENCFF379NOY,fill=col)) +  ylab("ASE ratio") + xlab("Compartment Level") +
    geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
    theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)



# my_comparisons2 <- list( c("AA", "AB"), c("AA", "BA"),c("AA", "BB"),c("AB", "BA"), c("AB", "BB"),c("BA", "BB"))
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp[which(comp_anno_total_ase_imp_merged_ENCODE_g$gene_type=="protein_coding"),], aes(x=combined_comp_raw, y=ASE_ratio_ENCFF675NTU,fill=col)) +  ylab("ASE ratio") + xlab("Compartment Level") +
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)
# 
# ggplot(comp_anno_total_ase_imp_merged_ENCODE_g_ase_exp, aes(x=combined_comp_raw, y=ASE_ratio_ENCFF379NOY,fill=col)) +  ylab("ASE ratio") + xlab("Compartment Level") +
#     geom_violin() + scale_fill_manual(values=c("lightskyblue2","limegreen","maroon")) + theme_classic() +
#     theme(legend.position="none") + geom_boxplot(width=0.1) + stat_summary(fun.y=mean, geom="point", size=2, color="black") + stat_compare_means(comparisons = my_comparisons2)
# 
# remove the NA genes, check how many genes are in each compartment category (table)

##### gene expression with SSIM 
# ASE with SSIM 
setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g = foverlaps(SSIM_1Mb_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 95145
setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(SSIM_1Mb_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 53199

# expressed ASE genes
par(mfrow=c(2,2))
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$log10_TPM_ENCFF675NTU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF675NTU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$log10_TPM_ENCFF379NOY,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF379NOY",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$log10_TPM_ENCFF910OBU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF910OBU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$log10_TPM_ENCFF174OMR,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF174OMR",
     ylab = "SSIM 1Mb")

# expressed genes
par(mfrow=c(2,2))
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$log10_TPM_ENCFF675NTU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF675NTU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$log10_TPM_ENCFF379NOY,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF379NOY",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU==TRUE),]$log10_TPM_ENCFF910OBU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF910OBU==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF910OBU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR==TRUE),]$log10_TPM_ENCFF174OMR,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF174OMR==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF174OMR",
     ylab = "SSIM 1Mb")


# all genes, expression SSIM
par(mfrow=c(2,2))
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF675NTU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF675NTU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF379NOY,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF379NOY",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF910OBU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF910OBU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF174OMR,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10_TPM_ENCFF174OMR",
     ylab = "SSIM 1Mb")

# try with Teresas data 
gene_exp_T=read.csv("/Users/jmarkow/Desktop/GAM/Co-Phasing/data/expression_data/genes_with_expression.bed",header = T, sep="\t")
colnames(gene_exp_T)[1] = "chrom"
dim(gene_exp_T) # 60725     9

red_chrom_set=paste("chr",1:22,sep='')
gene_exp_T=gene_exp_T[which(gene_exp_T$chrom%in%red_chrom_set),]
dim(gene_exp_T) # 57439     9
# sort by chromosome and position
gene_exp_T = gene_exp_T[order(gene_exp_T$start),]
gene_exp_T = gene_exp_T[order(as.numeric(substring(gene_exp_T$chrom, 4))),]

gene_exp_T_dt = as.data.table(gene_exp_T)

setkey(gene_exp_T_dt, chrom, start, end)
SSIM_1Mb_gene_exp_T = foverlaps(SSIM_1Mb_dt, gene_exp_T_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 95430

par(mfrow=c(2,2))
plot(x = SSIM_1Mb_gene_exp_T$TPM,
     y = SSIM_1Mb_gene_exp_T$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "TPM Teresa",
     ylab = "SSIM 1Mb")

plot(x = log10(SSIM_1Mb_gene_exp_T$TPM+1),
     y = SSIM_1Mb_gene_exp_T$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log10 TPM Teresa",
     ylab = "SSIM 1Mb")
##############################
### ASE, all genes
par(mfrow=c(2,2))
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "ASE_ratio_ENCFF675NTU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF379NOY,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "ASE_ratio_ENCFF379NOY",
     ylab = "SSIM 1Mb")

### SSIM LFC, all genes
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$log2foldchange_ENCFF675NTU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log2foldchange_ENCFF675NTU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$log2foldchange_ENCFF379NOY,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,
     xlab = "log2foldchange_ENCFF379NOY",
     ylab = "SSIM 1Mb")
##############################
### SSIM, ASE, expressed genes
par(mfrow=c(2,2))
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$ASE_ratio_ENCFF675NTU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, ylim = c(0.4,1),
     xlab = "ASE_ratio_ENCFF675NTU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$ASE_ratio_ENCFF379NOY,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, ylim = c(0.4,1),
     xlab = "ASE_ratio_ENCFF379NOY",
     ylab = "SSIM 1Mb")

### SSIM, LFC, expressed genes 
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$log2foldchange_ENCFF675NTU,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, ylim = c(0.4,1),
     xlab = "log2foldchange_ENCFF675NTU",
     ylab = "SSIM 1Mb")
plot(x = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$log2foldchange_ENCFF379NOY,
     y = SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(SSIM_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$SSIM_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, ylim = c(0.4,1),
     xlab = "log2foldchange_ENCFF379NOY",
     ylab = "SSIM 1Mb")


#### ASE at most common and most different contacts
# ssim_sliding_NPMI_10mb_noNA = ssim_sliding_10mb_noNA[[1]]
# sum_sliding_zmostdiff1_10mb_noNA = ssim_sliding_10mb_noNA[[2]]
# sum_sliding_zmostdiff2_10mb_noNA = ssim_sliding_10mb_noNA[[3]]
# sum_sliding_zdiff10mb_noNA = ssim_sliding_10mb_noNA[[4]]
# sum_sliding_zcom_10mb_noNA = ssim_sliding_10mb_noNA[[5]]
# nonNAvalues_10mb = ssim_sliding_10mb_noNA[[6]]


sum_sliding_zcom_1mb_noNA_df = lapply(seq_along(sum_sliding_zcom_1mb_noNA), function(chrom_idx){
    chrom = cbind(sum_sliding_zcom_1mb_noNA[[chrom_idx]],
                  paste("chr",chrom_idx,sep=""),
                  as.numeric(names(sum_sliding_zcom_1mb_noNA[[chrom_idx]])),
                  as.numeric(names(sum_sliding_zcom_1mb_noNA[[chrom_idx]]))+50000)
    colnames(chrom) = c("mostcom_1Mb","chrom","start","end")
    rownames(chrom) <- NULL
    return(chrom)
})
sum_sliding_zcom_1mb_noNA_df = do.call(rbind,sum_sliding_zcom_1mb_noNA_df)
mostcom_1Mb_dt = as.data.table(sum_sliding_zcom_1mb_noNA_df)
mostcom_1Mb_dt$mostcom_1Mb = as.numeric(mostcom_1Mb_dt$mostcom_1Mb)
mostcom_1Mb_dt$start = as.numeric(mostcom_1Mb_dt$start)
mostcom_1Mb_dt$start = mostcom_1Mb_dt$start+1
mostcom_1Mb_dt$end = as.numeric(mostcom_1Mb_dt$end)

# setkey(mostcom_1Mb_dt, chrom, start, end)
# comp_mostcom_1Mb_dt = foverlaps(comp_dt, mostcom_1Mb_dt, type="any", mult="all") 

setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g = foverlaps(mostcom_1Mb_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 103775
setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(mostcom_1Mb_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 57168

par(mfrow=c(3,2))
# all genes 
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 1Mb")
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF379NOY,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 1Mb")
# expressed genes 
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$ASE_ratio_ENCFF675NTU,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 1Mb")
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$ASE_ratio_ENCFF379NOY,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 1Mb")

### ase expressed genes 
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF675NTU,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 1Mb")
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF379NOY,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 1Mb")


### 4 MB


sum_sliding_zcom_4Mb_noNA_df = lapply(seq_along(sum_sliding_zcom_4mb_noNA), function(chrom_idx){
    chrom = cbind(sum_sliding_zcom_4mb_noNA[[chrom_idx]],
                  paste("chr",chrom_idx,sep=""),
                  as.numeric(names(sum_sliding_zcom_4mb_noNA[[chrom_idx]])),
                  as.numeric(names(sum_sliding_zcom_4mb_noNA[[chrom_idx]]))+50000)
    colnames(chrom) = c("mostcom_4Mb","chrom","start","end")
    rownames(chrom) <- NULL
    return(chrom)
})
sum_sliding_zcom_4Mb_noNA_df = do.call(rbind,sum_sliding_zcom_4Mb_noNA_df)
mostcom_4Mb_dt = as.data.table(sum_sliding_zcom_4Mb_noNA_df)
mostcom_4Mb_dt$mostcom_4Mb = as.numeric(mostcom_4Mb_dt$mostcom_4Mb)
mostcom_4Mb_dt$start = as.numeric(mostcom_4Mb_dt$start)
mostcom_4Mb_dt$start = mostcom_4Mb_dt$start+1
mostcom_4Mb_dt$end = as.numeric(mostcom_4Mb_dt$end)

# setkey(mostcom_4Mb_dt, chrom, start, end)
# comp_mostcom_4Mb_dt = foverlaps(comp_dt, mostcom_4Mb_dt, type="any", mult="all") 

setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g = foverlaps(mostcom_4Mb_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 100399
setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(mostcom_4Mb_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 55843

par(mfrow=c(3,2))
# all genes 
plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
     x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_4Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 4Mb")
plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF379NOY,
     x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_4Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 4Mb")
# expressed genes 
plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$ASE_ratio_ENCFF675NTU,
     x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$mostcom_4Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 4Mb")
plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$ASE_ratio_ENCFF379NOY,
     x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$mostcom_4Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 4Mb")

### ase expressed genes 
plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF675NTU,
     x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_4Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 4Mb")
plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF379NOY,
     x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_4Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 4Mb")

### 10 Mb


sum_sliding_zcom_10Mb_noNA_df = lapply(seq_along(sum_sliding_zcom_10mb_noNA), function(chrom_idx){
    chrom = cbind(sum_sliding_zcom_10mb_noNA[[chrom_idx]],
                  paste("chr",chrom_idx,sep=""),
                  as.numeric(names(sum_sliding_zcom_10mb_noNA[[chrom_idx]])),
                  as.numeric(names(sum_sliding_zcom_10mb_noNA[[chrom_idx]]))+50000)
    colnames(chrom) = c("mostcom_10Mb","chrom","start","end")
    rownames(chrom) <- NULL
    return(chrom)
})
sum_sliding_zcom_10Mb_noNA_df = do.call(rbind,sum_sliding_zcom_10Mb_noNA_df)
mostcom_10Mb_dt = as.data.table(sum_sliding_zcom_10Mb_noNA_df)
mostcom_10Mb_dt$mostcom_10Mb = as.numeric(mostcom_10Mb_dt$mostcom_10Mb)
mostcom_10Mb_dt$start = as.numeric(mostcom_10Mb_dt$start)
mostcom_10Mb_dt$start = mostcom_10Mb_dt$start+1
mostcom_10Mb_dt$end = as.numeric(mostcom_10Mb_dt$end)

# setkey(mostcom_10Mb_dt, chrom, start, end)
# comp_mostcom_10Mb_dt = foverlaps(comp_dt, mostcom_10Mb_dt, type="any", mult="all") 

setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g = foverlaps(mostcom_10Mb_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 100399
setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(mostcom_10Mb_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 55843

par(mfrow=c(3,2))
# all genes 
plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
     x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_10Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 10Mb")
plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF379NOY,
     x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_10Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 10Mb")
# expressed genes 
plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$ASE_ratio_ENCFF675NTU,
     x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$mostcom_10Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 10Mb")
plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$ASE_ratio_ENCFF379NOY,
     x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$mostcom_10Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 10Mb")

### ase expressed genes 
plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF675NTU,
     x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_10Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 10% most common contacts 10Mb")
plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF379NOY,
     x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_10Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF379NOY",
     xlab = "sum 10% most common contacts 10Mb")


###########################################################
###########################################################
# most differential contacts 

# sum_sliding_zmostdiff1_10mb_noNA = ssim_sliding_10mb_noNA[[2]]
# sum_sliding_zmostdiff2_10mb_noNA = ssim_sliding_10mb_noNA[[3]]

sum_sliding_zmostdiff1_1mb_noNA_df = lapply(seq_along(sum_sliding_zmostdiff1_1mb_noNA), function(chrom_idx){
    chrom = cbind(sum_sliding_zmostdiff1_1mb_noNA[[chrom_idx]],
                  paste("chr",chrom_idx,sep=""),
                  as.numeric(names(sum_sliding_zmostdiff1_1mb_noNA[[chrom_idx]])),
                  as.numeric(names(sum_sliding_zmostdiff1_1mb_noNA[[chrom_idx]]))+50000)
    colnames(chrom) = c("mostdiff1_1Mb","chrom","start","end")
    rownames(chrom) <- NULL
    return(chrom)
})
sum_sliding_zmostdiff1_1mb_noNA_df = do.call(rbind,sum_sliding_zmostdiff1_1mb_noNA_df)
sum_sliding_zmostdiff2_1mb_noNA_df = unlist(sum_sliding_zmostdiff2_1mb_noNA)
sum_sliding_zmostdiff_1mb_noNA_df = cbind(sum_sliding_zmostdiff1_1mb_noNA_df,sum_sliding_zmostdiff2_1mb_noNA_df)
colnames(sum_sliding_zmostdiff_1mb_noNA_df)[5] = "mostdiff2_1Mb"
mostdiff_1Mb_dt = as.data.table(sum_sliding_zmostdiff_1mb_noNA_df)
mostdiff_1Mb_dt$mostdiff1_1Mb = as.numeric(mostdiff_1Mb_dt$mostdiff1_1Mb)
mostdiff_1Mb_dt$mostdiff2_1Mb = as.numeric(mostdiff_1Mb_dt$mostdiff2_1Mb)
mostdiff_1Mb_dt$start = as.numeric(mostdiff_1Mb_dt$start)
mostdiff_1Mb_dt$start = mostdiff_1Mb_dt$start+1
mostdiff_1Mb_dt$end = as.numeric(mostdiff_1Mb_dt$end)
mostdiff_1Mb_dt$mostdiff_1Mb = mostdiff_1Mb_dt$mostdiff1_1Mb + mostdiff_1Mb_dt$mostdiff2_1Mb

# setkey(mostcom_1Mb_dt, chrom, start, end)
# comp_mostcom_1Mb_dt = foverlaps(comp_dt, mostcom_1Mb_dt, type="any", mult="all") 

setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g = foverlaps(mostdiff_1Mb_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 103775
setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(mostdiff_1Mb_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 57168

par(mfrow=c(3,2))
# all genes 
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 2")
# expressed genes 
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 2")

### ase expressed genes 
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 2")

# just chr 1
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$chrom=="chr1"),]$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$chrom=="chr1"),]$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$chrom=="chr1"),]$ASE_ratio_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$chrom=="chr1"),]$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, 
     ylab = "ASE_ratio_ENCFF675NTU",
     xlab = "sum 5% most differential contacts 1Mb - Genome 2")
### 4 MB
### down here not changed yet, just copied
# 
# sum_sliding_zcom_4Mb_noNA_df = lapply(seq_along(sum_sliding_zcom_4mb_noNA), function(chrom_idx){
#     chrom = cbind(sum_sliding_zcom_4mb_noNA[[chrom_idx]],
#                   paste("chr",chrom_idx,sep=""),
#                   as.numeric(names(sum_sliding_zcom_4mb_noNA[[chrom_idx]])),
#                   as.numeric(names(sum_sliding_zcom_4mb_noNA[[chrom_idx]]))+50000)
#     colnames(chrom) = c("mostcom_4Mb","chrom","start","end")
#     rownames(chrom) <- NULL
#     return(chrom)
# })
# sum_sliding_zcom_4Mb_noNA_df = do.call(rbind,sum_sliding_zcom_4Mb_noNA_df)
# mostcom_4Mb_dt = as.data.table(sum_sliding_zcom_4Mb_noNA_df)
# mostcom_4Mb_dt$mostcom_4Mb = as.numeric(mostcom_4Mb_dt$mostcom_4Mb)
# mostcom_4Mb_dt$start = as.numeric(mostcom_4Mb_dt$start)
# mostcom_4Mb_dt$start = mostcom_4Mb_dt$start+1
# mostcom_4Mb_dt$end = as.numeric(mostcom_4Mb_dt$end)
# 
# # setkey(mostcom_4Mb_dt, chrom, start, end)
# # comp_mostcom_4Mb_dt = foverlaps(comp_dt, mostcom_4Mb_dt, type="any", mult="all") 
# 
# setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
# mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g = foverlaps(mostcom_4Mb_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 100399
# setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
# mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(mostcom_4Mb_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 55843
# 
# par(mfrow=c(3,2))
# # all genes 
# plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
#      x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_4Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF675NTU",
#      xlab = "sum 10% most common contacts 4Mb")
# plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF379NOY,
#      x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_4Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF379NOY",
#      xlab = "sum 10% most common contacts 4Mb")
# # expressed genes 
# plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$ASE_ratio_ENCFF675NTU,
#      x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$mostcom_4Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF675NTU",
#      xlab = "sum 10% most common contacts 4Mb")
# plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$ASE_ratio_ENCFF379NOY,
#      x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$mostcom_4Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF379NOY",
#      xlab = "sum 10% most common contacts 4Mb")
# 
# ### ase expressed genes 
# plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF675NTU,
#      x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_4Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF675NTU",
#      xlab = "sum 10% most common contacts 4Mb")
# plot(y = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF379NOY,
#      x = mostcom_4Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_4Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF379NOY",
#      xlab = "sum 10% most common contacts 4Mb")
# 
# ### 10 Mb
# 
# 
# sum_sliding_zcom_10Mb_noNA_df = lapply(seq_along(sum_sliding_zcom_10mb_noNA), function(chrom_idx){
#     chrom = cbind(sum_sliding_zcom_10mb_noNA[[chrom_idx]],
#                   paste("chr",chrom_idx,sep=""),
#                   as.numeric(names(sum_sliding_zcom_10mb_noNA[[chrom_idx]])),
#                   as.numeric(names(sum_sliding_zcom_10mb_noNA[[chrom_idx]]))+50000)
#     colnames(chrom) = c("mostcom_10Mb","chrom","start","end")
#     rownames(chrom) <- NULL
#     return(chrom)
# })
# sum_sliding_zcom_10Mb_noNA_df = do.call(rbind,sum_sliding_zcom_10Mb_noNA_df)
# mostcom_10Mb_dt = as.data.table(sum_sliding_zcom_10Mb_noNA_df)
# mostcom_10Mb_dt$mostcom_10Mb = as.numeric(mostcom_10Mb_dt$mostcom_10Mb)
# mostcom_10Mb_dt$start = as.numeric(mostcom_10Mb_dt$start)
# mostcom_10Mb_dt$start = mostcom_10Mb_dt$start+1
# mostcom_10Mb_dt$end = as.numeric(mostcom_10Mb_dt$end)
# 
# # setkey(mostcom_10Mb_dt, chrom, start, end)
# # comp_mostcom_10Mb_dt = foverlaps(comp_dt, mostcom_10Mb_dt, type="any", mult="all") 
# 
# setkey(anno_total_ase_imp_merged_ENCODE_g_sort_dt, chrom, start_gene, end_gene)
# mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g = foverlaps(mostcom_10Mb_dt, anno_total_ase_imp_merged_ENCODE_g_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 100399
# setkey(anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, chrom, start_gene, end_gene)
# mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp = foverlaps(mostcom_10Mb_dt, anno_total_ase_imp_merged_ENCODE_g_ase_exp_sort_dt, type="any", mult="all",by.x = c("chrom","start","end")) # 55843
# 
# par(mfrow=c(3,2))
# # all genes 
# plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF675NTU,
#      x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_10Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF675NTU",
#      xlab = "sum 10% most common contacts 10Mb")
# plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$ASE_ratio_ENCFF379NOY,
#      x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_10Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF379NOY",
#      xlab = "sum 10% most common contacts 10Mb")
# # expressed genes 
# plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$ASE_ratio_ENCFF675NTU,
#      x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF675NTU==TRUE),]$mostcom_10Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF675NTU",
#      xlab = "sum 10% most common contacts 10Mb")
# plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$ASE_ratio_ENCFF379NOY,
#      x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g[which(mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g$expressed_ENCFF379NOY==TRUE),]$mostcom_10Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF379NOY",
#      xlab = "sum 10% most common contacts 10Mb")
# 
# ### ase expressed genes 
# plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF675NTU,
#      x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_10Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF675NTU",
#      xlab = "sum 10% most common contacts 10Mb")
# plot(y = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$ASE_ratio_ENCFF379NOY,
#      x = mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g_ase_exp$mostcom_10Mb,
#      pch = 19,col = adjustcolor("black", alpha.f = 0.5),
#      cex = 0.1, 
#      ylab = "ASE_ratio_ENCFF379NOY",
#      xlab = "sum 10% most common contacts 10Mb")
# 

##################################################################
# gene expression at most common, most different

#mostcom_10Mb_anno_total_ase_imp_merged_ENCODE_g
par(mfrow = c(2,2))
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF675NTU,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF675NTU",
     xlab = "sum zscore 10% most common contacts 1Mb")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF675NTU",
     xlab = "sum zscore 5% most differential contacts 1Mb - sum both Genomes")

plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF675NTU",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF675NTU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF675NTU",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 2")
###
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF379NOY,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF379NOY",
     xlab = "sum zscore 10% most common contacts 1Mb")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF379NOY,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF379NOY",
     xlab = "sum zscore 5% most differential contacts 1Mb - sum both Genomes")

plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF379NOY,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF379NOY",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF379NOY,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF379NOY",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 2")
###
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF910OBU,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF910OBU",
     xlab = "sum zscore 10% most common contacts 1Mb")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF910OBU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF910OBU",
     xlab = "sum zscore 5% most differential contacts 1Mb - sum both Genomes")

plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF910OBU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF910OBU",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF910OBU,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF910OBU",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 2")
###
plot(y = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF174OMR,
     x = mostcom_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostcom_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1, xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF174OMR",
     xlab = "sum zscore 10% most common contacts 1Mb")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF174OMR,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF174OMR",
     xlab = "sum zscore 5% most differential contacts 1Mb - sum both Genomes")

plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF174OMR,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff1_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF174OMR",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 1")
plot(y = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$log10_TPM_ENCFF174OMR,
     x = mostdiff_1Mb_anno_total_ase_imp_merged_ENCODE_g$mostdiff2_1Mb,
     pch = 19,col = adjustcolor("black", alpha.f = 0.5),
     cex = 0.1,  xlim = c(0,700), ylim = c(0,5),
     ylab = "log10_TPM_ENCFF174OMR",
     xlab = "sum zscore 5% most differential contacts 1Mb - Genome 2")
##################################################################
##################################################################
##################################################################
# all these RNA analysis based results are really depressing, 
# so I want to check in with Teresas results again, although
# those were not allele-specific, they might still give some 
# insight with the gene expression values

##################################################################
##################################################################
##################################################################
# plot big ass plot again including expressed genes, not only ase genes


# H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr 
# anno_total_ase_imp_merged_ENCODE_g_exp_chr 


for (chr in 1:22){
    #   chr=15
    print(chr)
    NPMI_1_chr = NPMI_1[[chr]]
    NPMI_2_chr = NPMI_2[[chr]]
    zscores_common_chr = zscores_common[[chr]] 
    zscores_diff_chr = zscores_diff[[chr]] 
    zscores_diff_1_chr = zscores_diff_1[[chr]]
    zscores_diff_2_chr = zscores_diff_2[[chr]]
    
    # let's calc a good ratio of width vs height. first thing: in order to try to keep the font size equal between chr, keep the height between chr plots, let's say 1000
    # now these 1000 height should fit 6 plots, so we got 166.6667 per plot.
    # these 166.6667 corresponds to 8000000 bp, that is 2.083333e-05 per bp
    p=(1000/6)/8000000
    width_plot = chromsize_hg38[chr,2]*p*1.5
    # I use 1666 height for 10 panels, so 166.6 per panel
    panel_height = 1000/6
    panels = 16
    if(chr==22){
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_expressed_genes_level_enhancers_SSIM_slidingsum_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*(panels-1))#,res=350) #1000 # # with 9 panels 1500
    }else{
        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_expressed_genes_level_enhancers_SSIM_slidingsum_chr",chr,".png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667
    } #        png(paste("/Users/jmarkow/Desktop/GAM/Co-Phasing/figures/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_chr_wise/phased_H1_NPMI_comparison_contacts_P1_P2_ase_expression_gene_level_enhancers_only_expr_ase_chr",chr,"larger.png",sep=""),width=2*width_plot,height=2*panel_height*panels)#,res=350) #1000 # with 10 panels 1667
    
    par(mar=c(3,4.1,0.1,0.1)) #bltr # margins
    par(mfrow=c(panels,1))
    if(chr==22){par(mfrow=c((panels-1),1))}
    plot_NPMI_w_both_t(NPMI_1_chr,NPMI_2_chr)
    abline(h=0,col="grey75")
    
    plot(y = ssim_sliding_NPMI_1mb_noNA[[chr]],
         x = as.numeric(names(ssim_sliding_NPMI_1mb_noNA[[chr]])),
         pch = 19,
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,1))
    points(y = ssim_sliding_NPMI_4mb_noNA[[chr]],
           x = as.numeric(names(ssim_sliding_NPMI_4mb_noNA[[chr]])), pch = 19, col = "hotpink3")
    points(y = ssim_sliding_NPMI_10mb_noNA[[chr]],
           x = as.numeric(names(ssim_sliding_NPMI_10mb_noNA[[chr]])), pch = 19, col = "deepskyblue3")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("black","hotpink3","deepskyblue3"), pch=19,legend = c("1Mb SSIM","4Mb SSIM","10Mb SSIM"),bty="n")
    
    plot_zscores_w(zscores_diff_chr,title="contact difference")
    abline(h=0,col="grey75")
    
    plot(y = sum_sliding_zdiff1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zdiff1mb_noNA[[chr]])),
         pch = 19,
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zdiff10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zdiff4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zdiff4mb_noNA[[chr]])), pch = 19, col = "hotpink3")
    points(y = sum_sliding_zdiff10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zdiff10mb_noNA[[chr]])), pch = 19, col = "deepskyblue3")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("black","hotpink3","deepskyblue3"), pch=19,legend = c("1Mb sum differences","4Mb  sum differences","10Mb  sum differences"),bty="n")
    
    plot_zscores_w(zscores_common_chr,title="most common contacts")
    abline(h=0,col="grey75")
    
    plot(y = sum_sliding_zcom_1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zcom_1mb_noNA[[chr]])),
         pch = 19,
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zcom_10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zcom_4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zcom_4mb_noNA[[chr]])), pch = 19, col = "hotpink3")
    points(y = sum_sliding_zcom_10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zcom_10mb_noNA[[chr]])), pch = 19, col = "deepskyblue3")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("black","hotpink3","deepskyblue3"), pch=19,legend = c("1Mb sum 10% most common","4Mb  sum 10% most common","10Mb  sum 10% most common"),bty="n")
    
    plot_zcores_w_both(zscores_diff_1_chr,zscores_diff_2_chr,title="most differential contacts")
    abline(h=0,col="grey75")
    
    plot(y = sum_sliding_zmostdiff1_1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zmostdiff1_1mb_noNA[[chr]])),
         pch = 19, col = "burlywood1",
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zmostdiff1_10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zmostdiff1_4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff1_4mb_noNA[[chr]])), pch = 19, col = "darkorange2")
    points(y = sum_sliding_zmostdiff1_10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff1_10mb_noNA[[chr]])), pch = 19, col = "darkorange4")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("burlywood1","darkorange2","darkorange4"), pch=19,legend = c("1Mb sum 5% most different Genome 1","4Mb  sum 5% most different Genome 1","10Mb  sum 5% most different Genome 1"),bty="n")
    
    plot(y = sum_sliding_zmostdiff2_1mb_noNA[[chr]],
         x = as.numeric(names(sum_sliding_zmostdiff2_1mb_noNA[[chr]])),
         pch = 19,col = "lavender",
         xaxt="none",xlab="",
         xlim=c(0,chromsize_hg38[chr,2]),
         ylim=c(0,max(sum_sliding_zmostdiff2_10mb_noNA[[chr]],na.rm=T)))
    points(y = sum_sliding_zmostdiff2_4mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff2_4mb_noNA[[chr]])), pch = 19, col = "mediumpurple2")
    points(y = sum_sliding_zmostdiff2_10mb_noNA[[chr]],
           x = as.numeric(names(sum_sliding_zmostdiff2_10mb_noNA[[chr]])), pch = 19, col = "purple4")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    legend("bottomleft",col = c("lavender","mediumpurple2","purple4"), pch=19,legend = c("1Mb sum 5% most different Genome 2","4Mb  sum 5% most different Genome 2","10Mb  sum 5% most different Genome 2"),bty="n")
    
    # enhancers and enhancer contact
    # H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr
    # H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr
    # enhancer contacts
    # max_y=max(abs(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mid_y))
    # plot(x =  c(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$enh_mid,H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$gene_mid_anno),
    #      y =  rep(0,times=2*nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]])),
    #      pch = 19,
    #      col = c(rep("orange",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]])),rep("cyan2",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]))),
    #      ylim = c((-0.2)*max_y,max_y),
    #      xlim=c(0,chromsize_hg38[chr,2]),
    #      xaxt="none",xlab="",
    #      ylab="H1 enhancers")
    # apply(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]],1,function(edge){
    #     x=c(edge['enh_mid'],edge['gene_mid_anno'],edge['mid_x'])
    #     y=c(0,0,edge['mid_y'])
    #     polygon(x,y,border=edge['col.x'])})
    
    max_y=max(abs(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mid_y))
    plot(x =  c(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$enh_mid,H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$gene_mid_anno),
         y =  rep(0,times=2*nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]])),
         pch = 19,
         col = c(rep("orange",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]])),rep("cyan2",times=nrow(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]))),
         ylim = c((-0.2)*max_y,max_y),
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="H1 enhancers")
    apply(H1_enh_gene_anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]],1,function(edge){
        x=c(edge['enh_mid'],edge['gene_mid_anno'],edge['mid_x'])
        y=c(0,0,edge['mid_y'])
        polygon(x,y,border=edge['col.x'])})
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    points(x=H1_enh_chr[[chr]]$enh_mid,
           y=rep((-0.1)*max_y,times=nrow(H1_enh_chr[[chr]])),
           pch=19, 
           xlim=c(0,chromsize_hg38[chr,2]),
           xaxt="none",xlab="",
           ylab="H1 enhancers")
    legend("bottomleft",col=c("black","orange","cyan2"), pch=19,legend = c("all enhancers","interacting enhancers","interacting genes"),bty="n")
    legend("left",col=c("black","grey75"),lty=1,lwd=3,legend=c("protein coding gene", "other gene type"),bty="n")
    
    # imprinted genes
    chrname = paste("chr",chr,sep="")
    if( chrname%in% names(imprinted_genes_chr)){
        rect(imprinted_genes_chr[[chrname]]$start,(-0.15)*max_y,imprinted_genes_chr[[chrname]]$end,(-0.2)*max_y,col=imprinted_genes_chr[[chrname]]$col)
    }
    
    # total expression gene level
    # anno_total_ase_imp_merged_ENCODE_g_exp_chr
    # plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
    #      y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$log10_TPM_ENCFF379NOY,
    #      pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
    #      xlim=c(0,chromsize_hg38[chr,2]),
    #      xaxt="none",xlab="",
    #      ylab="total gene expression")
    
    plot(x=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$log10_TPM_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="total gene expression")
    
    abline(h=0.3,col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("gene expression in log10(TPM+1)","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # ASE ratio
    # plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
    #      y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$ASE_ratio_ENCFF379NOY,
    #      pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
    #      xlim=c(0,chromsize_hg38[chr,2]),
    #      xaxt="none",xlab="",
    #      ylab="ASE ratio")
    plot(x=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$ASE_ratio_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="ASE ratio")
    abline(h=c(1/3,0.5,2/3),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("ASE ratio","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # LFC
    # plot(x=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$mids,
    #      y=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$log2foldchange_ENCFF379NOY,
    #      pch=19, col=anno_total_ase_imp_merged_ENCODE_g_ase_exp_chr[[chr]]$col,
    #      xlim=c(0,chromsize_hg38[chr,2]),
    #      xaxt="none",xlab="",
    #      ylab="L2FC")
    plot(x=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$mids,
         y=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$log2foldchange_ENCFF379NOY,
         pch=19, col=anno_total_ase_imp_merged_ENCODE_g_exp_chr[[chr]]$col,
         xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="L2FC")
    abline(h=c(-1,0,1),col="grey75")
    legend("bottomleft",col=c("black","red","black",adjustcolor("black",alpha.f = 0.3)), pch=19,legend = c("Log2foldchange","ASE gene","protein coding gene", "other gene type"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    
    
    
    # AB compartments    # somehow I only have this up until chr 21
    if(chr<22){
        plot(x=comp_chr[[chr]]$mid,
             y=comp_chr[[chr]]$Genome1_score1,
             type = "l",
             xlim=c(0,chromsize_hg38[chr,2]),
             #ylim = c(-7,7),
             xaxt="none",xlab="",
             ylab="AB compartments")#ylab="genome score")
        lines(x=comp_chr[[chr]]$mid,
              y=comp_chr[[chr]]$Genome2_score1,col="grey80")
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome1_score1,
               pch=19, col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome1_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=comp_chr[[chr]]$Genome2_score1,
               pch=19,  col=adjustcolor(c("limegreen","maroon")[comp_chr[[chr]]$Genome2_AB1],alpha.f=1))
        points(x=comp_chr[[chr]]$mid,
               y=rep(0,times=nrow(comp_chr[[chr]])),
               pch=19, col=adjustcolor(c("limegreen","lightskyblue2","maroon")[comp_chr[[chr]]$combined_comp],alpha.f=1))
        legend("bottomleft",col=c("limegreen","lightskyblue2","maroon"), pch=19,legend = c("A/AA","AB","B/BB"),bty="n")
        axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    }
    
    # SNP density
    plot(x=SNPs_H1_chr_binned_col[[chr]]$bin_mids,
         y=SNPs_H1_chr_binned_col[[chr]]$value,
         col=SNPs_H1_chr_binned_col[[chr]]$colour_100lim,
         pch=19, cex=2, xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="SNPs per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("SNP density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    # gene density # this is 38, let'S take 29 since all other results are based on this
    plot(x=annoV29_df_genes_chr_binned_col[[chr]]$bin_mids,
         y=annoV29_df_genes_chr_binned_col[[chr]]$value,
         col=annoV29_df_genes_chr_binned_col[[chr]]$colour_10lim,
         pch=19, cex=2, xlim=c(0,chromsize_hg38[chr,2]),
         xaxt="none",xlab="",
         ylab="genes per window")
    legend("bottomleft",col=c("black"), pch=19,legend = c("Gene density per 50kb bin"),bty="n")
    axis(1,at=seq(10e6,250e6,10e6),labels=paste(seq(10,250,10),"Mb",sep=""))
    
    dev.off()
}
###t ake care, the ase genes thing is not only protein coding genes probably, but everything there is

##################################################################
##################################################################
# count the number of expressed genes in sliding windows
# count number of ase genes in sliding windows?
# count number of ase / expressed genes in AA, AB, BA, BB compartments
# concentrate on protein coding genes

