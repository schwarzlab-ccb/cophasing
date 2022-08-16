### Comparison H1 haplotype-specific A/B compartments and unphased expression - 01.03.2021
# get a feeling for the compartment data from Sasha
# get a feeling forth expression data from Teresa
# figure out how to compare these, 
# I can imagine a boxplot with 4 categories A/A, A/B, B/B over expression
#
#
#-------------------------------------------------------------------------------
### load data
# compartments
comp=read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/compartments/210222.Compartments.AB.H1.phased.first.PCA.csv",header = T, sep="\t")

# expression
gene_exp=read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/expression/genes_with_expression.bed",header = T, sep="\t")
gro_50=read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/expression/groseq_intersected_50kb.bed",header = F, sep="\t")
gro_100=read.csv("/fast/groups/ag_schwarz/Projects/project-gam/data/Co-Phasing/expression/groseq_intersected_100kb.bed",header = F, sep="\t")

#-------------------------------------------------------------------------------
### clean data, give colnames, sort by chrom and pos, remove X chromosome, check all available entries in the columns
# compartments
comp = comp[order(comp$start),]
comp = comp[order(as.numeric(substring(comp$chrom, 4))),]
dim(comp) # 11749     7
unique(comp$chrom)
    # "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21"
comp$combined_comp=paste(comp$Genome1_AB1,comp$Genome2_AB1,sep='')
table(comp$combined_comp)/nrow(comp)*100
    #       AA       AB       BA       BB
    # 32.75172 18.59733 21.46566 27.18529

# expression
colnames(gro_50)=c("chrom","start","end","expression")
colnames(gro_100)=c("chrom","start","end","expression")

unique(gro_50$chrom) # unique(gro_100$chrom)
    #  "chr1"  "chr2"  "chr3"  "chr4"  "chr5"  "chr6"  "chr7"  "chr8"  "chr9"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX"

# dim(gro_50) # 60630     4
gro_50=gro_50[which(gro_50$chrom!="chrX"),]
# dim(gro_50) # 57509     4
# dim(gro_100) # 30321     4
gro_100=gro_100[which(gro_100$chrom!="chrX"),]
# dim(gro_100) # 28760     4

# boxplot all expression values within vs mean value within compartment windowsize
# check if all genes are protein coding or not, differentiate here

# intersect gro seq with compartments



require(data.table)
x = data.table(chr=c("Chr1", "Chr1", "Chr2", "Chr2", "Chr2"),
               start=c(5,10, 1, 25, 50), end=c(11,20,4,52,60))
y = data.table(chr=c("Chr1", "Chr1", "Chr2"), start=c(1, 15,1),
               end=c(4, 18, 55), geneid=letters[1:3])
setkey(y, chr, start, end)
foverlaps(x, y, type="any", which=TRUE)
foverlaps(x, y, type="any", which=TRUE, mult="all")
foverlaps(x, y, type="any")
foverlaps(x, y, type="any", nomatch=NULL)
foverlaps(x, y, type="within", which=TRUE)
foverlaps(x, y, type="within")
foverlaps(x, y, type="start")


comp_dt=as.data.table(comp)
# just so that intervals like 15-16 and 16-18 do not overlap
comp_dt$start = comp_dt$start + 1 
#gene_exp_dt=as.data.table(gene_exp)
gro_50_dt=as.data.table(gro_50)
gro_50_dt$start = gro_50_dt$start + 1 
setkey(gro_50_dt, chrom, start, end)
comp_gro50=foverlaps(comp_dt, gro_50_dt, type="any", mult="all") # , minoverlap=2 not yet implemented

gro_100_dt=as.data.table(gro_100)
gro_100_dt$start = gro_100_dt$start + 1 
setkey(gro_100_dt, chrom, start, end)
comp_gro100=foverlaps(comp_dt, gro_100_dt, type="any", mult="all") 

saveRDS(comp_gro50,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/comp_gro50.rds")
saveRDS(comp_gro100,"/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/comp_gro100.rds")

# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/comp_gro50.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/comp_gro50.rds
# scp bih:/fast/groups/ag_schwarz/Projects/project-gam/Co-Phasing/results/comp_gro100.rds /Users/jmarkow/Desktop/GAM/Co-Phasing/results/comp_gro100.rds

comp_gro50 = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/comp_gro50.rds")
comp_gro100 = readRDS("/Users/jmarkow/Desktop/GAM/Co-Phasing/results/comp_gro100.rds")
# class(comp_gro50$combined_comp) # "character"
comp_gro50$combined_comp=factor(comp_gro50$combined_comp, levels = c("AA","AB","BA","BB"))
# class(comp_gro50$combined_comp) # "factor"
class(comp_gro50$expression) # "integer"
#boxplot(comp_gro50,formula=expression~combined_comp, data=comp_gro50)
#Error in x[floor(d)] + x[ceiling(d)] : non-numeric argument to binary operator
comp_gro50_red=comp_gro50[,c("expression","combined_comp")]
boxplot(comp_gro50_red$expression~comp_gro50_red$combined_comp)
# looks promising
# without outliers
boxplot(comp_gro50_red$expression~comp_gro50_red$combined_comp, outline=F) 
# boojachacka

# significance!
# log?
boxplot(log(comp_gro50_red$expression)~comp_gro50_red$combined_comp)
# also ok I guess

comp_gro50_red_AB=comp_gro50_red
comp_gro50_red_AB$combined_comp[which(comp_gro50_red_AB$combined_comp=="BA")]="AB"
comp_gro50_red_AB$combined_comp=factor(comp_gro50_red_AB$combined_comp, levels = c("AA","AB","BB"))
boxplot(comp_gro50_red_AB$expression~comp_gro50_red_AB$combined_comp, outline=F,notch=T) 

### tukey test
# library
install.packages("multcompView")
library(multcompView)
# What is the effect of the treatment on the value ?
model=lm( comp_gro50_red_AB$expression~comp_gro50_red_AB$combined_comp )
ANOVA=aov(model)

# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=ANOVA, 'comp_gro50_red_AB$combined_comp', conf.level=0.95)
TUKEY
    # Tukey multiple comparisons of means
    # 95% family-wise confidence level
    # 
    # Fit: aov(formula = model)
    # 
    # $`comp_gro50_red_AB$combined_comp`
    # diff       lwr        upr     p adj
    # AB-AA -374.70226 -490.3785 -259.02599 0.0000000
    # BB-AA -456.06430 -583.4694 -328.65916 0.0000000
    # BB-AB  -81.36204 -203.3805   40.65638 0.2619061

# Tuckey test representation :
plot(TUKEY , las=1 , col="brown")

# include n, number of entries



