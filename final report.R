# AS.410.671.82.SU21 - Final Exam
# Author: Melissa Chua Wan Ying
# *** Primary Reference *** : https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS1962

# Title   :	
#	 - Neuronal and glioma-derived stem cell factor (SCF) expression on brain angiogenesis
# Summary :	
#	- Analysis of stem cell factor (SCF) overexpression by gliomas of different 
# - grades and its role in tumor angiogenesis in the brain.
# Organism:	
#	- Homo Sapiens
# Platform:	 
#	- GPL570: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# Citation:	
# - Sun L, Hui AM, Su Q, Vortmeyer A et al. 
# - "Neuronal and glioma-derived stem cell factor induces angiogenesis within the brain."
# - Cancer Cell 2006 Apr;9(4):287-300. 
# - PMID: 16616334
# - Reference series:	GSE4290


# Load packages
library(GEOquery)
library(tidyverse)
library(plyr)
library(gplots)
library(MASS)
library(outliers)
library(kernlab)
library(limma)
library(gplots)
library(ggplot2)
library(data.table)
library(pwr)
library(ssize)
library(gdata) 
library(e1071)
library(affyio)
library(viridis)
library(affy)
library(scatterplot3d) 
library(Rcpp)

# Set working directory path 
setwd("C:/Users/User/Dropbox/My PC (DESKTOP-VL7FAF0)/Documents/R/final")

# Source for svm 
source('msvmRFE.R')

# Open GDS1962 dataset 
gds <- getGEO("GDS1962")
eset <- GDS2eSet(gds, do.log2=F, getGPL=F)
dat <- exprs(eset)
dat <- as.data.frame(dat)

# Data information
dat.info <- data.frame(eset$sample, eset$disease.state, eset$tissue)

# Rearrange sample column positions based on glioma grades
dat <- cbind(dat[1:23], dat[24:30], dat[131:168], dat[31:49], dat[169:180], dat[50:130])
dat.info <- dat.info[order(match(dat.info$eset.sample, colnames(dat))),]
label <- as.factor(c(rep("NT", 23-1+1), rep("II", 68-24+1), rep("III", 99-69+1), rep("IV", 180-100+1)))
dat.info <- cbind(dat.info, label)

# Append glioma grades to corresponding samples
for(a in 1:23) {names(dat)[a] <- paste("NB", names(dat[a]), sep="_")}
for(b in 24:68) {names(dat)[b] <- paste("II", names(dat[b]), sep="_")}
for(c in 69:99) {names(dat)[c] <- paste("III", names(dat[c]), sep="_")}
for(d in 100:180) {names(dat)[d] <- paste("IV", names(dat[d]), sep="_")}

# Update dat.info sample names
dat.info$eset.sample <- colnames(dat)

# Corresponding annotation files
ann.file <- "GPL570.annot"
ann <- read.delim(ann.file, header=T, row.names=1, skip=27, sep="\t")
ann <- ann[1:nrow(dat), ]
gene.info <- data.frame(Description = ann$Gene.title, Symbol = ann$Gene.symbol)
rownames(gene.info) <- rownames(ann)
ann <- ann[!(is.na(ann$Gene.title) | ann$Gene.title==""),]
dat <- dat[rownames(dat) %in% rownames(ann),]

samp.matrix <- data.matrix(dat)
rownames(samp.matrix) <- rownames(ann)

samp.count <- ncol(samp.matrix)
id.count <- nrow(samp.matrix)

# Descriptive Satistics
stdev    <- apply(samp.matrix, 1, sd, na.rm = TRUE)
row.means <- rowMeans(samp.matrix, na.rm = TRUE)

# Histograms illustrating initial spread
png("Histogram_rowMeans.png")
hist(row.means, col="#97B9CE", xlab="Mean", ylab="Frequency", main = paste0("Histogram of mean expression values for",id.count,"profiles"))
dev.off()

png("Histogram_data.stdev.png")
hist(stdev, col="#7EF9FF", xlab="Standard Deviation", ylab="Frequency", main=paste0("Histogram of standard deviation expression values\n for ",id.count," gene profiles"))
dev.off()

########################
#   Outlier Analysis   #
########################

# Correlation Matrix
dat.cor  <- cor(samp.matrix, method="pearson", use="pairwise.complete.obs")

# Correlation plot - heatmap
cx <- colorRampPalette(viridis::mako(50))(50)
png("Heatmap_cor.matrix.png")
par(oma=c(0,0,1,0), mar=c(6,6,4,7), par(xpd=TRUE))
leg <- seq(from=0.1, to=1, length.out=10)
image(dat.cor, main="Correlation between Glioma vs Non-Tumor\n Gene Expression", col=cx, axes=F)
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]], cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]], cex.axis=0.9,las=2)
legend(1.1, 1.1, title="Correlation", legend=leg, fill=colorRampPalette(viridis::mako(50))(10))
dev.off()

#Hierarchical clustering dendrogram 
png("Dendrogram.png")
dat.dist <- dist(t(dat),method="euclidean") #calculate distance
dat.clust <- hclust(dat.dist,method="single") #calculate clusters
plot(dat.clust,labels=colnames(dat), cex=0.5, col="#132842",xlab="Distance", ylab="Height")
dev.off()

# CV-mean plot
dat.mean <- apply(log2(samp.matrix), 2, mean, na.rm=T) 
stdev.log <- sqrt(apply(log2(samp.matrix), 2, sd, na.rm=T))
cv <- stdev.log / dat.mean

png("Scatterplot_ColMeansCV.png")
plot(dat.mean, cv, xlab="log2(Mean)", ylab="log2(CV)", main="CV-Mean plot of Glioma Grades vs Non-Tumor", cex=1.5, type="n")
points(dat.mean, cv, col=c(rep("#c1a87d", length(dat[1:23])), rep("#7ef9ff", length(dat[24:68])), rep("#97b9ce", length(dat[69:99])), rep("#1f456e", length(dat[100:180]))), pch=c(rep(7, samp.count/4), rep(10, samp.count/4), rep(17, samp.count/4), rep(24, samp.count/4)))
legend("topright", c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), pch=c(7,10,17,24), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), bg="white", cex=0.6)
text(dat.mean, cv, labels=names(dat.mean), cex=0.5, offset=10, col="#132842")
dev.off()

# Average correlation plot
dat.avg <- apply(dat.cor, 1, mean, na.rm=T)
png("Scatterplot_CorMeans.png")
par(mfrow=c(2,1))
plot(c(1,length(dat.avg)), range(dat.avg), type="n", xlab="", ylab="Average correlation", main="Avg correlation for Glioma vs Non-Tumor samples", axes=F)
points(dat.avg, col=c(rep("#b58d3d", length(dat[1:23])), rep("#1f456e", length(dat[24:180]))), pch=c(rep(15, length(dat[1:23])), rep(21, length(dat[24:180]))))
axis(1, at=c(1:length(dat.avg)), labels=dimnames(samp.matrix)[[2]], las=2, cex.lab=0.4, cex.axis=0.6)
axis(2)
grid(nx=16, col="grey")
legend("bottomright", cex=0.8, c("Non-Tumor", "Glioma"), pch=c(15, 21), col=c("#b58d3d", "#1f456e"))

# Comparison between glioma grades
plot(c(1,length(dat.avg)), range(dat.avg), type="n", xlab="", ylab="Average correlation", main="Average Correlation for Various \n Glioma Grades vs Non-Tumor samples", axes=F)
points(dat.avg, col=c(rep("#c1a87d", length(dat[1:23])), rep("#7ef9ff", length(dat[24:68])), rep("#97b9ce", length(dat[69:99])), rep("#1f456e", length(dat[100:180]))), pch=c(rep(7, length(dat[1:23])), rep(10, length(dat[24:68])), rep(17, length(dat[69:99])), rep(24, length(dat[100:180]))))
axis(1, at=c(1:length(dat.avg)), labels=colnames(samp.matrix), las=2, cex.lab=0.4, cex.axis=0.6)
axis(2)
grid(nx=16, col="grey")
legend("bottomright", cex=0.4, c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), pch=c(7,10,17,24), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), bg="white")
dev.off()

# Identifying Outlier(s) via outlier()
o <- dat.avg <=  outlier(dat.avg)
outlier <- dat.avg[o]
cat(sprintf("%s Outlier(s) identified!\n", length(outlier)))
outlier
# IV_GSM97895 
#   0.7730638       

# Remove Outlier(s) 
dat <- dat[, !(colnames(dat) %in% c("II_GSM97878", "II_GSM97913", "IV_GSM97895", "IV_GSM97886"))]


########################
#     Filter Genes     #
########################

quantile(log2(rowMeans(samp.matrix)), na.rm=T)
# 0%         25%        50%       75%        100% 
# 0.9473545  6.5153297  8.0221536  9.5758638 15.0504017 

# Eliminate probes with rowMeans less than 0 on a log2 scale
dat.fil <- log2(subset(samp.matrix, rowMeans(samp.matrix) > 0))
removed <- nrow(samp.matrix) - nrow(dat.fil)
cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed))

# Row Means v. Column Variances on filtered data
fil.mean <- apply(dat.fil, 1, mean) 
fil.sd <- apply(dat.fil, 1, sd)  
f.cv <- fil.sd / fil.mean

# Plot probes with rowMeans more than 0 on a log2 scale
png("Scatterplot_RowMeansCV_1.png")
plot(fil.mean, f.cv, xlab="log2(RowMeans)", ylab="log2(CV)", main="Plot of Row Mean v. Row Variance\n for Filtered Non-Tumor vs Glioma samples", col=c(rep("#c1a87d", length(dat[1:23])), rep("#7ef9ff", length(dat[24:66])), rep("#97b9ce", length(dat[67:97])), rep("#1f456e", length(dat[98:176]))), pch=c(rep(7, length(dat[1:23])), rep(10, length(dat[24:66])), rep(17, length(dat[67:97])), rep(24, length(dat[98:176]))))
legend("topright", c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), pch=c(7,10,17,24), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"))
abline(v=3, col=2, lwd=2) # Threshold determined for stage 2/2 of the filtering process
dev.off()

# Eliminate probes with rowMeans less than 3 on a log2 scale
dat.filtered <- subset(dat.fil, rowMeans(dat.fil) > 3)
dat.filtered <- as.data.frame(dat.filtered)
removed  <- nrow(dat.fil) - nrow(dat.filtered)

cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed))

fil.mean.2 <- apply(dat.filtered, 1, mean) 
fil.sd.2 <- apply(dat.filtered, 1, sd)  
fil.cv.2   <- fil.sd.2 / fil.mean.2

png("Scatterplot_RowMeansCV_2.png")
plot(fil.mean.2, fil.cv.2, xlim=c(0,12), ylim=c(0,50), xlab="log2(RowMeans)", ylab="log2(CV)", main="Filtered CV-Mean plot for Glioma vs Non-Tumor samples", col=c(rep("#c1a87d", length(dat[1:23])), rep("#7ef9ff", length(dat[24:66])), rep("#97b9ce", length(dat[67:97])), rep("#1f456e", length(dat[98:176]))), pch=c(rep(7, length(dat[1:23])), rep(10, length(dat[24:66])), rep(17, length(dat[67:97])), rep(24, length(dat[98:176]))))
legend("topright", c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), pch=c(7,10,17,24), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"))
abline(v=3, col=2, lwd=2)
dev.off()

# Update gene.info data.frame
gene.info <- subset(gene.info, rownames(gene.info) %in% rownames(dat.filtered))

###############################
# SCF differential expression #
###############################
# SCF probe sets
scf.ann <- ann[ann$Gene.symbol=="KITLG",]
scf <- dat.filtered[rownames(dat.filtered) %in% rownames(scf.ann),]
scf.ctl.df <- scf[1:23]
scf.glioma.df <- scf[24:176]

#boxplot of KITLG probe 207029_at
dat.info <- dat.info[dat.info$eset.sample %in% colnames(dat.filtered),]
color=c("#7ef9ff", "#97b9ce", "#1f456e", "#c1a87d")
png("SCF_boxplot_1.png")
kitlg.207029_at <- as.numeric(dat.filtered["207029_at",])
boxplot(kitlg.207029_at~dat.info$label, dat.filtered, ylab='Expression',xlab='Glioma groups',main='Boxplot of KITLG probe 207029_at differential expresison \n between glioma groups',cex.main=0.9,col=color)
dev.off()

#boxplot of KITLG probe 211124_s_at
png("SCF_boxplot_2.png")
kitlg.211124_s_at <- as.numeric(as.numeric(dat.filtered["211124_s_at",]))
boxplot(kitlg.211124_s_at~dat.info$label, dat.filtered, ylab='Expression',xlab='Glioma groups',main='Boxplot of KITLG probe 211124_s_at differential expresison \n between glioma groups',cex.main=0.9,col=color)
dev.off()

#boxplot of KITLG probe 226534_at
png("SCF_boxplot_3.png")
kitlg.226534_at <- as.numeric(dat.filtered["226534_at",])
boxplot(kitlg.226534_at~dat.info$label, dat.filtered, ylab='Expression',xlab='Glioma groups',main='Boxplot of KITLG probe 226534_at differential expresison \n between glioma groups',cex.main=0.9,col=color)
dev.off()

# factors
scf.ctl <- grepl("^NB", colnames(scf))
scf.II <- grepl("^II", colnames(scf))
scf.III <- grepl("^III", colnames(scf))
scf.IV <- grepl("^IV", colnames(scf))

#1-factor ANOVA with 4 levels
aov.all.genes <- function(x,s1,s2,s3,s4) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])
  fac <- c(rep("A",length(x1)), rep("B",length(x2)), rep("C",length(x3)), rep("D",length(x4)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3,x4))
  names(a.dat) <- c("factor","express")
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]    # use to get F-statistic
  return(p.out)
}

# Raw p-value using One-way Anova for SCF differential expression across glioma grades
pv.scf <- apply(scf, 1, aov.all.genes,s1=scf.ctl,s2=scf.II,s3=scf.III,s4=scf.IV)
pv.scf <- sort(pv.scf)

# Plot the sorted raw P-values for SCF gene
png("SCF_Raw_PvaluePlot.png")
par(mfrow=c(1,2))
barplot(pv.scf, cex=0.9, col=c("#7ef9ff", "#97b9ce", "#1f456e"), las=2, cex.axis=0.6, xlab="", ylab="P-values", main="Raw P-values for KITLG gene probes\n", names=c("207029_at", "211124_s_at", "226534_at"))
abline(v=.01,col=2,lwd=2)
barplot(-log10(pv.scf), col=c("#7ef9ff", "#97b9ce", "#1f456e"), xlab="", ylab="log10(p-values)", main="-log10(P-values) for KITLG gene between\nall groups", cex.main=0.9, las=2, names=c("207029_at", "211124_s_at", "226534_at"))
abline(v= -log10(.01),col=2,lwd=2)
dev.off()

########################
#  Feature Selection   #
########################   

type <- lapply(colnames(dat.filtered), 
               function(x) {
                 if(regexpr("Ctl", x) < 1) {"Glioma"} else {"Ctl"}
               })

# Groups
dat.filtered <- as.data.frame((dat.filtered), na.rm=T)
is.ctl <- grepl("^NB", names(dat.filtered))
is.glioma <- grepl("^(?!NB)", names(dat.filtered), perl=T)

# Tissue grades
II <- grepl("^II", names(dat.filtered))
III <- grepl("^III", names(dat.filtered))
IV <- grepl("^IV", names(dat.filtered))

#1-factor ANOVA with 4 levels
aov.all.genes <- function(x,s1,s2,s3,s4) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])
  fac <- c(rep("A",length(x1)), rep("B",length(x2)), rep("C",length(x3)), rep("D",length(x4)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3,x4))
  names(a.dat) <- c("factor","express")
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]    # use to get F-statistic
  return(p.out)
}

#run AVOVA function
pv <- apply(dat.filtered, 1, aov.all.genes,s1=is.ctl,s2=II,s3=III,s4=IV)
pv <- as.data.frame(pv)

# Statistical significance of differential gene expression across glioma grades
png("expression_gliomagrades.png")
par(mfrow=c(1,2))
hist(as.numeric(unlist(pv)), col="#97B9CE", xlab="p-values", main="P-value dist'n between\nglioma grades", cex.main=0.9)
abline(v=.01,col=2,lwd=2)
hist(-log10(as.numeric(unlist(pv))),col="#97B9CE",xlab="log10(p-values)", main="-log10(pv) dist'n between\nglioma grades",cex.main=0.9)
abline(v= -log10(.01),col=2,lwd=2)
dev.off()

# Student's t-test function 
rawp <- c()
ptm <- proc.time()
t.test.all.genes <- function(x,s1,s2) {
  x <- as.numeric(x)
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

# Raw p-value for control and glioma samples
rawp <- apply(dat.filtered, 1, t.test.all.genes, s1=is.ctl, s2=is.glioma)
proc.time() - ptm
# user  system elapsed 
# 8.64    0.06   10.56 

# Output
cat(sprintf("Total number of genes with p-value < 0.05 is %s\n", sum(rawp < 0.05))) # 24946
cat(sprintf("Total number of genes with p-value < 0.01 is %s\n", sum(rawp < 0.01))) # 20351

# Bonferroni correction
Bonferroni <- 0.05/length(rawp)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni))
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(rawp < Bonferroni))) # 8481 

# Holm-Bonferroni correction
p.holm <- p.adjust(rawp, method="holm")
p.raw  <- sort(rawp)
p.holm <- sort(p.holm)
p1     <- as.matrix(p.raw)
p2     <- as.matrix(p.holm)
all.p  <- cbind(p1, p2)
colnames(all.p) <- c("Raw P-value", "Adjusted P-Value")

# Output
cat(sprintf("Total number of genes with p-value < 0.05 is %s\n", sum(p.holm < 0.05))) # 8663 
cat(sprintf("Total number of genes with p-value < 0.01 is %s\n", sum(p.holm < 0.01))) # 7401

Bonferroni.2 <- 0.05/length(p.holm)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni.2))
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(rawp < Bonferroni.2))) # 8481 

# Plot the sorted raw P-values
png("Raw_PvaluePlot.png")
plot(p.raw, type="b", main="Raw P-values for All Genes\n", xlab="Genes", ylab="P-values", pch=1, col="#7ef9ff")
dev.off()

# Plot Adjusted vs. Raw P-values
png("PvaluePlot.png")
matplot(all.p, type="b", pch=1, col=c("#1f456e", "#eee0b1"), xlab="Genes", ylab="P-values", main = "Adjusted Vs. Raw P-values for All Genes\n")
legend("bottomright", legend=colnames(all.p), pch=1, col=c("#1f456e", "#eee0b1"))
dev.off()

# Update gene.info data.frame with P-value information
gene.info$pvalue <- rawp

# Hypothesis Pass/Fail dataframe 
thresh <- Bonferroni.2
rnames <- rownames(dat.filtered)
p.test <- lapply(as.list(rawp), function(f){ if (f<thresh) TRUE else FALSE })
rawp.df <- as.data.frame(do.call(rbind, p.test), rname=rnames)
names(rawp.df) <- c("rawp.Test")
rownames(rawp.df) <- rownames(dat.filtered)

# Update gene.info data.frame with P-value test results (pass/fail)
gene.info$rawp.Test <- rawp.df$rawp.Test

cat(sprintf("Threshold: %s\n", thresh))
table(rawp.df$rawp.Test)
#   Number of True        : 8481        
#   Number of False       : 36201        
nrow(rawp.df)
#   Total                 : 44682 

# Retain only samples that passed the threshold
rawp.pass.df <- subset(rawp.df, rawp.Test=="TRUE", select=rawp.Test)
dat.filtered <- dat.filtered[rownames(dat.filtered) %in% rownames(rawp.pass.df),]
rawp <- as.data.frame(rawp)
rownames(rawp) <- rownames(gene.info)
rawp.pass <- rawp[rownames(rawp) %in% rownames(rawp.pass.df),]
rawp.pass <- as.numeric(unlist(rawp.pass))
rawp <- as.numeric(unlist(rawp))

# Histogram illustrating the distribution of P-values 
png("Histogram_rawp.png")
hist(rawp, col="#1f456e", xlab="P-Value", ylab="Frequency", main="Histogram of T-Test P-Values values for GDS1962 profiles\n prior to feature selection")
hist(p.holm, col="#7ef9ff", add=T)
legend("top", legend = colnames(all.p), pch = 16, col = c("#1f456e", "#7ef9ff"))
dev.off()

# Determining fold change between glioma and non-tumor samples in log2
ctl.matrix  <- as.matrix(dat.filtered[1:23])
ctl.mean    <- apply(ctl.matrix, 1, mean, na.rm=T)
glioma.matrix <- as.matrix(dat.filtered[23:180])
glioma.mean   <- apply(glioma.matrix, 1, mean, na.rm=T)

fold <- glioma.mean - ctl.mean
linear_fold <- 2^(fold) #to convert log scale to linear scale
max(linear_fold) #maximum fold change is 58.95058
min(linear_fold) #minimum fold change is 0.03848994

# Update gene.info data.frame with fold information
gene.info <- gene.info[rownames(gene.info) %in% rownames(dat.filtered),]
gene.info$fold <- fold	

# Determine probesets that demonstrate a 2x fold change
fold.test <- lapply(as.list(fold), function(x){ if (abs(x) > log2(2)) TRUE else FALSE })
fold.df   <- as.data.frame(do.call(rbind, fold.test))
names(fold.df) <- c("Fold.Test")

# Update gene.info data.frame with fold information
gene.info$Fold.Test <- fold.df$Fold.Test

# Number of genes with p-value <0.05 and linear fold >2
table(fold.df$Fold.Test)
# Number of True        : 3904           
# Number of False       : 4577           
nrow(fold.df)
# Total                 : 8481 

# New data.frame with genes that have passed both (Fold and rawp) tests
true.genes <- subset(gene.info, (Fold.Test & rawp.Test) == TRUE)
cat(sprintf("Total number of genes that pass both (rawp and Fold) tests: %s\n", nrow(true.genes))) # 5241

# Write these genes with their corresponding values to an output .txt file
write.table(true.genes, file="TrueGenes.csv", sep=",", col.names=NA, qmethod="double")

# Ordering the highest genes (by P-value) in the form of a data.frame
# Note: dat.filtered is still in log2 scale
best.genes    <- order(rawp.pass)[1:length(rawp.pass)]
best.genes.df <- data.frame(index=best.genes, rawp=rawp.pass[best.genes])
top.genes.matrix <- dat.filtered[best.genes, ]

# Feature Selection via svmRFE which utilizes the library e1071
t.dat  <- t(top.genes.matrix)
ptm <- proc.time()
svm.df <-data.frame(label, t.dat)
ranked.list <- svmRFE(svm.df, k=10, halve.above=100)
proc.time() - ptm
# Run-time:
# user  system elapsed 
# 47.06    0.33   50.25 

# Write the rankings to an output .txt file so that it can be read in later if needed
output <- data.frame(RankedOrder = ranked.list)
write.table(output, file = "RankedList.txt")
top.ranked.genes <- top.genes.matrix[ranked.list, ]
rownames(top.ranked.genes) <- rownames(top.genes.matrix[ranked.list, ])

# Create a new genes.info data.frame for the ranked genes
top.genes.info <- gene.info[rownames(top.ranked.genes ),]
tg <- top.genes.info$pvalue[top.genes.info$pvalue < thresh]

# Histogram
png("Histogram_ranked_rawp.png")
hist(tg, col="#7ef9ff", xlab="P-Value", ylab="Frequency", main="Histogram of T-Test P-Values values for ranked profiles")
abline(v=thresh, col=2, lwd=2)
dev.off()

############################
# Dimensionality reduction #
############################

# MDS Analysis via Kruskal's Non-metric Approach
dat.dist <- dist(t(top.ranked.genes))
dat.mds  <- isoMDS(dat.dist)
# initial  value 17.813889 
# final  value 17.811278 
# converged

# MDS plot
png("MDSplot_Kruskal.png")
par(mfrow=c(3,1))
plot(dat.mds$points, type="n", xlab="MDS Dimension 1", ylab="MDS Dimension 2")
points(dat.mds$points[, 1][as.numeric(label=="NT")==1], dat.mds$points[, 2][as.numeric(label=="NT")==1], col="#c1a87d", pch=7, cex=1.5)
points(dat.mds$points[, 1][as.numeric(label=="II")==1], dat.mds$points[, 2][as.numeric(label=="II")==1], col="#7ef9ff", pch=10, cex=1.5)
points(dat.mds$points[, 1][as.numeric(label=="III")==1], dat.mds$points[, 2][as.numeric(label=="III")==1], col="#97b9ce", pch=17, cex=1.5)
points(dat.mds$points[, 1][as.numeric(label=="IV")==1], dat.mds$points[, 2][as.numeric(label=="IV")==1], col="#1f456e", pch=24, cex=1.5)
title(main="Kruskal's Non-metric MDS map from distance for GDS1962 dataset")
legend("bottomleft", c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), col = c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), pch=c(7,10,17,24), cex=0.7)

# MDS rank order
dat.rank <- matrix((rank(as.matrix(dat.dist))), nrow=180)
dat.mds3 <- isoMDS(dat.rank)
# initial  value 20.826241 
# final  value 20.826241 
# converged

# MDS ranked plot
plot(dat.mds3$points, type="n", xlab="MDS Dimension 1", ylab="MDS Dimension 2")
points(dat.mds3$points[, 1][as.numeric(label=="NT")==1], dat.mds3$points[, 2][as.numeric(label=="NT")==1], col="#c1a87d", pch=7, cex=1.5)
points(dat.mds3$points[, 1][as.numeric(label=="II")==1], dat.mds3$points[, 2][as.numeric(label=="II")==1], col="#7ef9ff", pch=10, cex=1.5)
points(dat.mds3$points[, 1][as.numeric(label=="III")==1], dat.mds3$points[, 2][as.numeric(label=="III")==1], col="#97b9ce", pch=17, cex=1.5)
points(dat.mds3$points[, 1][as.numeric(label=="IV")==1], dat.mds3$points[, 2][as.numeric(label=="IV")==1], col="#1f456e", pch=24, cex=1.5)
title(main="Kruskal's MDS from rank distance for GDS1962 dataset")
legend("topleft", c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), col = c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), pch=c(7,10,17,24), cex=0.45)

# Determine the number of dimensions in Kruskal's non-metric scaling
mds1 <- isoMDS(dat.rank, k=1)  #non-metric MDS in 1 dimension 
mds2 <- isoMDS(dat.rank, k=2)  #non-metric MDS in 2 dimensions 
mds3 <- isoMDS(dat.rank, k=3)  #non-metric MDS in 3 dimensions 
mds4 <- isoMDS(dat.rank, k=4)  #non-metric MDS in 4 dimensions 
mds5 <- isoMDS(dat.rank, k=5)  #non-metric MDS in 5 dimensions 

stress = c(mds1$stress, mds2$stress, mds3$stress, mds4$stress, mds5$stress)
dimensions = 1:5

# Scree plot for the different dimensions for Kruskal's method
plot(dimensions, stress, type = "b", xlab = "Number of Dimensions", ylab = "Stress")
title(main="Determine the number of dimensions in Kruskal's Non-Metric approach")
dev.off()

# Weighted Laplacian Graph
k.speClust2 <- function (X, qnt=NULL) {
  dist2full <- function(dis) {
    n <- attr(dis, "Size")
    full <- matrix(0, n, n)
    full[lower.tri(full)] <- dis
    full + t(full)
  }
  dat.dis <- dist(t(X),"euc")^2
  if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
  if(is.null(qnt)) {eps <- min(dat.dis[dat.dis!=0])}
  kernal <- exp(-1 * dat.dis/(eps))
  K1 <- dist2full(kernal)
  diag(K1) <- 0
  D = matrix(0,ncol=ncol(K1),nrow=ncol(K1))
  tmpe <- apply(K1,1,sum)
  tmpe[tmpe>0] <- 1/sqrt(tmpe[tmpe>0])
  tmpe[tmpe<0] <- 0
  diag(D) <- tmpe
  L <- D%*% K1 %*% D
  X <- svd(L)$u
  Y <- X / sqrt(apply(X^2,1,sum))
}

temp <- t(top.ranked.genes)
temp <- scale(temp, center=T, scale=T)
phi <- k.speClust2(t(temp), qnt=NULL)

png("LaplacianPlot.png")
plot(range(phi[, 1]), range(phi[, 2]), xlab="Phi 1", ylab="Phi 2", main="Weighted Graph Laplacian Plot for GDS1962 dataset\nNon-Tumor vs. Glioma")
points(phi[, 1][as.numeric(label=="NT")==1], phi[, 2][as.numeric(label=="NT")==1], col="#c1a87d", pch=7, cex=1.5)
points(phi[, 1][as.numeric(label=="II")==1], phi[, 2][as.numeric(label=="II")==1], col="#7ef9ff", pch=10, cex=1.5)
points(phi[, 1][as.numeric(label=="III")==1], phi[, 2][as.numeric(label=="III")==1], col="#97b9ce", pch=17, cex=1.5)
points(phi[, 1][as.numeric(label=="IV")==1], phi[, 2][as.numeric(label=="IV")==1], col="#1f456e", pch=24, cex=1.5)
legend("bottomright", c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), col = c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), pch=c(7,10,17,24))
dev.off()


########################
#   Cluster Analysis   #
########################

# Hierarchical Clustering via Manhattan 
dat.dist <- dist(t(top.ranked.genes), method="manhattan") # calculate distance
dat.clust <- hclust(dat.dist, method="median") # calculate clusters
png("Dendrogram_TopGenes.png")
plot(dat.clust, labels=dat.info$label, xlab="Clustered Samples", ylab="Distance", main="Hierarchical Clustering Dendrogram\nRanked Classification for Glioma Grades")
dev.off()

# K-means clustering via PCA Analysis 
# Done though the kernlab library
# Extract out the first 5 component vectors and compute K-means with two centers 
dat.kpca <- kpca(t(top.ranked.genes), kernel="rbfdot", kpar=list(sigma=0.002), features=5)
pcv <- pcv(dat.kpca)
rot <- rotated(dat.kpca)
pcv.k2 <- kmeans(pcv, centers=2, iter.max=20)
rot.k2 <- kmeans(rot, centers=2, iter.max=20)
pcv.k3 <- kmeans(pcv, centers=3, iter.max=20)
rot.k3 <- kmeans(rot, centers=3, iter.max=20)

# 2D scatterplot of the first 5 eigenfeatures (from PCA)
png("Scatterplot_PCA.png")
par(mfrow=c(2, 2))
plot(pcv, col=pcv.k2$cluster, cex=1, main="PCA Scatter Plot with PC vectors \n(column-wise) nk=2", xlab="P1", ylab="P2")
points(pcv.k2$centers, col=pcv.k2$cluster, pch=1, cex=2.5)

plot(rot, col=rot.k2$cluster, cex=1, main="PCA Scatter Plot with PC \n Projected vectors k=2", xlab="P1", ylab="P2")
points(rot.k2$centers, col=rot.k2$cluster, pch=1, cex=2.5)

plot(pcv, col=pcv.k3$cluster, cex=1, main="PCA Scatter Plot with PC vectors \n (column-wise) k=3", xlab="P1", ylab="P2")
points(pcv.k3$centers, col=pcv.k3$cluster, pch=1, cex=2.5)

plot(rot, col=rot.k3$cluster, cex=1, main="PCA Scatter Plot with PC \n Projected vectors k=3", xlab="P1", ylab="P2")
points(rot.k3$centers, col=rot.k3$cluster, pch=1, cex=2.5)

dev.off()


###########################
#  Sample Classification  #
###########################  
# LDA Analysis
t.dat <- as.data.frame(t.dat)
# Train 15 control samples, 30 grade II glioma, 20 grade III glioma, 60 grade IV glioma
training <- rbind(t.dat[1:14,], t.dat[24:53,], t.dat[69:88,], t.dat[100:160,])
test <- t.dat[ !(rownames(t.dat) %in% rownames(training)), ]

te.names <- rownames(test)
test.names <- factor(gsub('_GSM[[:digit:]]+', '', te.names))
test$names <- NULL

tr.names <- rownames(training)
train.names <- factor(gsub('_GSM[[:digit:]]+', '', tr.names))
training$names <- NULL

ptm <- proc.time()
train.lda.2 <- lda(train.names~., data=training)
predicted.lda <- predict(train.lda.2, test)
proc.time() - ptm
# user  system elapsed 
# 31.71    0.63   32.80 
lda.cm <- table(predicted.lda$class, test.names)
# LDA Classification rate
lda.model <- lda.cm %>% prop.table() %>% round(3)

# Plot linear discriminant analysis
png("LDAPlot.png")
par(mfrow = c(2, 1))
plot(range(predicted.lda$x[, 1]), range(predicted.lda$x[, 2]), type="n", xlab="LD1", ylab="LD2", main="LDA Plot of all Genes\nDisciminant Scores on 2 Discriminant Variables")
points(predicted.lda$x[, 1], predicted.lda$x[, 2], col=c(rep("#c1a87d", 11), rep("#7ef9ff", 21), rep("#97b9ce", 3), rep("#1f456e", 20), pch=c(rep(7,11), rep(10,21), rep(17,3), rep(24,20))))
legend("topleft", cex=0.5, c("Non-Tumor", "Glioma II", "Glioma III", "Glioma IV"), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), pch=c(7,10,17,24))

x <- predicted.lda$x[, 1]
y <- predicted.lda$x[, 2]
z <- predicted.lda$x[, 3]
scatterplot3d(x, y, z, xlab="LD1", ylab="LD2", zlab="LD3", color="#97b9ce", pch=5, col.axis="#132842", col.grid="#7ef9ff", main="LDA Plot of all Genes\nDisciminant Scores on 3 Discriminant Variables")
dev.off()

# PCA analysis on the first 3 component vectors
pca       <- prcomp(t(top.ranked.genes))
dat.loadings <- pca$x[, 1:3] 

# Support vector machine
svp <- ksvm(dat.loadings, label, type="C-svc")
fit <- fitted(svp)

# error rates (incorrect classifications)
error.list <- NULL
label.NT <- rep("NT", length(label[label=="NT"]))
label.T <- rep("Glioma", length(label[!(label %in% label.NT)])) # Relabel all the Glioma samples as "glioma"
label.class <- c(label.NT, label.T)
er1 <- sum(fit[label.class=="NT"]=="Glioma") # number of incorrect normal classifications is 0
er2 <- sum(fit[label.class=="Glioma"]=="NT") # number of correct Glioma classifications is 4
er.total <- sum(er1,er2)/ncol(dat.loadings)	
er.total <- round(er.total*100,1) # 133.3 total misclassifications

# Error rate bassed on SVM fitting
svm.cm <- table(label, fit)

# SVM Classification rate
svm.model <- svm.cm %>% prop.table() %>% round(3)

# Scatter plot of PC1 vs PC2
nt <- dat.loadings[1:23,]
gII <- dat.loadings[24:68,]
gIII <- dat.loadings[69:99,]
gIV <- dat.loadings[100:180,]

png("Scatterplot_PCA.png")
par(mfrow = c(3,1))
plot(range(dat.loadings[, 1]), range(dat.loadings[, 2]), type="n", xlab="Principal Component 1", ylab="Principal Component 2", main="PCA Plot for GDS1962 data\n PC1 vs. PC2")
points(nt[, 1], nt[, 2], col="#c1a87d", pch=7)
points(gII[, 1], gII[, 2], col="#7ef9ff", pch=10)
points(gIII[, 1], gIII[, 2], col="#97b9ce", pch=17)
points(gIV[, 1], gIV[, 2], col="#1f456e", pch=24)
legend("topleft", c("Ctl", "II", "III", "IV"), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), pch=c(pch=c(7,10,17,24)))

# Scatter plot of PC1 vs PC3
plot(range(dat.loadings[, 1]), range(dat.loadings[, 3]), type="n", xlab="Principal Component 1", ylab="Principal Component 3", main="PCA Plot for GDS1962 data\n PC1 vs. PC3")
points(nt[, 1], nt[, 3], col="#c1a87d", pch=7)
points(gII[, 1], gII[, 3], col="#7ef9ff", pch=10)
points(gIII[, 1], gIII[, 3], col="#97b9ce", pch=17)
points(gIV[, 1], gIV[, 3], col="#1f456e", pch=24)
legend("bottomleft", c("Ctl", "II", "III", "IV"), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), pch=c(pch=c(7,10,17,24)))

# Scatter plot of PC2 vs PC3
plot(range(dat.loadings[, 2]), range(dat.loadings[, 3]), type="n", xlab="Principal Component 2", ylab="Principal Component 3", main="PCA Plot for GDS1962 data\n PC2 vs. PC3")
points(nt[, 2], nt[, 3], col="#c1a87d", pch=7)
points(gII[, 2], gII[, 3], col="#7ef9ff", pch=10)
points(gIII[, 2], gIII[, 3], col="#97b9ce", pch=17)
points(gIV[, 2], gIV[, 3], col="#1f456e", pch=24)
legend("bottomright", c("Ctl", "II", "III", "IV"), col=c("#c1a87d", "#7ef9ff", "#97b9ce", "#1f456e"), pch=c(pch=c(7,10,17,24)))
dev.off()

# pca.var
pca.var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
cat(sprintf("Note: Approximately %s variability is explained using only the first two eigenvalues.\n", 
            sum(pca.var[1:2])))

# Scree Plot illustrating Level of Variance
png("Screeplot.png")
plot(c(1:length(pca.var)), pca.var, type="b", xlab="Components", ylab="Percent Variance", bg="#132842", pch=21)
title("Scree Plot illustrating %-Variability Explained By Each Eigenvalue\n Non-Tumor/Glioma dataset")
dev.off()


##########################
# Top discriminant genes #
##########################
top.genes.info <- top.genes.info[rownames(top.genes.info) %in% rownames(ann),]
top.genes.info <- top.genes.info[order(top.genes.info$pvalue),]
top5    <- head(top.genes.info, n=5L, na.omit=T)
bottom5 <- tail(top.genes.info, n=5L, na.omit=T)
