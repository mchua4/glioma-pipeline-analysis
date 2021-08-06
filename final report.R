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
library(Biobase)
library(dplyr)
library(gplots)
library(MASS)
library(outliers)
library(kernlab)
library(limma)
library(gplots)
library(data.table)

# Set working directory 
setwd("C:/Users/User/Dropbox/My PC (DESKTOP-VL7FAF0)/Documents/R/final")

# Open GDS1962 dataset 
gds <- getGEO("GDS1962")
dat <- Table(gds)


# Corresponding annotation files
ann.file <- "GPL570.annot"
ann <- read.delim(ann.file, header=T, row.names=1, skip=27, sep="\t")
ann <- ann[1:nrow(dat), ]
gene.info <- data.frame(Description = ann$Gene.title, Symbol = ann$Gene.symbol)
rownames(gene.info) <- rownames(ann)

# Glioma grade data frame
astrocytomas <- cbind(dat[1:3], dat[,26:32], dat[,33:51])
glioblastomas <- cbind(dat[1:3], dat[,52:132])
oligodendrogliomas <- cbind(dat[1:3], dat[,133:170], dat[,171:182])

# Append tumor grades to corresponding samples
for(a in 26:32) {names(dat)[a] <- paste("II", names(dat[a]), sep="_")}
for(b in 33:51) {names(dat)[b] <- paste("III", names(dat[b]), sep="_")}
for(c in 52:132) {names(dat)[c] <- paste("IV", names(dat[c]), sep="_")}
for(d in 133:170) {names(dat)[d] <- paste("II", names(dat[d]), sep="_")}
for(e in 171:182) {names(dat)[e] <- paste("III", names(dat[e]), sep="_")}

# Append disease states to corresponding samples
for(i in 3:25) {names(dat)[i] <- paste("non-tumor", names(dat[i]), sep="_")}
for(j in 26:51) {names(dat)[j] <- paste("astrocytomas", names(dat[j]), sep="_")}
for(k in 52:132) {names(dat)[k] <- paste("glioblastomas", names(dat[k]), sep="_")}
for(m in 133:182) {names(dat)[m] <- paste("oligodendrogliomas", names(dat[m]), sep="_")}

# Convert dat to numeric matrix
samp.matrix <- data.matrix(dat[, (3:ncol(dat))])

# Obtain row and column info
rownames(samp.matrix) <- rownames(ann)
samp.count <- ncol(samp.matrix)
id.count  <- nrow(samp.matrix)

# Descriptive Satistics
stdev <- apply(samp.matrix, 1, sd, na.rm=T)
row.means <- rowMeans(samp.matrix, na.rm=T)

# Histograms illustrating initial spread
png("Histogram_data.rowMeans.png")
hist(row.means, col="#97B9CE", xlab="Mean", ylab="Frequency", main = paste0("Histogram of mean expression values for",id.count,"profiles"))

png("Histogram_data.stdev.png")
hist(stdev, col="#7EF9FF", xlab="Standard Deviation", ylab="Frequency", main=paste0("Histogram of standard deviation expression values for",id.count,"profiles"))
dev.off()


########################
#   Outlier Analysis   #
########################

# Correlation Matrix
dat.cor  <- cor(samp.matrix, method="pearson", use="pairwise.complete.obs")
cx <- c("#b58d3d", "#c1a87d", "#eee0b1", "#1f456e", "#132842", "#7ef9ff", "#97b9ce")
leg <- seq(min(dat.cor, na.rm=T), max(dat.cor, na.rm=T))

# Correlation plot - heatmap
png("Heatmap_dat.cor.png")	
par(oma=c(5,7,1,1))

image(dat.cor, main="Correlation between Glioma vs Non-Tumor Gene Expression")
axis(1, at=seq(0, 1, length=ncol(dat.cor)), label=dimnames(dat.cor)[[2]], cex.axis=0.9, las=2)
axis(2, at=seq(0, 1, length=ncol(dat.cor)), label=dimnames(dat.cor)[[2]], cex.axis=0.9, las=2)

image(as.matrix(leg), col=cx, axes=F)
tmp <- round(leg, 2)
axis(1, at=seq(0,1,length=length(leg)), labels=tmp, cex.axis=1)

dev.off()

# CV-mean plot
dat.mean <- apply(log2(samp.matrix), 2, mean) 
dat.var <- apply(log2(samp.matrix), 2, var)  
cv <- dat.var / dat.mean

png("Scatterplot_ColMeansCV.png")
plot(dat.mean, cv, xlab="Mean", ylab="CV", main="CV-Mean plot of Glioma vs Non-Tumor samples")
points(dat.mean, cv, col=c(rep("#1f456e", samp.count/2), rep("#eee0b1", samp.count/2)), pch=c(rep(23, samp.count/2), rep(49, samp.count/2)))
legend("topright", c("Non-Tumor", "Glioma"), pch=c(23, 49), col=c("#1f456e", "#eee0b1"))
text(dat.mean, cv, labels=names(dat.mean), cex=0.5, offset=10)

dev.off()

# Average correlation plot
dat.avg <- apply(dat.cor, 1, mean, na.rm=T)

png("Scatterplot_CorMeans_2groups.png")
plot(c(1,length(dat.avg)), range(dat.avg), type="n", xlab="", ylab="Average correlation", main="Avg correlation for Glioma vs Non-Tumor samples", axes=F)
points(dat.avg, col=c(rep("#1f456e", samp.count/2), rep("#eee0b1", samp.count/2)), pch=c(rep(23, samp.count/2), rep(49, samp.count/2)))

axis(1, at=c(1:length(dat.avg)), labels=colnames(samp.matrix), las=2, cex.lab=0.4, cex.axis=0.6)
axis(2)
grid(nx=16, col="grey")
legend("topright", c("Non-Tumor", "Glioma"), pch=c(23, 49), col=c("#1f456e", "#eee0b1"))

dev.off()

# png("Scatterplot_CorMeans_4groups.png")
# plot(
#   c(1,length(dat.avg)), 
#   range(dat.avg), 
#   type = "n", 
#   xlab = "",
#   ylab = "Average correlation",
#   main = "Avg correlation for Glioma vs Non-Tumor samples with SCF expression",
#   axes = FALSE
# )
# points(
#   dat.avg,
#   col = c(rep("#97b9ce", samp.count/4), rep("#1f456e", samp.count/4), rep("#c1a87d", samp.count/4), rep("#eee0b1", samp.count/4)),
#   pch = c(rep(1, samp.count/4), rep(23, samp.count/4), rep(24, samp.count/4), rep(49, samp.count/4))
# )
# axis(1, at=c(1:length(dat.avg)), labels = colnames(samp.matrix), las = 2, cex.lab = 0.4, cex.axis = 0.6)
# axis(2)
# grid(nx = 16, col = "grey")
# legend(
#   "bottomleft", 
#   c("Non-Tumor - 1", "Non-Tumor 2", "Astrocytomas - Tumor grade 2", "Astrocytomas - Tumor grade 3"), 
#   pch = c(1,23,24,49), col = c("#97b9ce", "#1f456e","#c1a87d", "#eee0b1"), bg = "white"
# )
# dev.off()


# Identifying Outlier(s) via outlier()
o <- dat.avg <=  outlier(dat.avg)
outlier <- dat.avg[o]
cat(sprintf("%s Outlier(s) identified!\n", length(outlier)))
outlier

# Remove Outlier(s) 
samp.matrix <- samp.matrix[, -(grep(names(outlier), colnames(samp.matrix)))]


########################
#     Filter Genes     #
########################

quantile(log2(rowMeans(samp.matrix)), na.rm=T)
#         0%        25%        50%        75%       100% 
#   0.9483085  6.2915759  7.7362583  9.2855611 15.0513555 

# Eliminating probes with rowMeans less than 0 on a log2 scale
dat.fil <- subset(samp.matrix, log2(rowMeans(samp.matrix)) > 0)
removed <- nrow(samp.matrix) - nrow(dat.fil)
cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed))

# Eliminate probes with rowMeans less than 3 on a log2 scale
dat.filtered <- subset(dat.fil, rowMeans(dat.fil) > 3)
dat.filtered <- as.data.frame(dat.filtered)
removed  <- nrow(dat.fil) - nrow(dat.filtered)
cat(sprintf("%s probes removed with rowMeans < 0 on a log2 scale\n", removed))

fil.mean <- apply(dat.filtered, 1, mean) 
fil.var  <- apply(dat.filtered, 1, var)  
fil.cv   <- fil.var / fil.mean

png("Scatterplot_RowMeansCV_Filtered.png")
plot(fil.mean, fil.cv, xlim=c(0,12), ylim=c(0,50), xlab="Means", ylab="CV", main="Filtered CV-Mean plot for Glioma vs Non-Tumor samples", col=c(rep("#1f456e", samp.count/2), rep("#eee0b1", samp.count/2)), pch=c(rep(23, samp.count/2), rep(49, samp.count/2)))
legend("topright", c("Non-Tumor", "Glioma"), pch = c(23, 49), col = c("#1f456e", "#eee0b1"))
abline(v=3, col=2, lwd=2)
dev.off()

# Update gene.info data.frame
gene.info <- subset(gene.info, rownames(gene.info) %in% rownames(dat.filtered))

########################
#  Feature Selection   #
########################   

# Any columns that contain KIT gene info
#KITLG <- gene.info[gene.info %like% "KIT",]
#KIT.dat <- subset(dat.filtered, dat$IDENTIFIER=="KIT" | dat$IDENTIFIER=="KITLG")

# Groups
is.ctl <- grepl("^non-tumor", names(dat.filtered))
is.glioma <- grepl("^(?!non-tumor)", names(dat.filtered), perl=T)


# Define types
type <- lapply(colnames(dat.filtered), 
               function(x) {
                 if(regexpr("Ctl", x) < 1) {"Glioma"} else {"Ctl"}
               })

# Student's t-test function
t.test.all.genes <- function(x,s1,s2) { 
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T) 
  out <- as.numeric(t.out$p.value)
  return(out)
}


# Run the Student's t-test function of KIT gene expression between control and glioma samples
rawp <- apply(dat.filtered, 1, t.test.all.genes, s1=is.ctl, s2=is.glioma)

# Bonferroni correction
Bonferroni <- 0.05/length(rawp)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni))
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(rawp < Bonferroni))) # 5953 <-- update this number

# Holm-Bonferroni correction
p.holm <- p.adjust(rawp, method="holm")
p.raw  <- sort(rawp)
p.holm <- sort(p.holm)
p1     <- as.matrix(p.raw)
p2     <- as.matrix(p.holm)
all.p  <- cbind(p1, p2)
colnames(all.p) <- c("Raw P-value", "Adjusted P-Value")

# Output
cat(sprintf("Total number of genes with p-value < 0.05 is %s\n", sum(p.holm < 0.05))) # 6111 <-- update this number
cat(sprintf("Total number of genes with p-value < 0.01 is %s\n", sum(p.holm < 0.01))) # 4878 <-- update this number

Bonferroni.2 <- 0.05/length(p.holm)
cat(sprintf("Bonferroni correction is %s\n", Bonferroni.2))
cat(sprintf("Total number of genes with p-value < the Bonferroni correction is %s\n", sum(rawp < Bonferroni.2))) # 5953 <-- update this number

# Plot the sorted raw P-values
png("Raw_rawpvaluePlot.png")
plot(p.raw, type="b", main="Raw P-values for All Genes\n", xlab="Genes", ylab="P-values", pch=1, col="#7ef9ff")
dev.off()

# Plot Adjusted vs. Raw P-values
png("Raw_PvaluePlot.png")
matplot(all.p, type="b", pch=1, col=c("#1f456e", "#eee0b1"), xlab="Genes", ylab="P-values", main = "Adjusted Vs. Raw P-values for All Genes\n")
legend("bottomright", legend=colnames(all.p), pch=1, col=c("#1f456e", "#eee0b1"))
dev.off()

# Update gene.info data.frame with P-value information
gene.info$pvalue <- rawp

# Hypothesis Pass/Fail dataframe
thresh  <- 0.05
rnames  <- rownames(dat.filtered)
p.test  <- lapply(as.list(rawp), function(f){ if (f < thresh) TRUE else FALSE })
rawp.df <- as.data.frame(do.call(rbind, p.test), rname = rnames)
names(rawp.df) <- c("rawp.Test")
rownames(rawp.df) <- rownames(dat.filtered) 

# Update gene.info data.frame with P-value test results (pass/fail)
gene.info$rawp.Test <- rawp.df$rawp.Test

cat(sprintf("Threshold: %s\n", thresh))
table(rawp.df$rawp.Test)
#   Number of True        : 28437   
#   Number of False       : 26169 
nrow(rawp.df)
#   Total                 : 54606 


# Histogram illustrating the distribution of P-values 
png("Histogram_rawp.png")
hist(rawp, col="#1f456e", xlab="P-Value", ylab="Frequency", main="Histogram of T-Test P-Values values for GDS1962 profiles")
hist(p.holm, col="#7ef9ff", add=T)
legend("top", legend = colnames(all.p), pch = 16, col = c("#1f456e", "#7ef9ff"))
dev.off()

# Determining fold change between glioma and non-tumor samples in log2
dat.log <- log2(dat.filtered)

ctl.mean <- apply(dat.log[3:25], 1, mean, na.rm=T)
glioma.mean <- apply(dat.log[26:51], 1, mean, na.rm=T)

fold <- ctl.mean - glioma.mean
linear <- 2^(fold)
max.fold <- max(linear) # maximum fold change is 180.2804 <-- update this number
min.fold <- min(linear) # minimum fold change is 0.0642937 <-- update this number


# Update gene.info data.frame with fold information
gene.info$fold <- fold	

# Via Lab_05, determine probesets that demonstrate a 2x fold change
fold.test <- lapply(as.list(fold), function(x){ if (abs(x) > log2(2)) TRUE else FALSE })
fold.df   <- as.data.frame(do.call(rbind, fold.test))
names(fold.df) <- c("Fold.Test")

# Update gene.info data.frame with fold information
gene.info$Fold.Test <- fold.df$Fold.Test

table(fold.df$Fold.Test)

# Number of True        : 4322   
# Number of False       : 50284 
(nrow(fold.df))
# Total                 : 54606 

# Histogram illustrating the distribution of log2(fold change)
png("Histogram_fold.png")
hist(fold, col="#c1a87d", xlab="Log2 Fold Change", ylab="Frequency", main=paste("Histogram of Fold Change values for GDS1962 profiles"))
abline(v=log2(2), col=2, lwd=2)
abline(v=-log2(2), col=2, lwd=2)
dev.off()

# Overall Volcano plot showing cutoffs and differences within the dataset
p.transformed <- (-1 * log10(rawp))
png("VolcanoPlot.png")
plot(range(p.transformed), range(fold), type="n", xlab="-1 * log10(P-Value)", ylab="Fold Change", main="Volcano Plot Illustrating  Glioma and Non-Tumor Differences")
points(p.transformed, fold, col=1, bg=1, pch=21)
points(p.transformed[(p.transformed > -log10(0.05) & fold > log2(2))], fold[(p.transformed > -log10(0.05) & fold > log2(2))], col=1, bg=2, pch=21)
points(p.transformed[(p.transformed > -log10(0.05) & fold < -log2(2))], fold[(p.transformed > -log10(0.05) & fold < -log2(2))], col=1, bg=3, pch=21)
abline(v=-log10(0.05))
abline(h=-log2(2))
abline(h=log2(2))
dev.off()

# New data.frame with genes that have passed both (Fold and rawp) tests
true.genes <- subset(gene.info, (Fold.Test & rawp.Test) == TRUE)
cat(sprintf("Total number of genes that pass both (rawp and Fold) tests: %s\n", nrow(true.genes))) # 9230 <-- update this number

# Write these genes with their corresponding values to an output .txt file
write.table(true.genes, file="TrueGenes.csv", sep=",", col.names=NA, qmethod="double")

# Ordering the highest genes (by P-value) in the form of a data.frame
# Note: dat.filtered is still in log2 scale
best.genes    <- order(rawp)[1:length(rawp)]
best.genes.df <- data.frame(index=best.genes, exp=2^dat.filtered[best.genes, ], rawp=rawp[best.genes])

# Expression matrix with the 'best.genes' in the original scale (based on P-value ranking)
top.genes.matrix <- 2^dat.filtered[best.genes, ]

# Feature Selection via svmRFE which utilizes the library e1071
t.dat  <- t(top.genes.matrix)
label  <- as.vector(unlist(type))
svm.df1 <- data.frame(type, t.dat)
dput(svm.df1[1:2,1:5])
svm.df2 <- gsub(",", "", svm.df1)
svm.df <- as.numeric(unlist(svm.df2))
dput(svm.df[1:100])


# svmRFE will automatically normalize/scale the dataframe
ptm <- proc.time()

svmRFE <- function(X, k=1, halve.above=5000) {
  
  X <- as.numeric(X)
  
  # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
  n = ncol(X) - 1
  
  # Scale data up front so it doesn't have to be redone each pass
  cat('Scaling data...')
  X[, -1] = scale(X[, -1])
  cat('Done!\n')
  flush.console()
  
  pb = txtProgressBar(1, n, 1, style=3)
  
  i.surviving = 1:n
  i.ranked    = n
  ranked.list = vector(length=n)
  
  # Recurse through all the features
  while(length(i.surviving) > 0) {
    if(k > 1) {
      # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
      folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
      folds = lapply(1:k, function(x) which(folds == x))
      
      # Obtain weights for each training set
      w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
      w = do.call(rbind, w)
      
      # Normalize each weights vector
      w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
      
      # Compute ranking criteria
      v    = w * w
      vbar = apply(v, 2, mean)
      vsd  = apply(v, 2, sd)
      c    = vbar / vsd
    } else {
      # Only do 1 pass (i.e. regular SVM-RFE)
      w = getWeights(NULL, X[, c(1, 1+i.surviving)])
      c = w * w
    }
    
    # Rank the features
    ranking = sort(c, index.return=T)$ix
    if(length(i.surviving) == 1) {
      ranking = 1
    }
    
    if(length(i.surviving) > halve.above) {
      # Cut features in half until less than halve.above
      nfeat = length(i.surviving)
      ncut  = round(nfeat / 2)
      n     = nfeat - ncut
      
      cat('Features halved from', nfeat, 'to', n, '\n')
      flush.console()
      
      pb = txtProgressBar(1, n, 1, style=3)
      
    } else ncut = 1
    
    # Update feature list
    ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
    i.ranked    = i.ranked - ncut
    i.surviving = i.surviving[-ranking[1:ncut]]
    
    setTxtProgressBar(pb, n-length(i.surviving))
    flush.console()
  }
  
  close(pb)
  
  return (ranked.list)
}


ranked.list <- svmRFE(svm.df, k=10, halve.above=1000)


proc.time() - ptm
# Run-time:
# user  system elapsed 
# 0.01    0.00    2.32

# Write the rankings to an output .txt file so that it can be read in later if needed
output <- data.frame(RankedOrder = ranked.list)
write.table(output, file = "RankedList.txt")

top.ranked.genes           <- top.genes.matrix[ranked.list, ]
rownames(top.ranked.genes) <- rownames(top.genes.matrix[ranked.list, ])

# Create a new genes.info data.frame for the ranked genes
top.genes.info             <- gene.info[rownames(top.ranked.genes), ]

tg <- top.genes.info$pvalue[top.genes.info$pvalue < thresh]

png("Histogram_ranked_rawp.png")
hist(tg, col="#7ef9ff", xlab="P-Value", ylab="Frequency", main="Histogram of T-Test P-Values values for ranked profiles")
abline(v=thresh, col=2, lwd=2)
dev.off()

top5    <- head(top.genes.info, n=5L)
bottom5 <- tail(top.genes.info, n=5L)


###########################
#  Sample Classification  #
###########################  

# PCA analysis on the first 3 component vectors
pca       <- prcomp(t(top.ranked.genes))
dat.loadings <- pca$x[, 1:3] 

# As per Lecture 10:
svp <- ksvm(dat.loadings,label,type="C-svc")
svp

# Get fitted values
fit <- fitted(svp)

# error rates (incorrect classifications)
er1 <- sum(fit[label == "Non-Tumor"] == "Glioma")	# 0 
er2 <- sum(fit[label == "Glioma"] == "Non-Tumor") # 0  


# Plot code adapted from Lab_07 and Lab_08
# Plotting: Component 1 vs. Component 2
#			Component 1 vs. Component 3
#			Component 2 vs. Component 3

png("Scatterplot_PCA_1v2.png")
plot(range(dat.loadings[, 1]), range(dat.loadings[, 2]), type="n", xlab="Principal Component 1", ylab="Principal Component 2", main="PCA Plot for GDS1962 data\n PC1 vs. PC2")
points(dat.loadings[, 1][as.numeric(type=="Non-Tumor")==1], dat.loadings[, 2][as.numeric(type=="Non-Tumor")==1], col="#132842", pch=15)
points(dat.loadings[, 1][as.numeric(type=="Glioma")==1], dat.loadings[, 2][as.numeric(type=="Glioma")==1], col="#c1a87d", pch=19)
legend("bottomleft", c("Non-Tumor", "Glioma"), col=c("#132842", "#c1a87d"), pch=c(15,19))
dev.off()

png("Scatterplot_PCA_1v3.png")
plot(range(dat.loadings[, 1]), range(dat.loadings[, 3]), type="n", xlab="Principal Component 1", ylab="Principal Component 3", main="PCA Plot for GDS1962 data\n PC1 vs. PC3")
points(dat.loadings[, 1][as.numeric(type=="Non-Tumor")==1], dat.loadings[, 3][as.numeric(type=="Non-Tumor")==1], col="#132842", pch=15)
points(dat.loadings[, 1][as.numeric(type=="Glioma")==1], dat.loadings[, 3][as.numeric(type=="Glioma")==1], col="#c1a87d", pch=19)
legend("bottomright", c("Non-Tumor", "Glioma"), col=c("#132842", "#c1a87d"), pch=c(15,19))
dev.off()

png("Scatterplot_PCA_2v3.png")
plot(range(dat.loadings[, 2]), range(dat.loadings[, 3]), type="n", xlab="Principal Component 2", ylab="Principal Component 3", main="PCA Plot for GDS1962 data\n PC2 vs. PC3")
points(dat.loadings[, 2][as.numeric(type=="Non-Tumor")==1], dat.loadings[, 3][as.numeric(type=="Non-Tumor")==1], col="#132842", pch=15)
points(dat.loadings[, 2][as.numeric(type=="Glioma")==1], dat.loadings[, 3][as.numeric(type=="Glioma")==1], col="#c1a87d", pch=19)
legend("bottomleft", c("Non-Tumor", "Glioma"), col=c("#132842", "#c1a87d"), pch=c(15,19))
dev.off()

pca.var <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
# pca.var
cat(sprintf("Note: Approximately %s variability is explained using the only the first two eigenvalues.\n", 
            sum(pca.var[1:2])))

# Scree Plot adapted from Lab_07 Illustrating Level of Variance
png("Screeplot.png")
plot(c(1:length(pca.var)), pca.var, type="b", xlab="Components", ylab="Percent Variance", bg="#132842", pch=21)
title("Scree Plot illustrating %-Variability Explained By Each Eigenvalue\n Non-Tumor/Glioma dataset")
dev.off()

# MDS Analysis via Kruskal's Non-metric Approach
dat.dist <- dist(t(top.ranked.genes))
dat.mds  <- isoMDS(dat.dist)

png("MDSplot_Kruskal.png")
plot(dat.mds$points, type = "n")
points(dat.mds$points[, 1][as.numeric(type=="Non-Tumor")==1], dat.mds$points[, 2][as.numeric(type=="Non-Tumor")==1], col="#c1a87d", pch=16, cex=1.5)
points(dat.mds$points[, 1][as.numeric(type=="Glioma")==1], dat.mds$points[, 2][as.numeric(type=="Glioma")==1], col="#132842", pch=16, cex=1.5)
title(main = "Kruskal's Non-metric MDS plot for GDS1962 dataset\nNon-Tumor vs. Glioma")
legend("bottomleft", c("Non-Tumor", "Glioma"), col = c("#c1a87d", "#132842"), fill = c("#c1a87d", "#132842"))
dev.off()

# MDS Analysis via the Classical Metric Approach
# Note: No stress value is provided
dat.loc  <- cmdscale(dat.dist) 

png("MDSplot_Classical.png")
plot(dat.loc, type = "n")
points(dat.loc[, 1][as.numeric(type=="Non-Tumor")==1], dat.loc[, 2][as.numeric(type=="Non-Tumor")==1], col="#c1a87d", pch=16, cex=1.5)
points(dat.loc[, 1][as.numeric(type=="Glioma")==1], dat.loc[, 2][as.numeric(type=="Glioma")==1], col="#132842", pch=16, cex=1.5)
title(main="Classical Metric MDS plot for GDS1962 dataset\nNon-Tumor vs. Glioma")
legend("bottomleft", c("Non-Tumor", "Glioma"), col=c("#c1a87d", "#132842"), fill=c("#c1a87d", "#132842"))
dev.off()

# Weighted Laplacian Graph
# Note: The function GeneratePhi was taken from Lab_07 and Lecture_08
GeneratePhi <- function (X, qnt = NULL) {
  dist2full <- function(dis) {
    n     <- attr(dis, "Size")
    full  <- matrix(0, n, n)
    full[lower.tri(full)] <- dis
    full + t(full)
  }
  dat.dis <- dist(t(X),"euc")^2
  if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
  if(is.null(qnt))  {eps <- min(dat.dis[dat.dis != 0])}
  kernal         <- exp(-1 * dat.dis / (eps))
  K1             <- dist2full(kernal)
  diag(K1)       <- 0
  D              <- matrix(0, ncol = ncol(K1), nrow = ncol(K1))
  tmpe           <- apply(K1, 1, sum)
  tmpe[tmpe > 0] <- 1/sqrt(tmpe[tmpe > 0])
  tmpe[tmpe < 0] <- 0
  diag(D)        <- tmpe
  L              <- D%*% K1 %*% D
  X              <- svd(L)$u
  Y              <- X / sqrt(apply(X^2, 1, sum))
}

temp <- t(top.ranked.genes)
temp <- scale(temp, center=T, scale=T)
phi <- GeneratePhi(t(temp), qnt=NULL)

png("LaplacianPlot.png")
plot(range(phi[, 1]), range(phi[, 2]), xlab="Phi 1", ylab="Phi 2", main="Weighted Graph Laplacian Plot for GDS1962 dataset\nNon-Tumor vs. Glioma")
points(phi[, 1][as.numeric(type=="Non-Tumor")==1], phi[, 2][as.numeric(type=="Non-Tumor")==1], col="#c1a87d", pch=16, cex=1.5)
points(phi[, 1][as.numeric(type=="Glioma")==1], phi[, 2][as.numeric(type=="Glioma")==1], col="#132842", pch=16, cex=1.5)
legend("top", c("Non-Tumor", "Glioma"), col=c("#c1a87d", "#132842"), fill=c("#c1a87d", "#132842"))
dev.off()





########################
#   Cluster Analysis   #
########################

# Hierarchical Clustering via Euclidean
top.ranked.genes <- t(top.ranked.genes) # transpose top.ranked.genes
dat.dist <- dist(top.ranked.genes, method="euclidean") # calculate distance
dat.clust <- hclust(dat.dist, method="centroid") # calculate clusters

png("Dendrogram_TopGenes.png")
plot(dat.clust, labels=colnames(top.ranked.genes), xlab="Clustered Samples", ylab="Distance", main="Hierarchical Clustering Dendrogram\nRanked Angiogenesis Classification")
dev.off()

# Heatmap of top.ranked.genes
# Note: Due to memory allocation, the top 50 ranked genes were plotted
png("Heatmap_Top50Genes.png")
heatmap(top.ranked.genes[1:50, ], col=color, xlab="Samples", ylab="Top Ranked Genes", main="Heatmap for the top 50 genes")
dev.off()

# Heatmap of top.ranked.genes
# Note: Due to memory allocation, 50 randomly ranked genes were plotted 
png("Heatmap_TopRandom50Genes.png")
heatmap(top.ranked.genes[sample(top.ranked.genes, 50), ], col=color, xlab="Samples", ylab="Randomly Ranked Genes", main="Heatmap for 50 randomly ranked genes")
dev.off()

# K-means clustering via PCA Analysis 
# Done though the kernlab library
# Extract out the first 5 component vectors and compute K-means with two centers 
dat.kpca <- kpca(t(top.ranked.genes), kernel="rbfdot", kpar=list(sigma=0.002), features=5)
pcv      <- pcv(dat.kpca)
rot      <- rotated(dat.kpca)
pcv.k    <- kmeans(pcv, centers=2, iter.max=20)
rot.k    <- kmeans(rot, centers=2, iter.max=20)

# 2D scatterplot of the first 5 eigenfeatures (from PCA)
png("Scatterplot_PCA.png")
par(mfrow=c(2, 1))
plot(pcv, col=pcv.k$cluster, cex=1, main="PCA Scatter Plot with PC vectors (column-wise)\nk=2", xlab="P1", ylab="P2")
points(pcv.k$centers, col=1:2, pch="*", cex=2.5)

plot(rot, col=rot.k$cluster, cex=1, main="PCA Scatter Plot with PC Projected vectors\nk = 2", xlab="P1", ylab="P2")
points(rot.k$centers, col=1:2, pch="*", cex=2.5)
dev.off()

# LDA Analysis
# Note: Code is adapted from Lab_09

training <- as.data.frame(rbind(t.dat[c(1, 2, 5, 6, 9, 10, 13, 14), ]))
test     <- as.data.frame(t.dat[ !(rownames(t.dat) %in% rownames(training)), ])

te.names <- rownames(test)
te.names[c(1, 2, 5, 6)] <- paste("66", te.names[c(1, 2, 5, 6)], sep="_")
te.names[c(3, 4, 7, 8)] <- paste("6", te.names[c(3, 4, 7, 8)], sep="_")
test.names <- factor(gsub('_GSM[[:digit:]]+', '', te.names))

tr.names <- rownames(training)
tr.names[c(1, 2, 5, 6)] <- paste("66", tr.names[c(1, 2, 5, 6)], sep="_")
tr.names[c(3, 4, 7, 8)] <- paste("6", tr.names[c(3, 4, 7, 8)], sep="_")
train.names <- factor(gsub('_GSM[[:digit:]]+', '', tr.names))

# Due to memory allocations LDA was run on the first 5000 genes
ptm <- proc.time()
train.lda.2      <- lda(train.names ~ ., data=training[, c(1:5000)])
train.pred.2.out <- predict(train.lda.2, test[, c(1:5000)])
proc.time() - ptm
table(train.pred.2.out$class, test.names)

png("LDAPlot_2D.png")
plot(range(train.pred.2.out$x[, 1]), range(train.pred.2.out$x[, 2]), type="n", xlab="LD1", ylab="LD2", main="LDA Plot of 5000 Genes\nDisciminant Scores on 2 Discriminant Variables",)
points(train.pred.2.out$x[, 1], train.pred.2.out$x[, 2], col=c(rep("#97b9ce", 2), rep("#1f456e", 2), rep("#c1a87d", 2), rep("#eee0b1", 2)), pch=c(rep(1, 23), rep(24, 49)))
legend("bottom", c("Non-Tumor - 1", "Non-Tumor 2", "Astrocytomas - Tumor grade 2", "Astrocytomas - Tumor grade 3"), col=c("#97b9ce", "#1f456e","#c1a87d", "#eee0b1"), pch=c(1,23,24,49))
dev.off()


col = c(rep("#97b9ce", 2), rep("#1f456e", 2), rep("#c1a87d", 2), rep("#eee0b1", 2))
png("LDAPlot_3D.png")
scatterplot3d(train.pred.2.out$x[, 1], train.pred.2.out$x[, 2], train.pred.2.out$x[, 3], xlab="LD1", ylab="LD2", zlab="LD3", col, pch=c(rep(16, 4), rep(17, 4)), col.axis="#132842", col.grid="#7ef9ff", main="LDA Plot of 5000 Genes\nDisciminant Scores on 3 Discriminant Variables")
dev.off()
