####### Main workflow in the paper
# If you only want to run part of the workflow, you can load the required data from the previous steps as indicated in each section below (load(...))
startTime <- Sys.time()
#######################################
###      Load required packages     ###
#######################################
library(ibb)
library(stringr)
library(BiocManager)
library(BiocParallel)
library(BioNet)
library(pheatmap)
library(countdata)
library(DESeq2)
library("vsn")
library('org.Hs.eg.db')
library(foreach)
library(doParallel)
library(clusterProfiler)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)
library(europepmc)
library(WGCNA)
options(stringsAsFactors=FALSE)
# Load required functions from the src directory
source('src/preprocessing.R')
source('src/differential_expression.R')
#source('src/add_uniprot_info.R')
source('src/clustering.R')
source('src/coexpression_analysis.R')
source('src/validation.R')


#######################################
###             Setup               ###
#######################################
### Set working directory
setwd('C:/Users/marle/Desktop/Y2/Internship/Project') #Replace with your working directory


### Create output directories
dir.create('output')
dir.create('figures')
dir.create('rdata')

### Load data
time1 <- Sys.time()
d <- read.table('./data/TCGA_rna_count_data.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)
group.data <- read.table('./data/non_silent_mutation_profile_crc.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)
time2 <- Sys.time()
print("Loading data:")
difftime(time2, time1, units="secs")

time1 <- Sys.time()
d <- d[,1:ceiling(ncol(d)/2)]
d <- d[,colnames(d) %in% rownames(group.data)]
group.data <- group.data[rownames(group.data) %in% colnames(d),]
new_order <- sort(colnames(d))
d <- d[, new_order]
new_order_rows <- sort(rownames(d))
d <- d[new_order_rows,]
d.raw <- d[,1:ncol(d)]
d.raw <- apply(d.raw, c(1,2), as.numeric)
d.raw <- d.raw[rowSums(d.raw) >= 10,]
samples <- colnames(d.raw)
groups <- group.data$KRAS
# Save input data
save(d.raw, group.data, groups, file='rdata/input_data.RData')
time2 <- Sys.time()
print("Processing data:")
difftime(time2, time1, units="secs")

group.KRAS <- as.data.frame(group.data[,"KRAS"])
rownames(group.KRAS) <- rownames(group.data)
colnames(group.KRAS) <- c("KRAS")

###PLOTS
png(file=paste0('figures/heatmap.png'), width=1920, height=1020)
pheatmap(d.raw, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, show_colnames = FALSE, annotation_col=group.KRAS, cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

colnames(d.raw) <- groups

png(file=paste0('figures/sampleCluster.png'), width=1920, height=1020)
plotClusterTreeSamples(t(d.raw), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

sampleDists <- dist(t(d.raw))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(d.raw)
colnames(sampleDistMatrix) <- colnames(d.raw)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(file=paste0('figures/sampleDistances.png'), width=1920, height=1020)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
dev.off()

par(mar = c(2, 1, 1, 1))
png(file=paste0('figures/meanExpression.png'), width=1920, height=1020)
barplot(apply(t(d.raw),1,mean, na.rm=T),
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

par(mar = c(1, 1, 1, 1))
png(file=paste0('figures/variance.png'), width=1920, height=1020)
meanSdPlot(d.raw, cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

###VST
d.adj <- varianceStabilizingTransformation(d.raw, blind = F)
#d.adj <- rlog(d.raw, blind=FALSE)

colnames(d.raw) <- samples

#ddds <- DESeqDataSetFromMatrix(countData = d.raw, colData = group.KRAS, design = ~ KRAS)
#ntd <- normTransform(ddds)
#d.adj <- assay(ddds)

###PLOTS AFTER VST
#dds <- DESeqDataSetFromMatrix(countData = d.raw, colData = group.KRAS, design = ~ KRAS)
#dds
#dds2 <- estimateSizeFactors(dds)
#d.norm <- counts(dds2, normalized = TRUE)

colnames(d.raw) <- groups

png(file=paste0('figures/clusterSamplesVST.png'), width=1920, height=1020)
sampleTree <- hclust(dist(t(d.adj)), method = "ave")
#png(file=paste0('figures/clusterSamplesNorm.png'), width=1920, height=1020)
#sampleTree <- hclust(dist(t(d.norm)), method = "ave")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 275, col = "red");
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 275, minSize = 10)
#clust <- cutreeStatic(sampleTree, cutHeight = 2250000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
dev.off()

colnames(d.raw) <- samples
d.raw <- d.raw[,keepSamples]
samples <- colnames(d.raw)
ids <- rownames(d.raw)
group.data <- group.data[keepSamples,]
groups <- group.data$KRAS
d.adj <- d.adj[,keepSamples]
colnames(d.adj) <- samples
#d.norm <- d.norm[,keepSamples]
#colnames(d.norm) <- samples

group.KRAS <- as.matrix(group.data[,"KRAS"])
rownames(group.KRAS) <- rownames(group.data)
colnames(group.KRAS) <- c("KRAS")

png(file=paste0('figures/heatmapNorm.png'), width=1920, height=1020)
#pheatmap(d.norm, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group.KRAS))
#dev.off()
#
#colnames(d.norm) <- groups
#
#sampleDists <- dist(t(d.norm))
#sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- colnames(d.norm)
#colnames(sampleDistMatrix) <- colnames(d.norm)
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#png(file=paste0('figures/sampleDistancesNorm.png'), width=1920, height=1020)
#pheatmap(sampleDistMatrix,
#         clustering_distance_rows=sampleDists,
#         clustering_distance_cols=sampleDists,
#         col=colors, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
#dev.off()
#
#par(mar = c(2, 1, 1, 1))
#png(file=paste0('figures/meanExpressionNorm.png'), width=1920, height=1020)
#barplot(apply(t(d.norm),1,mean, na.rm=T),
#        xlab = "Sample", ylab = "Mean expression",
#        main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
#dev.off()

png(file=paste0('figures/sampleClusterVST.png'), width=1920, height=1020)
plotClusterTreeSamples(t(d.adj), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

png(file=paste0('figures/heatmapVST.png'), width=1920, height=1020)
pheatmap(d.adj, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group.KRAS))
dev.off()

group.KRAS[,1] <- as.numeric(factor(group.KRAS[,1]))
colnames(d.adj) <- groups

par(mar = c(2, 1, 1, 1))
png(file=paste0('figures/meanExpressionVST.png'), width=1920, height=1020)
barplot(apply(t(d.adj),1,mean, na.rm=T),
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

par(mar = c(1, 1, 1, 1))
png(file=paste0('figures/varianceVST.png'), width=1920, height=1020)
meanSdPlot(d.adj, cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

sampleDists <- dist(t(d.adj))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(d.adj)
colnames(sampleDistMatrix) <- colnames(d.adj)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(file=paste0('figures/sampleDistancesVST.png'), width=1920, height=1020)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
dev.off()


#######################################
###         Preprocessing           ###
#######################################
#load(file='rdata/input_data.rdata')
dir.create('figures/differential expression')
#d.norm <- normalize.sample(d.raw)
#d.cs <- normalize.cs(d.norm)
#length <- length(groups)
## unique colnames in d.cs are necessary for clustering
#colnames(d.cs) <- paste(groups, 1:length, sep='_')
## save input data into a file
#save(d.norm, d.cs, ids, groups, file='rdata/normalized_data.rdata')

#DESeq2
time1 <- Sys.time()
dds <- DESeqDataSetFromMatrix(countData = d.raw, colData = group.KRAS, design = ~ KRAS)
dds
dds <- DESeq(dds)
res <- results(dds)
#cooks <- assays(dds)[["cooks"]]
#png(file=paste0('figures/differential expression/cooksCheck.png'), width=1920, height=1020)
#plot(x = rownames(cooks), y = cooks, type = "h")
#dev.off()
results.names <- resultsNames(dds)
res
results.names
summary(res)

colnames(d.raw) <- groups

d.summary <- data.frame(unlist(res$log2FoldChange), unlist(res$padj), row.names = rownames(res))
colnames(d.summary) <- c("log2FC", "Pvalue")
d.summary <- na.omit(d.summary)
d.adj <- d.adj[rownames(d.adj) %in% rownames(d.summary),]
ids <- rownames(d.adj)
#d.norm <- d.norm[rownames(d.norm) %in% rownames(d.summary),]
#ids <- rownames(d.norm)
time2<- Sys.time()
print("DESeq2:")
difftime(time2, time1, units="secs")

###DESEQ PLOTS
par(mar = c(2, 1, 1, 1))
png(file=paste0('figures/differential expression/MA.png'), width=1920, height=1020)
plotMA(res, ylim=c(-2,2), cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

png(file=paste0('figures/differential expression/counts.png'), width=1920, height=1020)
plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup="KRAS", cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

par(mar = c(1, 1, 1, 1))
png(file=paste0('figures/differential expression/dispersion.png'), width=1920, height=1020)
plotDispEsts(dds, cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

par(mar = c(2, 1, 1, 1))
png(file=paste0('figures/differential expression/rejections.png'), width=1920, height=1020)
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter", cex.lab = 2, cex.axis = 2, cex.main = 2)
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()

par(mar = c(1, 1, 1, 1))
png(file=paste0('figures/differential expression/outliers.png'), width=1920, height=1020)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, names = groups)
dev.off()

par(mar = c(2, 1, 1, 1))
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
png(file=paste0('figures/differential expression/Wald.png'), width=1920, height=1020)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3), cex.lab = 2, cex.axis = 2, cex.main = 2)
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))
dev.off()

use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
png(file=paste0('figures/differential expression/pass.png'), width=1920, height=1020)
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency", cex.lab = 2, cex.axis = 2, cex.main = 2)
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
dev.off()

#######################################
### Differential expression analysis###
#######################################
## Test significance of differential expression between the groups
#load(file='rdata/normalized_data.RData')
# dir.create('output/differential_expression')
# print(Sys.time())
# # Beta binomial test
# bb <- betaBinomial(d.norm, ids, groups, 'WT', 'SNV', 'two.sided')
# print(Sys.time())
# d.sign <- data.frame(log2FC=bb$table[ids,]$Log2ratio, Pvalue=bb$table[ids,]$Pvalue)
# rownames(d.sign) <- ids
# # Merge stats about P-value thresholds and number of significant proteins
# sign.table <- cbind(bb$FDRs)
# d.summary <- d.sign
# # Write output into files
# write.table(bb$table, file='output/differential_expression/sign_test.tsv', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
# write.table(sign.table, file='output/differential_expression/FDR_pval_summary.tsv', sep='\t', quote=FALSE)
# write.table(d.summary, file='output/differential_expression/log2_pval_uniprot.tsv', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
# save(bb, sign.table, d.summary, file='rdata/differential_expression.RData')

#######################################
###          Clustering             ###
#######################################
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
# dir.create('figures/clustering')
# # Set colors for the different significance levels
# ann_colors <- list(group=c(SNV='#F39C12', WT='#73B761'),
#                    SNV_WT=c(FDR1='#000066', FDR2='#0000CC', FDR3='#3399FF', FDR4='#99CCFF', insignificant='#CCCCCC'))
# names(ann_colors$SNV_WT) <- c('FDR = 0.01', 'FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')
# # Column group annotations (match to the column names of d.cs)
# ann_col <- data.frame(group=as.factor(groups))
# rownames(ann_col) <- paste(groups, 1:length(groups), sep='_')
# 
# 
# # Significance for each of the individual comparisons
# # TAU vs control
# thresholds <- sign.table[,'SNV_vs_WT_Pvalue']
# names(thresholds) <- rownames(sign.table)
# sign.SNVWT <- add.sign.labels(bb, thresholds)
# ann_row <- data.frame(SNV_WT=factor(sign.SNVWT, labels=c('FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')))
# 
# # Top 50 most significant in each differential expression analysis, combined
# sign.rows <- bb$table$HGCN[bb$table$Pvalue <= sort(bb$table$Pvalue)[50]]
# ids.sign.all <- ids[ids %in% sign.rows]
# ann_row4 <- data.frame(SNV_WT=factor(sign.SNVWT[ids.sign.all], labels=c('FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')))
# rownames(ann_row4) <- ids.sign.all
# 
# 
# # Create clustering heatmaps
# h.all <- Clustering(d.cs, ann_col=ann_col, anncol=ann_colors, f='figures/clustering/clustering_all.png')
# h.SNVWT <- Clustering(d.cs, ann_col=ann_col, ann_row=ann_row4, anncol=ann_colors, f='figures/clustering/clustering_SNVWT_all.png')
# h.SNVWT.sign <- Clustering(d.cs[sign.SNVWT[rownames(d.cs)]<5,], ann_col=ann_col, ann_row=data.frame(SNV_WT=ann_row[sign.SNVWT<5,], row.names=rownames(ann_row)[sign.SNVWT<5]),
#                             anncol=ann_colors, f='figures/clustering/clustering_SNVWT_sign.png')
# h.top50.sign <- Clustering(d.cs[ids.sign.all,], ann_col=ann_col, ann_row=ann_row4,
#                             anncol=ann_colors, f='figures/clustering/clustering_top50_sign.png')

#######################################
###   WGCNA coexpression analysis   ###
#######################################
#load(file='rdata/input_data.RData')
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#source('src/coexpression_analysis_prep.R')
dir.create('figures/coexpression')
dir.create('output/coexpression')

time1 <- Sys.time()
d.log2vals <- as.data.frame(d.summary[, c('log2FC')], row.names=row.names(d.summary))
colnames(d.log2vals) <- c('log2FC')
coexpression <- coexpression.analysis(t(d.adj), d.log2vals, 'output/coexpression', 'figures/coexpression')
wgcna.net <- coexpression[[1]]
module.significance <- coexpression[[2]]
power <- coexpression[[3]]
time2 <- Sys.time()
print("WGCNA:")
difftime(time2, time1, units="secs")
# Merge M18 (lightgreen) into M9 (magenta), since they were highly similar
#wgcna.net$colors <- replace(wgcna.net$colors, wgcna.net$colors==18, 9)
# Module labels and log2 values
moduleColors <- labels2colors(wgcna.net$colors)
modules <- levels(as.factor(moduleColors))
#gns<-mapIds(org.Hs.eg.db, keys = ids, column = "ENTREZID", keytype = "SYMBOL")
#GOenr.net <- GO.terms.modules(wgcna.net, ids, log2vals, gns, 'figures/coexpression', 'output/coexpression')
ids <- ids[moduleColors != "grey"]
log2vals <- d.summary[ids, "log2FC"]
pvals <- d.summary[ids, "Pvalue"]
names(log2vals) <- ids
names(pvals) <- ids
# Save data structures
# save(wgcna.net, GOenr.net, module.significance, ids, file='rdata/coexpression.RData')

### PLOTS
TOM <- TOMsimilarityFromExpr(t(d.adj[ids,]), power=power, TOMType="unsigned")
diss2 <- 1-TOM
hier2 <- hclust(as.dist(diss2), method="average")
colorDynamicTOM <- labels2colors(cutreeDynamic(hier2,method="tree"))
#restGenes <- (colorDynamicTOM != "grey")
diag(diss2) = NA;
moduleColors2 <- moduleColors[moduleColors != "grey"]

#restGenes <- (colorDynamicTOM != "grey")
#diss2 <- 1 - TOMsimilarityFromExpr(t(d.adj[restGenes,]), power = power)
#hier2 <- hclust(as.dist(diss2), method="average" )
#diag(diss2) = NA;

par(mar = c(1, 1, 1, 1))
png(file=paste0('figures/coexpression/TOMheatmap.png'), width=1920, height=1020)
TOMplot(diss2^4, hier2, main = "TOM heatmap plot, module genes", terrainColors = FALSE, colors = moduleColors2, cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

png(file=paste0('figures/coexpression/networkHeatmap.png'), width=1920, height=1020)
plotNetworkHeatmap(
  t(d.adj[ids,]),
  ids,
  useTOM = FALSE,
  power = power,
  networkType = "unsigned",
  main = "Heatmap of the network")
dev.off()

par(mar = c(2, 1, 1, 1))
png(file=paste0('figures/coexpression/moduleSignificance.png'), width=1920, height=1020)
plotModuleSignificance(
  d.summary[ids, "Pvalue"],
  moduleColors2,
  boxplot = FALSE,
  main = "Gene significance across modules,",
  ylab = "Gene Significance", cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

datME <- moduleEigengenes(t(d.adj[ids,]),colorDynamicTOM)$eigengenes
dissimME <- (1-t(cor(datME, use="p")))/2
hclustdatME <- hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
png(file=paste0('figures/coexpression/moduleEigengenes.png'), width=1920, height=1020)
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes", cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

#png(file=paste0('figures/coexpression/MEpairs.png'), width=1920, height=1020)
#MEpairs <- plotMEpairs(
#  datME,
#  main = "Relationship between module eigengenes",
#  clusterMEs = TRUE)
#dev.off()

GS1 <- as.numeric(cor(group.KRAS,t(d.adj[ids,])))
GeneSignificance <- abs(GS1)
ModuleSignificance <- tapply(GeneSignificance, colorDynamicTOM, mean, na.rm=T)
png(file=paste0('figures/coexpression/moduleSignificance2.png'), width=1920, height=1020)
plotModuleSignificance(GeneSignificance,colorDynamicTOM, cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()

disableWGCNAThreads()

MEs <- orderMEs(datME)
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(t(d.adj[ids,]), MEs, use = "p"));
names(geneModuleMembership) <- paste("MM", modNames, sep="");
GS1 <- as.data.frame(GS1)
names(GS1) = paste("GS.", names(groups), sep="")

corPVal <- corAndPvalue(t(d.adj[ids,]), use="pairwise.complete.obs")


########################################
####       Hierarchical HotNet       ###
########################################
#### Run hierarchical hotnet to obtain the most significant submodule within each of the identified modules
#### Note that the HotNet package was written in Python and needs to be intstalled separately on your machine
#### HotNet was run mostly with default parameters
#load(file='rdata/input_data.RData')
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#load(file='rdata/coexpression.RData')
dir.create('figures/hotnet')
dir.create('output/hotnet')
dir.create('output/hotnet/HotNet_input')

corN <- corPVal[["cor"]]
corN <- corN[corN < 0]
corP <- corPVal[["cor"]]
corP <- corP[corP > 0]

cutOff <- quantile(as.vector(TOM), probs=0.9)
print(paste("90th quantile: ", cutOff))

cutOffN <- quantile(as.vector(corN), probs=0.1)
cutOffP <- quantile(as.vector(corP), probs=0.9)
print(paste("90th negative quantile: ", cutOffN))
print(paste("90th positive quantile: ", cutOffP))

hotnetColors <- c()

#Setup backend to use many processors
totalCores = detectCores()
#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1)
registerDoParallel(cluster)

# Get table with interactions for each module
for (module in modules){
  if (module == "grey"){
    next
  } 
  else {
    print(module)
    par(mar = c(1, 1, 1, 1))
    png(file=paste0('figures/coexpression/', module, '_matrix.png'), width=1920, height=1020)
    plotMat(t(scale(t(d.adj[colorDynamicTOM==module,]))),rlabels=T,
            clabels=T,rcols=module,
            title=module, cex.lab = 2, cex.axis = 2, cex.main = 2)
    dev.off()
    
    ME=datME[, paste("ME",module, sep="")]
    par(mar = c(2, 1, 1, 1))
    png(file=paste0('figures/coexpression/', module, '_ME.png'), width=1920, height=1020)
    barplot(ME, col=module, main="", cex.main=2,
            ylab="eigengene expression",xlab="array sample", cex.lab = 2, cex.axis = 2, cex.main = 2)
    dev.off()
    # Select proteins in module
    inModule <- moduleColors2 == module
    sum(inModule)
    TOMmodule <- TOM[inModule, inModule]
    #Pmodule <- corPVal[["p"]][inModule, inModule]
    idsModule <- ids[inModule]
    log2Module <- log2vals[inModule]
    PModule <- pvals[inModule]
    
    geneNames <- mapIds(org.Hs.eg.db, keys = idsModule, column = "ENTREZID", keytype = "SYMBOL")
    goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
    if (!(is.null(goBP)) && nrow(goBP) > 0){
      print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust < 0.1,])))
      png(file=paste0('figures/coexpression/GO_', module, '.png'), width=1920, height=1020)
      print(dotplot(goBP, showCategory=10))
      dev.off()
    }
    
    par(mfrow = c(1,1));
    column <- match(module, modNames);
    png(file=paste0('figures/coexpression/', module, '_membershipVSsignficance.png'), width=1920, height=1020)
    verboseScatterplot(abs(geneModuleMembership[inModule, column]),
                       abs(GS1[inModule, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()
    
    for (i in 1:length(idsModule)){
      # If there is no gene name associated with the Uniprot ID, keep the uniprot ID
      if (is.na(idsModule[i])){
        idsModule[i] <- idsModule[i]
      }
    }
    node.frame <- data.frame(Symbol=idsModule, LogFC=log2Module, Pvalue=PModule)
    rownames(node.frame) <- 1:nrow(node.frame)
    #node.frame$InvPvalue <- 1 - PModuleProteins
    #write.table(node.frame, file=paste0('output/hotnet/nodes_', module, '.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
    # write.table(node.frame[, c('Symbol', 'InvPvalue')], file=paste0('output/hotnet/HotNet_input/g2s_Pval_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote     =FALSE)
    names(idsModule) <- idsModule
    
    write.table(node.frame[, c('Symbol', 'LogFC')], file=paste0('output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
    #nodeColorsModule <- node.colors[inModule]
    
    # Create empty dataframes
    time1 <- Sys.time()
    edges.matrix.1 <- data.frame(matrix(nrow=0, ncol=2))
    #  edges.matrix.2 <- data.frame(matrix(nrow=0, ncol=2))
    edges.matrix.1.num <- data.frame(matrix(nrow=0, ncol=2))
    edges.matrix.2.num <- data.frame(matrix(nrow=0, ncol=2))
    i2g.1 <- data.frame(matrix(nrow=0, ncol=2))
    i2g.2 <- data.frame(matrix(nrow=0, ncol=2))
    # Write tables: one with edges between all nodes, one with a treshold of 0.05 and one with custom thresholds
    #result = foreach(i = 1:(nrow(TOMmodule)-1), j = (i+1):nrow(TOMmodule)) %dopar% {
    combine <- function(x,y){
      #names(x[[3]]) <- names(y[[3]])
      w1 <- rbind(x[[1]], y[[1]])
      x1 <- rbind(x[[2]], y[[2]])
      y1 <- rbind(x[[3]], y[[3]])
      z1 <- rbind(x[[4]], y[[4]])
      y1 <- y1[!duplicated(y1), ]
      z1 <- z1[!duplicated(z1), ]
      return(list(w1,x1,y1,z1))
    }
    
    result <- foreach(i = 1:(nrow(TOMmodule)-1), .combine = 'combine', .inorder = TRUE) %dopar% {
      result1 <- data.frame(matrix(nrow=0,ncol=2))
      result2 <- data.frame(matrix(nrow=0,ncol=2))
      result3 <- data.frame(matrix(nrow=0,ncol=2))
      result4 <- data.frame(matrix(nrow=0,ncol=2))
      for (j in (i+1):nrow(TOMmodule)){
        if (TOMmodule[i,j] >= cutOff){
          #   # Add edge to list
          result1[nrow(result1)+1,] <- c(idsModule[i], idsModule[j])
          result2[nrow(result2)+1,] <- c(i,j)
          #   # Add node to node list
          if (!i %in% result3[,1]){
            result3[nrow(result3)+1,] <- c(i, idsModule[i])
            result4[nrow(result4)+1,] <- c(idsModule[i], log2Module[i])
          } 
          if (!j %in% result3[,1]){
            result3[nrow(result3)+1,] <- c(j, idsModule[j])
            result4[nrow(result4)+1,] <- c(idsModule[j], log2Module[j])
          }
        }
      }
      return(list(result1, result2, result3, result4))
    }
    time2 <- Sys.time()
    print("Slicing input:")
    print(difftime(time2, time1, units="secs"))
    time1 <- Sys.time()
    write.table(result[[1]], file=paste0('output/hotnet/HotNet_input/name_edges_expression_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
    write.table(result[[2]], file=paste0('output/hotnet/HotNet_input/edge_list_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
    write.table(result[[3]], file=paste0('output/hotnet/HotNet_input/i2g_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
    write.table(result[[4]], file=paste0('output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
    time2 <- Sys.time()
    print("Writing data:")
    print(difftime(time2, time1, units="secs"))
    if (file.info(paste('output/hotnet/HotNet_input/i2g_', module, '.tsv', sep=''))$size != 0){
      hotnetColors <- append(hotnetColors, module)
    }
  }
}
stopCluster(cluster)
#modulesHotNet <- c()
#for (module in modules){
#  if (file.info(paste0('output/hotnet/HotNet_input/edge_list_', module, '.tsv', sep=''))$size == 0){
#    next
#  }
#  else {
#    append(modulesHotNet, c(module))
#  }
#}
#
#print(modulesHotNet)


######### Run hierarchical hotnet to obtain the most significant submodule within each of the identified modules
######### Note that the HotNet package was written in Python and needs to be intstalled separately on your machine
time1 <- Sys.time()
system(paste('bash', paste(argv[1],'/src/run_hierarchicalHotnet_modules.sh "', sep=""), paste(hotnetColors, collapse=' '), '" ', cutOff))
system('module load R')
time2 <- Sys.time()
print("Hierarchical HotNet:")
difftime(time2, time1, units="secs")

time1 <- Sys.time()
dir.create('figures/GO')
dir.create('output/GO')

#GO
for (module in hotnetColors){
  if (file.exists(file=paste0('output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'))){
    p <- read.table(file=paste0('output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'), sep='\t', header= FALSE, comment.char="")[6, "V1"]
    p <- as.numeric(sub(".*: ", "", p))
    if (p <= 0.1){
      #Read in genes
      if (file.exists(paste('output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''))){
        if (file.info(paste('output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''))$size == 0){
          next
        }
        else {
          print(module)
          inSubnetwork <- as.vector(t(read.csv(paste('output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''), sep='\t', header = FALSE)[1,]))
          print(paste("In subnetwork: ", inSubnetwork))
          log2Subnetwork <- log2vals[inSubnetwork]
          print(paste("Log-fold changes: ", log2Subnetwork))
          PSubnetwork <- pvals[inSubnetwork]
          print(paste("P-values: ", PSubnetwork))
          geneNames <- mapIds(org.Hs.eg.db, keys = inSubnetwork, column = "ENTREZID", keytype = "SYMBOL")
          print(paste("Number of genes: ", length(geneNames)))
          if (is.null(geneNames)){
            next
          }
          else {
            goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
            goMF <- enrichGO(geneNames, ont="MF", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
            if (!(is.null(goBP)) && nrow(goBP) > 0){
              print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust < 0.1,])))
              write.table(goBP@result, file=paste0('output/GO/', module, 'BP.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
              png(file=paste0('figures/GO/GOBPBar_', module, '.png'), width=1920, height=1020)
              print(barplot(goBP, showCategory=10, cex.lab = 2, cex.axis = 2, cex.main = 2))
              dev.off()
              png(file=paste0('figures/GO/GOBPDot_', module, '.png'), width=1920, height=1020)
              print(dotplot(goBP, showCategory=10))
              dev.off()
              #              png(file=paste0('figures/GO/GOBPNetwork_', module, '.png'), width=1920, height=1020)
              #              print(goplot(goBP))
              #              dev.off()
              #              terms <- goBP$Description[1:10]
              #              png(file=paste0('figures/GO/GOBPPub_', module, '.png'), width=1920, height=1020)
              #              print(pmcplot(terms, 2010:2022))
              #              dev.off()
            }
            if (!(is.null(goMF)) && nrow(goMF) > 0){
              print(paste("Number of enriched molecular functions: ", nrow(goMF[goMF$p.adjust < 0.1,])))
              write.table(goMF@result, file=paste0('output/GO/', module, 'MF.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
              png(file=paste0('figures/GO/GOMFBar_', module, '.png'), width=1920, height=1020)
              print(barplot(goMF, showCategory=10, cex.lab = 2, cex.axis = 2, cex.main = 2))
              dev.off()
              png(file=paste0('figures/GO/GOMFDot_', module, '.png'), width=1920, height=1020)
              print(dotplot(goMF, showCategory=10))
              dev.off()
              #              png(file=paste0('figures/GO/GOMFNetwork_', module, '.png'), width=1920, height=1020)
              #              print(goplot(goMF))
              #              dev.off()
              #              terms <- goMF$Description[1:10]
              #              png(file=paste0('figures/GO/GOMFPub_', module, '.png'), width=1920, height=1020)
              #              print(pmcplot(terms, 2010:2022))
              #              dev.off()
            }
          }
        }
      }
    }  
  }
}

time2 <- Sys.time()
print("GO:")
difftime(time2, time1, units="secs")
endTime <- Sys.time()
print("Total time:")
difftime(endTime, startTime, units="secs")
#######################################
###          Validation             ###
#######################################
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#lo+ rownames(meta.val) <- meta.val$Gel.Code
# 
# ids1 <- d.val1$UniProt_accession
# ids1 <- as.character(unlist(lapply(ids1, function(x){strsplit(as.character(x), ';')[[1]][1]})))
# ids1 <- as.character(unlist(lapply(ids1, function(x){strsplit(as.character(x), '-')[[1]][1]})))
# ids2 <- d.val2$UniProt_accession
# ids2 <- as.character(unlist(lapply(ids2, function(x){strsplit(as.character(x), ';')[[1]][1]})))
# ids2 <- as.character(unlist(lapply(ids2, function(x){strsplit(as.character(x), '-')[[1]][1]}))) 
# 
# # Normalize data
# d.val1.norm <- normalize.sample(d.val1.val)
# d.val2.norm <- normalize.sample(d.val2.val)
# 
# # Beta binomial test on validation sets
# groups1 <- meta.val[colnames(d.val1.norm),"Condition"]
# groups2 <- meta.val[colnames(d.val2.norm),"Condition"]
# bb.1 <- betaBinomial(d.val1.norm, ids1, groups1, 'control', 'TAU', 'two.sided')
# bb.2 <- betaBinomial(d.val2.norm, ids2, groups2, 'control', 'TAU', 'two.sided')
# 
# # Write results
# results1.tab <- data.frame(Uniprot=ids1, GeneName=d.val1$UniProt_gene_symbol, Log2FC=bb.1$table[ids1, 'Log2ratio'], Pvalue=bb.1$table[ids1, 'Pvalue'])
# results2.tab <- data.frame(Uniprot=ids2, GeneName=d.val2$UniProt_gene_symbol, Log2FC=bb.2$table[ids2, 'Log2ratio'], Pvalue=bb.2$table[ids2, 'Pvalue'])
# rownames(results1.tab) <- ids1
# rownames(results2.tab) <- ids2
# 
# write.table(results1.tab, 'output/validation/log2_Pval_val_frontal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(results2.tab, 'output/validation/log2_Pval_val_temporal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
# write.table(bb.1$FDRs, 'output/validation/sign_test_val_frontal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
# write.table(bb.2$FDRs, 'output/validation/sign_test_val_temporal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
# 
# # Merge significance info with UniProt stats (if you run the workflow from the start the UniProt data is already loaded)
# d.uniprot <- read.table('data/uniprot_stats_9606.tab', header=TRUE, sep='\t', quote="", fill=TRUE, na.strings="")
# rownames(d.uniprot) <- d.uniprot$Entry
# d.uniprot.sub1 <- d.uniprot[ids1,]
# d.summary1 <- cbind(results1.tab, d.uniprot.sub1)
# d.uniprot.sub2 <- d.uniprot[ids2,]
# d.summary2 <- cbind(results2.tab, d.uniprot.sub2)
# 
# # Expression values
# ids.TAUcon <- rownames(d.summary)
# ids.valFrn <- rownames(d.summary1)
# ids.valTem <- rownames(d.summary2)
# log.TAUcon <- d.summary$log2FC_TAU_control
# P.TAUcon <- d.summary$Pvalue_TAU_control
# log.valFrn <- d.summary1$Log2FC
# P.valFrn <- d.summary1$Pvalue
# log.valTem <- d.summary2$Log2FC
# P.valTem <- d.summary2$Pvalue
# 
# names(log.TAUcon) <- ids.TAUcon
# names(P.TAUcon) <- ids.TAUcon
# names(log.valFrn) <- ids.valFrn
# names(P.valFrn) <- ids.valFrn
# names(log.valTem) <- ids.valTem
# names(P.valTem) <- ids.valTem
# names(ids.TAUcon) <- gene.names
# 
# # Get module IDs
# moduleColors <- labels2colors(wgcna.net$colors)
# modules <- levels(as.factor(moduleColors))
# names(moduleColors) <- ids.TAUcon
# 
# # Test against frontal dataset
# # Filter out proteins that are found in one of the two sets (Note that proteins in FC2/P2 that are not in FC1/P1 are
# # automatically ignored in fraction.correct, so there is no need to filter that dataset)
# ids.both.frn <- ids.TAUcon %in% ids.valFrn
# log.TAUcon.toTest.frn <- log.TAUcon[ids.both.frn]
# P.TAUcon.toTest.frn <- P.TAUcon[ids.both.frn]
# 
# # Permutation test on all proteins (frontal)
# outfolder <- 'output/validation/all_frontal/'
# figfolder <- 'figures/validation/all_frontal/'
# dir.create(outfolder, showWarnings=FALSE)
# dir.create(figfolder, showWarnings=FALSE)
# perm.allFrn <- permutation.test(log.TAUcon.toTest.frn, log.valFrn, P.TAUcon.toTest.frn, P.valFrn, outfolder, figfolder)
# save(perm.allFrn, file='rdata/validation/permutation_all_frontal.RData')
# 
# # Permutation test per module (frontal)
# moduleColorsBoth.frn <- moduleColors[ids.both.frn]
# for(m in modules){
#   print(m)
#   inModule <- moduleColorsBoth.frn == m
#   m.ids <- names(moduleColorsBoth.frn[inModule])
#   #log.TAUcon.module <- log.TAUcon.toTest[inModule]
#   #P.TAUcon.module <- P.TAUcon.toTest[inModule]
#   outfolder <- paste0('output/validation/module_', m, '_frontal/')
#   figfolder <- paste0('figures/validation/module_', m, '_frontal/')
#   dir.create(outfolder, showWarnings=FALSE)
#   dir.create(figfolder, showWarnings=FALSE)
#   perm.module <- permutation.test(log.TAUcon, log.valFrn, P.TAUcon, P.valFrn, outfolder, figfolder, m.ids=m.ids)
#   save(perm.module, file=paste0('rdata/validation/permutation_', m, '_frontal.RData'))
#   
#   # Do the same, but only for those ids that were significant in HotNet
#   #if (m == 'red'){
#   #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_001_', m, '.tsv'), header=FALSE, quote='', sep='\t')
#   #} else if (m == 'grey'){
#   #  next
#   #} else {
#   #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_003_', m, '.tsv'), header=FALSE, quote='', sep='\t')
#   #}
#   
#   #names.hotnet <- levels(as.factor(c(d.hotnet[,1], d.hotnet[,2])))
#   #ids.hotnet <- as.character(ids.TAUcon.named[names.hotnet])
#   #ids.both <- ids.TAUcon %in% ids.valFrn & ids.TAUcon %in% ids.hotnet
#   #inModule <- moduleColors[ids.both] == m
#   #m.ids <- names(moduleColors[ids.both][inModule])
#   #log.TAUcon.module <- log.TAUcon[ids.both]
#   #P.TAUcon.module <- P.TAUcon[ids.both]
#   #outfolder <- paste0('Figures/Validation4/Module_', m, '_frontal_hotnet/')
#   #dir.create(outfolder, showWarnings=FALSE)
#   #perm.module.hotnet <- permutation.test(log.TAUcon, log.valFrn, P.TAUcon, P.valFrn, outfolder, m.ids=m.ids)
#   #save(perm.module.hotnet, file=paste0('Output/Permutation_', m, '_frontal4_hotnet.RData'))
# }
# 
# 
# # Test against temporal dataset
# # Filter out proteins that are found in one of the two sets
# ids.both.tem <- ids.TAUcon %in% ids.valTem
# log.TAUcon.toTest.tem <- log.TAUcon[ids.both.tem]
# P.TAUcon.toTest.tem <- P.TAUcon[ids.both.tem]
# # Permutation test on all proteins
# outfolder <- 'output/validation/all_temporal/'
# figfolder <- 'figures/validation/all_temporal/'
# dir.create(outfolder, showWarnings=FALSE)
# dir.create(figfolder, showWarnings=FALSE)
# perm.allTem <- permutation.test(log.TAUcon.toTest.tem, log.valTem, P.TAUcon.toTest.tem, P.valTem, outfolder, figfolder)
# save(perm.allTem, file='rdata/validation/permutation_all_temporal.RData')
# 
# # Permutation test per module
# moduleColorsBoth.tem <- moduleColors[ids.both.tem]
# for(m in modules){
#   print(m)
#   inModule <- moduleColorsBoth.tem == m
#   #log.TAUcon.module <- log.TAUcon.toTest[inModule]
#   #P.TAUcon.module <- P.TAUcon.toTest[inModule]
#   m.ids <- names(moduleColorsBoth.tem[inModule])
#   outfolder <- paste0('output/validation/module_', m, '_temporal/')
#   figfolder <- paste0('figures/validation/module_', m, '_temporal/')
#   dir.create(outfolder, showWarnings=FALSE)
#   dir.create(figfolder, showWarnings=FALSE)
#   perm.module <- permutation.test(log.TAUcon, log.valTem, P.TAUcon, P.valTem, outfolder, figfolder, m.ids=m.ids)
#   save(perm.module, file=paste0('rdata/validation/permutation_', m, '_temporal.RData'))
#   
#   ## Do the same, but only for those ids that were significant in HotNet
#   #if (m == 'red'){
#   #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_001_', m, '.tsv'), header=FALSE, quote='', sep='\t')
#   #} else if (m == 'grey') {
#   #  # Skip empty files
#   #  next
#   #} else {
#   #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_003_', m, '.tsv'), header=FALSE, quote='', sep='\t')
#   #}
#   
#   #names.hotnet <- levels(as.factor(c(d.hotnet[,1], d.hotnet[,2])))
#   #ids.hotnet <- as.character(ids.TAUcon.named[names.hotnet])
#   #ids.both <- ids.TAUcon %in% ids.valTem & ids.TAUcon %in% ids.hotnet
#   #inModule <- moduleColors[ids.both] == m
#   ##log.TAUcon.module <- log.TAUcon[ids.both]
#   ##P.TAUcon.module <- P.TAUcon[ids.both]
#   #m.ids <- names(moduleColors[ids.both][inModule])
#   #outfolder <- paste0('Figures/Validation4/Module_', m, '_temporal_hotnet/')
#   #dir.create(outfolder, showWarnings=FALSE)
#   #perm.module.hotnet <- permutation.test(log.TAUcon, log.valTem, P.TAUcon, P.valTem, outfolder, m.ids=m.ids)
#   #save(perm.module.hotnet, file=paste0('Output/Permutation_', m, '_temporal4_hotnet.RData'))
# }
# 
# #######################################
# ###       Create module plots       ###
# #######################################
# #load(file='rdata/normalized_data.RData')
# #load(file='rdata/differential_expression.RData')
# #load(file='rdata/coexpression.RData')
# d.norm <- t(d.norm)
# moduleColors <- labels2colors(wgcna.net$colors)
# 
# d.array <- data.frame(matrix(nrow=5, ncol=18))
# colnames(d.array) <- labels2colors(0:17)
# rownames(d.array) <- c('FTLD-tau vs NHC', 'FTLD-TDP vs NHC', 'FTLD-tau vs FTLD-TDP', 'FTLD-tau vs NCH (validation frontal)', 'FTLD-tau vs NHC (validation temporal)')
# for (i in 0:17){
#   color <- labels2colors(i)
#   TAUcon <- d.summary$log2FC_TAU_control[moduleColors == color]
#   TDPcon <- d.summary$log2FC_TDP_control[moduleColors == color]
#   TAUTDP <- d.summary$log2FC_TAU_TDP[moduleColors == color]
#   
#   frn <- log.TAUcon.toTest.frn[moduleColorsBoth.frn==color]
#   tem <- log.TAUcon.toTest.tem[moduleColorsBoth.tem==color]
#   
#   mTAUcon <- median(TAUcon)
#   mTDPcon <- median(TDPcon)
#   mTAUTDP <- median(TAUTDP)
#   mfrn <- median(frn)
#   mtem <- median
#   d.array2 <- d.array
#   d.array[,color] <- c(mTAUcon, mTDPcon, mTAUTDP, mfrn, mtem)
#   
#   
#   TAUconP <- d.summary$Pvalue_TAU_control[moduleColors == color]
#   TDPconP <- d.summary$Pvalue_TDP_control[moduleColors == color]
#   TAUTDPP <- d.summary$Pvalue_TAU_TDP[moduleColors == color]
#   
#   frnP <- P.TAUcon.toTest.frn[moduleColorsBoth.frn==color]
#   temP <- P.TAUcon.toTest.tem[moduleColorsBoth.tem==color]
#   
#   
#   nTAUconP <- norm(TAUconP)
#   nTDPconP <- norm(TDPconP)
#   TAUTDPP <- norm(TAUTDPP)
#   nfrnP <- norm(frnP)
#   ntemP <- norm(tem)
#   d.array2[,color] <- c(nTAUconP, nTDPP, nTAUTDPP, nfrnP, ntemP)
#   
#   t.frn <- t.test(tem)
#   t.tem <- t.test(frn)
# }
# write.table(d.array, file='module_heatmap.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
# write.table(t(d.array2), file='module_pval.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)