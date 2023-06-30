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
source('src/similarity.R')


#######################################
###             Setup               ###
#######################################
### Set working directory
setwd('C:/Users/marle/Desktop/Y2/Internship/Project') #Replace with your working directory
dir.create('output')
#Setup backend to use many processors
totalCores = detectCores()
#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1)
registerDoParallel(cluster)

packages <- c('ibb','stringr', 'BiocManager', 'BiocParallel', 'BioNet', 'pheatmap', 'countdata', 'DESeq2', 'vsn', 'org.Hs.eg.db', 'clusterProfiler', 'ggplot2', 'pheatmap', 'RColorBrewer', 'foreach', 'doParallel', 'enrichplot', 'europepmc', 'WGCNA')

#d.norm <- d.norm[,keepSamples]
#colnames(d.norm) <- samples

###VST
#d.adj <- rlog(d.raw, blind=FALSE)

#colnames(d.raw) <- samples

#ddds <- DESeqDataSetFromMatrix(countData = d.raw, colData = group, design = ~ KRAS)
#ntd <- normTransform(ddds)
#d.adj <- assay(ddds)

###PLOTS AFTER VST
#dds <- DESeqDataSetFromMatrix(countData = d.raw, colData = group, design = ~ KRAS)
#dds
#dds2 <- estimateSizeFactors(dds)
#d.norm <- counts(dds2, normalized = TRUE)

#hits <- read.csv('./data/hitlist_snv_screen_coadread_tcga.csv', header = TRUE, sep=';', row.names=1)
hit <- "KRAS"
#foreach(hit = hits, .packages = packages) %dopar% {
### Create output directories
for (n in c("1","2","3")){
  #foreach(hit = hits, .packages = packages) %dopar% {
  ### Create output directories
  dir.create(paste0('output/', n))
  dir.create(paste0('output/', n, '/training'))
  dir.create(paste0('output/', n, '/training/output'))
  dir.create(paste0('output/', n, '/training/figures'))
  dir.create(paste0('output/', n, '/training/rdata'))
  
  ### Load data
  time1 <- Sys.time()
  d <- read.table('./data/TCGA_rna_count_data.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)
  group.data <- read.table('./data/non_silent_mutation_profile_crc.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)
  time2 <- Sys.time()
  print("Loading data:")
  difftime(time2, time1, units="secs")
  
  time1 <- Sys.time()
  d <- d[1:500,1:100]
  d <- d[,colnames(d) %in% rownames(group.data)]
  group.data <- group.data[rownames(group.data) %in% colnames(d),]
  groups <- group.data[,hit]
  new_order <- sort(colnames(d))
  d <- d[, new_order]
  new_order_rows <- sort(rownames(d))
  d <- d[new_order_rows,]
  d.raw <- d[,1:ncol(d)]
  d.raw <- apply(d.raw, c(1,2), as.numeric)
  d.raw <- d.raw[rowSums(d.raw) >= 10,]
  samples <- colnames(d.raw)
  
  d.adj <- varianceStabilizingTransformation(d.raw, blind = F)
  d.norm <- normalize.sample(d.raw)
  d.cs <- normalize.cs(d.norm)
  
  colnames(d.adj) <- groups
  
  png(file=paste0('output/', n, '/training/figures/clusterSamplesVST.png'), width=1920, height=1020)
  sampleTree <- hclust(dist(t(d.adj)), method = "ave")
  #png(file=paste0('figures/clusterSamplesNorm.png'), width=1920, height=1020)
  #sampleTree <- hclust(dist(t(d.norm)), method = "ave")
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  # Plot a line to show the cut
  abline(h = 225, col = "red");
  #abline(h = 260, col = "blue");
  #abline(h = 255, col = "green");
  # Determine cluster under the line
  clust <- cutreeStatic(sampleTree, cutHeight = 225, minSize = 10)
  #clust <- cutreeStatic(sampleTree, cutHeight = 2250000, minSize = 10)
  table(clust)
  # clust 1 contains the samples we want to keep.
  keepSamples <- (clust==1)
  dev.off()
  
  colnames(d.adj) <- samples
  d.raw <- d.raw[,keepSamples]
  group.data <- group.data[keepSamples,]
  #load(file='./data/sampleSelectionMore.RData')
  
  d.cs <- d.cs[,keepSamples]
  
  d.cs.WT <- d.cs[,rownames(group.data[group.data[,hit] == "WT",])]
  d.cs.SNV <- d.cs[,rownames(group.data[group.data[,hit] == "SNV",])]
  
  colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
  
  png(file=paste0('output/', n, '/training/figures/heatmapSNV.png'), width=1920, height=1020)
  pheatmap(d.cs.SNV,
           legend=TRUE,
           color=my_palette,
           breaks=colors,
           show_rownames=FALSE,
           show_colnames=FALSE
  )
  dev.off()
  
  png(file=paste0('output/', n, '/training/figures/heatmapWT.png'), width=1920, height=1020)
  pheatmap(d.cs.WT,
           legend=TRUE,
           color=my_palette,
           breaks=colors,
           show_rownames=FALSE,
           show_colnames=FALSE
  )
  dev.off()
  
  selection <- sample(colnames(d.raw), floor(ncol(d.raw)/2))
  save(selection, file=paste0('output/' , n, '/training/rdata/sampleSelection.RData'))
  d.training <- d.raw[,selection]
  samples <- colnames(d.training)
  ids <- rownames(d.training)
  group.data.training <- group.data[selection,]
  d.training.adj <- d.adj[,selection]
  colnames(d.training.adj) <- samples
  
  d.training.cs <- d.cs[,selection]
  
  d.training.cs.WT <- d.training.cs[,rownames(group.data.training[group.data.training[,hit] == "WT",])]
  d.training.cs.SNV <- d.training.cs[,rownames(group.data.training[group.data.training[,hit] == "SNV",])]
  #
  ##dr.WT <- dist(1-cor(t(d.training.cs.WT)))
  ##hr.WT <- hclust(dr.WT)
  ##dc.WT <- dist(1-cor(d.training.cs.WT))
  ##hc.WT <- hclust(dc.WT)
  ##
  ##dr.SNV <- dist(1-cor(t(d.training.cs.SNV)))
  ##hr.SNV <- hclust(dr.SNV)
  ##dc.SNV <- dist(1-cor(d.training.cs.SNV))
  ##hc.SNV <- hclust(dc.SNV)
  #
  colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
  
  png(file=paste0('output/', n, '/training/figures/heatmapSNVSubset.png'), width=1920, height=1020)
  pheatmap(d.training.cs.SNV,
           legend=TRUE,
           color=my_palette,
           breaks=colors,
           show_rownames=FALSE,
           show_colnames=FALSE
  )
  dev.off()
  
  png(file=paste0('output/', n, '/training/figures/heatmapWTSubset.png'), width=1920, height=1020)
  pheatmap(d.training.cs.WT,
           legend=TRUE,
           color=my_palette,
           breaks=colors,
           show_rownames=FALSE,
           show_colnames=FALSE
  )
  dev.off()
  ##
  #d.training.adj.WT <- d.training.adj[,rownames(group.data.training[group.data.training[,hit] == "WT",])]
  #d.training.adj.SNV <- d.training.adj[,rownames(group.data.training[group.data.training[,hit] == "SNV",])]
  
  #png(file=paste0('output/', n, '/training/figures/sampleClusterVSTWT.png'), width=1920, height=1020)
  #plotClusterTreeSamples(t(d.training.adj.WT), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #dev.off()
  #
  #png(file=paste0('output/', n, '/training/figures/sampleClusterVSTSNV.png'), width=1920, height=1020)
  #plotClusterTreeSamples(t(d.training.adj.SNV), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #dev.off()
  #
  #png(file=paste0('output/', n, '/training/figures/heatmapVSTWT.png'), width=1920, height=1020)
  #pheatmap(d.training.adj.WT, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE)
  #dev.off()
  #
  #png(file=paste0('output/', n, '/training/figures/heatmapVSTSNV.png'), width=1920, height=1020)
  #pheatmap(d.training.adj.SNV, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE)
  #dev.off()
  
  groups <- group.data.training[,hit]
  # Save input data
  save(d.training, group.data.training, groups, file=(paste0('output/', n, '/training/rdata/input_data.RData')))
  time2 <- Sys.time()
  print("Processing data:")
  difftime(time2, time1, units="secs")
  
  group <- as.matrix(group.data.training[,hit])
  rownames(group) <- rownames(group.data.training)
  colnames(group) <- c(hit)
  
  ###PLOTS
  png(file=paste0('output/', n, '/training/figures/heatmap.png'), width=1920, height=1020)
  pheatmap(d.training, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, show_colnames = FALSE, annotation_col=as.data.frame(group), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  colnames(d.training) <- groups
  
  png(file=paste0('output/', n, '/training/figures/sampleCluster.png'), width=1920, height=1020)
  plotClusterTreeSamples(t(d.training), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  sampleDists <- dist(t(d.training))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(d.training)
  colnames(sampleDistMatrix) <- colnames(d.training)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  png(file=paste0('output/', n, '/training/figures/sampleDistances.png'), width=1920, height=1020)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/meanExpression.png'), width=1920, height=1020)
  barplot(apply(t(d.training),1,mean, na.rm=T),
          xlab = "Sample", ylab = "Mean expression",
          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/variance.png'), width=1920, height=1020)
  meanSdPlot(d.training, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  #png(file=paste0('figures/heatmapNorm.png'), width=1920, height=1020)
  #pheatmap(d.training.norm, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group))
  #dev.off()
  #
  #colnames(d.training.norm) <- groups
  #
  #sampleDists <- dist(t(d.training.norm))
  #sampleDistMatrix <- as.matrix(sampleDists)
  #rownames(sampleDistMatrix) <- colnames(d.training.norm)
  #colnames(sampleDistMatrix) <- colnames(d.training.norm)
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
  #barplot(apply(t(d.training.norm),1,mean, na.rm=T),
  #        xlab = "Sample", ylab = "Mean expression",
  #        main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #dev.off()
  
  png(file=paste0('output/', n, '/training/figures/heatmapVST.png'), width=1920, height=1020)
  pheatmap(d.training.adj, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group))
  dev.off()
  
  group[,1] <- as.numeric(factor(group[,1]))
  colnames(d.training.adj) <- groups
  
  png(file=paste0('output/', n, '/training/figures/sampleClusterVST.png'), width=1920, height=1020)
  plotClusterTreeSamples(t(d.training.adj), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/meanExpressionVST.png'), width=1920, height=1020)
  barplot(apply(t(d.training.adj),1,mean, na.rm=T),
          xlab = "Sample", ylab = "Mean expression",
          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/varianceVST.png'), width=1920, height=1020)
  meanSdPlot(d.training.adj, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  sampleDists <- dist(t(d.training.adj))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(d.training.adj)
  colnames(sampleDistMatrix) <- colnames(d.training.adj)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  png(file=paste0('output/', n, '/training/figures/sampleDistancesVST.png'), width=1920, height=1020)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  dev.off()
  
  
  #######################################
  ###         Preprocessing           ###
  #######################################
  #load(file='rdata/input_data.rdata')
  dir.create(paste0('output/', n, '/training/figures/differential expression'))
  #d.training.norm <- normalize.sample(d.training)
  #d.training.cs <- normalize.cs(d.training.norm)
  #length <- length(groups)
  ## unique colnames in d.training.cs are necessary for clustering
  #colnames(d.training.cs) <- paste(groups, 1:length, sep='_')
  ## save input data into a file
  #save(d.training.norm, d.training.cs, ids, groups, file='rdata/normalized_data.rdata')
  
  #DESeq2
  time1 <- Sys.time()
  colnames(d.training) <- samples
  dds <- DESeqDataSetFromMatrix(countData = d.training, colData = group, design = formula(paste("~",hit)))
  dds
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  png(file=paste0('output/', n, '/training/figures/differential expression/PCA.png'), width=1920, height=1020)
  plotPCA(vsd, intgroup=c(hit))
  dev.off()
  
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
  
  colnames(d.training) <- groups
  
  d.training.summary <- data.frame(unlist(res$log2FoldChange), unlist(res$padj), row.names = rownames(res))
  colnames(d.training.summary) <- c("log2FC", "Pvalue")
  d.training.summary <- na.omit(d.training.summary)
  d.training.adj <- d.training.adj[rownames(d.training.adj) %in% rownames(d.training.summary),]
  ids <- rownames(d.training.adj)
  #d.training.norm <- d.training.norm[rownames(d.training.norm) %in% rownames(d.training.summary),]
  #ids <- rownames(d.training.norm)
  time2<- Sys.time()
  print("DESeq2:")
  difftime(time2, time1, units="secs")
  
  ###DESEQ PLOTS
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/differential expression/MA.png'), width=1920, height=1020)
  plotMA(res, ylim=c(-2,2), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  png(file=paste0('output/', n, '/training/figures/differential expression/counts.png'), width=1920, height=1020)
  plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup=hit, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/differential expression/dispersion.png'), width=1920, height=1020)
  plotDispEsts(dds, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/differential expression/rejections.png'), width=1920, height=1020)
  plot(metadata(res)$filterNumRej, 
       type="b", ylab="number of rejections",
       xlab="quantiles of filter", cex.lab = 2, cex.axis = 2, cex.main = 2)
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/differential expression/outliers.png'), width=1920, height=1020)
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, names = groups)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  W <- res$stat
  maxCooks <- apply(assays(dds)[["cooks"]],1,max)
  idx <- !is.na(W)
  png(file=paste0('output/', n, '/training/figures/differential expression/Wald.png'), width=1920, height=1020)
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
  png(file=paste0('output/', n, '/training/figures/differential expression/pass.png'), width=1920, height=1020)
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
  # bb <- betaBinomial(d.training.norm, ids, groups, 'WT', 'SNV', 'two.sided')
  # print(Sys.time())
  # d.training.sign <- data.frame(log2FC=bb$table[ids,]$Log2ratio, Pvalue=bb$table[ids,]$Pvalue)
  # rownames(d.training.sign) <- ids
  # # Merge stats about P-value thresholds and number of significant proteins
  # sign.table <- cbind(bb$FDRs)
  # d.training.summary <- d.training.sign
  # # Write output into files
  # write.table(bb$table, file='output/differential_expression/sign_test.tsv', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
  # write.table(sign.table, file='output/differential_expression/FDR_pval_summary.tsv', sep='\t', quote=FALSE)
  # write.table(d.training.summary, file='output/differential_expression/log2_pval_uniprot.tsv', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
  # save(bb, sign.table, d.training.summary, file='rdata/differential_expression.RData')
  
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
  # # Column group annotations (match to the column names of d.training.cs)
  # ann_col <- data.frame(group=as.factor(groups))
  # rownames(ann_col) <- paste(groups, 1:length(groups), sep='_')
  # 
  # 
  # # Significance for each of the individual comparisons
  # # TAU vs control
  # thresholds <- sign.table[,'SNV_vs_WT_Pvalue']
  # names(thresholds) <- rownames(sign.table)
  # sign.SNVWT <- add.training.sign.labels(bb, thresholds)
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
  # h.all <- Clustering(d.training.cs, ann_col=ann_col, anncol=ann_colors, f='figures/clustering/clustering_all.png')
  # h.SNVWT <- Clustering(d.training.cs, ann_col=ann_col, ann_row=ann_row4, anncol=ann_colors, f='figures/clustering/clustering_SNVWT_all.png')
  # h.SNVWT.sign <- Clustering(d.training.cs[sign.SNVWT[rownames(d.training.cs)]<5,], ann_col=ann_col, ann_row=data.frame(SNV_WT=ann_row[sign.SNVWT<5,], row.names=rownames(ann_row)[sign.SNVWT<5]),
  #                             anncol=ann_colors, f='figures/clustering/clustering_SNVWT_sign.png')
  # h.top50.sign <- Clustering(d.training.cs[ids.sign.all,], ann_col=ann_col, ann_row=ann_row4,
  #                             anncol=ann_colors, f='figures/clustering/clustering_top50_sign.png')
  
  #######################################
  ###   WGCNA coexpression analysis   ###
  #######################################
  #load(file='rdata/input_data.RData')
  #load(file='rdata/normalized_data.RData')
  #load(file='rdata/differential_expression.RData')
  #source('src/coexpression_analysis_prep.R')
  dir.create(paste0('output/', n, '/training/figures/coexpression'))
  dir.create(paste0('output/', n, '/training/output/coexpression'))
  
  time1 <- Sys.time()
  d.training.log2vals <- as.data.frame(d.training.summary[, c('log2FC')], row.names=row.names(d.training.summary))
  colnames(d.training.log2vals) <- c('log2FC')
  coexpression <- coexpression.analysis(t(d.training.adj), d.training.log2vals, paste0('output/', n, '/training/output/coexpression'), paste0('output/', n, '/training/figures/coexpression'), n, 'training')
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
  write.table(modules, file=paste0('output/', n, '/training/output/coexpression/modules.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  #gns<-mapIds(org.Hs.eg.db, keys = ids, column = "ENTREZID", keytype = "SYMBOL")
  #GOenr.net <- GO.terms.modules(wgcna.net, ids, log2vals, gns, 'figures/coexpression', 'output/coexpression')
  ids <- ids[moduleColors != "grey"]
  log2vals <- d.training.summary[ids, "log2FC"]
  pvals <- d.training.summary[ids, "Pvalue"]
  names(log2vals) <- ids
  names(pvals) <- ids
  # Save data structures
  #(wgcna.net, moduleColors, file='rdata/coexpression.RData')
  
  ### PLOTS
  TOM <- TOMsimilarityFromExpr(t(d.training.adj[ids,]), power=power, TOMType="unsigned")
  diss2 <- 1-TOM
  hier2 <- hclust(as.dist(diss2), method="average")
  colorDynamicTOM <- labels2colors(cutreeDynamic(hier2,method="tree"))
  diag(diss2) = NA;
  moduleColors2 <- moduleColors[moduleColors != "grey"]
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/coexpression/TOMheatmap.png'), width=1920, height=1020)
  TOMplot(diss2^4, hier2, main = "TOM heatmap plot, module genes", terrainColors = FALSE, colors = moduleColors2, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  png(file=paste0('output/', n, '/training/figures/coexpression/networkHeatmap.png'), width=1920, height=1020)
  plotNetworkHeatmap(
    t(d.training.adj[ids,]),
    ids,
    useTOM = FALSE,
    power = power,
    networkType = "unsigned",
    main = "Heatmap of the network")
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/training/figures/coexpression/moduleSignificance.png'), width=1920, height=1020)
  plotModuleSignificance(
    d.training.summary[ids, "Pvalue"],
    moduleColors2,
    boxplot = FALSE,
    main = "Gene significance across modules,",
    ylab = "Gene Significance", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  datME <- moduleEigengenes(t(d.training.adj[ids,]),colorDynamicTOM)$eigengenes
  dissimME <- (1-t(cor(datME, use="p")))/2
  hclustdatME <- hclust(as.dist(dissimME), method="average" )
  # Plot the eigengene dendrogram
  png(file=paste0('output/', n, '/training/figures/coexpression/moduleEigengenes.png'), width=1920, height=1020)
  par(mfrow=c(1,1))
  plot(hclustdatME, main="Clustering tree based of the module eigengenes", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  #png(file=paste0('figures/coexpression/MEpairs.png'), width=1920, height=1020)
  #MEpairs <- plotMEpairs(
  #  datME,
  #  main = "Relationship between module eigengenes",
  #  clusterMEs = TRUE)
  #dev.off()
  
  GS1 <- as.numeric(cor(group,t(d.training.adj[ids,])))
  GeneSignificance <- abs(GS1)
  ModuleSignificance <- tapply(GeneSignificance, colorDynamicTOM, mean, na.rm=T)
  png(file=paste0('output/', n, '/training/figures/coexpression/moduleSignificance2.png'), width=1920, height=1020)
  plotModuleSignificance(GeneSignificance,colorDynamicTOM, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  disableWGCNAThreads()
  
  MEs <- orderMEs(datME)
  modNames <- substring(names(MEs), 3)
  geneModuleMembership <- as.data.frame(cor(t(d.training.adj[ids,]), MEs, use = "p"));
  names(geneModuleMembership) <- paste("MM", modNames, sep="");
  GS1 <- as.data.frame(GS1)
  names(GS1) = paste("GS.", names(groups), sep="")
  
  corPVal <- corAndPvalue(t(d.training.adj[ids,]), use="pairwise.complete.obs")
  
  
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
  dir.create(paste0('output/', n, '/training/figures/hotnet'))
  dir.create(paste0('output/', n, '/training/output//hotnet'))
  dir.create(paste0('output/', n, '/training/output/hotnet/HotNet_input'))
  
  corN <- corPVal[["cor"]]
  corN <- corN[corN < 0]
  corP <- corPVal[["cor"]]
  corP <- corP[corP > 0]
  
  probs <- 0.9
  cutOff <- quantile(as.vector(TOM), probs=probs)
  print(paste("Quantile at ", probs, " :", cutOff))
  
  cutOffN <- quantile(as.vector(corN), probs=0.1)
  cutOffP <- quantile(as.vector(corP), probs=0.9)
  print(paste("90th negative quantile: ", cutOffN))
  print(paste("90th positive quantile: ", cutOffP))
  
  hotnetColors <- c()
  
  # Get table with interactions for each module
  for (module in modules){
    if (module == "grey"){
      next
    } 
    else {
      print(module)
      par(mar = c(1, 1, 1, 1))
      png(file=paste0('output/', n, '/training/figures/coexpression/', module, '_matrix.png'), width=1920, height=1020)
      plotMat(t(scale(t(d.training.adj[colorDynamicTOM==module,]))),rlabels=T,
              clabels=T,rcols=module,
              title=module, cex.lab = 2, cex.axis = 2, cex.main = 2)
      dev.off()
      
      ME=datME[, paste("ME",module, sep="")]
      par(mar = c(2, 1, 1, 1))
      png(file=paste0('output/', n, '/training/figures/coexpression/', module, '_ME.png'), width=1920, height=1020)
      barplot(ME, col=module, main="", cex.main=2,
              ylab="eigengene expression",xlab="array sample", cex.lab = 2, cex.axis = 2, cex.main = 2)
      dev.off()
      # Select proteins in module
      inModule <- moduleColors2 == module
      sum(inModule)
      TOMmodule <- TOM[inModule, inModule]
      corModule <- corPVal[["cor"]][inModule, inModule]
      #Pmodule <- corPVal[["p"]][inModule, inModule]
      idsModule <- ids[inModule]
      write.table(idsModule, file=paste0('output/', n, '/training/output/coexpression/', module, '.tsv'), sep='\t', row.names=FALSE, col.names=FALSE)
      log2Module <- log2vals[inModule]
      PModule <- pvals[inModule]
      
      geneNames <- mapIds(org.Hs.eg.db, keys = idsModule, column = "ENTREZID", keytype = "SYMBOL")
      goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
      if (!(is.null(goBP)) && nrow(goBP) > 0){
        print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust < 0.1,])))
        png(file=paste0('output/', n, '/training/figures/coexpression/GO_', module, '.png'), width=1920, height=1020)
        print(dotplot(goBP, showCategory=10))
        dev.off()
      }
      
      par(mfrow = c(1,1));
      column <- match(module, modNames);
      png(file=paste0('output/', n, '/training/figures/coexpression/', module, '_membershipVSsignficance.png'), width=1920, height=1020)
      verboseScatterplot(abs(geneModuleMembership[inModule, column]),
                         abs(GS1[inModule, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
      node.frame <- data.frame(Symbol=idsModule, LogFC=log2Module, Pvalue=PModule)
      rownames(node.frame) <- 1:nrow(node.frame)
      #node.frame$InvPvalue <- 1 - PModuleProteins
      #write.table(node.frame, file=paste0('output/hotnet/nodes_', module, '.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
      # write.table(node.frame[, c('Symbol', 'InvPvalue')], file=paste0('output/hotnet/HotNet_input/g2s_Pval_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote     =FALSE)
      names(idsModule) <- idsModule
      
      write.table(node.frame[, c('Symbol', 'LogFC')], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
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
      #result = foreach(n = 1:(nrow(TOMmodule)-1), j = (n+1):nrow(TOMmodule)) %dopar% {
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
          if (corModule[i,j] >= 0.5 || corModule[i,j] <= -0.5){
            #   # Add edge to list
            result1[nrow(result1)+1,] <- c(idsModule[i], idsModule[j])
            result2[nrow(result2)+1,] <- c(i,j)
            #   # Add node to node list
            if (!i %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(i, idsModule[i])
              result4[nrow(result4)+1,] <- c(idsModule[i], log2Module[i]*(1-PModule[i]))
            } 
            if (!j %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(j, idsModule[j])
              result4[nrow(result4)+1,] <- c(idsModule[j], log2Module[j]*(1-PModule[j]))
            }
          }
        }
        return(list(result1, result2, result3, result4))
      }
      time2 <- Sys.time()
      print("Slicing input:")
      print(difftime(time2, time1, units="secs"))
      time1 <- Sys.time()
      write.table(result[[1]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/name_edges_expression_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[2]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/edge_list_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[3]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/i2g_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[4]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      time2 <- Sys.time()
      print("Writing data:")
      print(difftime(time2, time1, units="secs"))
      if (file.info(paste('output/', n, '/training/output/hotnet/HotNet_input/i2g_', module, '.tsv', sep=''))$size != 0){
        hotnetColors <- append(hotnetColors, module)
      }
    }
  }
  
  
  ######### Run hierarchical hotnet to obtain the most significant submodule within each of the identified modules
  ######### Note that the HotNet package was written in Python and needs to be intstalled separately on your machine
  time1 <- Sys.time()
  system(paste('bash src/run_hierarchicalHotnet_modules.sh "', paste(modules, collapse=' '), '" ', cutOff, n, 'training'))
  time2 <- Sys.time()
  print("Hierarchical HotNet:")
  difftime(time2, time1, units="secs")
  
  time1 <- Sys.time()
  dir.create(paste0('output/', n, '/training/figures/GO'))
  dir.create(paste0('output/', n, '/training/output/GO'))
  
  BPterms <- c()
  MFterms <- c()
  genes <- c()
  subnetworks <- c()
  hotnetSubnetworks <- c()
  
  #GO
  for (module in hotnetColors){
    if (file.exists(file=paste0('output/', n, '/training/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'))){
      p <- read.table(file=paste0('output/', n, '/training/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'), sep='\t', header= FALSE, comment.char="")[6, "V1"]
      p <- as.numeric(sub(".*: ", "", p))
      if (file.exists(paste('output/', n, '/training/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''))){
        if (file.info(paste('output/', n, '/training/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''))$size == 0){
          next
        }
      }
      hotnetSubnetworks <- append(hotnetSubnetworks, module)
      if (p <= 0.1){
        #Read in genes
        print(module)
        subnetworks <- append(subnetworks, module)
        inSubnetwork <- as.vector(t(read.csv(paste('output/', n, '/training/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''), sep='\t', header = FALSE)[1,]))
        print(paste("In subnetwork: ", inSubnetwork))
        genes <- append(genes, inSubnetwork)
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
            print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust <= 0.1,])))
            write.table(goBP@result, file=paste0('output/', n, '/training/output/GO/', module, 'BP.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
            BPterms <- append(BPterms, rownames(goBP@result[goBP@result$p.adjust <= 0.1,]))
            png(file=paste0('output/', n, '/training/figures/GO/GOBPBar_', module, '.png'), width=1920, height=1020)
            print(barplot(goBP, showCategory=10, cex.lab = 2, cex.axis = 2, cex.main = 2))
            dev.off()
            png(file=paste0('output/', n, '/training/figures/GO/GOBPDot_', module, '.png'), width=1920, height=1020)
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
            print(paste("Number of enriched molecular functions: ", nrow(goMF[goMF$p.adjust <= 0.1,])))
            write.table(goMF@result, file=paste0('output/', n, '/training/output/GO/', module, 'MF.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
            MFterms <- append(MFterms, rownames(goMF@result[goMF@result$p.adjust <= 0.1,]))
            png(file=paste0('output/', n, '/training/figures/GO/GOMFBar_', module, '.png'), width=1920, height=1020)
            print(barplot(goMF, showCategory=10, cex.lab = 2, cex.axis = 2, cex.main = 2))
            dev.off()
            png(file=paste0('output/', n, '/training/figures/GO/GOMFDot_', module, '.png'), width=1920, height=1020)
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
  
  write.table(BPterms, file=paste0('output/', n, '/training/output/GO/BP.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(MFterms, file=paste0('output/', n, '/training/output/GO/MF.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(genes, file=paste0('output/', n, '/training/output/GO/genes.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(subnetworks, file=paste0('output/', n, '/training/output/GO/subnetworks.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(hotnetSubnetworks, file=paste0('output/', n, '/training/output/GO/hotnetSubnetworks.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  save(BPterms, file=paste0('output/', n, '/training/rdata/BP.RData'))
  save(MFterms, file=paste0('output/', n, '/training/rdata/MF.RData'))
  save(genes, file=paste0('output/', n, '/training/rdata/genes.RData'))
  save(subnetworks, file=paste0('output/', n, '/training/rdata/subnetworks.RData'))
  
  time2 <- Sys.time()
  print("GO:")
  difftime(time2, time1, units="secs")
  endTime <- Sys.time()
  print("Total time:")
  difftime(endTime, startTime, units="secs")
  ######################################
  ##       Internal validation       ###
  ######################################
  
  dir.create(paste0('output/', n, '/validation'))
  dir.create(paste0('output/', n, '/validation/output'))
  dir.create(paste0('output/', n, '/validation/figures'))
  dir.create(paste0('output/', n, '/validation/rdata'))
  
  #load(file='./data/sampleSelectionVal.RData')
  selection.val <- colnames(d.raw)[!(colnames(d.raw) %in% selection)]
  save(selection.val, file=paste0('output/', n, '/validation/rdata/sampleSelectionVal.RData'))
  
  d.validation <- d.raw[,selection.val]
  samples <- colnames(d.validation)
  ids <- rownames(d.validation)
  group.data.validation <- group.data[selection.val,]
  d.validation.adj <- d.adj[,selection.val]
  colnames(d.validation.adj) <- samples
  
  d.validation.cs <- d.cs[,selection.val]
  
  d.validation.cs.WT <- d.validation.cs[,rownames(group.data.validation[group.data.validation[,hit] == "WT",])]
  d.validation.cs.SNV <- d.validation.cs[,rownames(group.data.validation[group.data.validation[,hit] == "SNV",])]
  #
  ##dr.WT <- dist(1-cor(t(d.validation.cs.WT)))
  ##hr.WT <- hclust(dr.WT)
  ##dc.WT <- dist(1-cor(d.validation.cs.WT))
  ##hc.WT <- hclust(dc.WT)
  ##
  ##dr.SNV <- dist(1-cor(t(d.validation.cs.SNV)))
  ##hr.SNV <- hclust(dr.SNV)
  ##dc.SNV <- dist(1-cor(d.validation.cs.SNV))
  ##hc.SNV <- hclust(dc.SNV)
  #
  colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
  
  png(file=paste0('output/', n, '/validation/figures/heatmapSNVSubset.png'), width=1920, height=1020)
  pheatmap(d.validation.cs.SNV,
           legend=TRUE,
           color=my_palette,
           breaks=colors,
           show_rownames=FALSE,
           show_colnames=FALSE
  )
  dev.off()
  
  png(file=paste0('output/', n, '/validation/figures/heatmapWTSubset.png'), width=1920, height=1020)
  pheatmap(d.validation.cs.WT,
           legend=TRUE,
           color=my_palette,
           breaks=colors,
           show_rownames=FALSE,
           show_colnames=FALSE
  )
  dev.off()
  ##
  #d.validation.adj.WT <- d.validation.adj[,rownames(group.data.validation[group.data.validation[,hit] == "WT",])]
  #d.validation.adj.SNV <- d.validation.adj[,rownames(group.data.validation[group.data.validation[,hit] == "SNV",])]
  
  #png(file=paste0('output/', n, '/validation/figures/sampleClusterVSTWT.png'), width=1920, height=1020)
  #plotClusterTreeSamples(t(d.validation.adj.WT), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #dev.off()
  #
  #png(file=paste0('output/', n, '/validation/figures/sampleClusterVSTSNV.png'), width=1920, height=1020)
  #plotClusterTreeSamples(t(d.validation.adj.SNV), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #dev.off()
  #
  #png(file=paste0('output/', n, '/validation/figures/heatmapVSTWT.png'), width=1920, height=1020)
  #pheatmap(d.validation.adj.WT, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE)
  #dev.off()
  #
  #png(file=paste0('output/', n, '/validation/figures/heatmapVSTSNV.png'), width=1920, height=1020)
  #pheatmap(d.validation.adj.SNV, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE)
  #dev.off()
  
  groups <- group.data.validation[,hit]
  # Save input data
  save(d.validation, group.data.validation, groups, file=(paste0('output/', n, '/validation/rdata/input_data.RData')))
  time2 <- Sys.time()
  print("Processing data:")
  difftime(time2, time1, units="secs")
  
  group <- as.matrix(group.data.validation[,hit])
  rownames(group) <- rownames(group.data.validation)
  colnames(group) <- c(hit)
  
  ###PLOTS
  png(file=paste0('output/', n, '/validation/figures/heatmap.png'), width=1920, height=1020)
  pheatmap(d.validation, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, show_colnames = FALSE, annotation_col=as.data.frame(group), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  colnames(d.validation) <- groups
  
  png(file=paste0('output/', n, '/validation/figures/sampleCluster.png'), width=1920, height=1020)
  plotClusterTreeSamples(t(d.validation), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  sampleDists <- dist(t(d.validation))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(d.validation)
  colnames(sampleDistMatrix) <- colnames(d.validation)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  png(file=paste0('output/', n, '/validation/figures/sampleDistances.png'), width=1920, height=1020)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/meanExpression.png'), width=1920, height=1020)
  barplot(apply(t(d.validation),1,mean, na.rm=T),
          xlab = "Sample", ylab = "Mean expression",
          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/variance.png'), width=1920, height=1020)
  meanSdPlot(d.validation, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  #png(file=paste0('figures/heatmapNorm.png'), width=1920, height=1020)
  #pheatmap(d.validation.norm, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group))
  #dev.off()
  #
  #colnames(d.validation.norm) <- groups
  #
  #sampleDists <- dist(t(d.validation.norm))
  #sampleDistMatrix <- as.matrix(sampleDists)
  #rownames(sampleDistMatrix) <- colnames(d.validation.norm)
  #colnames(sampleDistMatrix) <- colnames(d.validation.norm)
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
  #barplot(apply(t(d.validation.norm),1,mean, na.rm=T),
  #        xlab = "Sample", ylab = "Mean expression",
  #        main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #dev.off()
  
  png(file=paste0('output/', n, '/validation/figures/heatmapVST.png'), width=1920, height=1020)
  pheatmap(d.validation.adj, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group))
  dev.off()
  
  group[,1] <- as.numeric(factor(group[,1]))
  colnames(d.validation.adj) <- groups
  
  png(file=paste0('output/', n, '/validation/figures/sampleClusterVST.png'), width=1920, height=1020)
  plotClusterTreeSamples(t(d.validation.adj), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/meanExpressionVST.png'), width=1920, height=1020)
  barplot(apply(t(d.validation.adj),1,mean, na.rm=T),
          xlab = "Sample", ylab = "Mean expression",
          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/varianceVST.png'), width=1920, height=1020)
  meanSdPlot(d.validation.adj, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  sampleDists <- dist(t(d.validation.adj))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(d.validation.adj)
  colnames(sampleDistMatrix) <- colnames(d.validation.adj)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  png(file=paste0('output/', n, '/validation/figures/sampleDistancesVST.png'), width=1920, height=1020)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  dev.off()
  
  
  #######################################
  ###         Preprocessing           ###
  #######################################
  #load(file='rdata/input_data.rdata')
  dir.create(paste0('output/', n, '/validation/figures/differential expression'))
  #d.validation.norm <- normalize.sample(d.validation)
  #d.validation.cs <- normalize.cs(d.validation.norm)
  #length <- length(groups)
  ## unique colnames in d.validation.cs are necessary for clustering
  #colnames(d.validation.cs) <- paste(groups, 1:length, sep='_')
  ## save input data into a file
  #save(d.validation.norm, d.validation.cs, ids, groups, file='rdata/normalized_data.rdata')
  
  #DESeq2
  time1 <- Sys.time()
  colnames(d.validation) <- samples
  dds <- DESeqDataSetFromMatrix(countData = d.validation, colData = group, design = formula(paste("~",hit)))
  dds
  
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  png(file=paste0('output/', n, '/validation/figures/differential expression/PCA.png'), width=1920, height=1020)
  plotPCA(vsd, intgroup=c(hit))
  dev.off()
  
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
  
  colnames(d.validation) <- groups
  
  d.validation.summary <- data.frame(unlist(res$log2FoldChange), unlist(res$padj), row.names = rownames(res))
  colnames(d.validation.summary) <- c("log2FC", "Pvalue")
  d.validation.summary <- na.omit(d.validation.summary)
  d.validation.adj <- d.validation.adj[rownames(d.validation.adj) %in% rownames(d.validation.summary),]
  ids <- rownames(d.validation.adj)
  #d.validation.norm <- d.validation.norm[rownames(d.validation.norm) %in% rownames(d.validation.summary),]
  #ids <- rownames(d.validation.norm)
  time2<- Sys.time()
  print("DESeq2:")
  difftime(time2, time1, units="secs")
  
  ###DESEQ PLOTS
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/differential expression/MA.png'), width=1920, height=1020)
  plotMA(res, ylim=c(-2,2), cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  png(file=paste0('output/', n, '/validation/figures/differential expression/counts.png'), width=1920, height=1020)
  plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup=hit, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/differential expression/dispersion.png'), width=1920, height=1020)
  plotDispEsts(dds, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/differential expression/rejections.png'), width=1920, height=1020)
  plot(metadata(res)$filterNumRej, 
       type="b", ylab="number of rejections",
       xlab="quantiles of filter", cex.lab = 2, cex.axis = 2, cex.main = 2)
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
  dev.off()
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/differential expression/outliers.png'), width=1920, height=1020)
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, names = groups)
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  W <- res$stat
  maxCooks <- apply(assays(dds)[["cooks"]],1,max)
  idx <- !is.na(W)
  png(file=paste0('output/', n, '/validation/figures/differential expression/Wald.png'), width=1920, height=1020)
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
  png(file=paste0('output/', n, '/validation/figures/differential expression/pass.png'), width=1920, height=1020)
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
  # bb <- betaBinomial(d.validation.norm, ids, groups, 'WT', 'SNV', 'two.sided')
  # print(Sys.time())
  # d.validation.sign <- data.frame(log2FC=bb$table[ids,]$Log2ratio, Pvalue=bb$table[ids,]$Pvalue)
  # rownames(d.validation.sign) <- ids
  # # Merge stats about P-value thresholds and number of significant proteins
  # sign.table <- cbind(bb$FDRs)
  # d.validation.summary <- d.validation.sign
  # # Write output into files
  # write.table(bb$table, file='output/differential_expression/sign_test.tsv', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
  # write.table(sign.table, file='output/differential_expression/FDR_pval_summary.tsv', sep='\t', quote=FALSE)
  # write.table(d.validation.summary, file='output/differential_expression/log2_pval_uniprot.tsv', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
  # save(bb, sign.table, d.validation.summary, file='rdata/differential_expression.RData')
  
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
  # # Column group annotations (match to the column names of d.validation.cs)
  # ann_col <- data.frame(group=as.factor(groups))
  # rownames(ann_col) <- paste(groups, 1:length(groups), sep='_')
  # 
  # 
  # # Significance for each of the individual comparisons
  # # TAU vs control
  # thresholds <- sign.table[,'SNV_vs_WT_Pvalue']
  # names(thresholds) <- rownames(sign.table)
  # sign.SNVWT <- add.validation.sign.labels(bb, thresholds)
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
  # h.all <- Clustering(d.validation.cs, ann_col=ann_col, anncol=ann_colors, f='figures/clustering/clustering_all.png')
  # h.SNVWT <- Clustering(d.validation.cs, ann_col=ann_col, ann_row=ann_row4, anncol=ann_colors, f='figures/clustering/clustering_SNVWT_all.png')
  # h.SNVWT.sign <- Clustering(d.validation.cs[sign.SNVWT[rownames(d.validation.cs)]<5,], ann_col=ann_col, ann_row=data.frame(SNV_WT=ann_row[sign.SNVWT<5,], row.names=rownames(ann_row)[sign.SNVWT<5]),
  #                             anncol=ann_colors, f='figures/clustering/clustering_SNVWT_sign.png')
  # h.top50.sign <- Clustering(d.validation.cs[ids.sign.all,], ann_col=ann_col, ann_row=ann_row4,
  #                             anncol=ann_colors, f='figures/clustering/clustering_top50_sign.png')
  
  #######################################
  ###   WGCNA coexpression analysis   ###
  #######################################
  #load(file='rdata/input_data.RData')
  #load(file='rdata/normalized_data.RData')
  #load(file='rdata/differential_expression.RData')
  #source('src/coexpression_analysis_prep.R')
  dir.create(paste0('output/', n, '/validation/figures/coexpression'))
  dir.create(paste0('output/', n, '/validation/output/coexpression'))
  
  time1 <- Sys.time()
  d.validation.log2vals <- as.data.frame(d.validation.summary[, c('log2FC')], row.names=row.names(d.validation.summary))
  colnames(d.validation.log2vals) <- c('log2FC')
  coexpression <- coexpression.analysis(t(d.validation.adj), d.validation.log2vals, paste0('output/', n, '/validation/output/coexpression'), paste0('output/', n, '/validation/figures/coexpression'), n, 'validation')
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
  write.table(modules, file=paste0('output/', n, '/validation/output/coexpression/modules.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  #gns<-mapIds(org.Hs.eg.db, keys = ids, column = "ENTREZID", keytype = "SYMBOL")
  #GOenr.net <- GO.terms.modules(wgcna.net, ids, log2vals, gns, 'figures/coexpression', 'output/coexpression')
  ids <- ids[moduleColors != "grey"]
  log2vals <- d.validation.summary[ids, "log2FC"]  
  pvals <- d.validation.summary[ids, "Pvalue"]
  names(log2vals) <- ids
  names(pvals) <- ids
  # Save data structures
  #(wgcna.net, moduleColors, file='rdata/coexpression.RData')
  
  ### PLOTS
  TOM <- TOMsimilarityFromExpr(t(d.validation.adj[ids,]), power=power, TOMType="unsigned")
  diss2 <- 1-TOM
  hier2 <- hclust(as.dist(diss2), method="average")
  colorDynamicTOM <- labels2colors(cutreeDynamic(hier2,method="tree"))
  diag(diss2) = NA;
  moduleColors2 <- moduleColors[moduleColors != "grey"]
  
  par(mar = c(1, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/coexpression/TOMheatmap.png'), width=1920, height=1020)
  TOMplot(diss2^4, hier2, main = "TOM heatmap plot, module genes", terrainColors = FALSE, colors = moduleColors2, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  png(file=paste0('output/', n, '/validation/figures/coexpression/networkHeatmap.png'), width=1920, height=1020)
  plotNetworkHeatmap(
    t(d.validation.adj[ids,]),
    ids,
    useTOM = FALSE,
    power = power,
    networkType = "unsigned",
    main = "Heatmap of the network")
  dev.off()
  
  par(mar = c(2, 1, 1, 1))
  png(file=paste0('output/', n, '/validation/figures/coexpression/moduleSignificance.png'), width=1920, height=1020)
  plotModuleSignificance(
    d.validation.summary[ids, "Pvalue"],
    moduleColors2,
    boxplot = FALSE,
    main = "Gene significance across modules,",
    ylab = "Gene Significance", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  datME <- moduleEigengenes(t(d.validation.adj[ids,]),colorDynamicTOM)$eigengenes
  dissimME <- (1-t(cor(datME, use="p")))/2
  hclustdatME <- hclust(as.dist(dissimME), method="average" )
  # Plot the eigengene dendrogram
  png(file=paste0('output/', n, '/validation/figures/coexpression/moduleEigengenes.png'), width=1920, height=1020)
  par(mfrow=c(1,1))
  plot(hclustdatME, main="Clustering tree based of the module eigengenes", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  #png(file=paste0('figures/coexpression/MEpairs.png'), width=1920, height=1020)
  #MEpairs <- plotMEpairs(
  #  datME,
  #  main = "Relationship between module eigengenes",
  #  clusterMEs = TRUE)
  #dev.off()
  
  GS1 <- as.numeric(cor(group,t(d.validation.adj[ids,])))
  GeneSignificance <- abs(GS1)
  ModuleSignificance <- tapply(GeneSignificance, colorDynamicTOM, mean, na.rm=T)
  png(file=paste0('output/', n, '/validation/figures/coexpression/moduleSignificance2.png'), width=1920, height=1020)
  plotModuleSignificance(GeneSignificance,colorDynamicTOM, cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  disableWGCNAThreads()
  
  MEs <- orderMEs(datME)
  modNames <- substring(names(MEs), 3)
  geneModuleMembership <- as.data.frame(cor(t(d.validation.adj[ids,]), MEs, use = "p"));
  names(geneModuleMembership) <- paste("MM", modNames, sep="");
  GS1 <- as.data.frame(GS1)
  names(GS1) = paste("GS.", names(groups), sep="")
  
  corPVal <- corAndPvalue(t(d.validation.adj[ids,]), use="pairwise.complete.obs")
  
  
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
  dir.create(paste0('output/', n, '/validation/figures/hotnet'))
  dir.create(paste0('output/', n, '/validation/output//hotnet'))
  dir.create(paste0('output/', n, '/validation/output/hotnet/HotNet_input'))
  
  corN <- corPVal[["cor"]]
  corN <- corN[corN < 0]
  corP <- corPVal[["cor"]]
  corP <- corP[corP > 0]
  
  probs <- 0.9
  cutOff <- quantile(as.vector(TOM), probs=probs)
  print(paste("Quantile at ", probs, " :", cutOff))
  
  cutOffN <- quantile(as.vector(corN), probs=0.1)
  cutOffP <- quantile(as.vector(corP), probs=0.9)
  print(paste("90th negative quantile: ", cutOffN))
  print(paste("90th positive quantile: ", cutOffP))
  
  hotnetColors <- c()
  
  # Get table with interactions for each module
  for (module in modules){
    if (module == "grey"){
      next
    } 
    else {
      print(module)
      par(mar = c(1, 1, 1, 1))
      png(file=paste0('output/', n, '/validation/figures/coexpression/', module, '_matrix.png'), width=1920, height=1020)
      plotMat(t(scale(t(d.validation.adj[colorDynamicTOM==module,]))),rlabels=T,
              clabels=T,rcols=module,
              title=module, cex.lab = 2, cex.axis = 2, cex.main = 2)
      dev.off()
      
      ME=datME[, paste("ME",module, sep="")]
      par(mar = c(2, 1, 1, 1))
      png(file=paste0('output/', n, '/validation/figures/coexpression/', module, '_ME.png'), width=1920, height=1020)
      barplot(ME, col=module, main="", cex.main=2,
              ylab="eigengene expression",xlab="array sample", cex.lab = 2, cex.axis = 2, cex.main = 2)
      dev.off()
      # Select proteins in module
      inModule <- moduleColors2 == module
      sum(inModule)
      TOMmodule <- TOM[inModule, inModule]
      corModule <- corPVal[["cor"]][inModule, inModule]
      #Pmodule <- corPVal[["p"]][inModule, inModule]
      idsModule <- ids[inModule]
      write.table(idsModule, file=paste0('output/', n, '/validation/output/coexpression/', module, '.tsv'), sep='\t', row.names=FALSE, col.names=FALSE)
      log2Module <- log2vals[inModule]
      PModule <- pvals[inModule]
      
      geneNames <- mapIds(org.Hs.eg.db, keys = idsModule, column = "ENTREZID", keytype = "SYMBOL")
      goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
      if (!(is.null(goBP)) && nrow(goBP) > 0){
        print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust < 0.1,])))
        png(file=paste0('output/', n, '/validation/figures/coexpression/GO_', module, '.png'), width=1920, height=1020)
        print(dotplot(goBP, showCategory=10))
        dev.off()
      }
      
      par(mfrow = c(1,1));
      column <- match(module, modNames);
      png(file=paste0('output/', n, '/validation/figures/coexpression/', module, '_membershipVSsignficance.png'), width=1920, height=1020)
      verboseScatterplot(abs(geneModuleMembership[inModule, column]),
                         abs(GS1[inModule, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
      
      node.frame <- data.frame(Symbol=idsModule, LogFC=log2Module, Pvalue=PModule)
      rownames(node.frame) <- 1:nrow(node.frame)
      #node.frame$InvPvalue <- 1 - PModuleProteins
      #write.table(node.frame, file=paste0('output/hotnet/nodes_', module, '.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
      # write.table(node.frame[, c('Symbol', 'InvPvalue')], file=paste0('output/hotnet/HotNet_input/g2s_Pval_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote     =FALSE)
      names(idsModule) <- idsModule
      
      write.table(node.frame[, c('Symbol', 'LogFC')], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
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
      #result = foreach(n = 1:(nrow(TOMmodule)-1), j = (n+1):nrow(TOMmodule)) %dopar% {
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
          if (corModule[i,j] >= 0.5 || corModule[i,j] <= -0.5){
            #   # Add edge to list
            result1[nrow(result1)+1,] <- c(idsModule[i], idsModule[j])
            result2[nrow(result2)+1,] <- c(i,j)
            #   # Add node to node list
            if (!i %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(i, idsModule[i])
              result4[nrow(result4)+1,] <- c(idsModule[i], log2Module[i]*(1-PModule[i]))
            } 
            if (!j %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(j, idsModule[j])
              result4[nrow(result4)+1,] <- c(idsModule[j], log2Module[j]*(1-PModule[j]))
            }
          }
        }
        return(list(result1, result2, result3, result4))
      }
      time2 <- Sys.time()
      print("Slicing input:")
      print(difftime(time2, time1, units="secs"))
      time1 <- Sys.time()
      write.table(result[[1]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/name_edges_expression_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[2]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/edge_list_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[3]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/i2g_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[4]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      time2 <- Sys.time()
      print("Writing data:")
      print(difftime(time2, time1, units="secs"))
      if (file.info(paste('output/', n, '/validation/output/hotnet/HotNet_input/i2g_', module, '.tsv', sep=''))$size != 0){
        hotnetColors <- append(hotnetColors, module)
      }
    }
  }
  
  
  ######### Run hierarchical hotnet to obtain the most significant submodule within each of the identified modules
  ######### Note that the HotNet package was written in Python and needs to be intstalled separately on your machine
  time1 <- Sys.time()
  system(paste('bash src/run_hierarchicalHotnet_modules.sh "', paste(modules, collapse=' '), '" ', cutOff, n, 'validation'))
  time2 <- Sys.time()
  print("Hierarchical HotNet:")
  difftime(time2, time1, units="secs")
  
  time1 <- Sys.time()
  dir.create(paste0('output/', n, '/validation/figures/GO'))
  dir.create(paste0('output/', n, '/validation/output/GO'))
  
  BPterms <- c()
  MFterms <- c()
  genes <- c()
  subnetworks <- c()
  hotnetSubnetworks <- c()
  
  #GO
  for (module in hotnetColors){
    if (file.exists(file=paste0('output/', n, '/validation/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'))){
      p <- read.table(file=paste0('output/', n, '/validation/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'), sep='\t', header= FALSE, comment.char="")[6, "V1"]
      p <- as.numeric(sub(".*: ", "", p))
      if (file.exists(paste('output/', n, '/validation/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''))){
        if (file.info(paste('output/', n, '/validation/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''))$size == 0){
          next
        }
      }
      hotnetSubnetworks <- append(hotnetSubnetworks, module)
      if (p <= 0.1){
        #Read in genes
        print(module)
        subnetworks <- append(subnetworks, module)
        inSubnetwork <- as.vector(t(read.csv(paste('output/', n, '/validation/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv', sep=''), sep='\t', header = FALSE)[1,]))
        print(paste("In subnetwork: ", inSubnetwork))
        genes <- append(genes, inSubnetwork)
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
            print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust <= 0.1,])))
            write.table(goBP@result, file=paste0('output/', n, '/validation/output/GO/', module, 'BP.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
            BPterms <- append(BPterms, rownames(goBP@result[goBP@result$p.adjust <= 0.1,]))
            png(file=paste0('output/', n, '/validation/figures/GO/GOBPBar_', module, '.png'), width=1920, height=1020)
            print(barplot(goBP, showCategory=10, cex.lab = 2, cex.axis = 2, cex.main = 2))
            dev.off()
            png(file=paste0('output/', n, '/validation/figures/GO/GOBPDot_', module, '.png'), width=1920, height=1020)
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
            print(paste("Number of enriched molecular functions: ", nrow(goMF[goMF$p.adjust <= 0.1,])))
            write.table(goMF@result, file=paste0('output/', n, '/validation/output/GO/', module, 'MF.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
            MFterms <- append(MFterms, rownames(goMF@result[goMF@result$p.adjust <= 0.1,]))
            png(file=paste0('output/', n, '/validation/figures/GO/GOMFBar_', module, '.png'), width=1920, height=1020)
            print(barplot(goMF, showCategory=10, cex.lab = 2, cex.axis = 2, cex.main = 2))
            dev.off()
            png(file=paste0('output/', n, '/validation/figures/GO/GOMFDot_', module, '.png'), width=1920, height=1020)
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
  
  write.table(BPterms, file=paste0('output/', n, '/validation/output/GO/BP.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(MFterms, file=paste0('output/', n, '/validation/output/GO/MF.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(genes, file=paste0('output/', n, '/validation/output/GO/genes.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(subnetworks, file=paste0('output/', n, '/validation/output/GO/subnetworks.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(hotnetSubnetworks, file=paste0('output/', n, '/validation/output/GO/hotnetSubnetworks.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  save(BPterms, file=paste0('output/', n, '/validation/rdata/BP.RData'))
  save(MFterms, file=paste0('output/', n, '/validation/rdata/MF.RData'))
  save(genes, file=paste0('output/', n, '/validation/rdata/genes.RData'))
  save(subnetworks, file=paste0('output/', n, '/validation/rdata/subnetworks.RData'))
  
  time2 <- Sys.time()
  print("GO:")
  difftime(time2, time1, units="secs")
  endTime <- Sys.time()
  print("Total time:")
  difftime(endTime, startTime, units="secs")
  
  dir.create(paste0('output/', n, '/validation/output/similarity'))
  
  similarity <- similarity.analysis(paste0('output/', n, '/training/output'), paste0('output/', n, '/validation/output'), paste0('output/', n, '/validation/output/similarity'))

print(paste('END OF RUN', n))
}

averageSimilarityModules1 <- mean(as.matrix(read.table(paste0('output/1/validation/output/similarity/moduleSimilarity.tsv'), sep='\t')))
averageSimilarityModules2 <- mean(as.matrix(read.table(paste0('output/2/validation/output/similarity/moduleSimilarity.tsv'), sep='\t')))
averageSimilarityModules3 <- mean(as.matrix(read.table(paste0('output/3/validation/output/similarity/moduleSimilarity.tsv'), sep='\t')))

averageSimilaritySubnetworks1 <- mean(as.matrix(read.table(paste0('output/1/validation/output/similarity/subnetworkSimilarity.tsv'), sep='\t')))
averageSimilaritySubnetworks2 <- mean(as.matrix(read.table(paste0('output/2/validation/output/similarity/subnetworkSimilarity.tsv'), sep='\t')))
averageSimilaritySubnetworks3 <- mean(as.matrix(read.table(paste0('output/3/validation/output/similarity/subnetworkSimilarity.tsv'), sep='\t')))

averageSimilaritySSubnetworks1 <- mean(as.matrix(read.table(paste0('output/1/validation/output/similarity/sSubnetworkSimilarity.tsv'), sep='\t')))
averageSimilaritySSubnetworks2 <- mean(as.matrix(read.table(paste0('output/2/validation/output/similarity/sSubnetworkSimilarity.tsv'), sep='\t')))
averageSimilaritySSubnetworks3 <- mean(as.matrix(read.table(paste0('output/3/validation/output/similarity/sSubnetworkSimilarity.tsv'), sep='\t')))

resultSimilarity1 <- read.table(paste0('output/1/validation/output/similarity/similarity.tsv'), sep='\t')
resultSimilarity2 <- read.table(paste0('output/2/validation/output/similarity/similarity.tsv'), sep='\t')
resultSimilarity3 <- read.table(paste0('output/3/validation/output/similarity/similarity.tsv'), sep='\t')

print(paste('Dissimilarity for modules across runs:', mean(c(averageSimilarityModules1, averageSimilarityModules2, averageSimilarityModules3))))
print(paste('Dissimilarity for subnetworks across runs:', mean(c(averageSimilaritySubnetworks1, averageSimilaritySubnetworks2, averageSimilaritySubnetworks3))))
print(paste('Dissimilarity for significant subnetworks across runs:', mean(c(averageSimilaritySSubnetworks1, averageSimilaritySSubnetworks2, averageSimilaritySSubnetworks3))))
print(paste('Similarity for significant genes in subnetworks across runs:', mean(resultSimilarity1['genes',1], resultSimilarity2['genes',1], resultSimilarity3['genes',1])))
print(paste('Similarity for significant BP terms in subnetworks across runs:', mean(resultSimilarity1['BP',1], resultSimilarity2['BP',1], resultSimilarity3['BP',1])))
print(paste('Similarity for significant MF terms in subnetworks across runs:', mean(resultSimilarity1['MF',1], resultSimilarity2['MF',1], resultSimilarity3['MF',1])))

######################################
##       External validation       ###
######################################
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#load(file='rdata/coexpression.RData')
#dir.create(paste0('output/', n, '/training/output/validation'))
#dir.create(paste0('output/', n, '/training/figures/validation'))
#dir.create(paste0('output/', n, '/training/rdata/validation'))
## 
## Load RIMOD data
#d.val <- read.table('./data/TCGA_rna_count_data.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)
#group.data.val <- read.table('./data/non_silent_mutation_profile_crc.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)
#d.norm.val <- normalize.sample(d.val)
#d.cs.val <- normalize.cs(d.norm.val)
#d.val.adj <- varianceStabilizingTransformation(d.val, blind = F)
##load(file='./data/sampleSelectionVal.RData')
#selection.val <- colnames(d.val)[!(colnames(d.val) %in% selection)]
#save(selection.val, file=paste0('output/', n, '/training/rdata/validation/sampleSelectionVal.RData'))
#d.val <- d.val[,selection.val]
#samples.val <- colnames(d.val)
#ids.val <- rownames(d.val)
#d.val.adj <- d.val.adj[,selection.val]
#group.data.val <- group.data.val[selection.val,]
#
#d.cs.val <- d.cs.val[,selection.val]
#
#d.cs.val.WT <- d.cs.val[,rownames(group.data.val[group.data.val[,hit] == "WT",])]
#d.cs.val.SNV <- d.cs.val[,rownames(group.data.val[group.data.val[,hit] == "SNV",])]
#
#colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
#my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
#
#png(file=paste0('output/', n, '/training/figures/validation/heatmapSNVSubset.png'), width=1920, height=1020)
#pheatmap(d.cs.val.SNV,
#                legend=TRUE,
#                color=my_palette,
#                breaks=colors,
#                show_rownames=FALSE,
#                show_colnames=FALSE
#)
#dev.off()
#
#png(file=paste0('output/', n, '/training/figures/validation/heatmapWTSubset.png'), width=1920, height=1020)
#pheatmap(d.cs.val.WT,
#                legend=TRUE,
#                color=my_palette,
#                breaks=colors,
#                show_rownames=FALSE,
#                show_colnames=FALSE
#)
#dev.off()
#
#group.val <- as.matrix(group.data.val[,hit])
#rownames(group.val) <- rownames(group.data.val)
#colnames(group.val) <- c(hit)
#
#group.val[,1] <- as.numeric(factor(group.val[,1]))
#colnames(d.val.adj) <- groups.val
#
#dds.val <- DESeqDataSetFromMatrix(countData = d.val, colData = group.val, design = formula(paste("~",hit)))
#dds.val
#dds.val <- DESeq(dds.val)
#res.val <- results(dds.val)
#
#colnames(d.val) <- groups.val
#
#d.val.summary <- data.frame(unlist(res.val$log2FoldChange), unlist(res.val$padj), row.names = rownames(res.val))
#colnames(d.val.summary) <- c("log2FC", "Pvalue")
#d.val.summary <- na.omit(d.val.summary)
#d.val.adj <- d.val.adj[rownames(d.val.adj) %in% rownames(d.val.summary),]
#ids.val <- rownames(d.val.adj)
#
#log2vals.val <- d.val.summary[ids.val, "log2FC"]
#pvals.val <- d.val.summary[ids.val, "Pvalue"]
#names(log2vals.val) <- ids.val
#names(pvals.val) <- ids.val
#
## write.table(bb.1$FDRs, 'output/validation/sign_test_val_frontal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
## write.table(bb.2$FDRs, 'output/validation/sign_test_val_temporal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
#
## 
## # Get module IDs
## moduleColors <- labels2colors(wgcna.net$colors)
## modules <- levels(as.factor(moduleColors))
## names(moduleColors) <- ids.TAUcon
## 
## # Test against frontal dataset
## # Filter out proteins that are found in one of the two sets (Note that proteins in FC2/P2 that are not in FC1/P1 are
## # automatically ignored in fraction.correct, so there is no need to filter that dataset)
#ids.both <- ids %in% ids.val
#log.toTest<- log2vals[ids.both]
#P.toTest <- pvals[ids.both]
## 
## # Permutation test on all proteins (frontal)
#outfolder <- paste0('output/', n, '/training/output/validation/')
#figfolder <- paste0('output/', n, '/training/figures/validation/')
#dir.create(outfolder, showWarnings=FALSE)
#dir.create(figfolder, showWarnings=FALSE)
##perm <- permutation.test(log.toTest, log2vals.val, P.toTest, pvals.val, outfolder, figfolder)
## save(perm.allFrn, file='rdata/validation/permutation_all_frontal.RData')
## 
## # Permutation test per module
#moduleColorsBoth <- moduleColors[ids.both]
#for(m in modules){
#  if (m == 'grey'){
#    next
#  }
#  else{
#   print(m)
#   inModule <- moduleColorsBoth == m
#   m.ids <- names(moduleColorsBoth[inModule])
#   #log.TAUcon.module <- log.TAUcon.toTest[inModule]
#   #P.TAUcon.module <- P.TAUcon.toTest[inModule]
#   outfolder <- paste0('output/', n, '/training/output/validation/module_', m)
#   figfolder <- paste0('output/', n, '/training/figures/validation/module_', m)
#   dir.create(outfolder, showWarnings=FALSE)
#   dir.create(figfolder, showWarnings=FALSE)
#   perm.module <- permutation.test(log2vals, log2vals.val, pvals, pvals.val, outfolder, figfolder, m.ids=m.ids)
#   #save(perm.module, file=paste0('rdata/validation/permutation_', m, '_frontal.RData'))
#   
#   
#   ids.hotnet <- as.vector(t(read.csv(paste0(hit, '/output/hotnet/HotNet_results/consensus_nodes_log2_', m, '.tsv', sep=''), sep='\t', header = FALSE)[1,]))
#   
#   ids.both <- ids %in% ids.val & ids %in% ids.hotnet
#   inModule <- moduleColors[ids.both] == m
#   m.ids <- names(moduleColors[ids.both][inModule])
#   log.module <- log2vals[ids.both]
#   P.module <- pvals[ids.both]
#   outfolder <- paste0('output/', n, '/training/figures/validation/Module_', m, '_hotnet/')
#   dir.create(outfolder, showWarnings=FALSE)
#   perm.module.hotnet <- permutation.test(log2vals, log2vals.val, pvals, pvals.val, outfolder, m.ids=m.ids)
#   #save(perm.module.hotnet, file=paste0('Output/Permutation_', m, '_frontal4_hotnet.RData'))
#  }
#}
stopCluster(cluster)
# 
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
# for (n in 0:17){
#   color <- labels2colors(n)
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