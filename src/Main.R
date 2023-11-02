argv = commandArgs(trailingOnly=TRUE)
startTime <- Sys.time()
#######################################
#####          Packages           #####
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
library(clusterProfiler)
library(ggplot2)
library(ggprism)
library(ggrepel)
library(ggraph)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(enrichplot)
library(data.table)
library(reshape2)
library(tidyr)
library(tidygraph)
library(msigdbr)
library(dplyr, exclude = c("collapse", "combine", "count", "desc", "filter", "first", "lag", "rename", "select", "slice"))
library(WGCNA)
options(stringsAsFactors=FALSE)
source(paste0(argv[1],'/src/coexpression_analysis.R'))
source(paste0(argv[1],'/src/similarity.R'))

#######################################
######          Setup            ######
#######################################
setwd(argv[1])
dir.create('output')

totalCores = detectCores()
cluster <- makeCluster(31)
registerDoParallel(cluster)

#hits <- read.csv('./data/hitlist_snv_screen_coadread_tcga.csv', header = TRUE, sep=';', row.names=1)

hit <- argv[2]

time1 <- Sys.time()
d <- data.frame(fread('./data/TCGA_rna_count_data.csv', header=TRUE), check.names = FALSE, row.names = 1)
group.data <- data.frame(fread('./data/non_silent_mutation_profile_crc.csv', header=TRUE), row.names=1)
if (argv[4] != ""){
  mut <- read.csv('./data/TCGA_mutation_definition_data.csv', check.names=FALSE)
  variant <- subset(mut, grepl(argv[4], HGVSp_Short)  &  Hugo_Symbol == hit)
  group.data[substr(variant$Tumor_Sample_Barcode,1,16), hit] <- "SNV"
  group.data[!(rownames(group.data) %in% substr(variant$Tumor_Sample_Barcode,1,16)), hit] <- "WT"
  rm(variant)
  rm(mut)
}
time2 <- Sys.time()
print(paste("Loading data:", difftime(time2, time1, units="secs")))

time1 <- Sys.time()
d <- d[,colnames(d) %in% rownames(group.data)]
group.data <- group.data[rownames(group.data) %in% colnames(d),]
groups <- group.data[,hit]
new_order <- sort(colnames(d))
d <- d[, new_order]
new_order_rows <- sort(rownames(d))
d <- d[new_order_rows,]
d.raw <- apply(d, c(1,2), as.numeric)
d.raw <- d.raw[rowSums(d.raw) >= 10,]
samples <- colnames(d.raw)

group <- as.matrix(group.data[,hit])
rownames(group) <- rownames(group.data)
colnames(group) <- c(hit)

d.adj <- varianceStabilizingTransformation(d.raw, blind = F)

pca <- prcomp(t(d.adj))
d.f <- cbind(group.data, pca$x)

png(file='output/PCA.png', width=1920, height=1020)
print(ggplot(d.f) + 
        geom_point(aes(x=PC1, y=PC2,colour = KRAS),size = 3) +
        theme_prism(base_size = 35) +
        theme(legend.position = c(0.9, 0.8),
              legend.key.height = unit(12, "pt"),
              legend.title = element_text(size=35),
              legend.text = element_text(size=35)))
dev.off()

colnames(d.adj) <- groups

sampleTree <- hclust(dist(t(d.adj)), method = "ave")

png(file='output/clusterSamplesVST.png', width=1920, height=1020)
plot(sampleTree, main = "Sample clustering to detect outliers", cex = 0.5, lwd = 2, axes=FALSE, labels=FALSE)
abline(h = 225, col = "red",lwd = 2);
axis(side = 2, at = seq(0, 400, 100),
     labels = TRUE, lwd = 2)
clust <- cutreeStatic(sampleTree, cutHeight = 225, minSize = 10)
table(clust)
keepSamples <- (clust==1)
dev.off()

colnames(d.adj) <- samples

d.raw <- d.raw[,keepSamples]
group.data <- group.data[keepSamples,]

count_sd_data <- data.frame(
  count = rowMeans(d.raw),
  deviation = apply(d.raw, 1, sd)
)

png(file='output/variance.png', width=1920, height=1020)
print(ggplot(count_sd_data, aes(x = count, y = deviation)) +
        geom_point(size=2) +
        geom_smooth(method = "lm", se = FALSE, size=2) + 
        labs(x = "Mean of counts", y = "Standard Deviation") +
        ggtitle("Standard Deviation vs. Mean Plot") +
        theme_prism(border=TRUE,base_size = 30))
dev.off()

count_sd_data <- data.frame(
  count = rowMeans(d.adj),
  deviation = apply(d.adj, 1, sd)
)

png(file='output/varianceVST.png', width=1920, height=1020)
print(ggplot(count_sd_data, aes(x = count, y = deviation)) +
        geom_point(size=2) +
        geom_smooth(method = "lm", se = FALSE,size=2) +
        labs(x = "Mean of counts", y = "Standard Deviation") +
        ggtitle("Standard Deviation vs. Mean Plot") +
        theme_prism(border=TRUE,base_size = 30))
dev.off()

means_data <- data.frame(
  means = colMeans(d.raw),
  samples = colnames(d.raw)
)

png(file='output/meanExpression.png', width=1920, height=1020)
print(ggplot(means_data, aes(x = samples, y = means, fill=samples)) +
        geom_bar(stat = "identity", colour='black',width = 1) +
        labs(
          title = "Mean Counts per Sample",
          x = "Sample",
          y = "Mean Counts"
        ) +
        theme_prism(base_size = 30) +
        theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(expand = c(0, 0)) +
        guides(fill = FALSE))
dev.off()

group <- as.matrix(group.data[,hit])
rownames(group) <- rownames(group.data)
colnames(group) <- c(hit)

dd <- DESeqDataSetFromMatrix(countData = d.raw, colData = group, design = formula(paste("~",hit)))
dd.norm <- estimateSizeFactors(dd)
normalized_counts <- counts(dd.norm, normalized=TRUE)

means_data <- data.frame(
  means = colMeans(normalized_counts),
  samples = colnames(normalized_counts)
)

png(file='output/meanExpressionNorm.png', width=1920, height=1020)
print(ggplot(means_data, aes(x = samples, y = means, fill=samples)) +
        geom_bar(stat = "identity", colour='black',width = 1) +
        labs(
          title = "Mean Counts per Sample",
          x = "Sample",
          y = "Mean Counts"
        ) +
        theme_prism(base_size = 30) +
        theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(expand = c(0, 0)) +
        guides(fill = FALSE))
dev.off()
  
for (n in c("1","2","3")){
  #foreach(hit = hits, .packages = packages) %dopar% {
  dir.create(paste0('output/', n))
  dir.create(paste0('output/', n, '/training'))
  dir.create(paste0('output/', n, '/training/output'))
  dir.create(paste0('output/', n, '/training/figures'))
  dir.create(paste0('output/', n, '/training/rdata'))
  
  selection <- sample(colnames(d.raw), floor(ncol(d.raw)/2))
  selection <- sort(selection)
  save(selection, file=paste0('output/' , n, '/training/rdata/sampleSelection.RData'))
  
  d.training <- d.raw[,selection]
  samples <- colnames(d.training)
  ids <- rownames(d.training)
  group.data.training <- group.data[selection,]
  groups <- group.data.training[,hit]
  d.training.adj <- d.adj[,selection]
  colnames(d.training.adj) <- samples
  
  time2 <- Sys.time()
  print(paste("Processing data:", difftime(time2, time1, units="secs")))
  
  group <- as.matrix(group.data.training[,hit])
  rownames(group) <- rownames(group.data.training)
  colnames(group) <- c(hit)
  
  #  d.norm <- normalize.sample(d.raw)
  #  d.cs <- normalize.cs(d.norm)
  #
  ### PLOTS
  #
  #  d.cs <- d.cs[,keepSamples]
  #  d.cs.WT <- d.cs[,rownames(group.data[group.data[,hit] == "WT",])]
  #  d.cs.SNV <- d.cs[,rownames(group.data[group.data[,hit] == "SNV",])]
  #  
  #  colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
  #  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
  #  
  #  png(file=paste0('output/heatmap.png'), width=1920, height=1020)
  #    pheatmap(d.cs, color=my_palette, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, show_colnames = FALSE, annotation_col=as.data.frame(group), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #  
  #  png(file=paste0('output/heatmapSNV.png'), width=1920, height=1020)
  #    pheatmap(d.cs.SNV,
  #             legend=TRUE,
  #             color=my_palette,
  #             breaks=colors,
  #             show_rownames=FALSE,
  #             show_colnames=FALSE
  #    )
  #  dev.off()
  #  
  #  png(file=paste0('output/heatmapWT.png'), width=1920, height=1020)
  #    pheatmap(d.cs.WT,
  #             legend=TRUE,
  #             color=my_palette,
  #             breaks=colors,
  #             show_rownames=FALSE,
  #             show_colnames=FALSE
  #    )
  #  dev.off()
  #
  #  d.training.cs <- d.cs[,selection]
  #  d.training.cs.WT <- d.training.cs[,rownames(group.data.training[group.data.training[,hit] == "WT",])]
  #  d.training.cs.SNV <- d.training.cs[,rownames(group.data.training[group.data.training[,hit] == "SNV",])]
  #  
  #  png(file=paste0('output/', n, '/training/figures/heatmapSNVSubset.png'), width=1920, height=1020)
  #    pheatmap(d.training.cs.SNV,
  #             legend=TRUE,
  #             color=my_palette,
  #             breaks=colors,
  #             show_rownames=FALSE,
  #             show_colnames=FALSE
  #    )
  #  dev.off()
  #  
  #  png(file=paste0('output/', n, '/training/figures/heatmapWTSubset.png'), width=1920, height=1020)
  #    pheatmap(d.training.cs.WT,
  #             legend=TRUE,
  #             color=my_palette,
  #             breaks=colors,
  #             show_rownames=FALSE,
  #             show_colnames=FALSE
  #    )
  #  dev.off()   
  #
  #  png(file=paste0('output/', n, '/training/figures/heatmap.png'), width=1920, height=1020)
  #    pheatmap(d.training.cs, color=my_palette, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, show_colnames = FALSE, annotation_col=as.data.frame(group), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  
  colnames(d.training) <- groups
  
  #  png(file=paste0('output/', n, '/training/figures/sampleCluster.png'), width=1920, height=1020)
  #    plotClusterTreeSamples(t(d.training), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #  
  #  sampleDists <- dist(t(d.training))
  #  sampleDistMatrix <- as.matrix(sampleDists)
  #  rownames(sampleDistMatrix) <- colnames(d.training)
  #  colnames(sampleDistMatrix) <- colnames(d.training)
  #  png(file=paste0('output/', n, '/training/figures/sampleDistances.png'), width=1920, height=1020)
  #    pheatmap(sampleDistMatrix,
  #             clustering_distance_rows=sampleDists,
  #             clustering_distance_cols=sampleDists,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #  
  #  geneDists <- dist(d.training)
  #  geneDistMatrix <- as.matrix(geneDists)
  #  rownames(geneDistMatrix) <- rownames(d.training)
  #  colnames(geneDistMatrix) <- rownames(d.training)
  #  png(file=paste0('output/', n, '/training/figures/geneDistances.png'), width=1920, height=1020)
  #    pheatmap(geneDistMatrix,
  #             clustering_distance_rows=geneDists,
  #             clustering_distance_cols=geneDists,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #
  #  png(file=paste0('output/', n, '/training/figures/meanExpressionSelection.png'), width=1920, height=1020)
  #  barplot(apply(t(d.training),1,mean, na.rm=T),
  #          xlab = "Sample", ylab = "Mean expression",
  #          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #  
  #  png(file=paste0('output/', n, '/training/figures/heatmapVST.png'), width=1920, height=1020)
  #    pheatmap(d.training.adj, color=my_palette, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group))
  #  dev.off()
  
  colnames(d.training.adj) <- groups
  
  #  png(file=paste0('output/', n, '/training/figures/sampleClusterVST.png'), width=1920, height=1020)
  #    plotClusterTreeSamples(t(d.training.adj), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #  
  #  sampleDists.adj <- dist(t(d.training.adj))
  #  sampleDistMatrix.adj <- as.matrix(sampleDists.adj)
  #  rownames(sampleDistMatrix.adj) <- colnames(d.training.adj)
  #  colnames(sampleDistMatrix.adj) <- colnames(d.training.adj)
  #  png(file=paste0('output/', n, '/training/figures/sampleDistancesVST.png'), width=1920, height=1020)
  #    pheatmap(sampleDistMatrix.adj,
  #             clustering_distance_rows=sampleDists.adj,
  #             clustering_distance_cols=sampleDists.adj,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #  
  #  geneDists.adj <- dist(d.training.adj)
  #  geneDistMatrix.adj <- as.matrix(geneDists.adj)
  #  rownames(geneDistMatrix.adj) <- rownames(d.training.adj)
  #  colnames(geneDistMatrix.adj) <- rownames(d.training.adj)
  #  png(file=paste0('output/', n, '/training/figures/geneDistancesVST.png'), width=1920, height=1020)
  #    pheatmap(geneDistMatrix.adj,
  #             clustering_distance_rows=geneDists.adj,
  #             clustering_distance_cols=geneDists.adj,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #
  #  png(file=paste0('output/', n, '/training/figures/meanExpressionSelectionVST.png'), width=1920, height=1020)
  #  barplot(apply(t(d.training.adj),1,mean, na.rm=T),
  #          xlab = "Sample", ylab = "Mean expression",
  #          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  
  ### Save input data
  save(d.training, d.training.adj, group.data.training, group, file=(paste0('output/', n, '/training/rdata/input_data.RData')))
  
  #######################################
  #######         DESeq2           ######
  #######################################
  dir.create(paste0('output/', n, '/training/figures/differential expression'))
  dir.create(paste0('output/', n, '/training/output/differential expression'))
  
  time1 <- Sys.time()
  colnames(d.training) <- samples
  group[,1] <- factor(group[,1], levels = c("WT","SNV"))
  dds <- DESeqDataSetFromMatrix(countData = d.training, colData = group, design = formula(paste("~",hit)))
  dds
  dds <- DESeq(dds)
  res <- results(dds)
  res
  summary(res)
  
  colnames(d.training) <- groups
  
  d.training.summary <- data.frame(unlist(res$log2FoldChange), unlist(res$padj), row.names = rownames(res))
  colnames(d.training.summary) <- c("log2FC", "Pvalue")
  d.training.summary <- na.omit(d.training.summary)
  d.training.adj <- d.training.adj[rownames(d.training.adj) %in% rownames(d.training.summary),]
  ids <- rownames(d.training.adj)
  
  time2<- Sys.time()
  print(paste("DESeq2:", difftime(time2, time1, units="secs")))
  
  write.table(d.training.summary[,"log2FC", drop = FALSE], file=paste0('output/', n, '/training/output/differential expression/DE_LFC.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  ###DESEQ PLOTS
  d.training.summary2 <- d.training.summary %>% 
    mutate(
      Expression = case_when(log2FC >= 0.1 & Pvalue <= 0.1 ~ "Up-regulated",
                             log2FC <= -0.1 & Pvalue <= 0.1 ~ "Down-regulated",
                             TRUE ~ "Unchanged")
    ) %>%
    arrange(desc(log2FC))
  
  d.training.summary2 <- cbind(gene=rownames(d.training.summary2), d.training.summary2)
  
  png(file=paste0('output/', n, '/training/figures/differential expression/volcano.png'), width=1920, height=1020)
  print(ggplot(d.training.summary2, aes(log2FC, -log(Pvalue,10))) +
          geom_point(aes(color = Expression), size = 2) +
          xlab(expression("log"[2]*"FC")) + 
          ylab(expression("-log"[10]*"adjusted p-value")) +
          scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Unchanged" = "grey")) +
          scale_x_continuous(
            limits = c(-4, 4),
            breaks = seq(-4, 4, by = 0.5),
            guide = "prism_minor") +
          scale_y_continuous(
            limits = c(0, 15),
            breaks = seq(0, 15, by = 3)) +
          geom_label_repel(data=head(filter(d.training.summary2, Pvalue<=0.1),10), aes(label=gene,fontface='bold',segment.size=1), box.padding = 0.5, size=6, xlim  = c(3,4), min.segment.length = 0) +
          geom_label_repel(data=tail(filter(d.training.summary2, Pvalue<=0.1),10), aes(label=gene,fontface='bold',segment.size=1), box.padding = 0.5, size=6, xlim  = c(-4,-3), min.segment.length = 0) +
          guides(colour = guide_legend(override.aes = list(size=3))) +
          theme_prism(border=TRUE,base_size = 25) +
          theme(legend.key.height = unit(30, "pt"),
                legend.title = element_text(size=25),
                legend.text = element_text(size=25),
                legend.margin = margin(0, 30, 0, 0)) +
          labs(title = "DESeq2 log-fold changes and p-values"))
  dev.off()
  
  dds2 <- as.data.frame(mcols(dds)) %>% 
    select(baseMean, dispGeneEst, dispFit, dispersion) %>% 
    melt(id.vars="baseMean") %>% 
    filter(baseMean>0) %>%
    arrange(factor(variable, levels=c("dispGeneEst","dispersion","dispFit")))
  
  png(file=paste0('output/', n, '/training/figures/differential expression/dispersion.png'), width=1920, height=1020)
  print(ggplot(dds2, aes(x=baseMean, y=value, colour=variable)) + 
          geom_point(size=1.5) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_bw() + 
          ylab("Dispersion") + 
          xlab("BaseMean") +
          scale_colour_manual(
            values=c("Black", "#377eb8", "Red"), 
            breaks=c("dispGeneEst", "dispersion", "dispFit"), 
            labels=c("Estimate", "Final", "Fit"),
            name=""
          ) +
          guides(colour = guide_legend(override.aes = list(size=2))) +
          theme_prism(base_size = 30) +
          theme(legend.key.height = unit(30, "pt"),
                legend.title = element_text(size=25),
                legend.text = element_text(size=25)) +
          labs(title = "DESeq2 dispersion"))
  dev.off()
  
  #  use <- res$baseMean > metadata(res)$filterThreshold
  #  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  #  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  #  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  #  png(file=paste0('output/', n, '/training/figures/differential expression/pass.png'), width=1920, height=1020)
  #  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
  #          col = colori, space = 0, main = "", ylab="frequency", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
  #       adj = c(0.5,1.7), xpd=NA)
  #  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  #  dev.off()
  
  #######################################
  #######          GSEA           #######
  #######################################
  dir.create(paste0('output/', n, '/training/figures/GO'))
  dir.create(paste0('output/', n, '/training/output/GO'))
  
  geneList.training <- d.training.summary[,1]
  names(geneList.training) <- rownames(d.training.summary)
  geneList.training <- sort(geneList.training, decreasing = TRUE)
  write.table(names(geneList.training), file=paste0('output/', n, '/training/output/differential expression/DE_genes.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  rows <- names(geneList.training)
  try({rows[substr(rows,1,3) == "LOC"] <- mapIds(org.Hs.eg.db, keys = gsub('LOC', '', rows[substr(rows,1,3) == "LOC"]), column = "SYMBOL", keytype="ENTREZID")}, silent=TRUE)
  try({rows[which(is.na(mapIds(org.Hs.eg.db, keys = rows, column = "ENTREZID", keytype="SYMBOL")))] <- mapIds(org.Hs.eg.db, keys = rows[which(is.na(mapIds(org.Hs.eg.db, keys = rows, column = "ENTREZID", keytype="SYMBOL")))],     column = "SYMBOL", keytype="ALIAS")}, silent=TRUE)
  rows[which(rows == "LOC101929759")] <- "SLCO5A1-AS1"
  rows[which(rows == "LOC107984285")] <- "LARP4B"
  rows[which(rows == "UNQ6494")] <- "LINC03062"
  rows[which(rows == "FLJ22447")] <- "LINC03033"
  rows[which(rows == "FLJ37453")] <- "SPEN-AS1"
  rows[which(rows == "MPP6")] <- "PALS2"
  rows[which(rows == "TIAF1")] <- "MYO18A"
  rows[which(rows == "H4-16")] <- "H4C16"
  rows[which(rows == "CBWD3")] <- "ZNG1C"
  rows[which(rows == "LCA10")] <- "L1CAM-AS1"  
  names(geneList.training) <- mapIds(org.Hs.eg.db, keys = rows, column = "ENTREZID", keytype = "SYMBOL")
  
  gsea.training.BP <- gseGO(geneList     = geneList.training,
                            OrgDb        = org.Hs.eg.db,
                            ont          = "BP",
                            pvalueCutoff = 0.1,
                            eps=0)
  
  gsea.training.MF <- gseGO(geneList     = geneList.training,
                            OrgDb        = org.Hs.eg.db,
                            ont          = "MF",
                            pvalueCutoff = 0.1,
                            eps=0)
  
  
  #gsea.training.KEGG <- GSEA(geneList.training, pvalueCutoff = 0.1, TERM2GENE = msigdbr(species = "Homo sapiens", category = "C2", subcategory="KEGG") %>% 
  #dplyr::select(gs_name, entrez_gene))
  
  
  top_BP <- head(gsea.training.BP@result[order(gsea.training.BP@result$p.adjust), ], 5)
  top_MF <- head(gsea.training.MF@result[order(gsea.training.MF@result$p.adjust), ], 5)
  #top_KEGG <- head(gsea.training.KEGG[order(gsea.training.KEGG@result$p.adjust), ], 5)
  
  top_BP <- top_BP %>% mutate(Description = paste("GO:BP", Description))
  top_MF <- top_MF %>% mutate(Description = paste("GO:MF", Description))
  #top_KEGG <- top_KEGG %>% mutate(Description = paste("KEGG", Description))
  
  #top <- rbind(top_BP,top_MF,top_KEGG)
  top <- rbind(top_BP,top_MF)
  top <- top[order(top$p.adjust), ]
  
  png(file=paste0('output/', n, '/training/figures/GO/GSEA.png'), width=1920, height=1020)
  print(ggplot(top, aes(x = NES, y = reorder(Description, -p.adjust), fill = p.adjust)) +
          geom_bar(stat = "identity") +
          scale_fill_gradient(low = "red", high = "blue") +
          labs(title = "Gene Ontology Enrichment Analysis",
               x = "NES",
               fill = "Adjusted p-value",
               y="") +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
          theme_prism(base_size = 28) +
          theme(plot.title = element_text(hjust = 0.75)) +
          theme(legend.text = element_text(size=24),
                legend.title = element_text(size=26)))
  dev.off()
  
  write.table(rownames(gsea.training.BP@result[gsea.training.BP@result$p.adjust <= 0.1,]), file=paste0('output/', n, '/training/output/GO/GSEA_BP.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(rownames(gsea.training.MF@result[gsea.training.MF@result$p.adjust <= 0.1,]), file=paste0('output/', n, '/training/output/GO/GSEA_MF.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  #write.table(rownames(gsea.training.KEGG@result[gsea.training.KEGG@result$p.adjust <= 0.1,]), file=paste0('output/', n, '/training/output/GO/GSEA_KEGG.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(names(geneList.training), file=paste0('output/', n, '/training/output/differential expression/DE_genes.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  
  save(geneList.training, d.training.summary, dds, top, file=paste0('output/', n, '/training/rdata/DE_data.RData'))
  
  #######################################
  ######           WGCNA           ######
  #######################################
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
  print(paste("WGCNA:", difftime(time2, time1, units="secs")))
  moduleColors <- labels2colors(wgcna.net$colors)
  modules <- levels(as.factor(moduleColors))
  write.table(modules, file=paste0('output/', n, '/training/output/coexpression/modules.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  
  ids <- ids[moduleColors != "grey"]
  log2vals <- d.training.summary[ids, "log2FC"]
  pvals <- d.training.summary[ids, "Pvalue"]
  names(log2vals) <- ids
  names(pvals) <- ids
  
  ### PLOTS
  TOM <- TOMsimilarityFromExpr(t(d.training.adj[ids,]), power=power, TOMType="unsigned")
  diss2 <- 1-TOM
  hier2 <- hclust(as.dist(diss2), method="average")
  colorDynamicTOM <- labels2colors(cutreeDynamic(hier2,method="tree"))
  diag(diss2) = NA;
  moduleColors2 <- moduleColors[moduleColors != "grey"]
  
  datME <- moduleEigengenes(t(d.training.adj[ids,]),colorDynamicTOM)$eigengenes
  dissimME <- (1-t(cor(datME, use="p")))/2
  hclustdatME <- hclust(as.dist(dissimME), method="average")
  par(mfrow=c(1,1))
  png(file=paste0('output/', n, '/training/figures/coexpression/moduleEigengenes.png'), width=1920, height=1020)
  plot(hclustdatME, main="Clustering tree based of the module eigengenes", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  #  GS1 <- as.numeric(cor(group,t(d.training.adj[ids,])))
  #  GeneSignificance <- abs(GS1)
  #  ModuleSignificance <- tapply(GeneSignificance, colorDynamicTOM, mean, na.rm=T)
  #  png(file=paste0('output/', n, '/training/figures/coexpression/moduleSignificance2.png'), width=1920, height=1020)
  #    plotModuleSignificance(GeneSignificance,colorDynamicTOM, cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #
  #  MEs <- orderMEs(datME)
  #  modNames <- substring(names(MEs), 3)
  #  geneModuleMembership <- as.data.frame(cor(t(d.training.adj[ids,]), MEs, use = "p"));
  #  names(geneModuleMembership) <- paste0("MM", modNames);
  #  GS1 <- as.data.frame(GS1)
  #  names(GS1) = paste0("GS.", names(groups))
  
  disableWGCNAThreads()
  
  save(hclustdatME, file=(paste0('output/', n, '/training/rdata/coexpression_data.RData')))
  
  ########################################
  ######     Hierarchical HotNet    ######
  ########################################
  dir.create(paste0('output/', n, '/training/figures/hotnet'))
  dir.create(paste0('output/', n, '/training/output/hotnet'))
  dir.create(paste0('output/', n, '/training/output/hotnet/HotNet_input'))
  
  corPVal <- corAndPvalue(t(d.training.adj[ids,]), use="pairwise.complete.obs")
  
  cutOff <- 0.5
  print(paste("Cut off:", cutOff))
  
  hotnetColors <- c()
  
  for (module in modules){
    if (module == "grey"){
      next
    } 
    else {
      print(module)
      inModule <- moduleColors2 == module
      sum(inModule)
      TOMmodule <- TOM[inModule, inModule]
      corModule <- corPVal[["cor"]][inModule, inModule]
      idsModule <- ids[inModule]
      write.table(idsModule, file=paste0('output/', n, '/training/output/coexpression/', module, '.tsv'), sep='\t', row.names=FALSE, col.names=FALSE)
      log2Module <- log2vals[inModule]
      PModule <- pvals[inModule]
      
      node.frame <- data.frame(Symbol=idsModule, LogFC=abs(log2Module), Pvalue=PModule)
      rownames(node.frame) <- 1:nrow(node.frame)
      #node.frame$InvPvalue <- 1 - PModuleProteins
      # write.table(node.frame[, c('Symbol', 'InvPvalue')], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/g2s_Pval_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      names(idsModule) <- idsModule
      
      write.table(node.frame[, c('Symbol', 'LogFC')], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      
      time1 <- Sys.time()
      combine <- function(x,y){
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
            result1[nrow(result1)+1,] <- c(idsModule[i], idsModule[j])
            result2[nrow(result2)+1,] <- c(i,j)
            if (!i %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(i, idsModule[i])
              result4[nrow(result4)+1,] <- c(idsModule[i], abs(log2Module[i]))
            } 
            if (!j %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(j, idsModule[j])
              result4[nrow(result4)+1,] <- c(idsModule[j], abs(log2Module[j]))
            }
          }
        }
        return(list(result1, result2, result3, result4))
      }
      time2 <- Sys.time()
      print(paste("Slicing input:", difftime(time2, time1, units="secs")))
      time1 <- Sys.time()
      write.table(result[[1]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/name_edges_expression_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[2]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/edge_list_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[3]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/i2g_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[4]], file=paste0('output/', n, '/training/output/hotnet/HotNet_input/g2s_log2', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      time2 <- Sys.time()
      print(paste("Writing data:", difftime(time2, time1, units="secs")))
      if (file.info(paste0('output/', n, '/training/output/hotnet/HotNet_input/i2g_', module, '.tsv'))$size != 0){
        hotnetColors <- append(hotnetColors, module)
      }
    }
    
    idsModule2 <- idsModule
    
    try({idsModule2[substr(idsModule2,1,3) == "LOC"] <- mapIds(org.Hs.eg.db, keys = gsub('LOC', '', idsModule2[substr(idsModule2,1,3) == "LOC"]), column = "SYMBOL", keytype="ENTREZID")}, silent=TRUE)
    
    try({idsModule2[which(is.na(mapIds(org.Hs.eg.db, keys = idsModule2, column = "ENTREZID", keytype="SYMBOL")))] <- mapIds(org.Hs.eg.db, keys = idsModule2[which(is.na(mapIds(org.Hs.eg.db, keys = idsModule2, column = "ENTREZID", keytype="SYMBOL")))], column = "SYMBOL", keytype="ALIAS")}, silent=TRUE)
    
    idsModule2[which(idsModule2 == "LOC101929759")] <- "SLCO5A1-AS1"
    idsModule2[which(idsModule2 == "LOC107984285")] <- "LARP4B"
    idsModule2[which(idsModule2 == "UNQ6494")] <- "LINC03062"
    idsModule2[which(idsModule2 == "FLJ22447")] <- "LINC03033"
    idsModule2[which(idsModule2 == "FLJ37453")] <- "SPEN-AS1"
    idsModule2[which(idsModule2 == "MPP6")] <- "PALS2"
    idsModule2[which(idsModule2 == "TIAF1")] <- "MYO18A"
    idsModule2[which(idsModule2 == "H4-16")] <- "H4C16"
    idsModule2[which(idsModule2 == "CBWD3")] <- "ZNG1C"
    idsModule2[which(idsModule2 == "LCA10")] <- "L1CAM-AS1"  
    
    geneNames <- c()
    try({geneNames <- mapIds(org.Hs.eg.db, keys = idsModule2, column = "ENTREZID", keytype = "SYMBOL")}, silent=TRUE)
    
    if (is.null(geneNames)){next}
    
    goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
    goMF <- enrichGO(geneNames, ont="MF", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
    #goKEGG <- enricher(geneNames, pvalueCutoff = 0.1, TERM2GENE = msigdbr(species = "Homo sapiens", category = "C2", subcategory="KEGG") %>% 
    #dplyr::select(gs_name, entrez_gene))
    
    top_BP <- data.frame()
    if (!is.null(goBP)){
      top_BP <- head(goBP@result[order(goBP@result$p.adjust), ], 5)
      top_BP <- top_BP %>% mutate(Description = paste("GO:BP", Description))  
    }
    top_MF <- data.frame()
    if (!is.null(goMF)){
      top_MF <- head(goMF@result[order(goMF@result$p.adjust), ], 5)
      top_MF <- top_MF %>% mutate(Description = paste("GO:MF", Description))        
    }
    #      top_KEGG <- data.frame()
    #      if (!is.null(goKEGG)){
    #        top_KEGG <- head(goKEGG@result[order(goKEGG@result$p.adjust), ], 5)
    #        top_KEGG <- top_KEGG %>% mutate(Description = paste("KEGG", Description))        
    #      }
    
    #top <- rbind(top_BP,top_MF,top_KEGG)
    top <- rbind(top_BP,top_MF)
    
    if (nrow(top) > 0){
      top <- top[order(top$p.adjust), ]   
      png(file=paste0('output/', n, '/training/figures/coexpression/GO_', module, '.png'), width=1920, height=1020)
      print(ggplot(top, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
              geom_bar(stat = "identity") +
              scale_fill_gradient(low = "red", high = "blue") +
              labs(title = "Gene Ontology Enrichment Analysis",
                   x = "Count",
                   fill = "Adjusted p-value",
                   y="") +
              scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
              theme_prism(base_size = 28) +
              theme(plot.title = element_text(hjust = 0.75)) +
              theme(legend.text = element_text(size=24),
                    legend.title = element_text(size=26)))
      dev.off()
    }
    
    write.table(goBP@result[goBP@result$p.adjust <= 0.1,], file=paste0('output/', n, '/training/output/coexpression/', module, 'BP.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
    write.table(goMF@result[goMF@result$p.adjust <= 0.1,], file=paste0('output/', n, '/training/output/coexpression/', module, 'MF.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
    #   write.table(goKEGG@result[goKEGG@result$p.adjust <= 0.1,], file=paste0('output/', n, '/training/output/coexpression/', module, 'KEGG.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
    write.table(top, file=paste0('output/', n, '/training/output/coexpression/GO_', module, '.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
    
    #    
    #      par(mar = c(1, 1, 1, 1))
    #      png(file=paste0('output/', n, '/training/figures/coexpression/', module, '_matrix.png'), width=1920, height=1020)
    #        plotMat(t(scale(t(d.training.adj[colorDynamicTOM==module,]))),rlabels=T,
    #                clabels=T,rcols=module,
    #                title=module, cex.lab = 2, cex.axis = 2, cex.main = 2)
    #      dev.off()
    #      
    #      ME=datME[, paste0("ME",module)]
    #      par(mar = c(2, 1, 1, 1))
    #      png(file=paste0('output/', n, '/training/figures/coexpression/', module, '_ME.png'), width=1920, height=1020)
    #      barplot(ME, col=module, main="", cex.main=2,
    #              ylab="eigengene expression",xlab="array sample", cex.lab = 2, cex.axis = 2, cex.main = 2)
    #      dev.off()
    #      
    #      par(mfrow = c(1,1));
    #      column <- match(module, modNames);
    #      png(file=paste0('output/', n, '/training/figures/coexpression/', module, '_membershipVSsignficance.png'), width=1920, height=1020)
    #        verboseScatterplot(abs(geneModuleMembership[inModule, column]),
    #                           abs(GS1[inModule, 1]),
    #                           xlab = paste("Module Membership in", module, "module"),
    #                           ylab = "Gene significance",
    #                           main = paste("Module membership vs. gene significance\n"),
    #                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    #      dev.off()
    #      
    #      }
  }
  
  time1 <- Sys.time()
  system(paste('bash', paste(argv[1],'/src/run_hierarchicalHotnet_modules.sh "', sep=""), paste(hotnetColors, collapse=' '), '" ', cutOff, n, 'training'))
  time2 <- Sys.time()
  print(paste("Hierarchical HotNet:", difftime(time2, time1, units="secs")))
  
  time1 <- Sys.time()
  
  BPterms <- c()
  MFterms <- c()
  genes <- c()
  subnetworks <- c()
  hotnetSubnetworks <- c()
  
  for (module in hotnetColors){
    if (file.exists(file=paste0('output/', n, '/training/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'))){
      p <- read.table(file=paste0('output/', n, '/training/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'), sep='\t', header= FALSE, comment.char="")[6, "V1"]
      p <- as.numeric(sub(".*: ", "", p))
      if (file.exists(paste0('output/', n, '/training/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv')) &&
          file.info(paste0('output/', n, '/training/output/hotnet//HotNet_results/consensus_nodes_log2_', module, '.tsv'))$size == 0){
        next
      }
      if (!(module %in% hotnetSubnetworks)){hotnetSubnetworks <- append(hotnetSubnetworks, module)}
      if (p <= 0.1){
        print(module)
        if (!(module %in% subnetworks)){subnetworks <- append(subnetworks, module)}
        nodes <- read.csv(paste0('output/', n, '/training/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv'), sep='\t', header = FALSE)
        edges <- read.csv(paste0('output/', n, '/training/output/hotnet/HotNet_results/consensus_edges_log2_', module, '.tsv'), sep='\t', header = FALSE)
        tidy_basic <- edges %>% 
          as_tbl_graph()
        scores <- d.training.summary[unlist(nodes),"log2FC"]
        tidy_basic_manipulated <- tidy_basic %>% 
          activate(nodes) %>%
          mutate(colours = scores)  
        
        png(file=paste0('output/', n, '/training/figures/GO_', module, '_subnetwork.png'), width=1920, height=1020)
        print(ggraph(tidy_basic_manipulated, layout = 'circle')+
                geom_edge_link(width = 1.2,alpha=0.2) +
                geom_node_point(aes(color = scores), size=15) +
                coord_fixed() +
                geom_node_text(aes(label = name,fontface='bold'),size=5.5, repel=TRUE) +
                scale_color_gradient2(low = "blue", mid="white", high = "red", limits = c(3,3), breaks = c(-2.5,-2,-1,0,1,2,2.5)) +
                scale_x_continuous(limits=c(-2,2)) +
                scale_y_continuous(limits=c(-2,2)) +
                labs(colour="Log-fold change", title=paste(module, "subnetwork")) +
                theme_prism(base_size = 25,base_fontface = "bold") +
                theme(legend.key.height = unit(30, "pt"),
                      legend.title = element_text(size=20),
                      legend.key.width = unit(15, "pt"),
                      legend.text = element_text(size=20),
                      legend.position = c(1.1, 0.5)))
        dev.off()
        
        for (rown in nrow(nodes)){
          nodes2 <- unlist(nodes[rown,])
          print(paste("In subnetwork: ", nodes2))
          genes <- append(genes, nodes2)
          
          log2Subnetwork <- log2vals[nodes2]
          print(paste("Log-fold changes: ", log2Subnetwork))
          
          PSubnetwork <- pvals[nodes2]
          print(paste("P-values: ", PSubnetwork))
          
          try({nodes2[substr(nodes2,1,3) == "LOC"] <- mapIds(org.Hs.eg.db, keys = gsub('LOC', '', nodes2[substr(nodes2,1,3) == "LOC"]), column = "SYMBOL", keytype="ENTREZID")}, silent=TRUE)
          
          try({nodes2[which(is.na(mapIds(org.Hs.eg.db, keys = nodes2, column = "ENTREZID", keytype="SYMBOL")))] <- mapIds(org.Hs.eg.db, keys = nodes2[which(is.na(mapIds(org.Hs.eg.db, keys = nodes2, column = "ENTREZID", keytype="SYMBOL")))], column = "SYMBOL", keytype="ALIAS")}, silent=TRUE)
          
          nodes2[which(nodes2 == "LOC101929759")] <- "SLCO5A1-AS1"
          nodes2[which(nodes2 == "LOC107984285")] <- "LARP4B"
          nodes2[which(nodes2 == "UNQ6494")] <- "LINC03062"
          nodes2[which(nodes2 == "FLJ22447")] <- "LINC03033"
          nodes2[which(nodes2 == "FLJ37453")] <- "SPEN-AS1"
          nodes2[which(nodes2 == "MPP6")] <- "PALS2"
          nodes2[which(nodes2 == "TIAF1")] <- "MYO18A"
          nodes2[which(nodes2 == "H4-16")] <- "H4C16"
          nodes2[which(nodes2 == "CBWD3")] <- "ZNG1C"
          nodes2[which(nodes2 == "LCA10")] <- "L1CAM-AS1"  
          
          geneNames <- c()
          try({geneNames <- mapIds(org.Hs.eg.db, keys = nodes2, column = "ENTREZID", keytype = "SYMBOL")}, silent=TRUE)
          
          print(paste("Number of genes: ", length(geneNames)))
          if (is.null(geneNames)){next}
          else {
            goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
            goMF <- enrichGO(geneNames, ont="MF", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
            top_BP <- data.frame()
            top_MF <- data.frame()
            if (!(is.null(goBP))){
              if (nrow(goBP) > 0){
                print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust <= 0.1,])))
                
                write.table(goBP@result[goBP@result$p.adjust <= 0.1,], file=paste0('output/', n, '/training/output/GO/', module, rown, 'BP.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
                BPterms <- append(BPterms, rownames(goBP@result[goBP@result$p.adjust <= 0.1,]))
                
                top_BP <- head(goBP@result[order(goBP@result$p.adjust), ], 5)
                top_BP <- top_BP %>% mutate(Description = paste("GO:BP", Description))  
              }
            }
            if (!(is.null(goMF))){
              if (nrow(goMF) > 0){
                print(paste("Number of enriched molecular functions: ", nrow(goMF[goMF$p.adjust <= 0.1,])))
                
                write.table(goMF@result[goMF@result$p.adjust <= 0.1,], file=paste0('output/', n, '/training/output/GO/', module, rown, 'MF.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
                MFterms <- append(MFterms, rownames(goMF@result[goMF@result$p.adjust <= 0.1,]))
                
                top_MF <- head(goMF@result[order(goMF@result$p.adjust), ], 5)
                top_MF <- top_MF %>% mutate(Description = paste("GO:MF", Description))
              }
            }
            top <- rbind(top_BP,top_MF)
            if (nrow(top) > 0){
              top <- top[order(top$p.adjust), ]
              png(file=paste0('output/', n, '/training/figures/GO_', module, rown, '.png'), width=1920, height=1020)
              print(ggplot(top, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
                      geom_bar(stat = "identity") +
                      scale_fill_gradient(low = "red", high = "blue") +
                      labs(title = "Gene Ontology Enrichment Analysis",
                           x = "Count",
                           fill = "Adjusted p-value",
                           y="") +
                      scale_y_discrete(labels = function(x) str_wrap(x, width = 65)) +
                      theme_prism(base_size = 18) +
                      theme(legend.text = element_text(size=25),
                            legend.title = element_text(size=25)))
              dev.off()
              write.table(top, file=paste0('output/', n, '/training/output/GO_', module, rown, '.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
            }
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
  print(paste("GO:", difftime(time2, time1, units="secs")))
  endTime <- Sys.time()
  print(paste("Total time:", difftime(endTime, startTime, units="secs")))
  
  ######################################
  ######    Internal validation   ######
  ######################################
  dir.create(paste0('output/', n, '/validation'))
  dir.create(paste0('output/', n, '/validation/output'))
  dir.create(paste0('output/', n, '/validation/figures'))
  dir.create(paste0('output/', n, '/validation/rdata'))
  
  selection.val <- colnames(d.raw)[!(colnames(d.raw) %in% selection)]
  selection.val <- sort(selection.val)
  save(selection.val, file=paste0('output/', n, '/validation/rdata/sampleSelectionVal.RData'))
  d.validation <- d.raw[,selection.val]
  samples <- colnames(d.validation)
  ids <- rownames(d.validation)
  group.data.validation <- group.data[selection.val,]
  d.validation.adj <- d.adj[,selection.val]
  colnames(d.validation.adj) <- samples
  groups <- group.data.validation[,hit]
  
  time2 <- Sys.time()
  print(paste("Processing data:", difftime(time2, time1, units="secs")))
  
  group <- as.matrix(group.data.validation[,hit])
  rownames(group) <- rownames(group.data.validation)
  colnames(group) <- c(hit)
  
  ### PLOTS 
  #
  #  d.validation.cs <- d.cs[,selection.val]
  #  d.validation.cs.WT <- d.validation.cs[,rownames(group.data.validation[group.data.validation[,hit] == "WT",])]
  #  d.validation.cs.SNV <- d.validation.cs[,rownames(group.data.validation[group.data.validation[,hit] == "SNV",])]
  #
  #  colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
  #  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
  #  png(file=paste0('output/', n, '/validation/figures/heatmapSNVSubset.png'), width=1920, height=1020)
  #    pheatmap(d.validation.cs.SNV,
  #             legend=TRUE,
  #             color=my_palette,
  #             breaks=colors,
  #             show_rownames=FALSE,
  #             show_colnames=FALSE
  #    )
  #  dev.off()
  #  
  #  png(file=paste0('output/', n, '/validation/figures/heatmapWTSubset.png'), width=1920, height=1020)
  #    pheatmap(d.validation.cs.WT,
  #             legend=TRUE,
  #             color=my_palette,
  #             breaks=colors,
  #             show_rownames=FALSE,
  #             show_colnames=FALSE
  #    )
  #  dev.off()
  #
  #  png(file=paste0('output/', n, '/validation/figures/heatmap.png'), width=1920, height=1020)
  #    pheatmap(d.validation.cs, color=my_palette, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, show_colnames = FALSE, annotation_col=as.data.frame(group), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  
  colnames(d.validation) <- groups
  
  #  png(file=paste0('output/', n, '/validation/figures/sampleCluster.png'), width=1920, height=1020)
  #  plotClusterTreeSamples(t(d.validation), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #  
  #  sampleDists <- dist(t(d.validation))
  #  sampleDistMatrix <- as.matrix(sampleDists)
  #  rownames(sampleDistMatrix) <- colnames(d.validation)
  #  colnames(sampleDistMatrix) <- colnames(d.validation)
  #  png(file=paste0('output/', n, '/validation/figures/sampleDistances.png'), width=1920, height=1020)
  #    pheatmap(sampleDistMatrix,
  #             clustering_distance_rows=sampleDists,
  #             clustering_distance_cols=sampleDists,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #  
  #  geneDists <- dist(d.validation)
  #  geneDistMatrix <- as.matrix(geneDists)
  #  rownames(geneDistMatrix) <- rownames(d.validation)
  #  colnames(geneDistMatrix) <- rownames(d.validation)
  #  png(file=paste0('output/', n, '/validation/figures/geneDistances.png'), width=1920, height=1020)
  #    pheatmap(geneDistMatrix,
  #             clustering_distance_rows=geneDists,
  #             clustering_distance_cols=geneDists,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #  
  #  png(file=paste0('output/', n, '/validation/figures/meanExpressionSelection.png'), width=1920, height=1020)
  #  barplot(apply(t(d.validation),1,mean, na.rm=T),
  #          xlab = "Sample", ylab = "Mean expression",
  #          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #
  #  png(file=paste0('output/', n, '/validation/figures/heatmapVST.png'), width=1920, height=1020)
  #    pheatmap(d.validation.adj, color=my_palette, cluster_rows=TRUE, show_colnames=FALSE, cluster_cols=TRUE, show_rownames = FALSE, annotation_col=as.data.frame(group))
  #  dev.off()
  
  colnames(d.validation.adj) <- groups
  
  #  png(file=paste0('output/', n, '/validation/figures/sampleClusterVST.png'), width=1920, height=1020)
  #    plotClusterTreeSamples(t(d.validation.adj), as.numeric(factor(groups)), cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #  
  #  sampleDists.adj <- dist(t(d.validation.adj))
  #  sampleDistMatrix.adj <- as.matrix(sampleDists.adj)
  #  rownames(sampleDistMatrix.adj) <- colnames(d.validation.adj)
  #  colnames(sampleDistMatrix.adj) <- colnames(d.validation.adj)
  #  png(file=paste0('output/', n, '/validation/figures/sampleDistancesVST.png'), width=1920, height=1020)
  #    pheatmap(sampleDistMatrix.adj,
  #             clustering_distance_rows=sampleDists.adj,
  #             clustering_distance_cols=sampleDists.adj,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #  
  #  geneDists.adj <- dist(d.validation.adj)
  #  geneDistMatrix.adj <- as.matrix(geneDists.adj)
  #  rownames(geneDistMatrix.adj) <- rownames(d.validation.adj)
  #  colnames(geneDistMatrix.adj) <- rownames(d.validation.adj)
  #  png(file=paste0('output/', n, '/validation/figures/geneDistancesVST.png'), width=1920, height=1020)
  #    pheatmap(geneDistMatrix.adj,
  #             clustering_distance_rows=geneDists.adj,
  #             clustering_distance_cols=geneDists.adj,
  #             col=my_palette, cex.lab = 2, cex.axis = 2, cex.main = 2, show_colnames = FALSE, show_rownames = FALSE)
  #  dev.off()
  #
  #  png(file=paste0('output/', n, '/validation/figures/meanExpressionSelectionVST.png'), width=1920, height=1020)
  #  barplot(apply(t(d.validation.adj),1,mean, na.rm=T),
  #          xlab = "Sample", ylab = "Mean expression",
  #          main ="Mean expression across samples", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  
  ### Save input data
  save(d.validation, d.validation.adj, group.data.validation, group, file=(paste0('output/', n, '/validation/rdata/input_data.RData')))
  
  #######################################
  #######         DESeq2           ######
  #######################################
  dir.create(paste0('output/', n, '/validation/figures/differential expression'))
  dir.create(paste0('output/', n, '/validation/output/differential expression'))
  
  time1 <- Sys.time()
  colnames(d.validation) <- samples
  group[,1] <- factor(group[,1], levels = c("WT","SNV"))
  dds <- DESeqDataSetFromMatrix(countData = d.validation, colData = group, design = formula(paste("~",hit)))
  dds
  
  dds <- DESeq(dds)
  res <- results(dds)
  res
  summary(res)
  
  colnames(d.validation) <- groups
  
  d.validation.summary <- data.frame(unlist(res$log2FoldChange), unlist(res$padj), row.names = rownames(res))
  colnames(d.validation.summary) <- c("log2FC", "Pvalue")
  d.validation.summary <- na.omit(d.validation.summary)
  d.validation.adj <- d.validation.adj[rownames(d.validation.adj) %in% rownames(d.validation.summary),]
  ids <- rownames(d.validation.adj)
  time2<- Sys.time()
  print(paste("DESeq2:", difftime(time2, time1, units="secs")))
  
  write.table(d.validation.summary[,"log2FC", drop = FALSE], file=paste0('output/', n, '/validation/output/differential expression/DE_LFC.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  ###DESEQ PLOTS
  d.validation.summary2 <- d.validation.summary %>% 
    mutate(
      Expression = case_when(log2FC >= 0.1 & Pvalue <= 0.1 ~ "Up-regulated",
                             log2FC <= -0.1 & Pvalue <= 0.1 ~ "Down-regulated",
                             TRUE ~ "Unchanged")
    )  %>%
    arrange(desc(log2FC))
  
  d.validation.summary2 <- cbind(gene=rownames(d.validation.summary2), d.validation.summary2)
  
  png(file=paste0('output/', n, '/validation/figures/differential expression/volcano.png'), width=1920, height=1020)
  print(ggplot(d.validation.summary2, aes(log2FC, -log(Pvalue,10))) +
          geom_point(aes(color = Expression), size = 2) +
          xlab(expression("log"[2]*"FC")) + 
          ylab(expression("-log"[10]*"adjusted p-value")) +
          scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Unchanged" = "grey")) +
          scale_x_continuous(
            limits = c(-4, 4),
            breaks = seq(-4, 4, by = 0.5),
            guide = "prism_minor") +
          scale_y_continuous(
            limits = c(0, 15),
            breaks = seq(0, 15, by = 3)) +
          geom_label_repel(data=head(filter(d.validation.summary2, Pvalue<=0.1),10), aes(label=gene,fontface='bold',segment.size=1), box.padding = 0.5, size=6, xlim  = c(3,4), min.segment.length = 0) +
          geom_label_repel(data=tail(filter(d.validation.summary2, Pvalue<=0.1),10), aes(label=gene,fontface='bold',segment.size=1), box.padding = 0.5, size=6, xlim  = c(-4,-3), min.segment.length = 0) +
          guides(colour = guide_legend(override.aes = list(size=3))) +
          theme_prism(border=TRUE,base_size = 25) +
          theme(legend.key.height = unit(30, "pt"),
                legend.title = element_text(size=25),
                legend.text = element_text(size=25),
                legend.margin = margin(0, 30, 0, 0)) +
          labs(title = "DESeq2 log-fold changes and p-values"))
  dev.off()
  
  dds2 <- as.data.frame(mcols(dds)) %>% 
    select(baseMean, dispGeneEst, dispFit, dispersion) %>% 
    melt(id.vars="baseMean") %>% 
    filter(baseMean>0) %>%
    arrange(factor(variable, levels=c("dispGeneEst","dispersion","dispFit")))
  
  png(file=paste0('output/', n, '/validation/figures/differential expression/dispersion.png'), width=1920, height=1020)
  print(ggplot(dds2, aes(x=baseMean, y=value, colour=variable)) + 
          geom_point(size=1.5) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_bw() + 
          ylab("Dispersion") + 
          xlab("BaseMean") +
          scale_colour_manual(
            values=c("Black", "#377eb8", "Red"), 
            breaks=c("dispGeneEst", "dispersion", "dispFit"), 
            labels=c("Estimate", "Final", "Fit"),
            name=""
          ) +
          guides(colour = guide_legend(override.aes = list(size=2))) +
          theme_prism(base_size = 30) +
          theme(legend.key.height = unit(30, "pt"),
                legend.title = element_text(size=25),
                legend.text = element_text(size=25)) +
          labs(title = "DESeq2 dispersion"))
  dev.off()
  
  #  use <- res$baseMean > metadata(res)$filterThreshold
  #  h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
  #  h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
  #  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  #  png(file=paste0('output/', n, '/validation/figures/differential expression/pass.png'), width=1920, height=1020)
  #    barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
  #            col = colori, space = 0, main = "", ylab="frequency", cex.lab = 2, cex.axis = 2, cex.main = 2)
  #    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
  #         adj = c(0.5,1.7), xpd=NA)
  #    legend("topright", fill=rev(colori), legend=rev(names(colori)))
  #  dev.off()
  
  #######################################
  #######          GSEA           #######
  #######################################
  dir.create(paste0('output/', n, '/validation/figures/GO'))
  dir.create(paste0('output/', n, '/validation/output/GO'))
  
  geneList.validation <- d.validation.summary[,1]
  names(geneList.validation) <- rownames(d.validation.summary)
  geneList.validation <- sort(geneList.validation, decreasing = TRUE)
  write.table(names(geneList.validation), file=paste0('output/', n, '/validation/output/differential expression/DE_genes.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  rows <- names(geneList.validation)
  try({rows[substr(rows,1,3) == "LOC"] <- mapIds(org.Hs.eg.db, keys = gsub('LOC', '', rows[substr(rows,1,3) == "LOC"]), column = "SYMBOL", keytype="ENTREZID")}, silent=TRUE)
  try({rows[which(is.na(mapIds(org.Hs.eg.db, keys = rows, column = "ENTREZID", keytype="SYMBOL")))] <- mapIds(org.Hs.eg.db, keys = rows[which(is.na(mapIds(org.Hs.eg.db, keys = rows, column = "ENTREZID", keytype="SYMBOL")))],     column = "SYMBOL", keytype="ALIAS")}, silent=TRUE)
  rows[which(rows == "LOC101929759")] <- "SLCO5A1-AS1"
  rows[which(rows == "LOC107984285")] <- "LARP4B"
  rows[which(rows == "UNQ6494")] <- "LINC03062"
  rows[which(rows == "FLJ22447")] <- "LINC03033"
  rows[which(rows == "FLJ37453")] <- "SPEN-AS1"
  rows[which(rows == "MPP6")] <- "PALS2"
  rows[which(rows == "TIAF1")] <- "MYO18A"
  rows[which(rows == "H4-16")] <- "H4C16"
  rows[which(rows == "CBWD3")] <- "ZNG1C"
  rows[which(rows == "LCA10")] <- "L1CAM-AS1" 
  names(geneList.validation) <- mapIds(org.Hs.eg.db, keys = rows, column = "ENTREZID", keytype = "SYMBOL")
  
  gsea.validation.BP <- gseGO(geneList     = geneList.validation,
                              OrgDb        = org.Hs.eg.db,
                              ont          = "BP",
                              pvalueCutoff = 0.1,
                              eps = 0)
  
  gsea.validation.MF <- gseGO(geneList     = geneList.validation,
                              OrgDb        = org.Hs.eg.db,
                              ont          = "MF",
                              pvalueCutoff = 0.1,
                              eps = 0)
  
  
  #gsea.validation.KEGG <- GSEA(geneList.validation, pvalueCutoff = 0.1, TERM2GENE = msigdbr(species = "Homo sapiens", category = "C2", subcategory="KEGG") %>% 
  #dplyr::select(gs_name, entrez_gene))
  
  top_BP <- head(gsea.validation.BP@result[order(gsea.validation.BP@result$p.adjust), ], 5)
  top_MF <- head(gsea.validation.MF@result[order(gsea.validation.MF@result$p.adjust), ], 5)
  #top_KEGG <- head(gsea.validation.KEGG[order(gsea.validation.KEGG@result$p.adjust), ], 5)
  
  top_BP <- top_BP %>% mutate(Description = paste("GO:BP", Description))
  top_MF <- top_MF %>% mutate(Description = paste("GO:MF", Description))
  #top_KEGG <- top_KEGG %>% mutate(Description = paste("KEGG", Description))
  
  #top <- rbind(top_BP,top_MF,top_KEGG)
  top <- rbind(top_BP,top_MF)
  top <- top[order(top$p.adjust), ]
  
  png(file=paste0('output/', n, '/validation/figures/GO/GSEA.png'), width=1920, height=1020)
  print(ggplot(top, aes(x = NES, y = reorder(Description, -p.adjust), fill = p.adjust)) +
          geom_bar(stat = "identity") +
          scale_fill_gradient(low = "red", high = "blue") +  # You can adjust colors as needed
          labs(title = "Gene Ontology Enrichment Analysis",
               x = "NES",
               fill = "Adjusted p-value",
               y="") +
          scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
          theme_prism(base_size = 28) +
          theme(plot.title = element_text(hjust = 0.75)) +
          theme(legend.text = element_text(size=24),
                legend.title = element_text(size=26)))
  dev.off()
  
  write.table(rownames(gsea.validation.BP@result[gsea.validation.BP@result$p.adjust <= 0.1,]), file=paste0('output/', n, '/validation/output/GO/GSEA_BP.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  write.table(rownames(gsea.validation.MF@result[gsea.validation.MF@result$p.adjust <= 0.1,]), file=paste0('output/', n, '/validation/output/GO/GSEA_MF.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  #write.table(rownames(gsea.validation.KEGG@result[gsea.validation.KEGG@result$p.adjust <= 0.1,]), file=paste0('output/', n, '/validation/output/GO/GSEA_KEGG.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  save(geneList.validation, d.validation.summary, dds, top, file=paste0('output/', n, '/validation/rdata/DE_data.RData'))
  
  #######################################
  ######           WGCNA           ######
  #######################################
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
  print(paste("WGCNA:", difftime(time2, time1, units="secs")))
  moduleColors <- labels2colors(wgcna.net$colors)
  modules <- levels(as.factor(moduleColors))
  write.table(modules, file=paste0('output/', n, '/validation/output/coexpression/modules.tsv'), row.names=FALSE, col.names=FALSE, sep='\t')
  
  ids <- ids[moduleColors != "grey"]
  log2vals <- d.validation.summary[ids, "log2FC"]  
  pvals <- d.validation.summary[ids, "Pvalue"]
  names(log2vals) <- ids
  names(pvals) <- ids
  
  ### PLOTS
  TOM <- TOMsimilarityFromExpr(t(d.validation.adj[ids,]), power=power, TOMType="unsigned")
  diss2 <- 1-TOM
  hier2 <- hclust(as.dist(diss2), method="average")
  colorDynamicTOM <- labels2colors(cutreeDynamic(hier2,method="tree"))
  diag(diss2) = NA;
  moduleColors2 <- moduleColors[moduleColors != "grey"]
  
  datME <- moduleEigengenes(t(d.validation.adj[ids,]),colorDynamicTOM)$eigengenes
  dissimME <- (1-t(cor(datME, use="p")))/2
  hclustdatME <- hclust(as.dist(dissimME), method="average")
  par(mfrow=c(1,1))
  png(file=paste0('output/', n, '/validation/figures/coexpression/moduleEigengenes.png'), width=1920, height=1020)
  plot(hclustdatME, main="Clustering tree based of the module eigengenes", cex.lab = 2, cex.axis = 2, cex.main = 2)
  dev.off()
  
  #  GS1 <- as.numeric(cor(group,t(d.validation.adj[ids,])))
  #  GeneSignificance <- abs(GS1)
  #  ModuleSignificance <- tapply(GeneSignificance, colorDynamicTOM, mean, na.rm=T)
  #  png(file=paste0('output/', n, '/validation/figures/coexpression/moduleSignificance2.png'), width=1920, height=1020)
  #    plotModuleSignificance(GeneSignificance,colorDynamicTOM, cex.lab = 2, cex.axis = 2, cex.main = 2)
  #  dev.off()
  #  
  #  MEs <- orderMEs(datME)
  #  modNames <- substring(names(MEs), 3)
  #  geneModuleMembership <- as.data.frame(cor(t(d.validation.adj[ids,]), MEs, use = "p"));
  #  names(geneModuleMembership) <- paste0("MM", modNames);
  #  GS1 <- as.data.frame(GS1)
  #  names(GS1) = paste0("GS.", names(groups))
  
  disableWGCNAThreads()
  
  save(hclustdatME, file=(paste0('output/', n, '/validation/rdata/coexpression_data.RData')))
  
  ########################################
  ######     Hierarchical HotNet    ######
  ########################################
  dir.create(paste0('output/', n, '/validation/figures/hotnet'))
  dir.create(paste0('output/', n, '/validation/output//hotnet'))
  dir.create(paste0('output/', n, '/validation/output/hotnet/HotNet_input'))
  
  corPVal <- corAndPvalue(t(d.validation.adj[ids,]), use="pairwise.complete.obs")
  
  cutOff <- 0.5
  print(paste("Cut off:", cutOff))
  
  hotnetColors <- c()
  
  for (module in modules){
    if (module == "grey"){
      next
    } 
    else {
      print(module)
      
      inModule <- moduleColors2 == module
      sum(inModule)
      TOMmodule <- TOM[inModule, inModule]
      corModule <- corPVal[["cor"]][inModule, inModule]
      idsModule <- ids[inModule]
      write.table(idsModule, file=paste0('output/', n, '/validation/output/coexpression/', module, '.tsv'), sep='\t', row.names=FALSE, col.names=FALSE)
      log2Module <- log2vals[inModule]
      PModule <- pvals[inModule]
      
      node.frame <- data.frame(Symbol=idsModule, LogFC=abs(log2Module), Pvalue=PModule)
      rownames(node.frame) <- 1:nrow(node.frame)
      #node.frame$InvPvalue <- 1 - PModuleProteins
      # write.table(node.frame[, c('Symbol', 'InvPvalue')], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/g2s_Pval_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      names(idsModule) <- idsModule
      
      write.table(node.frame[, c('Symbol', 'LogFC')], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      
      time1 <- Sys.time()
      
      combine <- function(x,y){
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
            result1[nrow(result1)+1,] <- c(idsModule[i], idsModule[j])
            result2[nrow(result2)+1,] <- c(i,j)
            if (!i %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(i, idsModule[i])
              result4[nrow(result4)+1,] <- c(idsModule[i], abs(log2Module[i]))
            } 
            if (!j %in% result3[,1]){
              result3[nrow(result3)+1,] <- c(j, idsModule[j])
              result4[nrow(result4)+1,] <- c(idsModule[j], abs(log2Module[j]))
            }
          }
        }
        return(list(result1, result2, result3, result4))
      }
      time2 <- Sys.time()
      print(paste("Slicing input:", difftime(time2, time1, units="secs")))
      time1 <- Sys.time()
      write.table(result[[1]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/name_edges_expression_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[2]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/edge_list_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[3]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/i2g_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      write.table(result[[4]], file=paste0('output/', n, '/validation/output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
      time2 <- Sys.time()
      print(paste("Writing data:", difftime(time2, time1, units="secs")))
      if (file.info(paste0('output/', n, '/validation/output/hotnet/HotNet_input/i2g_', module, '.tsv'))$size != 0){
        hotnetColors <- append(hotnetColors, module)
      }
    }
    
    idsModule2 <- idsModule
    
    try({idsModule2[substr(idsModule2,1,3) == "LOC"] <- mapIds(org.Hs.eg.db, keys = gsub('LOC', '', idsModule2[substr(idsModule2,1,3) == "LOC"]), column = "SYMBOL", keytype="ENTREZID")},silent=TRUE)
    
    try({idsModule2[which(is.na(mapIds(org.Hs.eg.db, keys = idsModule2, column = "ENTREZID", keytype="SYMBOL")))] <- mapIds(org.Hs.eg.db, keys = idsModule2[which(is.na(mapIds(org.Hs.eg.db, keys = idsModule2, column = "ENTREZID", keytype="SYMBOL")))], column = "SYMBOL", keytype="ALIAS")},silent=TRUE)
    
    idsModule2[which(idsModule2 == "LOC101929759")] <- "SLCO5A1-AS1"
    idsModule2[which(idsModule2 == "LOC107984285")] <- "LARP4B"
    idsModule2[which(idsModule2 == "UNQ6494")] <- "LINC03062"
    idsModule2[which(idsModule2 == "FLJ22447")] <- "LINC03033"
    idsModule2[which(idsModule2 == "FLJ37453")] <- "SPEN-AS1"
    idsModule2[which(idsModule2 == "MPP6")] <- "PALS2"
    idsModule2[which(idsModule2 == "TIAF1")] <- "MYO18A"
    idsModule2[which(idsModule2 == "H4-16")] <- "H4C16"
    idsModule2[which(idsModule2 == "CBWD3")] <- "ZNG1C"
    idsModule2[which(idsModule2 == "LCA10")] <- "L1CAM-AS1"  
    
    geneNames <- c()
    try({geneNames <- mapIds(org.Hs.eg.db, keys = idsModule2, column = "ENTREZID", keytype = "SYMBOL")}, silent=TRUE)
    
    if (is.null(geneNames)){next}
    
    goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
    goMF <- enrichGO(geneNames, ont="MF", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
    
    top_BP <- data.frame()
    if (!is.null(goBP)){
      top_BP <- head(goBP@result[order(goBP@result$p.adjust), ], 5)
      top_BP <- top_BP %>% mutate(Description = paste("GO:BP", Description))  
    }
    top_MF <- data.frame()
    if (!is.null(goMF)){
      top_MF <- head(goMF@result[order(goMF@result$p.adjust), ], 5)
      top_MF <- top_MF %>% mutate(Description = paste("GO:MF", Description))        
    }
    #      top_KEGG <- data.frame()
    #      if (!is.null(goKEGG)){
    #        top_KEGG <- head(goKEGG@result[order(goKEGG@result$p.adjust), ], 5)
    #        top_KEGG <- top_KEGG %>% mutate(Description = paste("KEGG", Description))        
    #      }
    
    #top <- rbind(top_BP,top_MF,top_KEGG)
    top <- rbind(top_BP,top_MF)
    
    if (nrow(top) > 0){
      top <- top[order(top$p.adjust), ]
      png(file=paste0('output/', n, '/validation/figures/coexpression/GO_', module, '.png'), width=1920, height=1020)
      print(ggplot(top, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
              geom_bar(stat = "identity") +
              scale_fill_gradient(low = "red", high = "blue") +
              labs(title = "Gene Ontology Enrichment Analysis",
                   x = "Count",
                   fill = "Adjusted p-value",
                   y="") +
              scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
              theme_prism(base_size = 28) +
              theme(plot.title = element_text(hjust = 0.75)) +
              theme(legend.text = element_text(size=24),
                    legend.title = element_text(size=26)))
      dev.off()
    }
    
    write.table(goBP@result[goBP@result$p.adjust <= 0.1,], file=paste0('output/', n, '/validation/output/coexpression/', module, 'BP.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
    write.table(goMF@result[goMF@result$p.adjust <= 0.1,], file=paste0('output/', n, '/validation/output/coexpression/', module, 'MF.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
    #   write.table(goKEGG@result[goKEGG@result$p.adjust <= 0.1,], file=paste0('output/', n, '/validation/output/coexpression/', module, 'KEGG.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
    write.table(top, file=paste0('output/', n, '/validation/output/coexpression/GO_', module, '.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
    
    #      par(mar = c(1, 1, 1, 1))
    #      png(file=paste0('output/', n, '/validation/figures/coexpression/', module, '_matrix.png'), width=1920, height=1020)
    #        plotMat(t(scale(t(d.validation.adj[colorDynamicTOM==module,]))),rlabels=T,
    #                clabels=T,rcols=module,
    #                title=module, cex.lab = 2, cex.axis = 2, cex.main = 2)
    #      dev.off()
    #      
    #      ME=datME[, paste0("ME",module)]
    #      par(mar = c(2, 1, 1, 1))
    #      png(file=paste0('output/', n, '/validation/figures/coexpression/', module, '_ME.png'), width=1920, height=1020)
    #        barplot(ME, col=module, main="", cex.main=2,
    #                ylab="eigengene expression",xlab="array sample", cex.lab = 2, cex.axis = 2, cex.main = 2)
    #      dev.off()
    #
    #      par(mfrow = c(1,1));
    #      column <- match(module, modNames);
    #      png(file=paste0('output/', n, '/validation/figures/coexpression/', module, '_membershipVSsignficance.png'), width=1920, height=1020)
    #        verboseScatterplot(abs(geneModuleMembership[inModule, column]),
    #                           abs(GS1[inModule, 1]),
    #                           xlab = paste("Module Membership in", module, "module"),
    #                           ylab = "Gene significance",
    #                           main = paste("Module membership vs. gene significance\n"),
    #                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    #      dev.off()
    #      }
  }
  
  
  ######### Run hierarchical hotnet to obtain the most significant submodule within each of the identified modules
  ######### Note that the HotNet package was written in Python and needs to be intstalled separately on your machine
  time1 <- Sys.time()
  system(paste('bash', paste(argv[1],'/src/run_hierarchicalHotnet_modules.sh "', sep=""), paste(hotnetColors, collapse=' '), '" ', cutOff, n, 'validation'))
  time2 <- Sys.time()
  print(paste("Hierarchical HotNet:", difftime(time2, time1, units="secs")))
  
  time1 <- Sys.time()
  
  BPterms <- c()
  MFterms <- c()
  genes <- c()
  subnetworks <- c()
  hotnetSubnetworks <- c()
  
  for (module in hotnetColors){
    if (file.exists(file=paste0('output/', n, '/validation/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'))){
      p <- read.table(file=paste0('output/', n, '/validation/output/hotnet/HotNet_results/clusters_hierarchies_log2_', module, '.tsv'), sep='\t', header= FALSE, comment.char="")[6, "V1"]
      p <- as.numeric(sub(".*: ", "", p))
      if (file.exists(paste0('output/', n, '/validation/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv')) &&
          file.info(paste0('output/', n, '/validation/output/hotnet//HotNet_results/consensus_nodes_log2_', module, '.tsv'))$size == 0){
        next
      }
      if (!(module %in% hotnetSubnetworks)){hotnetSubnetworks <- append(hotnetSubnetworks, module)}
      if (p <= 0.1){
        print(module)
        if (!(module %in% subnetworks)){subnetworks <- append(subnetworks, module)}
        nodes <- read.csv(paste0('output/', n, '/validation/output/hotnet/HotNet_results/consensus_nodes_log2_', module, '.tsv'), sep='\t', header = FALSE)
        edges <- read.csv(paste0('output/', n, '/validation/output/hotnet/HotNet_results/consensus_edges_log2_', module, '.tsv'), sep='\t', header = FALSE)
        tidy_basic <- edges %>% 
          as_tbl_graph()
        scores <- d.validation.summary[unlist(nodes),"log2FC"]
        tidy_basic_manipulated <- tidy_basic %>% 
          activate(nodes) %>%
          mutate(colours = scores)  
        
        png(file=paste0('output/', n, '/validation/figures/GO_', module, '_subnetwork.png'), width=1920, height=1020)
        print(ggraph(tidy_basic_manipulated, layout = 'circle')+
                geom_edge_link(width = 1.2,alpha=0.2) +
                geom_node_point(aes(color = scores), size=15) +
                coord_fixed() +
                geom_node_text(aes(label = name,fontface='bold'),size=5.5, repel=TRUE) +
                scale_color_gradient2(low = "blue", mid="white", high = "red", limits = c(3,3), breaks = c(-2.5,-2,-1,0,1,2,2.5)) +
                scale_x_continuous(limits=c(-2,2)) +
                scale_y_continuous(limits=c(-2,2)) +
                labs(colour="Log-fold change", title=paste(module, "subnetwork")) +
                theme_prism(base_size = 25,base_fontface = "bold") +
                theme(legend.key.height = unit(30, "pt"),
                      legend.title = element_text(size=20),
                      legend.key.width = unit(15, "pt"),
                      legend.text = element_text(size=20),
                      legend.position = c(1.1, 0.5)))
        dev.off()
        
        for (rown in nrow(nodes)){
          nodes2 <- unlist(nodes[rown,])
          print(paste("In subnetwork: ", nodes2))
          genes <- append(genes, nodes2)
          
          log2Subnetwork <- log2vals[nodes2]
          print(paste("Log-fold changes: ", log2Subnetwork))
          
          PSubnetwork <- pvals[nodes2]
          print(paste("P-values: ", PSubnetwork))
          
          try({nodes2[substr(nodes2,1,3) == "LOC"] <- mapIds(org.Hs.eg.db, keys = gsub('LOC', '', nodes2[substr(nodes2,1,3) == "LOC"]), column = "SYMBOL", keytype="ENTREZID")}, silent=TRUE)
          
          try({nodes2[which(is.na(mapIds(org.Hs.eg.db, keys = nodes2, column = "ENTREZID", keytype="SYMBOL")))] <- mapIds(org.Hs.eg.db, keys = nodes2[which(is.na(mapIds(org.Hs.eg.db, keys = nodes2, column = "ENTREZID", keytype="SYMBOL")))], column = "SYMBOL", keytype="ALIAS")}, silent=TRUE)
          
          nodes2[which(nodes2 == "LOC101929759")] <- "SLCO5A1-AS1"
          nodes2[which(nodes2 == "LOC107984285")] <- "LARP4B"
          nodes2[which(nodes2 == "UNQ6494")] <- "LINC03062"
          nodes2[which(nodes2 == "FLJ22447")] <- "LINC03033"
          nodes2[which(nodes2 == "FLJ37453")] <- "SPEN-AS1"
          nodes2[which(nodes2 == "MPP6")] <- "PALS2"
          nodes2[which(nodes2 == "TIAF1")] <- "MYO18A"
          nodes2[which(nodes2 == "H4-16")] <- "H4C16"
          nodes2[which(nodes2 == "CBWD3")] <- "ZNG1C"
          nodes2[which(nodes2 == "LCA10")] <- "L1CAM-AS1"  
          
          geneNames <- c()
          try({geneNames <- mapIds(org.Hs.eg.db, keys = nodes2, column = "ENTREZID", keytype = "SYMBOL")}, silent=TRUE)
          
          print(paste("Number of genes: ", length(geneNames)))
          if (is.null(geneNames)){next}
          else {
            goBP <- enrichGO(geneNames, ont="BP", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
            goMF <- enrichGO(geneNames, ont="MF", keyType = "ENTREZID", pvalueCutoff = 0.1, OrgDb=org.Hs.eg.db)
            top_BP <- data.frame()
            top_MF <- data.frame()
            if (!(is.null(goBP))){
              if (nrow(goBP) > 0){
                print(paste("Number of enriched biological processes: ", nrow(goBP[goBP$p.adjust <= 0.1,])))
                
                write.table(goBP@result[goBP@result$p.adjust <= 0.1,], file=paste0('output/', n, '/validation/output/GO/', module, rown, 'BP.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
                BPterms <- append(BPterms, rownames(goBP@result[goBP@result$p.adjust <= 0.1,]))
                
                top_BP <- head(goBP@result[order(goBP@result$p.adjust), ], 5)
                top_BP <- top_BP %>% mutate(Description = paste("GO:BP", Description))  
              }
            }
            if (!(is.null(goMF))){
              if (nrow(goMF) > 0){
                print(paste("Number of enriched molecular functions: ", nrow(goMF[goMF$p.adjust <= 0.1,])))
                
                write.table(goMF@result[goMF@result$p.adjust <= 0.1,], file=paste0('output/', n, '/validation/output/GO/', module, rown, 'MF.tsv'), col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)
                MFterms <- append(MFterms, rownames(goMF@result[goMF@result$p.adjust <= 0.1,]))
                
                top_MF <- head(goMF@result[order(goMF@result$p.adjust), ], 5)
                top_MF <- top_MF %>% mutate(Description = paste("GO:MF", Description))
              }
            }
            top <- rbind(top_BP,top_MF)
            if (nrow(top) > 0){
              top <- top[order(top$p.adjust), ]
              png(file=paste0('output/', n, '/validation/figures/GO_', module, rown, '.png'), width=1920, height=1020)
              print(ggplot(top, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
                      geom_bar(stat = "identity") +
                      scale_fill_gradient(low = "red", high = "blue") +
                      labs(title = "Gene Ontology Enrichment Analysis",
                           x = "Count",
                           fill = "Adjusted p-value",
                           y="") +
                      scale_y_discrete(labels = function(x) str_wrap(x, width = 65)) +
                      theme_prism(base_size = 18) +
                      theme(legend.text = element_text(size=25),
                            legend.title = element_text(size=25)))
              dev.off()
              write.table(top, file=paste0('output/', n, '/validation/output/GO_', module, rown, '.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
            }
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
  print(paste("GO:", difftime(time2, time1, units="secs")))
  endTime <- Sys.time()
  print(paste("Total time:", difftime(endTime, startTime, units="secs")))
  
  dir.create(paste0('output/', n, '/validation/output/similarity'))
  dir.create(paste0('output/', n, '/validation/figures/similarity'))
  
  similarity <- similarity.analysis(paste0('output/', n, '/training/output'), paste0('output/', n, '/validation/output'), paste0('output/', n, '/validation/output/similarity'))
  
  print(paste('END OF RUN', n))
}

similarityModules1 <- as.matrix(read.table(paste0('output/1/validation/output/similarity/moduleSimilarity.tsv'), sep='\t'))
similarityModules2 <- as.matrix(read.table(paste0('output/2/validation/output/similarity/moduleSimilarity.tsv'), sep='\t'))
similarityModules3 <- as.matrix(read.table(paste0('output/3/validation/output/similarity/moduleSimilarity.tsv'), sep='\t'))

averageSimilarityModules1 <- mean(apply(similarityModules1,2,max))
averageSimilarityModules2 <- mean(apply(similarityModules2,2,max))
averageSimilarityModules3 <- mean(apply(similarityModules3,2,max))

similaritySubnetworks1 <- as.matrix(read.table(paste0('output/1/validation/output/similarity/subnetworkSimilarity.tsv'), sep='\t'))
similaritySubnetworks2 <- as.matrix(read.table(paste0('output/2/validation/output/similarity/subnetworkSimilarity.tsv'), sep='\t'))
similaritySubnetworks3 <- as.matrix(read.table(paste0('output/3/validation/output/similarity/subnetworkSimilarity.tsv'), sep='\t'))

averageSimilaritySubnetworks1 <- mean(apply(similaritySubnetworks1,2,max))
averageSimilaritySubnetworks2 <- mean(apply(similaritySubnetworks2,2,max))
averageSimilaritySubnetworks3 <- mean(apply(similaritySubnetworks3,2,max))

similaritySSubnetworks1 <- as.matrix(read.table(paste0('output/1/validation/output/similarity/sSubnetworkSimilarity.tsv'), sep='\t'))
similaritySSubnetworks2 <- as.matrix(read.table(paste0('output/2/validation/output/similarity/sSubnetworkSimilarity.tsv'), sep='\t'))
similaritySSubnetworks3 <- as.matrix(read.table(paste0('output/3/validation/output/similarity/sSubnetworkSimilarity.tsv'), sep='\t'))

averageSimilaritySSubnetworks1 <- mean(apply(similaritySSubnetworks1,2,max))
averageSimilaritySSubnetworks2 <- mean(apply(similaritySSubnetworks2,2,max))
averageSimilaritySSubnetworks3 <- mean(apply(similaritySSubnetworks3,2,max))

resultSimilarity1 <- as.matrix(read.table(paste0('output/1/validation/output/similarity/similarity.tsv'), sep='\t'))
resultSimilarity2 <- as.matrix(read.table(paste0('output/2/validation/output/similarity/similarity.tsv'), sep='\t'))
resultSimilarity3 <- as.matrix(read.table(paste0('output/3/validation/output/similarity/similarity.tsv'), sep='\t'))

GSEASimilarity1 <- as.matrix(read.table(paste0('output/1/validation/output/similarity/gsea.tsv'), sep='\t'))
GSEASimilarity2 <- as.matrix(read.table(paste0('output/2/validation/output/similarity/gsea.tsv'), sep='\t'))
GSEASimilarity3 <- as.matrix(read.table(paste0('output/3/validation/output/similarity/gsea.tsv'), sep='\t'))

DESimilarity1 <- as.matrix(read.table(paste0('output/1/validation/output/similarity/de.tsv'), sep='\t'))
DESimilarity2 <- as.matrix(read.table(paste0('output/2/validation/output/similarity/de.tsv'), sep='\t'))
DESimilarity3 <- as.matrix(read.table(paste0('output/3/validation/output/similarity/de.tsv'), sep='\t'))

end <- data.frame(Type=character(0), Similarity=numeric(0))

for (l in apply(similarityModules1,2,max)){
  end[nrow(end) + 1,] <-  c('Module', as.numeric(l))
}
for (l in apply(similaritySubnetworks1,2,max)){
  end[nrow(end) + 1,] <- c('Subnetwork', as.numeric(l))
}
for (l in apply(similaritySSubnetworks1,2,max)){
  end[nrow(end) + 1,] <- c('Significant subnetwork', as.numeric(l))
}

end[nrow(end) + 1,] <- c('GO genes', as.numeric(resultSimilarity1['genes',1]))
end[nrow(end) + 1,] <- c('GO BP', as.numeric(resultSimilarity1['BP',1]))
end[nrow(end) + 1,] <- c('GO MF', as.numeric(resultSimilarity1['MF',1]))
end[nrow(end) + 1,] <- c('GSEA BP', as.numeric(GSEASimilarity1['BP',1]))
end[nrow(end) + 1,] <- c('GSEA MF', as.numeric(GSEASimilarity1['MF',1]))
end[nrow(end) + 1,] <- c('DE', as.numeric(DESimilarity1[1,1]))

for (l in apply(similarityModules2,2,max)){
  end[nrow(end) + 1,] <-  c('Module', as.numeric(l))
}
for (l in apply(similaritySubnetworks2,2,max)){
  end[nrow(end) + 1,] <- c('Subnetwork', as.numeric(l))
}
for (l in apply(similaritySSubnetworks2,2,max)){
  end[nrow(end) + 1,] <- c('Significant subnetwork', as.numeric(l))
}

end[nrow(end) + 1,] <- c('GO genes', as.numeric(resultSimilarity2['genes',1]))
end[nrow(end) + 1,] <- c('GO BP', as.numeric(resultSimilarity2['BP',1]))
end[nrow(end) + 1,] <- c('GO MF', as.numeric(resultSimilarity2['MF',1]))
end[nrow(end) + 1,] <- c('GSEA BP', as.numeric(GSEASimilarity2['BP',1]))
end[nrow(end) + 1,] <- c('GSEA MF', as.numeric(GSEASimilarity2['MF',1]))
end[nrow(end) + 1,] <- c('DE', as.numeric(DESimilarity2[1,1]))

for (l in apply(similarityModules3,2,max)){
  end[nrow(end) + 1,] <-  c('Module', as.numeric(l))
}
for (l in apply(similaritySubnetworks3,2,max)){
  end[nrow(end) + 1,] <- c('Subnetwork', as.numeric(l))
}
for (l in apply(similaritySSubnetworks3,2,max)){
  end[nrow(end) + 1,] <- c('Significant subnetwork', as.numeric(l))
}

end[nrow(end) + 1,] <- c('GO genes', as.numeric(resultSimilarity3['genes',1]))
end[nrow(end) + 1,] <- c('GO BP', as.numeric(resultSimilarity3['BP',1]))
end[nrow(end) + 1,] <- c('GO MF', as.numeric(resultSimilarity3['MF',1]))
end[nrow(end) + 1,] <- c('GSEA BP', as.numeric(GSEASimilarity3['BP',1]))
end[nrow(end) + 1,] <- c('GSEA MF', as.numeric(GSEASimilarity3['MF',1]))
end[nrow(end) + 1,] <- c('DE', as.numeric(DESimilarity3[1,1]))

end$Similarity <- as.numeric(end$Similarity)

png(file=paste0('output/similarity.png'), width=1920, height=1020)
print(ggplot(end, aes(x=Type, y=Similarity, fill=Type)) +
        geom_boxplot(width=0.5, color="black",show.legend = FALSE, size = 0.8) +
        theme_prism(base_size = 16) +
        scale_y_continuous(name = "Jaccard score", limits=c(0,1), breaks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1)) +
        scale_x_discrete(limits=c('DE', 'GSEA BP', 'GSEA MF',"Module", "Subnetwork", "Significant subnetwork", 'GO genes', 'GO BP', 'GO MF'), labels = str_wrap(c('DE','GSEA BP', 'GSEA MF', 'Modules', 'Subnetworks',         'Significant subnetworks', 'GO genes', 'GO BP', 'GO MF'), width = 10)) +
        labs(x = "Pipeline stage",
             title = "Jaccard scores through pipeline", y="Jaccard score"))
dev.off()

modules1 <- data.frame(matrix(nrow = length(colnames(similarityModules1)), ncol = 6))
rownames(modules1) <- colnames(similarityModules1)
colnames(modules1) <- c('Name1', 'Name2', 'Similarity', 'Size1', 'Size2', 'Run')
modules1[, 'Name1'] <- rownames(modules1)
modules1[, 'Name2'] <- rownames(similarityModules1)[apply(similarityModules1,2,which.max)]
modules1[,'Run'] <- 1
modules1[, 'Similarity'] <- apply(similarityModules1,2,max)
for (i in 1:nrow(modules1)){
  modules1[i, 'Size1'] <- length(as.vector(t(read.table(paste0('output/1/training/output/coexpression/', row.names(modules1[i,]), '.tsv'), sep='\t', header = FALSE))))
  modules1[i, 'Size2'] <- length(as.vector(t(read.table(paste0('output/1/validation/output/coexpression/', modules1[i, 'Name2'], '.tsv'), sep='\t', header = FALSE))))
}

modules2 <- data.frame(matrix(nrow = length(colnames(similarityModules2)), ncol = 6))
rownames(modules2) <- colnames(similarityModules2)
colnames(modules2) <- c('Name1', 'Name2', 'Similarity', 'Size1', 'Size2', 'Run')
modules2[, 'Name1'] <- rownames(modules2)
modules2[, 'Name2'] <- rownames(similarityModules2)[apply(similarityModules2,2,which.max)]
modules2[,'Run'] <- 2
modules2[, 'Similarity'] <- apply(similarityModules2,2,max)
for (i in 1:nrow(modules2)){
  modules2[i, 'Size1'] <- length(as.vector(t(read.table(paste0('output/2/training/output/coexpression/', row.names(modules2[i,]), '.tsv'), sep='\t', header = FALSE))))
  modules2[i, 'Size2'] <- length(as.vector(t(read.table(paste0('output/2/validation/output/coexpression/', modules2[i, 'Name2'], '.tsv'), sep='\t', header = FALSE))))
}

modules3 <- data.frame(matrix(nrow = length(colnames(similarityModules3)), ncol = 6))
rownames(modules3) <- colnames(similarityModules3)
colnames(modules3) <- c('Name1', 'Name2', 'Similarity', 'Size1', 'Size2','Run')
modules3[, 'Name1'] <- rownames(modules3)
modules3[, 'Name2'] <- rownames(similarityModules3)[apply(similarityModules3,2,which.max)]
modules3[, 'Similarity'] <- apply(similarityModules3,2,max)
modules3[,'Run'] <- 3
for (i in 1:nrow(modules3)){
  modules3[i, 'Size1'] <- length(as.vector(t(read.table(paste0('output/3/training/output/coexpression/', row.names(modules3[i,]), '.tsv'), sep='\t', header = FALSE))))
  modules3[i, 'Size2'] <- length(as.vector(t(read.table(paste0('output/3/validation/output/coexpression/', modules3[i, 'Name2'], '.tsv'), sep='\t', header = FALSE))))
}

combined_data <- rbind(modules1, modules2, modules3)

combined_data$Run <- as.character(combined_data$Run)

png(file='output/moduleSize.png', width=1920, height=1020)
print(ggplot(combined_data, aes(x = Size1, y = Similarity, color = Run)) +
        geom_point(size = 3) +
        labs(x = "Number of genes in module (exploration)", y = "Jaccard score", title = "Jaccard scores for different module sizes (exploration)") +
        scale_color_manual(values = c("1" = "magenta", "2" = "purple", "3" = "blue"), name = "Run") +
        theme_prism(border=TRUE,base_size = 35) +
        theme(legend.key.height = unit(5, "pt"),
              legend.title = element_text(size=35)))
dev.off()

print(paste('Similarity for DE genes:', mean(DESimilarity1, DESimilarity2, DESimilarity3)))
print(paste('Similarity for GSEA BP across runs:', mean(GSEASimilarity1['BP',1], GSEASimilarity2['BP',1], GSEASimilarity3['BP',1])))
print(paste('Similarity for GSEA MF across runs:', mean(GSEASimilarity1['MF',1], GSEASimilarity2['MF',1], GSEASimilarity3['MF',1])))
print(paste('Similarity for modules across runs:', mean(c(averageSimilarityModules1, averageSimilarityModules2, averageSimilarityModules3))))
print(paste('Similarity for subnetworks across runs:', mean(c(averageSimilaritySubnetworks1, averageSimilaritySubnetworks2, averageSimilaritySubnetworks3))))
print(paste('Similarity for significant subnetworks across runs:', mean(c(averageSimilaritySSubnetworks1, averageSimilaritySSubnetworks2, averageSimilaritySSubnetworks3))))
print(paste('Similarity for significant genes in subnetworks across runs:', mean(resultSimilarity1['genes',1], resultSimilarity2['genes',1], resultSimilarity3['genes',1])))
print(paste('Similarity for significant BP terms in subnetworks across runs:', mean(resultSimilarity1['BP',1], resultSimilarity2['BP',1], resultSimilarity3['BP',1])))
print(paste('Similarity for significant MF terms in subnetworks across runs:', mean(resultSimilarity1['MF',1], resultSimilarity2['MF',1], resultSimilarity3['MF',1])))

stopCluster(cluster)