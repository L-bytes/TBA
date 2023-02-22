####### Main workflow in the paper
# If you only want to run part of the workflow, you can load the required data from the previous steps as indicated in each section below (load(...))
### Set working directory to scratch on LISA
argv = commandArgs(trailingOnly=TRUE)
#######################################
###      Load required packages     ###
#######################################
library(ibb)
library(stringr)
library(BiocManager)
library(BioNet)
library(pheatmap)
library(countdata)
library(DESeq2)
library('org.Hs.eg.db')
library(WGCNA)
options(stringsAsFactors=FALSE)
# Load required functions from the src directory
source(paste(argv[2],'/src/preprocessing.R', sep=""))
source(paste(argv[2],'/src/differential_expression.R', sep=""))
#source('src/add_uniprot_info.R')
source(paste(argv[2],'/src/clustering.R', sep=""))
source(paste(argv[2],'/src/coexpression_analysis.R', sep=""))
source(paste(argv[2],'/src/validation.R', sep=""))


#######################################
###             Setup               ###
#######################################
### Set working directory to scratch on LISA
setwd(argv[1]) #Replace with your working directory

### Create output directories
dir.create('output')
dir.create('figures')
dir.create('rdata')

### Load data
d <- read.table('./data/TCGA_rna_count_data.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)
group.data <- read.table('./data/non_silent_mutation_profile_crc.txt', header=TRUE, sep='\t', quote="", row.names = 1, check.names=FALSE)

# Get only the columns with expression values

d <- d[1800:3600,1:10]
d <- d[,colnames(d) %in% rownames(group.data)]
group.data <- group.data[rownames(group.data) %in% colnames(d),]

new_order <- sort(colnames(d))
d <- d[, new_order]

new_order_rows <- sort(rownames(d))
d <- d[new_order_rows,]

d.raw <- d[,1:ncol(d)]
d.raw <- apply(d.raw, c(1,2), as.numeric)
# Get identifiers of genes
ids <- rownames(d.raw)
groups <- group.data$KRAS
# colnames(d.raw) <- groups
count.groups <- d.raw
colnames(count.groups) <- groups
# Save input data
save(d.raw, group.data, groups, ids, file='rdata/input_data.RData')


#######################################
###         Preprocessing           ###
#######################################
#load(file='rdata/input_data.rdata')
d.norm <- normalize.sample(d.raw)
d.cs <- normalize.cs(d.norm)
length <- length(groups)
# unique colnames in d.cs are necessary for clustering
colnames(d.cs) <- paste(groups, 1:length, sep='_')
# save input data into a file
save(d.norm, d.cs, ids, groups, file='rdata/normalized_data.rdata')

#DE-SEQ
group.KRAS <- as.matrix(group.data[,"KRAS"])
rownames(group.KRAS) <- rownames(group.data)
colnames(group.KRAS) <- c("KRAS")
group.KRAS[,1] <- as.numeric(factor(group.KRAS[,1]))
dds <- DESeqDataSetFromMatrix(countData = d.raw, colData = group.KRAS, design = ~ KRAS)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds)
results.names <- resultsNames(dds)
res
results.names
colnames(d.raw) <- groups

#Normalization
dds2 <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds2, blind=F)
d.adj <- counts(dds2, normalized=TRUE)
colnames(d.adj) <- groups
ids <- rownames(d.adj)

d.summary <- data.frame(unlist(res$log2FoldChange), unlist(res$padj), row.names = row.names(res))
colnames(d.summary) <- c("log2FC", "Pvalue")

# resLFC <- lfcShrink(dds, coef="KRAS_2_vs_1", type="apeglm")
# resLFC

#Some plots
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]
# plotCounts(dds, gene=which.min(res$padj), intgroup="KRAS")
# ntd <- normTransform
# BiocManager::install("vsn")
# library("vsn")
# meanSdPlot(assay(ntd))
# library("pheatmap")
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)
# df <- as.data.frame(colData(dds)[,c("KRAS")])
# rownames(df) <- colnames(assay(ntd)[select,])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
# plotDispEsts(dds)

#######################################
### Differential expression analysis###
#######################################
### Test significance of differential expression between the groups
# load(file='rdata/normalized_data.RData')
# dir.create('output/differential_expression')
# # Beta binomial test
# bb <- betaBinomial(d.norm, ids, groups, 'WT', 'SNV', 'two.sided')
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

d.log2vals <- as.data.frame(d.summary[, c('log2FC')], row.names=row.names(d.summary))
colnames(d.log2vals) <- c('log2FC')
coexpression <- coexpression.analysis(t(d.adj), d.log2vals, 'output/coexpression', 'figures/coexpression')
wgcna.net <- coexpression[[1]]
module.significance <- coexpression[[2]]
# Merge M18 (lightgreen) into M9 (magenta), since they were highly similar
#wgcna.net$colors <- replace(wgcna.net$colors, wgcna.net$colors==18, 9)
log2vals <- d.summary$log2FC
names(log2vals) <- ids
gns<-mapIds(org.Hs.eg.db, keys = ids, column = "ENTREZID", keytype = "SYMBOL")
GOenr.net <- GO.terms.modules(wgcna.net, ids, log2vals, gns, 'figures/coexpression', 'output/coexpression')
# Save data structures
save(wgcna.net, GOenr.net, module.significance, ids, file='rdata/coexpression.RData')

#######################################
###       Hierarchical HotNet       ###
#######################################
### Run hierarchical hotnet to obtain the most significant submodule within each of the identified modules
### Note that the HotNet package was written in Python and needs to be intstalled separately on your machine
### HotNet was run mostly with default parameters
#load(file='rdata/input_data.RData')
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#load(file='rdata/coexpression.RData')
dir.create('figures/hotnet')
dir.create('output/hotnet')
dir.create('output/hotnet/HotNet_input')

# Get TOM matrix
TOM <- TOMsimilarityFromExpr(t(d.adj), power=2)

# Module labels and log2 values
moduleColors <- labels2colors(wgcna.net$colors)
modules <- levels(as.factor(moduleColors))
log2.SNV.WT <- d.summary$log2FC
P.SNV.WT <- d.summary$Pvalue

# Get table with interactions for each module
for (module in modules){
  print(module)
  
  # Select proteins in module
  inModule <- moduleColors == module
  sum(inModule)
  TOMmodule <- TOM[inModule, inModule]
  idsModule <- ids[inModule]
  log2ModuleProteins <- log2.SNV.WT[inModule]
  PModuleProteins <- P.SNV.WT[inModule]
  
  # Convert Uniprot IDs to Gene symbols
  namesModule <- idsModule
  for (i in 1:length(namesModule)){
    # If there is no gene name associated with the Uniprot ID, keep the uniprot ID
    if (is.na(namesModule[i])){
      namesModule[i] <- idsModule[i]
    }
  }
  node.frame <- data.frame(Symbol=idsModule, LogFC=log2ModuleProteins, Pvalue=PModuleProteins)
  rownames(node.frame) <- 1:nrow(node.frame)
  node.frame$InvPvalue <- 1 - PModuleProteins
  write.table(node.frame, file=paste0('output/hotnet/nodes_', module, '.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
  write.table(node.frame[, c('Symbol', 'InvPvalue')], file=paste0('output/hotnet/HotNet_input/g2s_Pval_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
  names(namesModule) <- idsModule
  
  write.table(node.frame[, c('Symbol', 'LogFC')], file=paste0('output/hotnet/HotNet_input/g2s_log2_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
  
  #nodeColorsModule <- node.colors[inModule]
  
  # Create empty dataframes
  edges.matrix.1 <- data.frame(matrix(nrow=0, ncol=2))
  edges.matrix.2 <- data.frame(matrix(nrow=0, ncol=2))
  edges.matrix.1.num <- data.frame(matrix(nrow=0, ncol=2))
  edges.matrix.2.num <- data.frame(matrix(nrow=0, ncol=2))
  i2g.1 <- data.frame(matrix(nrow=0, ncol=2))
  i2g.2 <- data.frame(matrix(nrow=0, ncol=2))
  
  # Write tables: one with edges between all nodes, one with a treshold of 0.05 and one with custom thresholds
  for (i in 1:(nrow(TOMmodule)-1)){
    #print(i)
    for (j in (i+1):nrow(TOMmodule)){
      if (TOMmodule[i,j] > 0.03){
        # Add edge to list
        edges.matrix.1[nrow(edges.matrix.1)+1,] <- c(namesModule[i],namesModule[j])
        edges.matrix.1.num[nrow(edges.matrix.1.num)+1,] <- c(i,j)
        # Add node to node list
        if (!i %in% i2g.1[,1]){
          i2g.1[nrow(i2g.1)+1,] <- c(i, namesModule[i])
        }
        if (!j %in% i2g.1[,1]){
          i2g.1[nrow(i2g.1)+1,] <- c(j, namesModule[j])
        }
      }
      if (TOMmodule[i,j] > 0.01){
        # Add edge to list
        edges.matrix.2[nrow(edges.matrix.2)+1,] <- c(namesModule[i],namesModule[j])
        edges.matrix.2.num[nrow(edges.matrix.2.num)+1,] <- c(i,j)
        # Add node to node list
        if (!i %in% i2g.2[,1]){
          i2g.2[nrow(i2g.2)+1,] <- c(i, namesModule[i])
        }
        if (!j %in% i2g.2[,1]){
          i2g.2[nrow(i2g.2)+1,] <- c(j, namesModule[j])
        }
      }
    }
  }
  write.table(edges.matrix.1, file=paste0('output/hotnet/HotNet_input/name_edges_expression_003_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
  write.table(edges.matrix.2, file=paste0('output/hotnet/HotNet_input/name_edges_expression_001_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
  write.table(edges.matrix.1.num, file=paste0('output/hotnet/HotNet_input/edge_list_003_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
  write.table(edges.matrix.2.num, file=paste0('output/hotnet/HotNet_input/edge_list_001_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
  write.table(i2g.1, file=paste0('output/hotnet/HotNet_input/i2g_003_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
  write.table(i2g.2, file=paste0('output/hotnet/HotNet_input/i2g_001_', module, '.tsv'), col.names=FALSE, row.names=FALSE, sep='\t', quote=FALSE)
}


######### Run hierarchical hotnet to obtain the most significant submodule within each of the identified modules
######### Note that the HotNet package was written in Python and needs to be intstalled separately on your machine
system(paste('bash', paste(argv[2],'/src/run_hierarchicalHotnet_modules.sh', sep="")))
system('module load R')

#######################################
###          Validation             ###
#######################################
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#load(file='rdata/coexpression.RData')
# dir.create('output/validation')
# dir.create('figures/validation')
# dir.create('rdata/validation')
# 
# # Load RIMOD data
# d.val1 <- read.table('data/frontal_TAU_con.tsv', sep='\t', header=TRUE, quote="")
# d.val1.val <- apply(d.val1[,4:ncol(d.val1)], c(1,2), as.numeric)
# d.val2 <- read.table('data/temporal_TAU_con.tsv', sep='\t', header=TRUE, quote="")
# d.val2.val <- apply(d.val2[,4:ncol(d.val2)], c(1,2), as.numeric)
# meta.val <- read.table('data/metadata_validation.txt', sep='\t', header=TRUE, quote="")
# rownames(meta.val) <- meta.val$Gel.Code
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