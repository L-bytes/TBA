####### Main workflow in the paper
# If you only want to run part of the workflow, you can load the required data from the previous steps as indicated in each section below (load(...))

#######################################
###      Load required packages     ###
#######################################
library(ibb)
library(stringr)
library(BiocManager)
library(BioNet)
library(pheatmap)
library(WGCNA)
library(countdata)
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
setwd('C:/Users/marle/Desktop/Y2/Internship/FTD_PRODIA-main2/FTD_PRODIA-main') #Replace with your working directory

### Create output directories
dir.create('output')
dir.create('figures')
dir.create('rdata')

### Load data
d <- read.table('data/proteins.tsv', header=TRUE, sep='\t', quote='')
meta.data <- read.table('data/metadata.tsv', header=TRUE, sep='\t', quote='')

# Get only the columns with expression values
d.raw <- d[,5:ncol(d)]
d.raw <- apply(d.raw, c(1,2), as.numeric)
# Get uniprot identifiers of the peptides
ids <- d$uniprot.identifier
ids <- as.character(unlist(lapply(ids, function(x){strsplit(as.character(x), ';')[[1]][1]})))
ids <- as.character(unlist(lapply(ids, function(x){strsplit(as.character(x), '-')[[1]][1]})))
rownames(d.raw) <- ids
# Meta data rename groups (TAU, TDP and control)
meta.data$group <- unlist(lapply(meta.data$group, function(x){strsplit(as.character(x), ' ')[[1]][1]}))
groups <- meta.data$group
groups <- replace(groups, groups=='NDC', 'control')
groups <- replace(groups, groups=='FTD_TAU', 'TAU')
groups <- replace(groups, groups=='FTD_TDP', 'TDP')
colnames(d.raw) <- groups
#colnames(d.raw) <- paste(groups, meta.data$sample, sep='_')
# Save input data
save(d.raw, meta.data, groups, ids, file='rdata/input_data.RData')


#######################################
###         Preprocessing           ###
#######################################
#load(file='rdata/input_data.RData')
d.norm <- normalize.sample(d.raw)
d.cs <- normalize.cs(d.norm)
# Unique colnames in d.cs are necessary for clustering
colnames(d.cs) <- paste(groups, 1:26, sep='_')
# Save input data into a file
save(d.norm, d.cs, ids, groups, file='rdata/normalized_data.RData')


#######################################
### Differential expression analysis###
#######################################
### Test significance of differential expression between the groups
#load(file='rdata/normalized_data.RData')
dir.create('output/differential_expression')
# Beta binomial test
bb.TAU.con <- betaBinomial(d.norm, ids, groups, 'control', 'TAU', 'two.sided')
bb.TDP.con <- betaBinomial(d.norm, ids, groups, 'control', 'TDP', 'two.sided')
bb.TAU.TDP <- betaBinomial(d.norm, ids, groups, 'TDP', 'TAU', 'two.sided')
d.sign <- data.frame(log2FC_TAU_control=bb.TAU.con$table[ids,]$Log2ratio, Pvalue_TAU_control=bb.TAU.con$table[ids,]$Pvalue,
                     log2FC_TDP_control=bb.TDP.con$table[ids,]$Log2ratio, Pvalue_TDP_control=bb.TDP.con$table[ids,]$Pvalue,
                     log2FC_TAU_TDP=bb.TAU.TDP$table[ids,]$Log2ratio, Pvalue_TAU_TDP=bb.TAU.TDP$table[ids,]$Pvalue)
rownames(d.sign) <- ids
# Merge stats about P-value thresholds and number of significant proteins
sign.table <- cbind(bb.TAU.con$FDRs, bb.TDP.con$FDRs, bb.TAU.TDP$FDRs)
sign.table <- sign.table[,c(1,3,5,2,4,6)]
# Add additional info from the UniProt database to the log2FC and P-values
d.uniprot <- read.table('data/uniprot_stats_9606.tab', header=TRUE, sep='\t', quote="", fill=TRUE, na.strings="")
rownames(d.uniprot) <- d.uniprot$Entry
d.uniprot.sub <- d.uniprot[ids,]
#s.colnames <- c('Uniprot', 'Entrez', 'log2FC_TAU_control', 'Pvalue_TAU_control', 'log2FC_TAU_TDP', 'Pvalue_TAU_TDP', 'log2FC_TDP_control',
#                'Pvalue_TDP_control', 'Entry.name', 'Gene.names', 'Protein.names', 'Status', 'Organism', 'Length', 'Keywords',
#                'Keyword.ID', 'Annotation', 'Features', 'Pathway', 'Gene.ontology.IDs', 'Gene.ontology..biological.process.',
#                'Gene.ontology..cellular.component.', 'Gene.ontology..molecular.function.', 'Function..CC')
#d.summary <- add_uniprot(d.sign, d.uniprot, s.colnames)
d.summary <- cbind(d.sign, d.uniprot.sub)
# Write output into files
write.table(bb.TAU.con$table, file='output/differential_expression/sign_test_TAU_control.tsv',
            row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
write.table(bb.TAU.TDP$table, file='output/differential_expression/sign_test_TAU_TDP.tsv',
            row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
write.table(bb.TDP.con$table, file='output/differential_expression/sign_test_TDP_control.tsv',
            row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
write.table(sign.table, file='output/differential_expression/FDR_pval_summary.tsv', sep='\t', quote=FALSE)
write.table(d.summary, file='output/differential_expression/log2_pval_uniprot.tsv', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
save(bb.TAU.con, bb.TAU.TDP, bb.TDP.con, sign.table, d.summary, file='rdata/differential_expression.RData')
# Write input format for GSEA (Note GSEA is run as a separate tool, this is not part of the R workflow)
#d.TAU.con <- data.frame(Uniprot=bb.TAU.con$table$Uniprot, diff=bb.TAU.con$table$Log2ratio)
#write.table(d.TAU.con[order(d.TAU.con$diff, decreasing=TRUE),], file='Output/GSEA/diff_TAU_control_log2.rnk',
#            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
#d.TAU.TDP <- data.frame(Uniprot=bb.TAU.TDP$table$Uniprot, diff=bb.TAU.TDP$table$Log2ratio)
#write.table(d.TAU.TDP[order(d.TAU.TDP$diff, decreasing=TRUE),], file='Output/GSEA/diff_TAU_TDP_log2.rnk',
#            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
#d.TDP.con <- data.frame(Uniprot=bb.TDP.con$table$Uniprot, diff=bb.TDP.con$table$Log2ratio)
#write.table(d.TDP.con[order(d.TDP.con$diff, decreasing=TRUE),], file='Output/GSEA/diff_TDP_control_log2.rnk',
#            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

#######################################
###          Clustering             ###
#######################################
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
dir.create('figures/clustering')
# Set colors for the different significance levels
ann_colors <- list(group=c(TAU='#F39C12', TDP='#8E44AD', control='#73B761'),
                   TAU_control=c(FDR1='#000066', FDR2='#0000CC', FDR3='#3399FF', FDR4='#99CCFF', insignificant='#CCCCCC'),
                   TDP_control=c(FDR1='#000066', FDR2='#0000CC', FDR3='#3399FF', FDR4='#99CCFF', insignificant='#CCCCCC'),
                   TAU_TDP=c(FDR1='#000066', FDR2='#0000CC', FDR3='#3399FF', FDR4='#99CCFF', insignificant='#CCCCCC'))
names(ann_colors$TAU_control) <- c('FDR = 0.01', 'FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')
names(ann_colors$TDP_control) <- c('FDR = 0.01', 'FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')
names(ann_colors$TAU_TDP) <- c('FDR = 0.01', 'FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')
# Column group annotations (match to the column names of d.cs)
ann_col <- data.frame(group=as.factor(groups))
rownames(ann_col) <- paste(groups, 1:26, sep='_')


# Significance for each of the individual comparisons
# TAU vs control
thresholds1 <- sign.table[, 'TAU_vs_control_Pvalue']
names(thresholds1) <- rownames(sign.table)
sign.TAUcon <- add.sign.labels(bb.TAU.con, thresholds1)
ann_row1 <- data.frame(TAU_control=factor(sign.TAUcon, labels=c('FDR = 0.01', 'FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')))

# TDP vs control
thresholds2 <- sign.table[, 'TDP_vs_control_Pvalue']
names(thresholds2) <- rownames(sign.table)
sign.TDPcon <- add.sign.labels(bb.TDP.con, thresholds2)
ann_row2 <- data.frame(TDP_control=factor(sign.TDPcon, labels=c('FDR = 0.1', 'FDR = 0.2', 'insignificant')))

# TAU vs TDP
thresholds3 <- sign.table[, 'TAU_vs_TDP_Pvalue']
names(thresholds3) <- rownames(sign.table)
sign.TAUTDP <- add.sign.labels(bb.TAU.TDP, thresholds3)
ann_row3 <- data.frame(TAU_TDP=factor(sign.TAUTDP, labels=c('FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')))

# Top 50 most significant in each differential expression analysis, combined
sign.rows1 <- bb.TAU.con$table$Uniprot[bb.TAU.con$table$Pvalue <= sort(bb.TAU.con$table$Pvalue)[50]]
sign.rows2 <- bb.TDP.con$table$Uniprot[bb.TDP.con$table$Pvalue <= sort(bb.TDP.con$table$Pvalue)[50]]
sign.rows3 <- bb.TAU.TDP$table$Uniprot[bb.TAU.TDP$table$Pvalue <= sort(bb.TAU.TDP$table$Pvalue)[50]]
ids.sign.all <- ids[ids %in% sign.rows1 | ids %in% sign.rows2 | ids %in% sign.rows3]
ann_row4 <- data.frame(TAU_control=factor(sign.TAUcon[ids.sign.all], labels=c('FDR = 0.01', 'FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')),
                      TDP_control=factor(sign.TDPcon[ids.sign.all], labels=c('FDR = 0.1', 'FDR = 0.2', 'insignificant')),
                      TAU_TDP=factor(sign.TAUTDP[ids.sign.all], labels=c('FDR = 0.05', 'FDR = 0.1', 'FDR = 0.2', 'insignificant')))
rownames(ann_row4) <- ids.sign.all


# Create clustering heatmaps
h.all <- Clustering(d.cs, ann_col=ann_col, anncol=ann_colors, f='figures/clustering/clustering_all.png')
h.TAUcon <- Clustering(d.cs, ann_col=ann_col, ann_row=ann_row1, anncol=ann_colors, f='figures/clustering/clustering_TAUcon_all.png')
h.TAUcon.sign <- Clustering(d.cs[sign.TAUcon[rownames(d.cs)]<5,], ann_col=ann_col, ann_row=data.frame(TAU_control=ann_row1[sign.TAUcon<5,], row.names=rownames(ann_row1)[sign.TAUcon<5]),
                            anncol=ann_colors, f='figures/clustering/clustering_TAUcon_sign.png')
h.TDPcon <- Clustering(d.cs, ann_col=ann_col, ann_row=ann_row2, anncol=ann_colors, f='figures/clustering/clustering_TDPcon_all.png')
h.TDPcon.sign <- Clustering(d.cs[sign.TDPcon[rownames(d.cs)]<5,], ann_col=ann_col, ann_row=data.frame(TDP_control=ann_row2[sign.TAUcon<5,], row.names=rownames(ann_row2)[sign.TAUcon<5]),
                            anncol=ann_colors, f='figures/clustering/clustering_TDPcon_sign.png')
h.TAUTDP <- Clustering(d.cs, ann_col=ann_col, ann_row=ann_row3, anncol=ann_colors, f='figures/clustering/clustering_TAUTDP_all.png')
h.TAUTDP.sign <- Clustering(d.cs[sign.TAUTDP[rownames(d.cs)]<5,], ann_col=ann_col, ann_row=data.frame(TAU_TDP=ann_row3[sign.TAUcon<5,], row.names=rownames(ann_row3)[sign.TAUcon<5]),
                            anncol=ann_colors, f='figures/clustering/clustering_TAUTDP_sign.png')
h.top50.sign <- Clustering(d.cs[ids.sign.all,], ann_col=ann_col, ann_row=ann_row4,
                            anncol=ann_colors, f='figures/clustering/clustering_top50_sign.png')

#######################################
###   WGCNA coexpression analysis   ###
#######################################
#load(file='rdata/input_data.RData')
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#source('src/coexpression_analysis_prep.R')
dir.create('figures/coexpression')
dir.create('output/coexpression')

d.log2vals <- d.summary[, c('log2FC_TAU_control', 'log2FC_TDP_control', 'log2FC_TAU_TDP')]
coexpression <- coexpression.analysis(t(d.norm), d.log2vals, 'output/coexpression', 'figures/coexpression')
wgcna.net <- coexpression[[1]]
module.significance <- coexpression[[2]]
# Merge M18 (lightgreen) into M9 (magenta), since they were highly similar
#wgcna.net$colors <- replace(wgcna.net$colors, wgcna.net$colors==18, 9)
gene.names <- unlist(lapply(d.summary$Gene.names, function(x){strsplit(x, split=' ')[[1]][1]}))
names(gene.names) <- ids
entrez.ids <- unlist(lapply(d.summary$Cross.reference..GeneID., function(x){strsplit(x, split=';')[[1]][1]}))
names(entrez.ids) <- ids
log2vals <- d.summary$log2FC_TAU_control
names(log2vals) <- ids
GOenr.net <- GO.terms.modules(wgcna.net, ids, log2vals, gene.names, entrez.ids, 'figures/coexpression', 'output/coexpression')
# Save data structures
save(wgcna.net, GOenr.net, module.significance, gene.names, file='rdata/coexpression.RData')

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

# Get TOM matrix
TOM <- TOMsimilarityFromExpr(t(d.norm), power=7)

# Module labels and log2 values
moduleColors <- labels2colors(wgcna.net$colors)
modules <- levels(as.factor(moduleColors))
log2.TAU.con <- d.summary$log2FC_TAU_control
P.TAU.con <- d.summary$Pvalue_TAU_control

# Get table with interactions for each module
for (module in modules){
  print(module)
  
  # Select proteins in module
  inModule <- moduleColors == module
  sum(inModule)
  TOMmodule <- TOM[inModule, inModule]
  idsModule <- ids[inModule]
  log2ModuleProteins <- log2.TAU.con[inModule]
  PModuleProteins <- P.TAU.con[inModule]
  
  # Convert Uniprot IDs to Gene symbols
  namesModule <- unlist(lapply(d.summary[idsModule,'Gene.names'], function(x){strsplit(x, split=' ')[[1]][1]}))
  for (i in 1:length(namesModule)){
    # If there is no gene name associated with the Uniprot ID, keep the uniprot ID
    if (is.na(namesModule[i])){
      namesModule[i] <- idsModule[i]
    }
  }
  node.frame <- data.frame(ID=idsModule, LogFC=log2ModuleProteins, Pvalue=PModuleProteins, Symbol=namesModule)
  rownames(node.frame) <- 1:nrow(node.frame)
  write.table(node.frame, file=paste0('output/hotnet/nodes_', module, '.tsv'), col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
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

system('bash src/run_hierarchicalHotnet_modules.sh')


#######################################
###          Validation             ###
#######################################
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#load(file='rdata/coexpression.RData')
dir.create('output/validation')
dir.create('figures/validation')
dir.create('rdata/validation')

# Load RIMOD data
d.val1 <- read.table('data/frontal_TAU_con.tsv', sep='\t', header=TRUE, quote="")
d.val1.val <- apply(d.val1[,4:ncol(d.val1)], c(1,2), as.numeric)
d.val2 <- read.table('data/temporal_TAU_con.tsv', sep='\t', header=TRUE, quote="")
d.val2.val <- apply(d.val2[,4:ncol(d.val2)], c(1,2), as.numeric)
meta.val <- read.table('data/metadata_validation.txt', sep='\t', header=TRUE, quote="")
rownames(meta.val) <- meta.val$Gel.Code

ids1 <- d.val1$UniProt_accession
ids1 <- as.character(unlist(lapply(ids1, function(x){strsplit(as.character(x), ';')[[1]][1]})))
ids1 <- as.character(unlist(lapply(ids1, function(x){strsplit(as.character(x), '-')[[1]][1]})))
ids2 <- d.val2$UniProt_accession
ids2 <- as.character(unlist(lapply(ids2, function(x){strsplit(as.character(x), ';')[[1]][1]})))
ids2 <- as.character(unlist(lapply(ids2, function(x){strsplit(as.character(x), '-')[[1]][1]}))) 

# Normalize data
d.val1.norm <- normalize.sample(d.val1.val)
d.val2.norm <- normalize.sample(d.val2.val)

# Beta binomial test on validation sets
groups1 <- meta.val[colnames(d.val1.norm),"Condition"]
groups2 <- meta.val[colnames(d.val2.norm),"Condition"]
bb.1 <- betaBinomial(d.val1.norm, ids1, groups1, 'control', 'TAU', 'two.sided')
bb.2 <- betaBinomial(d.val2.norm, ids2, groups2, 'control', 'TAU', 'two.sided')

# Write results
results1.tab <- data.frame(Uniprot=ids1, GeneName=d.val1$UniProt_gene_symbol, Log2FC=bb.1$table[ids1, 'Log2ratio'], Pvalue=bb.1$table[ids1, 'Pvalue'])
results2.tab <- data.frame(Uniprot=ids2, GeneName=d.val2$UniProt_gene_symbol, Log2FC=bb.2$table[ids2, 'Log2ratio'], Pvalue=bb.2$table[ids2, 'Pvalue'])
rownames(results1.tab) <- ids1
rownames(results2.tab) <- ids2

write.table(results1.tab, 'output/validation/log2_Pval_val_frontal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(results2.tab, 'output/validation/log2_Pval_val_temporal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(bb.1$FDRs, 'output/validation/sign_test_val_frontal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
write.table(bb.2$FDRs, 'output/validation/sign_test_val_temporal.tsv', sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)

# Merge significance info with UniProt stats (if you run the workflow from the start the UniProt data is already loaded)
d.uniprot <- read.table('data/uniprot_stats_9606.tab', header=TRUE, sep='\t', quote="", fill=TRUE, na.strings="")
rownames(d.uniprot) <- d.uniprot$Entry
d.uniprot.sub1 <- d.uniprot[ids1,]
d.summary1 <- cbind(results1.tab, d.uniprot.sub1)
d.uniprot.sub2 <- d.uniprot[ids2,]
d.summary2 <- cbind(results2.tab, d.uniprot.sub2)

# Expression values
ids.TAUcon <- rownames(d.summary)
ids.valFrn <- rownames(d.summary1)
ids.valTem <- rownames(d.summary2)
log.TAUcon <- d.summary$log2FC_TAU_control
P.TAUcon <- d.summary$Pvalue_TAU_control
log.valFrn <- d.summary1$Log2FC
P.valFrn <- d.summary1$Pvalue
log.valTem <- d.summary2$Log2FC
P.valTem <- d.summary2$Pvalue

names(log.TAUcon) <- ids.TAUcon
names(P.TAUcon) <- ids.TAUcon
names(log.valFrn) <- ids.valFrn
names(P.valFrn) <- ids.valFrn
names(log.valTem) <- ids.valTem
names(P.valTem) <- ids.valTem
names(ids.TAUcon) <- gene.names

# Get module IDs
moduleColors <- labels2colors(wgcna.net$colors)
modules <- levels(as.factor(moduleColors))
names(moduleColors) <- ids.TAUcon

# Test against frontal dataset
# Filter out proteins that are found in one of the two sets (Note that proteins in FC2/P2 that are not in FC1/P1 are
# automatically ignored in fraction.correct, so there is no need to filter that dataset)
ids.both.frn <- ids.TAUcon %in% ids.valFrn
log.TAUcon.toTest.frn <- log.TAUcon[ids.both.frn]
P.TAUcon.toTest.frn <- P.TAUcon[ids.both.frn]

# Permutation test on all proteins (frontal)
outfolder <- 'output/validation/all_frontal/'
figfolder <- 'figures/validation/all_frontal/'
dir.create(outfolder, showWarnings=FALSE)
dir.create(figfolder, showWarnings=FALSE)
perm.allFrn <- permutation.test(log.TAUcon.toTest.frn, log.valFrn, P.TAUcon.toTest.frn, P.valFrn, outfolder, figfolder)
save(perm.allFrn, file='rdata/validation/permutation_all_frontal.RData')

# Permutation test per module (frontal)
moduleColorsBoth.frn <- moduleColors[ids.both.frn]
for(m in modules){
  print(m)
  inModule <- moduleColorsBoth.frn == m
  m.ids <- names(moduleColorsBoth.frn[inModule])
  #log.TAUcon.module <- log.TAUcon.toTest[inModule]
  #P.TAUcon.module <- P.TAUcon.toTest[inModule]
  outfolder <- paste0('output/validation/module_', m, '_frontal/')
  figfolder <- paste0('figures/validation/module_', m, '_frontal/')
  dir.create(outfolder, showWarnings=FALSE)
  dir.create(figfolder, showWarnings=FALSE)
  perm.module <- permutation.test(log.TAUcon, log.valFrn, P.TAUcon, P.valFrn, outfolder, figfolder, m.ids=m.ids)
  save(perm.module, file=paste0('rdata/validation/permutation_', m, '_frontal.RData'))
  
  # Do the same, but only for those ids that were significant in HotNet
  #if (m == 'red'){
  #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_001_', m, '.tsv'), header=FALSE, quote='', sep='\t')
  #} else if (m == 'grey'){
  #  next
  #} else {
  #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_003_', m, '.tsv'), header=FALSE, quote='', sep='\t')
  #}
  
  #names.hotnet <- levels(as.factor(c(d.hotnet[,1], d.hotnet[,2])))
  #ids.hotnet <- as.character(ids.TAUcon.named[names.hotnet])
  #ids.both <- ids.TAUcon %in% ids.valFrn & ids.TAUcon %in% ids.hotnet
  #inModule <- moduleColors[ids.both] == m
  #m.ids <- names(moduleColors[ids.both][inModule])
  #log.TAUcon.module <- log.TAUcon[ids.both]
  #P.TAUcon.module <- P.TAUcon[ids.both]
  #outfolder <- paste0('Figures/Validation4/Module_', m, '_frontal_hotnet/')
  #dir.create(outfolder, showWarnings=FALSE)
  #perm.module.hotnet <- permutation.test(log.TAUcon, log.valFrn, P.TAUcon, P.valFrn, outfolder, m.ids=m.ids)
  #save(perm.module.hotnet, file=paste0('Output/Permutation_', m, '_frontal4_hotnet.RData'))
}


# Test against temporal dataset
# Filter out proteins that are found in one of the two sets
ids.both.tem <- ids.TAUcon %in% ids.valTem
log.TAUcon.toTest.tem <- log.TAUcon[ids.both.tem]
P.TAUcon.toTest.tem <- P.TAUcon[ids.both.tem]
# Permutation test on all proteins
outfolder <- 'output/validation/all_temporal/'
figfolder <- 'figures/validation/all_temporal/'
dir.create(outfolder, showWarnings=FALSE)
dir.create(figfolder, showWarnings=FALSE)
perm.allTem <- permutation.test(log.TAUcon.toTest.tem, log.valTem, P.TAUcon.toTest.tem, P.valTem, outfolder, figfolder)
save(perm.allTem, file='rdata/validation/permutation_all_temporal.RData')

# Permutation test per module
moduleColorsBoth.tem <- moduleColors[ids.both.tem]
for(m in modules){
  print(m)
  inModule <- moduleColorsBoth.tem == m
  #log.TAUcon.module <- log.TAUcon.toTest[inModule]
  #P.TAUcon.module <- P.TAUcon.toTest[inModule]
  m.ids <- names(moduleColorsBoth.tem[inModule])
  outfolder <- paste0('output/validation/module_', m, '_temporal/')
  figfolder <- paste0('figures/validation/module_', m, '_temporal/')
  dir.create(outfolder, showWarnings=FALSE)
  dir.create(figfolder, showWarnings=FALSE)
  perm.module <- permutation.test(log.TAUcon, log.valTem, P.TAUcon, P.valTem, outfolder, figfolder, m.ids=m.ids)
  save(perm.module, file=paste0('rdata/validation/permutation_', m, '_temporal.RData'))
  
  ## Do the same, but only for those ids that were significant in HotNet
  #if (m == 'red'){
  #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_001_', m, '.tsv'), header=FALSE, quote='', sep='\t')
  #} else if (m == 'grey') {
  #  # Skip empty files
  #  next
  #} else {
  #  d.hotnet <- read.table(paste0('C:/Users/jhmva/surfdrive/PhD/ftdproject_related_stuff/HotNet/Modules/HotNet_results_100/consensus_edges_Pval_003_', m, '.tsv'), header=FALSE, quote='', sep='\t')
  #}
  
  #names.hotnet <- levels(as.factor(c(d.hotnet[,1], d.hotnet[,2])))
  #ids.hotnet <- as.character(ids.TAUcon.named[names.hotnet])
  #ids.both <- ids.TAUcon %in% ids.valTem & ids.TAUcon %in% ids.hotnet
  #inModule <- moduleColors[ids.both] == m
  ##log.TAUcon.module <- log.TAUcon[ids.both]
  ##P.TAUcon.module <- P.TAUcon[ids.both]
  #m.ids <- names(moduleColors[ids.both][inModule])
  #outfolder <- paste0('Figures/Validation4/Module_', m, '_temporal_hotnet/')
  #dir.create(outfolder, showWarnings=FALSE)
  #perm.module.hotnet <- permutation.test(log.TAUcon, log.valTem, P.TAUcon, P.valTem, outfolder, m.ids=m.ids)
  #save(perm.module.hotnet, file=paste0('Output/Permutation_', m, '_temporal4_hotnet.RData'))
}

#######################################
###       Create module plots       ###
#######################################
#load(file='rdata/normalized_data.RData')
#load(file='rdata/differential_expression.RData')
#load(file='rdata/coexpression.RData')
d.norm <- t(d.norm)
moduleColors <- labels2colors(wgcna.net$colors)

d.array <- data.frame(matrix(nrow=5, ncol=18))
colnames(d.array) <- labels2colors(0:17)
rownames(d.array) <- c('FTLD-tau vs NHC', 'FTLD-TDP vs NHC', 'FTLD-tau va FTLD-TDP', 'FTLD-tau vs NCH (validation frontal)', 'FTLD-tau vs NHC (validation temporal)')
for (i in 0:17){
  color <- labels2colors(i)
  TAUcon <- d.summary$log2FC_TAU_control[moduleColors == color]
  TDPcon <- d.summary$log2FC_TDP_control[moduleColors == color]
  TAUTDP <- d.summary$log2FC_TAU_TDP[moduleColors == color]
  
  frn <- log.TAUcon.toTest.frn[moduleColorsBoth.frn==color]
  tem <- log.TAUcon.toTest.tem[moduleColorsBoth.tem==color]
  
  mTAUcon <- median(TAUcon)
  mTDPcon <- median(TDPcon)
  mTAUTDP <- median(TAUTDP)
  mfrn <- median(frn)
  mtem <- median(tem)
  d.array[,color] <- c(mTAUcon, mTDPcon, mTAUTDP, mfrn, mtem)
  
  
  t.frn <- t.test(tem)
  t.tem <- t.test(frn)
}
write.table(d.array, file='module_heatmap.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(t(pval.modules.norm), file='module_pval.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)