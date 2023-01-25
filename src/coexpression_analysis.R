### Coexpression analysis workflow using the WGCNA package. The WGCNA package has a tutorial itself, which also includes some
### preparatory steps to check for sample outliers etc. As these are very specific to the input data, we have not included them
### as part of the workflow here but instead assume these have been followed already. We recommend going back to the original
### documentation for the quality check (and also the general workflow, as we generally followed their guidelines):
### (https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/)
### The quality check for our dataset can be found in "coexpression_analysis_prep.R"

#' Coexpression analysis workflow using the WGCNA package
#' 
#' @param d Input a dataframe with expression values, using the Protein IDs as rows, samples as columns.
#' @param d.log2fc log2FC values of the groups
#' @param outfolder Output folder for figures
#' @return Coexpression network of \code{d} and significance of log2FC value distribution(s) within each of the modules

coexpression.analysis <- function(d, d.log2fc, outfolder, figfolder, power=FALSE){
  # Pick acceptable value for beta
  # Choose a set of soft-thresholding powers (beta)
  powers <- c(c(1:11), seq(from=12, to=26, by=2))
  # Call the network topology analysis function
  soft.tresh <- pickSoftThreshold(d, powerVector=powers, verbose=5)
  # Plot scale-free topology parameters as function of beta
  png(filename=paste0(figfolder, '/Scale_free_threshold.png'), width=960, height=480)
  par(mfrow = c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index as a function of beta
  plot(soft.tresh $fitIndices[,1], -sign(soft.tresh $fitIndices[,3])*soft.tresh $fitIndices[,2], xlab="Beta",ylab="Scale Free Topology Model Fit, R^2",type="n", main = paste("Scale independence"))
  text(soft.tresh $fitIndices[,1], -sign(soft.tresh $fitIndices[,3])*soft.tresh $fitIndices[,2], labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.80,col="red")
  # Mean connectivity as a function of beta
  plot(soft.tresh $fitIndices[,1], soft.tresh $fitIndices[,5], xlab="Beta", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
  text(soft.tresh $fitIndices[,1], soft.tresh $fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  # If power not provided, pick from the scale independence analysis
  if (!power){
    # Pick first value where scale independence is higher than 0.8
    sf.values <- -sign(soft.tresh $fitIndices[,3])*soft.tresh $fitIndices[,2]
    sf.bool <- sf.values > 0.8
    power <- soft.tresh $fitIndices[,1][sf.bool][1]
  }
  # Create network from chosen value of beta
  rownames(d) <- paste(rownames(d), 1:26, sep='_')
  coex.net <- blockwiseModules(d, power=power,
                               TOMType="unsigned", minModuleSize=15,
                               reassignThreshold=0, mergeCutHeight=0.25,
                               numericLabels=TRUE, pamRespectsDendro=FALSE,
                               saveTOMs=TRUE, saveTOMFileBase="rdata/coexpression_discovery", verbose=3)
  # Merge M18 (lightgreen) into M9 (magenta), since they were highly similar
  coex.net$colors <- replace(coex.net$colors, coex.net$colors==18, 9)
  # Plot modules and gene hierarchical clustering
  moduleColors <- labels2colors(coex.net$colors)
  modules <- levels(as.factor(moduleColors))
  n.modules <- length(modules)

  png(file=paste0(figfolder, '/module_dendrogram.png'), width=960, height=480)
  plotDendroAndColors(coex.net$dendrograms[[1]], moduleColors[coex.net$blockGenes[[1]]],
                      "Module colors", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)#,
                      #abHeight=0.85, abCol='red')
  dev.off()
  # Plot distribution of eigengene expression values for each module, and calculate significance of differential expression within the modules
  MEs <- coex.net$MEs
  ME.control <- MEs[str_detect(rownames(MEs), 'control'),]
  ME.TAU <- MEs[str_detect(rownames(MEs), 'TAU'),]
  ME.TDP <- MEs[str_detect(rownames(MEs), 'TDP'),]
  
  # Table with P-values for distribution of log2 values
  pvals.modules <- data.frame(matrix(nrow=n.modules, ncol=ncol(d.log2fc)))
  colnames(pvals.modules) <- colnames(d.log2fc)
  rownames(pvals.modules) <- modules
  
  for (i in 0:(n.modules-1)){
    module <- labels2colors(i)
    # Plot module eigengene expression
    png(paste0(figfolder, '/ME_distribution_', i, '_', module, '.png'))
    boxplot(ME.TAU[,paste0('ME', i)], ME.TDP[,paste0('ME', i)], ME.control[,paste0('ME', i)], names=c('TAU', 'TDP', 'Control'),
            col=module, xlab='group', ylab='ME value')
    dev.off()
    # Calculate if distribution of log2FC values between the groups are significantly different from zero in each module
    #pval.modules <- data.frame()
    for (j in 1:ncol(d.log2fc)){
      #logvals <- d.log2fc[names(coex.net$colors[coex.net$colors==i]), j]
      logvals <- d.log2fc[names(coex.net$colors[coex.net$colors==i]), j]
      t.stat <- t.test(logvals)
      pvals.modules[module, j] <- t.stat$p.value
      print(c(module, median(logvals), sd(logvals), length(logvals)))
    }
  }
  # Multiple testing correction
  p.adj <- p.adjust(unlist(pvals.modules))
  pval.modules.adj <- pvals.modules
  for (j in 1:ncol(pvals.modules)){
    pvals.modules[,j] <- p.adj[((j-1)*n.modules+1):(j*n.modules)]
  }
  write.table(pval.modules.adj, file=paste0(outfolder, '/sign_log2fc_modules.tsv'), sep='\t',
              row.names=TRUE, col.names=TRUE, quote=FALSE)
  return(list(coex.net, pval.modules.adj))
}

#' GO term enrichment analysis of the modules identified by WGCNA
#' 
#' @param coex.net Coexpression network generated by WGCNA
#' @param IDs Protein identifier
#' @param log2values Log2FC values for each of the proteins
#' @param names Gene/protein names
#' @param entrez.ids Entrez identifiers for the proteins
#' @param fig.folder Directory to save figure files
#' @param out.folder Output folder for data files
#' @return Object with GO enrichment statistics per module and a heatmap showing the log2FC values of the top 10 enriched GO terms

GO.terms.modules <- function(coex.net, IDs, log2values, names, entrez.ids, fig.folder, out.folder){
  # GO enrichment analysis of the modules
  #expr2ID <- match(ids, d.summary$Entry)
  allEntrezIDs <- unlist(lapply(d.summary$Cross.reference..GeneID., function(x){strsplit(x, split=';')[[1]][1]}))
  names(allEntrezIDs) <- IDs
  moduleColors <- labels2colors(coex.net$colors)
  modules <- levels(as.factor(moduleColors))
  n.modules <- length(modules)
  for (i in 0:(n.modules-1)) {
    # Save Entrez IDs per module
    module <- labels2colors(i)
    entrez.module <- entrez.ids[moduleColors==module]
    filename <- paste0(out.folder, '/EntrezIDs_', i, '_', module, '.txt')
    write.table(as.data.frame(entrez.module), file=filename, row.names=FALSE, col.names=FALSE)
  }
  # As background in the enrichment analysis, we will use all probes in the analysis.
  filename <- paste0(out.folder, '/EntrezIDs_all.txt')
  write.table(as.data.frame(allEntrezIDs), file=filename, row.names=FALSE, col.names=FALSE)
  
  # GO enrichment
  # Save 100 most significant GO terms
  GO.enrichment <- GOenrichmentAnalysis(moduleColors, entrez.ids, organism='human', nBestP=100, includeOffspring=TRUE)
  enrichment.tab <- GO.enrichment$bestPTerms[[4]]$enrichment
  write.table(enrichment.tab, file=paste0(out.folder, '/GOEnrichmentTable100.tsv'),
              sep='\t', quote=FALSE, row.names=FALSE)
  # Summary with most relevant columns only
  rel.cols <- c(1, 2, 5, 6, 7, 11, 12, 13)
  sub.enr.tab <- enrichment.tab[,rel.cols]
  write.table(sub.enr.tab, file=paste0(out.folder, '/GOEnrichmentTable100_summary.tsv'),
              sep='\t', quote=FALSE, row.names=FALSE)
  
  ### Plot heatmaps of the top ten enriched GO terms per module
  signGOTerms <- data.frame(modules=GO.enrichment$bestPTerms$`BP, CC, MF`$enrichment$module, GOTerms=GO.enrichment$bestPTerms$`BP, CC, MF`$enrichment$termID)
  
  # Color palette for the heatmap
  colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
  
  for (i in 0:(n.modules-1)){
    module <- labels2colors(i)
    #print(module)
    # Select proteins in the module
    inModule <- moduleColors == module
    # Original IDs
    idsModule <- IDs[inModule]
    # Gene names
    namesModule <- names[idsModule]
    for (j in 1:length(namesModule)){
      # If there is no gene name associated with the Uniprot ID, keep the uniprot ID
      if (is.na(namesModule[j])){
        namesModule[j] <- idsModule[j]
      }
    }
    # Log2FC values of proteins in module
    log2ModuleProteins <- log2values[idsModule]
    #names(log2ModuleProteins) <- idsModule
    # GO IDs significantly enriched in module
    GOTermsModule <- signGOTerms[signGOTerms$modules == module, 'GOTerms']
    # Entrez IDs of proteins in module
    entrezModule <- entrez.ids[idsModule]
    entrezModule.noNA <- unique(entrezModule[!is.na(entrezModule)])
    # GO mapping
    egGO <- org.Hs.egGO
    # Include the offspring for each GO term
    GO.offspr <- list(bp=as.list(GOBPOFFSPRING), cc=as.list(GOCCOFFSPRING), mf=as.list(GOMFOFFSPRING))
    # All Entrez IDs that have at least one GO term
    mappedProteins <- mappedkeys(egGO)
    # Get the indeces of the Entrez IDs in our list
    indeces <- match(entrezModule.noNA, mappedProteins)
    if (sum(is.na(indeces)) != 0) {
      entrezModule.noNA <- entrezModule.noNA[!is.na(indeces)]
      indeces <- indeces[!is.na(indeces)]
    }
    annotations <- as.list(egGO[mappedProteins[indeces]])
    
    # For each GO term, count the number of Entrez IDs with that GO term and create a matrix
    # to be plotted as a heatmap
    d.GO <- data.frame(matrix(data=NA, ncol=sum(inModule), nrow=length(GOTermsModule)))
    rownames(d.GO) <- unlist(lapply(GOTermsModule, function(x){GOTERM[[x]]@Term}))
    colnames(d.GO) <- namesModule
    GOcounts <- list()
    
    for (term in GOTermsModule) {
      GOcounts[[term]] <- 0
      # Check whether the current term or one of the direct children of the term is associated with the Entrez ID
      children <- c(GO.offspr$bp[[term]], GO.offspr$cc[[term]], GO.offspr$mf[[term]])
      for (entr in names(annotations)) {
        if (is.na(entr)){
          print('Entrez ID unknown')
        }
        if (term %in% names(annotations[[entr]]) | sum(children %in% names(annotations[[entr]])) > 0){
          GOcounts[[term]] <- GOcounts[[term]] + 1
          # Match column number
          colnum <- match(entr, entrezModule)
          d.GO[GOTERM[[term]]@Term, colnum] <- log2ModuleProteins[colnum]
        }
      }
    }
    # Sort columns by log2FC values
    d.GO <- d.GO[,order(log2ModuleProteins)]
    
    # Check for very large overlap in GO terms (i.e. GO terms that are probably related)
    #print(head(rownames(d.GO)))
    terms.mapping <- select(GO.db, keys=rownames(d.GO), columns='GOID', keytype='TERM')
    terms.mapping$print <- rep(TRUE, nrow(terms.mapping))
    for (j in 1:nrow(terms.mapping)) {
      ID <- terms.mapping$GOID[j]
      term <- terms.mapping$TERM[j]
      children <- c(GO.offspr$bp[[ID]], GO.offspr$cc[[ID]], GO.offspr$mf[[ID]])
      # Don't put a term in the heatmap if one of its direct descendants is also in the list
      if (sum(children %in% terms.mapping$GOID > 0)){
        terms.mapping$print[j] <- FALSE
      }
    }
    #print(sum(terms.mapping$print))
    #print(terms.mapping[terms.mapping$print,])
    n.include <- sum(terms.mapping$print)
    #print(n.include)
    if (n.include > 10){
      d.GO <- d.GO[terms.mapping$print,][1:10,]
    } else {
      d.GO <- d.GO[terms.mapping$print,]
    }
    #print(dim(d.GO))
    
    # Remove empty columns
    d.GO <- d.GO[,apply(d.GO, 2, function(x){sum(!is.na(x)) > 0})]
    
    # Create heatmap of clusters
    plotfile <- paste0(fig.folder, '/GO_heatmap_', i, '_', module, '.png')
    ylab.size <- max(unlist(lapply(rownames(d.GO), nchar)))/15
    h <- pheatmap(d.GO,
                  cluster_rows=FALSE,
                  cluster_cols=FALSE,
                  legend=TRUE,
                  color=my_palette,
                  breaks=colors,
                  show_rownames=TRUE,
                  show_colnames=TRUE,
                  filename=plotfile,
                  width=ncol(d.GO)/5+ylab.size,
                  height=nrow(d.GO)/4,
                  na_col='#E3E3E3'
    )
    write.table(d.GO, file=paste0(out.folder, '/GO_heatmapvalues_', i, '_', module, '.tsv'),
                sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
  }
  return(list(GO.enrichment, d.GO))
}