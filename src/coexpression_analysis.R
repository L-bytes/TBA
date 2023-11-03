#' Coexpression analysis workflow using the WGCNA package
#' 
#' @param d Input a dataframe with expression values, using the gene symbols as rows, samples as columns.
#' @param d.log2fc log2FC values
#' @param outfolder Output folder for tables
#' @param figfolder Output folder for figures
#' @param i Number of run
#' @param stage Exploration or validation round
#' @return Coexpression network of \code{d} and significance of log2FC value distribution(s) within each of the modules

coexpression.analysis <- function(d, d.log2fc, outfolder, figfolder, i, stage, power=FALSE){
  powers <- c(c(1:11), seq(from=12, to=26, by=2))
  soft.tresh <- pickSoftThreshold(d, powerVector=powers, verbose=5)
  png(filename=paste0(figfolder, '/Scale_free_threshold.png'), width=960, height=480)
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(soft.tresh $fitIndices[,1], -sign(soft.tresh $fitIndices[,3])*soft.tresh $fitIndices[,2], xlab="Beta",ylab="Scale Free Topology Model Fit, R^2",type="n", main = paste("Scale independence"))
  text(soft.tresh $fitIndices[,1], -sign(soft.tresh $fitIndices[,3])*soft.tresh $fitIndices[,2], labels=powers,cex=cex1,col="red")
  abline(h=0.80,col="red")
  plot(soft.tresh $fitIndices[,1], soft.tresh $fitIndices[,5], xlab="Beta", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
  text(soft.tresh $fitIndices[,1], soft.tresh $fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=150,col="red")
  dev.off()
  if (!power){
    sf.values <- -sign(soft.tresh $fitIndices[,3])*soft.tresh $fitIndices[,2]
    sf.values2 <- -sign(soft.tresh $fitIndices[,5])
    threshold <- 0.8
    print(paste0("Threshold: ", threshold))
    sf.bool <- sf.values > threshold
    sf.bool <- sf.bool & (sf.values2 < 100)
    power <- soft.tresh $fitIndices[,1][sf.bool][1]
  }
  print(paste0("Power : ", power))
  rownames(d) <- paste(rownames(d), length(groups), sep='_')
  dSplit <- 2
  mergeCut <- 0.25
  print(paste0("deepSplit: ", dSplit))
  print(paste0("mergeCutHeight: ", mergeCut))
  coex.net <- blockwiseModules(d, power=power,
                                TOMType="unsigned", minModuleSize=30,
                               reassignThreshold=0, mergeCutHeight=mergeCut,
                               numericLabels=TRUE, pamRespectsDendro=FALSE, deepSplit = dSplit,
                               saveTOMs=TRUE, saveTOMFileBase=paste0("output/", i, "/", stage, "/rdata/coexpression_discovery"), verbose=3, maxBlockSize=ncol(d))
  moduleColors <- labels2colors(coex.net$colors)
  modules <- levels(as.factor(moduleColors))
  n.modules <- length(modules)
  par(mar = c(1, 1, 1, 1))
  png(file=paste0(figfolder, '/module_dendrogram.png'), width=1920, height=1080)
  plotDendroAndColors(coex.net$dendrograms[[1]], moduleColors,
                      "Module colors", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, cex.lab = 2, cex.axis = 2, cex.main = 2, cex.colorLabels = 2)
  dev.off()
  MEs <- coex.net$MEs
  ME.WT <- MEs[str_detect(rownames(MEs), 'WT'),]
  ME.SNV <- MEs[str_detect(rownames(MEs), 'SNV'),]
  
  pvals.modules <- data.frame(matrix(nrow=n.modules, ncol=ncol(d.log2fc)))
  colnames(pvals.modules) <- colnames(d.log2fc)
  rownames(pvals.modules) <- modules
  
  for (i in 1:(n.modules-1)){
    module <- labels2colors(i)
    logvals <- d.log2fc[names(coex.net$colors[coex.net$colors==i]),]
    t.stat <- t.test(logvals)
    pvals.modules[module] <- t.stat$p.value
    print(c(module, median(logvals), sd(logvals), length(logvals)))
#    par(mar = c(2, 1, 1, 1))
#    png(paste0(figfolder, '/ME_distribution_', i, '_', module, '.png'), width=1920, height=1080)
#      boxplot(ME.SNV[,paste0('ME', i)], ME.WT[,paste0('ME', i)], names=c('SNV', 'WT'),
#              col=module, xlab='group', ylab='ME value', cex.lab = 2, cex.axis = 2, cex.main = 2, cex.colorLabels = 2)
#    dev.off()
  }
  p.adj <- p.adjust(unlist(pvals.modules))
  pval.modules.adj <- pvals.modules
  for (j in 1:ncol(pvals.modules)){
    pvals.modules[,j] <- p.adj[((j-1)*n.modules+1):(j*n.modules)]
  }
  write.table(pval.modules.adj, file=paste0(outfolder, '/sign_log2fc_modules.tsv'), sep='\t',
              row.names=TRUE, col.names=TRUE, quote=FALSE)
  return(list(coex.net, pval.modules.adj, power))
}