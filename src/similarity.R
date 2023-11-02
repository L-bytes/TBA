#' Similarity analysis
#' 
#' @param folder1 Folder of exploration dataset output
#' @param folder2 Folder of validation dataset output
#' @param outfolder Output folder for similarity output

similarity.analysis <- function(folder1, folder2, outfolder){

  jaccard_similarity <- function(A, B) {
    intersection = length(intersect(A, B))
    union = length(A) + length(B) - intersection
    return (intersection/union)
  }
  
  #GSEA
  dir1 <- paste0(folder1, '/GO/')
  dir2 <- paste0(folder2, '/GO/')
  gsea <- matrix(0, ncol = 1, nrow = 2)
  rownames(gsea) <- c('BP', 'MF')
  if (file.exists(paste0(dir1, '/GSEA_BP.tsv')) && file.exists(paste0(dir2, '/GSEA_BP.tsv'))){
    if (file.info(paste0(dir1, '/GSEA_BP.tsv'))$size != 0 && file.info(paste0(dir2, '/GSEA_BP.tsv'))$size != 0){
      SetA <- as.vector(t(read.table(paste0(dir1, 'GSEA_BP.tsv'), sep='\t', header = FALSE)))
      SetB <- as.vector(t(read.table(paste0(dir2, 'GSEA_BP.tsv'), sep='\t', header = FALSE)))
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      gsea['BP',1] <- jaccard_similarity(SetA,SetB)
      print(paste('Similarity for GSEA BP:', Jaccard_Similarity))
    }
  }
  
  if (file.exists(paste0(dir1, '/GSEA_MF.tsv')) && file.exists(paste0(dir2, '/GSEA_MF.tsv'))){
    if (file.info(paste0(dir1, '/GSEA_MF.tsv'))$size != 0 && file.info(paste0(dir2, '/GSEA_MF.tsv'))$size != 0){
      SetA <- as.vector(t(read.table(paste0(dir1, 'GSEA_MF.tsv'), sep='\t', header = FALSE)))
      SetB <- as.vector(t(read.table(paste0(dir2, 'GSEA_MF.tsv'), sep='\t', header = FALSE)))
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      gsea['MF', 1] <- Jaccard_Similarity
      print(paste('Similarity for GSEA MF:', Jaccard_Similarity))
    }
  }
  
  write.table(gsea, file=paste0(outfolder, '/gsea.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  #DESEQ2
  dir1 <- paste0(folder1, '/differential expression/')
  dir2 <- paste0(folder2, '/differential expression/')
  
  DE <- matrix(0, ncol = 1, nrow = 1)

  DF1 <- read.table(paste0(dir1, 'DE_LFC.tsv'), sep='\t', header = TRUE, row.names = 1)
  DF2 <- read.table(paste0(dir2, 'DE_LFC.tsv'), sep='\t', header = TRUE, row.names = 1)
  
  R1 <- rownames(DF1)
  R2 <- rownames(DF2)
  
  DF2 <- DF2[rownames(DF2) %in% R1,]
  DF1 <- DF1[rownames(DF1) %in% R2,]
  
  DF3 <- DF1-DF2
  LFC <- abs(DF3)
  print(paste('Similarity for DE genes', ':',  length(LFC[LFC <= 0.25])/length(DF3)))
  
  DE[1,1] <- length(LFC[LFC <= 0.25])/length(DF3)
  
  write.table(DE, file=paste0(outfolder, '/de.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  DE2 <- matrix(0, ncol = 2, nrow = 5)
  colnames(DE2) <- c('Positive', 'Negative')
  rownames(DE2) <- c('Top 10', 'Top 50', 'Top 100', 'Top 500', 'Top 1000')

  geneList1 <- as.vector(t(read.table(paste0(dir1, 'DE_genes.tsv'), sep='\t', header=FALSE)))
  geneList2 <- as.vector(t(read.table(paste0(dir1, 'DE_genes.tsv'), sep='\t', header=FALSE)))

  SetA <- geneList1[1:10]
  SetB <- geneList2[1:10]
  DE2['Top 10','Positive'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- geneList1[1:50]
  SetB <- geneList2[1:50]
  DE2['Top 50','Positive'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- geneList1[1:100]
  SetB <- geneList2[1:100]
  DE2['Top 100','Positive'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- geneList1[1:500]
  SetB <- geneList2[1:500]
  DE2['Top 500','Positive'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- geneList1[1:1000]
  SetB <- geneList2[1:1000]
  DE2['Top 1000','Positive'] <- jaccard_similarity(SetA,SetB)

  SetA <- tail(geneList1,10)
  SetB < tail(geneList2,10)
  DE2['Top 10','Negative'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- tail(geneList1,50)
  SetB <- tail(geneList2,50)
  DE2['Top 50','Negative'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- tail(geneList1,100)
  SetB <- tail(geneList2,100)
  DE2['Top 100','Negative'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- tail(geneList1,500)
  SetB <- tail(geneList2,500)
  DE2['Top 500','Negative'] <- jaccard_similarity(SetA,SetB)
  
  SetA <- tail(geneList1,1000)
  SetB <- tail(geneList2,1000)
  DE2['Top 1000','Negative'] <- jaccard_similarity(SetA,SetB)

  write.table(DE2, file=paste0(outfolder, '/de2.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  #MODULES
  dir1 <- paste0(folder1, '/coexpression/')
  dir2 <- paste0(folder2, '/coexpression/')
  colorsA <- as.vector(t(read.table(paste0(dir1, 'modules.tsv'), sep='\t', header = FALSE)))
  colorsB <- as.vector(t(read.table(paste0(dir2, 'modules.tsv'), sep='\t', header = FALSE)))
  exclusiveA <- colorsA[!(colorsA == 'grey')]
  exclusiveB <- colorsB[!(colorsB == 'grey')]
  modules <- matrix(0, ncol = length(exclusiveA), nrow = length(exclusiveB))
  colnames(modules) <- exclusiveA
  rownames(modules) <- exclusiveB
  for (colorA in exclusiveA){
    best <- 0
    bestColor <- ""
    for (colorB in exclusiveB){
      SetA <- as.vector(t(read.table(paste0(dir1, colorA, '.tsv'), sep='\t', header = FALSE)))
      
      SetB <- as.vector(t(read.table(paste0(dir2, colorB, '.tsv'), sep='\t', header = FALSE)))
      
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      modules[colorB, colorA] <- Jaccard_Similarity
      if (Jaccard_Similarity > best){
        best <- Jaccard_Similarity
        bestColor <- colorB
      }
    }
    print(paste('Similarity for', colorA, 'and', bestColor, 'modules:', best))
  }
  
  write.table(modules, file=paste0(outfolder, '/moduleSimilarity.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  #SUBNETWORKS
  dir1 <- paste0(folder1, '/hotnet/HotNet_results/consensus_nodes_log2_')
  dir2 <- paste0(folder2, '/hotnet/HotNet_results/consensus_nodes_log2_')
  
  colorsA <- as.vector(t(read.table(paste0(folder1, '/GO/hotnetSubnetworks.tsv'), sep='\t', header = FALSE)))
  colorsB <- as.vector(t(read.table(paste0(folder2, '/GO/hotnetSubnetworks.tsv'), sep='\t', header = FALSE)))
  exclusiveA <- colorsA[!(colorsA == 'grey')]
  exclusiveB <- colorsB[!(colorsB == 'grey')]
  subnetworks <- matrix(0, ncol = length(exclusiveA), nrow = length(exclusiveB))
  colnames(subnetworks) <- exclusiveA
  rownames(subnetworks) <- exclusiveB
  for (colorA in exclusiveA){
    best <- 0
    bestColor <- ""
    for (colorB in exclusiveB){
      if (file.info(paste0(dir1, colorA, '.tsv'))$size == 0 || file.info(paste0(dir2, colorB, '.tsv'))$size == 0){
        next
      }
      SetA <- unlist(read.csv(paste0(dir1, colorA, '.tsv'), sep='\t', header = FALSE))
      SetB <- unlist(read.csv(paste0(dir2, colorB, '.tsv'), sep='\t', header = FALSE))
      
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      subnetworks[colorB, colorA] <- Jaccard_Similarity
      if (Jaccard_Similarity > best){
        best <- Jaccard_Similarity
        bestColor <- colorB
      }
    }
    print(paste('Similarity for', colorA, 'and', bestColor, 'subnetworks:', best))
  }
  
  write.table(subnetworks, file=paste0(outfolder, '/subnetworkSimilarity.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  
  sSubnetworks <- matrix(0, ncol = length(exclusiveA), nrow = length(exclusiveB))
  colnames(sSubnetworks) <- exclusiveA
  rownames(sSubnetworks) <- exclusiveB
  
  overallSimilarity <- matrix(0, nrow=3, ncol=1)
  rownames(overallSimilarity) <- c('genes', 'BP', 'MF')
  
  write.table(sSubnetworks, file=paste0(outfolder, '/sSubnetworkSimilarity.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  write.table(overallSimilarity, file=paste0(outfolder, '/similarity.tsv'), row.names=TRUE, sep='\t')
  
  if (file.info(paste0(folder1, '/GO/subnetworks.tsv'))$size != 0 && file.info(paste0(folder2, '/GO/subnetworks.tsv'))$size != 0) {
    #SIGNFICANT SUBNETWORKS
    colorsA <- as.vector(t(read.table(paste0(folder1, '/GO/subnetworks.tsv'), sep='\t', header = FALSE)))
    colorsB <- as.vector(t(read.table(paste0(folder2, '/GO/subnetworks.tsv'), sep='\t', header = FALSE)))
    exclusiveA <- colorsA[!(colorsA == 'grey')]
    exclusiveB <- colorsB[!(colorsB == 'grey')]
    sSubnetworks <- matrix(0, ncol = length(exclusiveA), nrow = length(exclusiveB))
    colnames(sSubnetworks) <- exclusiveA
    rownames(sSubnetworks) <- exclusiveB
    for (colorA in exclusiveA){
      best <- 0
      bestColor <- ""
      for (colorB in exclusiveB){
        if (file.info(paste0(dir1, colorA, '.tsv'))$size == 0 || file.info(paste0(dir2, colorB, '.tsv'))$size == 0){
          next
        }
        SetA <- unlist(read.csv(paste0(dir1, colorA, '.tsv'), sep='\t', header = FALSE))
        SetB <- unlist(read.csv(paste0(dir2, colorB, '.tsv'), sep='\t', header = FALSE))
        
        Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
        sSubnetworks[colorB, colorA] <- Jaccard_Similarity
        if (Jaccard_Similarity > best){
          best <- Jaccard_Similarity
          bestColor <- colorB
        }
      }
      print(paste('Similarity for', colorA, 'and', bestColor, 'signficant subnetworks:', best))
    }
    
    dir1 <- paste0(folder1, '/GO/')
    dir2 <- paste0(folder2, '/GO/') 
    for (colorA in exclusiveA){
      for (colorB in exclusiveB){
        for (n in nrow(read.csv(paste0(folder1, '/hotnet/HotNet_results/consensus_nodes_log2_', colorA, '.tsv'), sep='\t', header = FALSE))){
          for (m in nrow(read.csv(paste0(folder2, '/hotnet/HotNet_results/consensus_nodes_log2_', colorB, '.tsv'), sep='\t', header = FALSE))){
            try({SetA <- as.vector(t(read.table(paste0(dir1, colorA, 'BP', n, '.tsv'), sep='\t', header = FALSE)))}, silent=TRUE)
            try({SetB <- as.vector(t(read.csv(paste0(dir2, colorB, 'BP', n, '.tsv'), sep='\t', header = FALSE)))}, silent=TRUE)
            print(paste('Similarity for', colorA, 'and', colorB, 'BP:', jaccard_similarity(SetA,SetB)))
            
            try({SetA <- as.vector(t(read.table(paste0(dir1, colorA, 'MF', n, '.tsv'), sep='\t', header = FALSE)))}, silent=TRUE)
            try({SetB <- as.vector(t(read.csv(paste0(dir2, colorB, 'MF', n, '.tsv'), sep='\t', header = FALSE)))}, silent=TRUE)
            print(paste('Similarity for', colorA, 'and', colorB, 'MF:', jaccard_similarity(SetA,SetB)))
          }
        }
      }
    }
    
    #GENES
    SetA <- as.vector(t(read.csv(paste0(dir1, '/genes.tsv'), sep='\t', header = FALSE)))
    SetB <- as.vector(t(read.table(paste0(dir2, '/genes.tsv'), sep='\t', header = FALSE)))
    Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
    print(paste('Similarity for significant genes:', Jaccard_Similarity))
    overallSimilarity['genes',1] <- Jaccard_Similarity
    
    #BP 
    if (file.info(paste0(dir1, '/BP.tsv'))$size != 0 && file.info(paste0(dir2, '/BP.tsv'))$size != 0){
      SetA <- as.vector(t(read.table(paste0(dir1, '/BP.tsv'), sep='\t', header = FALSE)))
      SetB <- as.vector(t(read.table(paste0(dir2, '/BP.tsv'), sep='\t', header = FALSE)))
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      print(paste('Similarity for biological processes:', Jaccard_Similarity))
      overallSimilarity['BP',1] <- Jaccard_Similarity
    }
    
    #MF
    if (file.info(paste0(dir1, '/MF.tsv'))$size != 0 && file.info(paste0(dir2, '/MF.tsv'))$size != 0){
      SetA <- as.vector(t(read.table(paste0(dir1, '/MF.tsv'), sep='\t', header = FALSE)))
      SetB <- as.vector(t(read.table(paste0(dir2, '/MF.tsv'), sep='\t', header = FALSE)))
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      print(paste('Similarity for molecular functions:', Jaccard_Similarity))
      overallSimilarity['MF',1] <- Jaccard_Similarity
    }
    write.table(sSubnetworks, file=paste0(outfolder, '/sSubnetworkSimilarity.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
    write.table(overallSimilarity, file=paste0(outfolder, '/similarity.tsv'), row.names=TRUE, sep='\t')
  }
}