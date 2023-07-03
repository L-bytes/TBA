jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  union = length(A) + length(B) - intersection
  return (intersection/union)
}

similarity.analysis <- function(folder1, folder2, outfolder){
  #MODULES
  dir1 <- paste0(folder1, '/coexpression/')
  dir2 <- paste0(folder2, '/coexpression/')
  colorsA <- as.vector(t(read.table(paste0(dir1, 'modules.tsv'), sep='\t', header = FALSE)))
  colorsB <- as.vector(t(read.table(paste0(dir2, 'modules.tsv'), sep='\t', header = FALSE)))
  exclusiveA <- colorsA[!(colorsA == 'grey')]
  exclusiveB <- colorsB[!(colorsB == 'grey')]
  modules <- matrix(ncol = length(exclusiveA), nrow = length(exclusiveB))
  colnames(modules) <- exclusiveA
  rownames(modules) <- exclusiveB
  for (colorA in exclusiveA){
    best <- 2
    bestColor <- ""
    for (colorB in exclusiveB){
      SetA <- as.vector(t(read.table(paste0(dir1, colorA, '.tsv'), sep='\t', header = FALSE)))
      
      SetB <- as.vector(t(read.table(paste0(dir2, colorB, '.tsv'), sep='\t', header = FALSE)))
      
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      Jaccard_Distance = 1 - Jaccard_Similarity
      modules[colorB, colorA] <- Jaccard_Distance
      if (Jaccard_Distance < best){
        best <- Jaccard_Distance
        bestColor <- colorB
      }
    }
    print(paste('Dissimilarity for', colorA, 'and', bestColor, 'modules', ':', best))
  }
  
  write.table(modules, file=paste0(outfolder, '/moduleSimilarity.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  #SIGNFICANT SUBNETWORKS
  dir1 <- paste0(folder1, '/hotnet/HotNet_results/consensus_nodes_log2_')
  dir2 <- paste0(folder2, '/hotnet/HotNet_results/consensus_nodes_log2_')
  
  #SUBNETWORKS
  colorsA <- as.vector(t(read.table(paste0(folder1, '/GO/hotnetSubnetworks.tsv'), sep='\t', header = FALSE)))
  colorsB <- as.vector(t(read.table(paste0(folder2, '/GO/hotnetSubnetworks.tsv'), sep='\t', header = FALSE)))
  exclusiveA <- colorsA[!(colorsA == 'grey')]
  exclusiveB <- colorsB[!(colorsB == 'grey')]
  subnetworks <- matrix(ncol = length(exclusiveA), nrow = length(exclusiveB))
  colnames(subnetworks) <- exclusiveA
  rownames(subnetworks) <- exclusiveB
  for (colorA in exclusiveA){
    best <- 2
    bestColor <- ""
    for (colorB in exclusiveB){
      if (file.info(paste0(dir1, colorA, '.tsv'))$size == 0 || file.info(paste0(dir2, colorB, '.tsv'))$size == 0){
        next
      }
      SetA <- as.vector(t(read.csv(paste0(dir1, colorA, '.tsv'), sep='\t', header = FALSE)[1,]))
      
      SetB <- as.vector(t(read.csv(paste0(dir2, colorB, '.tsv'), sep='\t', header = FALSE)[1,]))
      
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      Jaccard_Distance = 1 - Jaccard_Similarity
      subnetworks[colorB, colorA] <- Jaccard_Distance
      if (Jaccard_Distance < best){
        best <- Jaccard_Distance
        bestColor <- colorB
      }
    }
    print(paste('Dissimilarity for', colorA, 'and', bestColor, 'subnetworks', ':', best))
  }
  
  write.table(subnetworks, file=paste0(outfolder, '/subnetworkSimilarity.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  
  
  sSubnetworks <- matrix(1, ncol = length(exclusiveA), nrow = length(exclusiveB))
  colnames(sSubnetworks) <- exclusiveA
  rownames(sSubnetworks) <- exclusiveB
  
  overallSimilarity <- matrix(0, nrow=3, ncol=1)
  rownames(overallSimilarity) <- c('genes', 'BP', 'MF')
  
  if (file.info(paste0(folder1, '/GO/subnetworks.tsv'))$size != 0 && file.info(paste0(folder2, '/GO/subnetworks.tsv'))$size != 0) {
    colorsA <- as.vector(t(read.table(paste0(folder1, '/GO/subnetworks.tsv'), sep='\t', header = FALSE)))
    colorsB <- as.vector(t(read.table(paste0(folder2, '/GO/subnetworks.tsv'), sep='\t', header = FALSE)))
    exclusiveA <- colorsA[!(colorsA == 'grey')]
    exclusiveB <- colorsB[!(colorsB == 'grey')]
    for (colorA in exclusiveA){
      best <- 2
      bestColor <- ""
      for (colorB in exclusiveB){
        if (file.info(paste0(dir1, colorA, '.tsv'))$size == 0 || file.info(paste0(dir2, colorB, '.tsv'))$size == 0){
          next
        }
        SetA <- as.vector(t(read.csv(paste0(dir1, colorA, '.tsv'), sep='\t', header = FALSE)[1,]))
        
        SetB <- as.vector(t(read.csv(paste0(dir2, colorB, '.tsv'), sep='\t', header = FALSE)[1,]))
        
        Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
        Jaccard_Distance = 1 - Jaccard_Similarity
        sSubnetworks[colorB, colorA] <- Jaccard_Distance
        if (Jaccard_Distance < best){
          best <- Jaccard_Distance
          bestColor <- colorB
        }
      }
      print(paste('Dissimilarity for', colorA, 'and', bestColor, 'signficant subnetworks', ':', best))
    }
    
    #GENES
    dir1 <- paste0(folder1, '/GO/')
    dir2 <- paste0(folder2, '/GO/') 
    SetA <- as.vector(t(read.csv(paste0(dir1, '/genes.tsv'), sep='\t', header = FALSE)))
    SetB <- as.vector(t(read.table(paste0(dir2, '/genes.tsv'), sep='\t', header = FALSE)))
    Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
    print(paste('Similarity for significant genes:', Jaccard_Similarity))
    overallSimilarity['genes',1] <- Jaccard_Similarity
    
    #BP 
    if (file.info(paste0(dir1, '/BP.tsv'))$size != 0){
      SetA <- as.vector(t(read.table(paste0(dir1, '/BP.tsv'), sep='\t', header = FALSE)))
      SetB <- as.vector(t(read.table(paste0(dir2, '/BP.tsv'), sep='\t', header = FALSE)))
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      print(paste('Similarity for biological processes:', Jaccard_Similarity))
      overallSimilarity['BP',1] <- Jaccard_Similarity
    }
    
    #MF
    if (file.info(paste0(dir1, '/MF.tsv'))$size != 0){
      SetA <- as.vector(t(read.table(paste0(dir1, '/MF.tsv'), sep='\t', header = FALSE)))
      SetB <- as.vector(t(read.table(paste0(dir2, '/MF.tsv'), sep='\t', header = FALSE)))
      Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
      print(paste('Similarity for molecular functions:', Jaccard_Similarity))
      overallSimilarity['MF',1] <- Jaccard_Similarity
    }
  }
  write.table(sSubnetworks, file=paste0(outfolder, '/sSubnetworkSimilarity.tsv'), row.names=TRUE, col.names=TRUE, sep='\t')
  write.table(overallSimilarity, file=paste0(outfolder, '/similarity.tsv'), row.names=TRUE, sep='\t')
}