argv = commandArgs(trailingOnly=TRUE)

type <- argv[1]
if (type == 'subnetworks' || type== 'modules'){
  colorsA <- as.vector(t(read.table(argv[1], sep='\t', header = FALSE)[1,]))
  colorsB <- as.vector(t(read.table(argv[1], sep='\t', header = FALSE)[1,]))
  colorsBoth <- intersect(colorsA,colorsB)
  
  dir1 <- ""
  dir2 <- ""
  
  if (type == 'modules'){
    dir1 <- paste0(argv[2], '/coexpression/')
    dir2 <- paste0(argv[3], '/coexpression/')
  } else {
    dir1 <- paste0(argv[2], '/hotnet/HotNet_results/consensus_nodes_log2_')
    dir2 <- paste0(argv[3], '/hotnet/HotNet_results/consensus_nodes_log2_')
  }
  for (colors in colorsBoth){
    # Set A - numeric vector    
    SetA <- as.vector(t(read.table(paste0(dir1, colors, '.tsv'), row.names = FALSE, col.names = FALSE, sep='\t')[1,]))
    
    # Set B - numeric vector
    SetB <- as.vector(t(read.table(paste0(dir2, colors, '.tsv'), row.names = FALSE, col.names = FALSE, sep='\t')[1,]))
    
    # Function for computing Jaccard Similarity
    jaccard_similarity <- function(A, B) {
      intersection = length(intersect(A, B))
      union = length(A) + length(B) - intersection
      return (intersection/union)
    }
    
    # Jaccard Similarity between sets, A and B
    Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
    print(paste('Similarity for', colors, type, ':', Jaccard_Similarity))
    
    # Jaccard Dissimilarity/Distance between sets, A and B 
    Jaccard_Distance = 1 - Jaccard_Similarity
    Jaccard_Distance
  }
} else {
  # Set A - numeric vector    
  SetA <- as.vector(t(read.table(argv[2], row.names = FALSE, col.names = FALSE)[1,]))
  
  # Set B - numeric vector
  SetB <- as.vector(t(read.table(argv[3], row.names = FALSE, col.names = FALSE)[1,]))
  
  # Function for computing Jaccard Similarity
  jaccard_similarity <- function(A, B) {
    intersection = length(intersect(A, B))
    union = length(A) + length(B) - intersection
    return (intersection/union)
  }
  
  # Jaccard Similarity between sets, A and B
  Jaccard_Similarity <- jaccard_similarity(SetA,SetB)
  print(paste('Similarity for', type, ':', Jaccard_Similarity))
  
  # Jaccard Dissimilarity/Distance between sets, A and B 
  Jaccard_Distance = 1 - Jaccard_Similarity
  Jaccard_Distance
}