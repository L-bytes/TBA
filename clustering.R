### Data clustering

#' Create clustered heatmaps of the samples
#' 
#' @param d Input a dataframe with expression values, using the Protein IDs as rows, samples as columns. Works best if values are centered and scaled.
#' @param ann_col Dataframe with groupings per sample. Rownames are the sample names, single column with factor for group values
#' @param ann_row Dataframe with groupings per protein. Rownames are the protein IDs, single column with factor for group values
#' @param anncol Color codes with colors for each of the groups listed in \code{ann_col} and \code{ann_row}
#' @param f Filename for the file to be saved to
#' @return Dataframe with the normalized protein expression values from \code{d}.

Clustering <- function(d, ann_col=NA, anncol=NA, ann_row=NA, f=NA) {
  # Create distance matrices and clustering hierarchy
  dr <- dist(1-cor(t(d)))
  hr <- hclust(dr)
  dc <- dist(1-cor(d))
  hc <- hclust(dc)
  
  # Color palette
  colors = c(seq(-5,-1,length=1000),seq(-.999999,.999999,length=1000),seq(1, 5,length=1000))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 2999)
  # Create heatmap of clusters
  h <- pheatmap(d,
                cluster_rows=hr,
                cluster_cols=hc,
                legend=TRUE,
                color=my_palette,
                breaks=colors,
                show_rownames=FALSE,
                show_colnames=FALSE,
                annotation_col=ann_col,
                annotation_row=ann_row,
                annotation_colors=anncol,
                filename=f
  )
  return(h)
}


#' Add row labels for clustering
#' 
#' @param bb.test.obj Output of the differential expression analysis from the bb.test function
#' @param p.thresholds Named vector of P values to be used for the different colors. Names are the corresponding FDR values
#' @return Significance labels for each protein in \code{bb.test.obj} that can be used for coloring in the \code{Clustering} function
#' 
add.sign.labels <- function(bb.test.obj, p.thresholds){
  uniprot.ids <- rownames(bb.test.obj$table)
  FDRs <- names(p.thresholds)
  n.thr <- length(FDRs)
  sign.label <- c()
  for (i in 1:nrow(bb.test.obj$table)){
    uniprot <- uniprot.ids[i]
    l <- 0
    j <- 1
    while (j <= n.thr){
      # If P value is lower than the threshold, add label
      if (bb.test.obj$table[i, 'Pvalue'] < p.thresholds[FDRs[j]]){
        l <- j
        break
      } else {
        # Else move on to the next threshold
        j <- j+1
        if (j == n.thr+1){
          l <- n.thr+1
        }
      }
    }
    # Add label to the vector
    sign.label <- c(sign.label, l)
  }
  names(sign.label) <- uniprot.ids
  return(sign.label)
}