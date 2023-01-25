### Perform the beta binomial test to assess significance of differential expression

#' Perform the beta binomial test to assess significance of differential expression
#' Subsequently fits a beta uniform mixture model to match P-values with FDR values
#' 
#' @param d Input a dataframe with the Protein IDs as rows, samples as columns
#' @param ids Protein identifiers (e.g. Uniprot IDs)
#' @param groups Names of the sample groups. The length must match the number of columns in \code{d}
#' @param g1 Name of the first sample group to compare
#' @param g2 Name of the second sample group to compare
#' @param test.type Alternative hypothesis. Should be one of: "two.sided", "less", "greater"
#' @return List with two items. The first is a dataframe containing the protein IDs (\code{ids}), Log2 fold change values, and P-values. The second is a table stating which P-value threshold should be used to obtain a given FDR (0.01, 0.05, 0.1, 0.2), and how many proteins would be significant under this criterion.
#' @examples
#' betaBinomial(data.frame(matrix(data=rnorm(50), ncol=5)), c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), c("A", "A", "A", "B", "B"), "A", "B", "two.sided")


betaBinomial <- function(d, ids, groups, g1, g2, test.type='two.sided'){
  # Get groups to compare
  d.sub <- d[, groups %in% g1 | groups %in% g2]
  groups.sub <- groups[groups %in% g1 | groups %in% g2]
  # Calculate mean difference between the groups (for output only, not used in the test)
  mean1 <- apply(d[, groups %in% g1], 1, mean)
  mean2 <- apply(d[, groups %in% g2], 1, mean)
  mean.diff <- log2(mean2 / mean1)
  names(mean.diff) <- ids
  
  # Beta-binomial test
  totals <- apply(d.sub, 2, sum)
  new.groups <- c()
  for (i in 1:length(groups)){
    if (groups[i] %in% g1){
      new.groups <- c(new.groups, 'g1')
    } else if (groups[i] %in% g2){
      new.groups <- c(new.groups, 'g2')
    }
  }
  result <- bb.test(d.sub, totals, group=new.groups, alternative=test.type, n.threads=-1)
  #p.fdr <- p.adjust(result$p.value, method='BH')
  # Data frame with Uniprot IDs and corresponding P-values
  sign <- data.frame(Uniprot=ids, Log2ratio=mean.diff, Pvalue=result$p.value, row.names=ids)
  sign <- sign[order(sign$Pvalue, decreasing=FALSE),]
  # Fit BUM model
  bum <- fitBumModel(result$p.value, plot=FALSE)
  # Pick P-value threshold with given FDR value
  thresholds <- data.frame(Pvalues=c(0.01, 0.05, 0.1, 0.2), counts=rep(NA, 4))
  colnames(thresholds) <- c(paste(g2, 'vs', g1, 'Pvalue', sep='_'), paste(g2, 'vs', g1, 'sign', sep='_'))
  rownames(thresholds) <- c(0.01, 0.05, 0.1, 0.2)
  thresholds[,1] <- apply(thresholds, 1, function(x){fdrThreshold(x[1], bum)})
  thresholds[,2] <- unlist(lapply(thresholds[,1], function(x){sum(sign$Pvalue < x)}))
  
  #p.sign <- fdrThreshold(FDR, bum)
  # Return data frame of Uniprots with P-value lower than threshold
  return(list(table=sign, FDRs=thresholds))
}


