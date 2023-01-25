### Merge protein differential expression values with additional info from Uniprot
#' Merge protein differential expression values with additional info from Uniprot
#' Change the colnames of d.summary (line 14-25) according to the features you would like to add
#' 
#' @param d Input a dataframe with the differential expression (or other values of interest), using Protein IDs as rows, features as columns
#' @param d.uniprot Database info to be added to the \code{d}. Row names of \code{d} and \code{d.uniprot} must be matched. "Entry" must be one of the columns downloaded from UniProt to be able to match IDs
#' @param col.names A vector with the desired column names of the output table
#' @return Dataframe with Uniprot info added to \code{d}.

add_uniprot <- function(d, d.uniprot, col.names){
  rownames(d.uniprot) <- d.uniprot$Entry
  ids <- rownames(d)
  
  nrows <- nrow(d)
  ncols <- ncol(d) + ncol(d.uniprot)
  d.summary <- data.frame(matrix(ncol=ncols, nrow=nrows))
  colnames(d.summary) <- col.names
  rownames(d.summary) <- ids
  
  # Load obsolete UniProt IDs
  uniprot.obsolete <- read.table('Data/old_uniprot_IDs.txt', sep='\t', header=FALSE, row.names=1)
  for (i in 1:nrows) {
    
    # Take care of some obsolete entries, or entries that have received a new UniProt ID
    if (ids[i] %in% rownames(d.uniprot)){
      uniprot.info <- d.uniprot[ids[i],]
    } else if (uniprot.obsolete[ids[i], 1] == 'Obsolete') {
      uniprot.info <- c(ids[i], rep('', 23))
      names(uniprot.info) <- colnames(d.uniprot)
    } else if (uniprot.obsolete[ids[i], 1] %in% rownames(d.uniprot)){
      new.ID <- uniprot.obsolete[ids[i], 1]
      uniprot.info <- d.uniprot[new.ID,]
      uniprot.info[1] <- ids[i]
    } else {
      print('Something is wrong')
      print(i)
      print(ids[i])
    }
    # Add values to the dataframe
    d.summary[ids[i], colnames(d)] <- d[ids[i],]
    d.summary[ids[i], -colnames(d)] <- uniprot.info
    
    
    #d.summary[ids[i], 1] <- uniprot.info['Entry']
    #d.summary[ids[i], 2] <- uniprot.info['Cross.reference..GeneID.']
    #d.summary[ids[i], 3:8] <- c(d$log2TAUcon, d$pvalTAUcon, d$log2TAUTDP, d$pvalTAUTDP, d$log2TDPcon, d$pvalTDPcon)
    #d.summary[ids[i], 9] <- uniprot.info['Entry.name']
    #d.summary[ids[i], 10] <- uniprot.info['Gene.names']
    #d.summary[ids[i], 11] <- uniprot.info['Protein.names']
    #d.summary[ids[i], 12] <- uniprot.info['Status']
    #d.summary[ids[i], 13] <- uniprot.info['Organism']
    #d.summary[ids[i], 14] <- uniprot.info['Length']
    #d.summary[ids[i], 15] <- uniprot.info['Keywords']
    #d.summary[ids[i], 16] <- uniprot.info['Keyword.ID']
    #d.summary[ids[i], 17] <- uniprot.info['Annotation']
    #d.summary[ids[i], 18] <- uniprot.info['Features']
    #d.summary[ids[i], 19] <- uniprot.info['Pathway']
    #d.summary[ids[i], 20] <- uniprot.info['Gene.ontology.IDs']
    #d.summary[ids[i], 21] <- uniprot.info['Gene.ontology..biological.process.']
    #d.summary[ids[i], 22] <- uniprot.info['Gene.ontology..cellular.component.']
    #d.summary[ids[i], 23] <- uniprot.info['Gene.ontology..molecular.function.']
    #d.summary[ids[i], 24] <- uniprot.info['Function..CC.']
  }
  return(d.summary)
}