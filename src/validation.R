### Functions to perform the permutation test to validate the identified modules

#' Shuffle protein labels
#' 
#' @param FC Named vector with log2FC values. Names should be protein identifiers (e.g. UniProt IDs)
#' @param P Named vector with P-values. Names should be protein identifiers (e.g. UniProt IDs). Should be identical to the names of FC
#' @return FC and P with shuffled labels
#' 
randomize <- function(FC, P){
  " Randomly shuffles the labels of the protein values "
  ids <- names(FC)
  ids.random <- sample(ids)
  names(FC) <- ids.random
  names(P) <- ids.random
  return(list(FC, P))
}


#' Calculate which fraction of the expression values is different in the same direction (positive/negative log2FC), and
#' has a P-value below a certain threshold in both sets
#' 
#' @param FC1 Named vector with log2FC values of the discovery set. Names should be protein identifiers (e.g. UniProt IDs)
#' @param P1 Named vector with P-values of the discovery set. Names should be protein identifiers (e.g. UniProt IDs). Should be identical to the names of FC1
#' @param FC2 Named vector with log2FC values of the validation set. Names should be protein identifiers (e.g. UniProt IDs). IDs must be of the same type as FC1 to allow matching
#' @param P2 Named vector with P-values of the validation set. Names should be protein identifiers (e.g. UniProt IDs). Should be identical to the names of FC2
#' @param t P-value threshold for significance
#' @param randomize Boolean. Do the labels need to be shuffled for permutation?
#' @return fraction of proteins with log2FC values in the same direction and P-values below \code{t} in both sets

fraction.correct <- function(FC1, FC2, P1, P2, t, randomize=FALSE){
  " Calculates what fraction is of expression values is
  significantly different (at threshold t) in the same
  direction in both groups "
  # Get IDs
  ids1 <- names(FC1)
  ids2 <- names(FC2)
  
  # Shuffle labels for permutation test
  if(randomize){
    shuffled <- randomize(FC1, P1)
    FC1 <- shuffled[[1]]
    P1 <- shuffled[[2]]
  }
  
  # Calculate the fraction of proteins significantly different in both groups
  count <- 0
  total <- length(ids1)
  #count <- sum((P1 < t) & (P2 < t) & ((FC1*FC2) > 0))
  for(id in ids1){
    # Check if P-values are lower than t and log2FC expression goes in the same direction
    if(P1[id]<t && P2[id]<t && FC1[id]*FC2[id]>0){
      count <- count + 1
    }
  }
  return(count/total)
}

#' Calculate which fraction of the expression values is different in the same direction (positive/negative log2FC), and
#' has a P-value below a certain threshold in both sets
#' 
#' @param FC1 Named vector with log2FC values of the discovery set. Names should be protein identifiers (e.g. UniProt IDs)
#' @param P1 Named vector with P-values of the discovery set. Names should be protein identifiers (e.g. UniProt IDs). Should be identical to the names of FC1
#' @param FC2 Named vector with log2FC values of the validation set. Names should be protein identifiers (e.g. UniProt IDs). IDs must be of the same type as FC1 to allow matching
#' @param P2 Named vector with P-values of the validation set. Names should be protein identifiers (e.g. UniProt IDs). Should be identical to the names of FC2
#' @param t.range vector with P-value thresholds for significance
#' @param randomize Boolean. Do the labels need to be shuffled for permutation?
#' @return vector with fraction of proteins with log2FC values in the same direction and P-values below \code{t} in both sets for each threshold \code{t} in \code{t.range}

calculate.curve <- function(FC1, FC2, P1, P2, t.range, randomize=FALSE){
  values <- c()
  for(t in t.range){
    values <- c(values, fraction.correct(FC1, FC2, P1, P2, t, randomize=randomize))
  }
  return(values)
}

#' Perform permutation test
#' 
#' @param FC1 Named vector with log2FC values of the discovery set. Names should be protein identifiers (e.g. UniProt IDs)
#' @param P1 Named vector with P-values of the discovery set. Names should be protein identifiers (e.g. UniProt IDs). Should be identical to the names of FC1
#' @param FC2 Named vector with log2FC values of the validation set. Names should be protein identifiers (e.g. UniProt IDs). IDs must be of the same type as FC1 to allow matching
#' @param P2 Named vector with P-values of the validation set. Names should be protein identifiers (e.g. UniProt IDs). Should be identical to the names of FC2
#' @param out.folder Folder for output data files
#' @param fig.folder Folder for output figure files
#' @param m.ids \code{FALSE} or vector with protein identifiers of proteins within a certain module. All proteins in \code{m.ids} must be in \code{FC1}, \code{P1}, \code{FC2} and \code{P2}
#' @param tmin Lowest P-value threshold for significance to be tested
#' @param tmax Highest P-value threshold for significance to be tested
#' @param steps Number of P-value thresholds to be checked
#' @param checks P-value thresholds over which the GESS score is calculated
#' @param nPerm Number of permutations
#' @return vector with fraction of proteins with log2FC values in the same direction and P-values below \code{t} in both sets for each threshold \code{t} in \code{t.range}

permutation.test <- function(FC1, FC2, P1, P2, out.folder, fig.folder, m.ids=NULL, tmin=0, tmax=1, steps=31, checks=c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8), nPerm=1000){
  if(max(checks)>tmax){stop("Maximum value in checks is larger than tmax")}
  if(min(checks)<tmin){stop("Minimum value in checks is smaller than tmin")}
  # Add checkpoints for which GESS score is calculated
  #for(i in c(0.01, 0.1, 0.2, 0.5)){
  for(i in c(0.1, 0.2, 0.5, 0.8)){
    if(!i %in% checks){
      checks <- c(checks, i)
    }
  }
  checks <- sort(checks)
  
  t.values <- seq(tmin,tmax,length.out=steps)
  for(x in checks){
    if(!x %in% t.values){
      t.values <- c(t.values, x)
    }
  }
  t.values <- sort(t.values)
  
  # Calculate permuted values
  if(nPerm > 0){
    d.perm <- data.frame(matrix(nrow=nPerm, ncol=length(t.values)))
    colnames(d.perm) <- t.values
    
    # Get module proteins for the validation set
    if(!is.null(m.ids)){
      m.FC2 <- FC2[m.ids]
      m.P2 <- P2[m.ids]
    } else {
      m.FC2 <- FC2
      m.P2 <- P2
    }
    
    for(i in 1:nPerm){
      # Shuffle labels of discovery set
      shuffled <- randomize(FC1, P1)
      random.FC <- shuffled[[1]]
      random.P <- shuffled[[2]]
      
      # Filter for module proteins
      if(!is.null(m.ids)){
        random.FC <- random.FC[m.ids]
        random.P <- random.P[m.ids]
      }
      
      values <- calculate.curve(random.FC, m.FC2, random.P, m.P2, t.values, randomize=FALSE)
      d.perm[i,] <- values
    }
  }
  # Calculate actual values
  if(!is.null(m.ids)){
    m.FC1 <- FC1[m.ids]
    m.P1 <- P1[m.ids]
  } else {
    m.FC1 <- FC1
    m.P1 <- P1
  }
  values.real <- calculate.curve(m.FC1, m.FC2, m.P1, m.P2, t.values)
  names(values.real) <- t.values
  
  # Calculate GESS score
  # Method 1
  to.check <- unlist(lapply(c(0.1, 0.2, 0.5, 0.8), as.character))
  #real.gess <- values.real[to.check]
  #real.gess <- sum(real.gess)/length(to.check)
  #values.real['GESS'] <- real.gess
  # Method 2
  #real.gess2 <- weighted.gess(values.real[to.check], counts.real[to.check])
  #values.real['GESS2'] <- real.gess2
  
  # Method 1
  #d.gess <- d.perm[,to.check]
  #perm.gess <- apply(d.gess, 1, sum)
  #perm.gess <- unlist(lapply(perm.gess, function(x){x/length(to.check)}))
  #d.perm[,'GESS'] <- perm.gess
  #p.gess <- pnorm(real.gess, mean=mean(perm.gess), sd=sd(perm.gess), lower.tail=FALSE)
  # Method 2
  #perm.gess2 <- c()
  #for (i in 1:nrow(d.perm)){
  #  perm.gess2[i] <- weighted.gess(values, counts)
  #}
  #d.perm[,'GESS2'] <- perm.gess2
  #p.gess2 <- pnorm(real.gess2, mean=mean(perm.gess2), sd=sd(perm.gess2), lower.tail=FALSE)
  
  # Plot histograms with the score distributions
  #checks <- c(checks, 'GESS')
  pvals <- rep(NA, length(checks))
  names(pvals) <- checks
  for(x in checks){
    x <- as.character(x)
    png(file=paste0(fig.folder, 'permutation_density_', x, '.png'), width=720, height=480)
    hist(d.perm[,x], breaks=50, freq=TRUE, xlim=c(0,max(c(max(d.perm[,x])+0.05), values.real[x]+0.05)))
    #lines(density(d.perm[,x]), lwd=2)
    points(values.real[x], 0, pch=4, lwd=2)
    dev.off()
    #p.value <- t.test(d.perm[,x], mu=values.real[x])$p.value
    p.value <- pnorm(values.real[x], mean=mean(d.perm[,x]), sd=sd(d.perm[,x]), lower.tail=FALSE)
    pvals[x] <- p.value
  }
  p.gess <- mean(pvals[to.check])
  #pvals['GESS'] <- p.gess
  #pvals['GESS2'] <- p.gess2
  pvals['GESS'] <- p.gess
  
  # Plot curves
  png(file=paste0(fig.folder, 'permutation_error_curve.png'), width=480, height=480)
  for(i in 1:nrow(d.perm)){
    if(i==1){
      plot(t.values, d.perm[i,], type='l', lwd=1, col='grey', ylim=c(0,1), xlab='P-value threshold', ylab='Fraction validated', cex.lab=1.5, cex.axis=1.5, las=1)#xaxt='n', yaxt='n')
      #axis(las=1)
    } else{
      lines(t.values, d.perm[i,], type='l', lwd=1, col='grey')
    }
    
  }
  lines(t.values, values.real[], type='l', lwd=2, col='blue')
  # Gess 3 score
  text(0.2, 1, paste0('GESS = ', formatC(p.gess, format="e", digits=2)), cex=1.5)
  #text(0.2, 0.9, paste0('P-GESS2 = ', formatC(p.gess2, format="e", digits=2)), cex=1.5)
  #text(0.2, 1, paste0('P-GESS = ', formatC(p.gess, format="e", digits=2)), cex=1.5)
  dev.off()
  
  write.table(data.frame(Fraction_correct=values.real), paste0(out.folder, '/error_curve_real.tsv'), row.names=TRUE, col.names=TRUE, sep='\t', quote=FALSE)
  write.table(data.frame(d.perm), paste0(out.folder, '/error_curve_permutations.tsv'), row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
  write.table(data.frame(Pvalue=pvals), paste0(out.folder, '/error_curve_pvals.tsv'), row.names=TRUE, col.names=TRUE, sep='\t', quote=FALSE)
  return(list(values.real, d.perm, pvals))
}

