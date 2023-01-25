### Preprocess the data files

#' Normalize raw data based on total protein count per sample.
#' 
#' @param d Input a dataframe with expression values, using the Protein IDs as rows, samples as columns
#' @return Dataframe with the normalized protein expression values from \code{d}.
#' @examples
#' normalize.sample(data.frame(matrix(data=rnorm(50), ncol=5)))

normalize.sample <- function(d){
  # Normalize raw data based on total protein count per sample. input the Protein IDs as rows, samples as columns
  m <- apply(d, 2, mean)
  d.norm <- t(t(d)/m) * mean(d)
  return(d.norm)
}

#' Center and scale the data.
#' 
#' @param d Input a dataframe with the Protein IDs as rows, samples as columns
#' @return Dataframe with the normalized protein expression values from \code{d}.
#' @examples
#' normalize.cs(data.frame(matrix(data=rnorm(50), ncol=5)))

normalize.cs <- function(d){
  m <- apply(d, 1, mean)
  sd <- apply(d, 1, sd)
  d.cs <- (d - m) / sd
  return(d.cs)
}
