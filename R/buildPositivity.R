library(dplyr)
library(ggplot2)

buildPositivity <- function(mat, metadata, gene1, gene2){
  
  # build an object that specifies which cells are positive/negative
  # for gene1 and gene2 - outputs proportions
  exprs.col1 <- paste0(gene1, ".exprs")
  exprs.col2 <- paste0(gene2, ".exprs")
  metadata[,exprs.col1] <- unlist(mat[gene1,])
  metadata[,exprs.col2] <- unlist(mat[gene2,])
  metadata[,gene1] <- ifelse(metadata[,exprs.col1] > 0, paste0(gene1, "pos"), paste0(gene1, "neg"))
  metadata[,gene2] <- ifelse(metadata[,exprs.col2] > 0, paste0(gene2, "pos"), paste0(gene2, "neg"))
  metadata$status <- paste0(metadata[,gene1], ":", metadata[,gene2])
  return(metadata)
}


plotPositivity <- function(mat, metadata, gene1, gene2){
  
  positivity <- buildPositivity(mat, metadata, gene1, gene2)
  
  # labels
  labels <- positivity %>% dplyr::count(get(gene1), get(gene2)) %>% mutate(prop=n/sum(n)*100)
  return(labels)
  # have given up on this for the time being
   
}