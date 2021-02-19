library(broom)
library(dplyr)

# Functions for doing Wilcoxon rank sum test on
# microbiome and host modules

codeZeros <- function(abundances){
  categories <- ifelse(is.na(abundances), "Zero", "Non-zero")
  return(categories)
}

#######################################
#######################################
#######################################

wilcoxEigengenes <- function(clr.matrix, mes){

  # make sure that clr.matrix and mes have same
  # columns
  stopifnot(all.equal(colnames(clr.matrix), colnames(mes)))
  
  # make sure they are in the same order
  clr.matrix <- clr.matrix[, colnames(mes)]
  
  results <- list()
  c <- 0
  for (i in 1:nrow(clr.matrix)){
    feature <- rownames(clr.matrix[i,])
    categories <- codeZeros(clr.matrix[feature,])
    if (!("Zero" %in% categories)){next}
    for (j in 1:nrow(mes)){
      c <- c+1
      module <- rownames(mes[j,])
      exprs <- unlist(mes[j,])
      res <- wilcox.test(exprs ~ as.factor(categories))
      res <- data.frame(broom::tidy(res))
      res$feature <- feature
      res$Module <- module
      results[[c]] <- res
    }
  }
results <- bind_rows(results)
results$padj <- p.adjust(results$p.value, method="BH")
return(results)
}  



  