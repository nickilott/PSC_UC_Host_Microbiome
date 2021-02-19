###############################################
###############################################
###############################################
# run wilcoxon rank sum test on a multi-feature
# matrix and output a data frame of stats
###############################################
###############################################
###############################################


run_wilcox <- function(feature.matrix, group){

    # NB make sure feature matrix columns are
    # in the same order as the group variable

    # create empty result matrix
    result <- matrix(ncol=2, nrow=nrow(feature.matrix))

    # iterate over input matrix and perform wilcox.test
    # on each feature. Add results to result matrix
    for (i in 1:nrow(feature.matrix)){
        values <- unlist(feature.matrix[i,], use.names=FALSE)
        test.results <- wilcox.test(values~group)
        w <- test.results$statistic[[1]]
	p <- test.results$p.value
	result[i,1] <- w
	result[i,2] <- p
    }
    result <- data.frame(result)
    colnames(result) <- c("W", "p.value")
    result$padj <- p.adjust(result$p, method="BH")
    rownames(result) <- rownames(feature.matrix)
    return(result)
    }

