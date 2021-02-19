
filterRows <- function(mat, k=10, a=1){

    # filter matrix based on counts being over a
    # in at least k samples
    mat <- mat[rowSums(mat > a) >= k,]
    return(mat)
}

relab <- function(counts){

    relab <- (sweep(counts, 2, colSums(counts), "/"))*100
    return(relab)
    }

pprev <- function(counts){

    nsamples <- ncol(counts)
    prev <- (rowSums(counts > 0)/nsamples)*100
    prev.df <- data.frame(Taxon=rownames(counts), Prevalence=prev)
    return(prev.df)
    }
