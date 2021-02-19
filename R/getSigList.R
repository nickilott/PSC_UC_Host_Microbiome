

getSigList <- function(result.table, level="ASV", padj=0.05){

    # return a vactor of significant feratures from
    # a multiDE analysis
    asv.df <- result.table[result.table$level == level,]
    sig <- asv.df$test_id[asv.df$padj < padj & !(is.na(asv.df$padj))]
    return(sig)
    }