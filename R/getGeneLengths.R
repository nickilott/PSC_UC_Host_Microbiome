#' Get gene lengths from a GTF file for input to GOSeq
#'
#' GOSeq requires length information to calculate biases in enrichments for pathways
#' analysis. I have found that the easiest way to deal with this is to first create
#' a data frame of these from a gtf file of interest.
#' @param gtf gtf file as input
#' @return data frame that maps gene ids to median transcript length for each gene in gtf.
#' @import GenomicFeatures plyr
#' @examples
#' getGeneLengths(gtf="geneset_all.gtf.gz")
#' @export

getGeneLengths <- function(gtf){

    txdata <- makeTxDbFromGFF(gtf)
    lengths <- transcriptLengths(txdata)
    lengths <- data.frame(gene_id=lengths$gene_id, length=lengths$tx_len)
    lengths <- ddply(lengths, "gene_id", summarize, median(length))

    colnames(lengths) <- c("gene_id", "length")
    rownames(lengths) <- lengths$gene_id
    lengths$gene_id <- as.character(lengths$gene_id)

    return(lengths)
    }