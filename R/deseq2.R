runDESeq2Pair <- function(mat, metadata, a=0, k=1, geneid2gene_name=NULL, tissue="Ileum", column="Tissue.location", results.file=NULL, c1="HEALTHY", c2="PSC/UC"){

    metadata <- metadata[metadata$Disease %in% c(c1, c2),]
    metadata <- metadata[metadata[,column] == tissue,]
    mat <- mat[,rownames(metadata)]

    dds <- DESeqDataSetFromMatrix(countData = mat,
                                  colData = metadata,
			                      design = ~ Disease)


    # get gene names
    featureData <- makeFeatureData(geneid2gene_name, mat)

    # filtering > a in at least k samples (i.e. all of one condition)
    keep <- rowSums(counts(dds) >= a) >= k
    dds <- dds[keep,]

    # make sure conditions are factors
    dds$Disease <- factor(dds$Disease, levels=c(c1, c2))

    # run analysis
    dds.wald <- DESeq(dds, test="Wald", fitType="local")
    res <- results(dds.wald)

    # get transformed counts
    # rld <- rlog(dds, fitType="local", blind=TRUE)
    # rldf <- data.frame(assay(rld))
    # rldf$gene_id <- rownames(rldf)

    # write table out
    # write.table(rldf, file="export.dir/genes_rlog.tsv", row.names=F, quote=F, sep="\t")

    # plot mean-variance relationship and dispersion
    # estimates

    # rld <- read.csv("export.dir/genes_rlog.tsv", header=T, stringsAsFactors=F, sep="\t")
    # rownames(rld) <- rld$gene_id
    # rld <- rld[,c(1:ncol(rld)-1)]

    # par(mfrow=c(1,2))
    # plotMeanSd(rld)
    # plotDispEsts(dds.lrt)

    results.table <- getResultsTable(res, featureData)
    write.table(results.table, file=results.file, sep="\t", row.names=FALSE, quote=FALSE)
    return(results.table)
    }

runDESeq2PairWithBatch <- function(mat, metadata, geneid2gene_name=NULL, a=0, k=1, tissue="Ileum", results.file=NULL, c1="HEALTHY", c2="PSC/UC"){
    
    metadata <- metadata[metadata$Disease %in% c(c1, c2),]
    metadata <- metadata[metadata$Tissue.location == tissue,]
    mat <- mat[,rownames(metadata)]
    
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                  colData = metadata,
                                  design = ~ Disease + batch)
    
    
    # get gene names
    featureData <- makeFeatureData(geneid2gene_name, mat)
    
    # filtering > a in at least k samples (i.e. all of one condition)
    keep <- rowSums(counts(dds) >= a) >= k
    dds <- dds[keep,]
    
    # make sure conditions are factors
    dds$Disease <- factor(dds$Disease, levels=c(c1, c2))
    
    # run analysis
    dds.lrt <- DESeq(dds, test="LRT", reduced = ~ batch)
    res <- results(dds.lrt)
    
    results.table <- getResultsTable(res, featureData)
    write.table(results.table, file=results.file, sep="\t", row.names=FALSE, quote=FALSE)
    return(results.table)
}

addSigColumn <- function(resultset, l2fold=1){
    
    resultset$sig <- ifelse(resultset$padj < 0.05 & resultset$log2FoldChange > l2fold & !(is.na(resultset$padj)) & !(is.na(resultset$log2FoldChange)), "Up", NA)
    resultset$sig <- ifelse(resultset$padj < 0.05 & resultset$log2FoldChange < -(l2fold) & !(is.na(resultset$padj)) & !(is.na(resultset$log2FoldChange)), "Down", resultset$sig)
    return(resultset)
    }

scaleWithinTissue <- function(mat, metadata, tissue.var="Tissue.location", tissues=c()){
    
    tissues.scaled <- list()
    for(i in 1:length(tissues)){
        m <- mat[metadata[,tissue.var] == tissues[i]]
        m.s <- data.frame(t(apply(m, 1, scale)))
        colnames(m.s) <- colnames(m)
        tissues.scaled[[i]] <- m.s}
    result <- bind_cols(tissues.scaled)
}

