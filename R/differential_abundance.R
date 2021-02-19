###########################################
###########################################
###########################################
# DESeq2 analysis at multiple levels
###########################################
###########################################
###########################################

library(dplyr)
#source("/gfs/devel/nilott/NGSKit/R/deseq2_helper.R")


#####################
#####################
#####################

multiDE <- function(counts.files, metadata, a=10, k=10, feature.column=1, model.formula=~condition, reduced.model=~1){

    # metadata is a df
    # counts.files is a list of filenames
    # covariates must be present in metadata file

    # at the moment uses hard-coded parameters
    # so won't neccesarily work in all settings

    # NOTE: metadata should have sample names as rownames
    # NOTE: counts files should have samples as colnames
    # NOTE: colnames(counts) == rownames(metadata)

    # NOTE: the level must be specified in the counts file
    # as level_abundance.tsv

    # Returns a dataframe of results at each level

    result.set <- list()
    for (i in 1:length(counts.files)){
        countData <- read.csv(counts.files[i], header=T, stringsAsFactors=F, sep="\t", row.names=feature.column, quote="")

        # filter
	countData <- countData[rowSums(countData >= k) >= a,]

        # in case files are in subdirectory
        level <- tail(unlist(strsplit(counts.files[i], "/")), n=1)
        level <- unlist(strsplit(level, "_"))[1] 

        samples <- intersect(colnames(countData), rownames(metadata))

        # make sure order of counts and metadata are the same
        metadata <- metadata[samples,]
	countData <- countData[,samples]

        dds <- DESeqDataSetFromMatrix(countData = countData,
                                      colData = metadata,
			              design = model.formula)

        dds.lrt <- DESeq(dds, test="LRT", fitType="local", reduced=reduced.model) 
        res <- results(dds.lrt)
        res2 <- data.frame(res@listData)
	res2$level <- level
        res2$test_id <- rownames(res)
        result.set[[i]] <- res2
    }
    results <- bind_rows(result.set)
}

#####################
#####################
#####################

# multiDEMetagenomeSeq <- function(counts.files, metadata, a=10, k=10, feature.column=1, model.formula=~condition, reduced.model=~1){

#     # metadata is a df
#     # counts.files is a list of filenames
#     # covariates must be present in metadata file

#     # at the moment uses hard-coded parameters
#     # so won't neccesarily work in all settings

#     # NOTE: metadata should have sample names as rownames
#     # NOTE: metadata should only have the relevant covariates
#     # NOTE: counts files should have samples as colnames
#     # NOTE: colnames(counts) == rownames(metadata)

#     # NOTE: the level must be specified in the counts file
#     # as level_abundance.tsv

#     # Returns a dataframe of results at each level

#     result.set <- list()
#     for (i in 1:length(counts.files)){
#         countData <- read.csv(counts.files[i], header=T, stringsAsFactors=F, sep="\t", row.names=feature.column, quote="")

#         # filter
# 	countData <- countData[rowSums(countData >= k) >= a,]

#         # in case files are in subdirectory
#         level <- tail(unlist(strsplit(counts.files[i], "/")), n=1)
#         level <- unlist(strsplit(level, "_"))[1] 

#         samples <- intersect(colnames(countData), rownames(metadata))

#         # make sure order of counts and metadata are the same
#         metadata <- metadata[samples,]
# 	countData <- countData[,samples]

#         phenotypeData <- AnnotatedDataFrame(metadata)
#         annoCountData <- AnnotatedDataFrame(countData)

#         mrobj <- newMRexperiment(countData, phenoData=phenotypeData, featureData=annoCountData)

#         # normalisation
#         p <- cumNormStatFast(mrobj)
#         mrobj.norm <- cumNorm(mrobj, p=p)

#         # Use fitZig function for including covariates

#         mod <- model.matrix(model.formula, data=pdata(mrobj))
#         res <- fitZig(mrobj.norm, mod)
#         results.table <- MRcoefs(res, number = nrow(countData))




#         dds <- DESeqDataSetFromMatrix(countData = countData,
#                                       colData = metadata,
# 			              design = model.formula)

#         dds.lrt <- DESeq(dds, test="LRT", fitType="local", reduced=reduced.model) 
#         res <- results(dds.lrt)
#         res2 <- data.frame(res@listData)
# 	res2$level <- level
#         res2$test_id <- rownames(res)
#         result.set[[i]] <- res2
#     }
#     results <- bind_rows(result.set)
# }


#####################
#####################
#####################

getNumberDiff <- function(multiDE.result, padj=0.05){

    levels <- unique(multiDE.result$level)
    ndiff.all <- data.frame(matrix(nrow=1, ncol=length(levels)))
    colnames(ndiff.all) <- levels
    
    for (i in 1:length(levels)){
        dat <- multiDE.result[multiDE.result$level == levels[i],]
	ndiff <- nrow(dat[dat$padj < padj,])
	ndiff.all[,levels[i]] <- ndiff
	}
    return(ndiff.all)
}