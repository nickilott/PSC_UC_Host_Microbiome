##############################################
##############################################
##############################################
# calculate within individual dissimilarity
##############################################
##############################################
##############################################


library(phyloseq)
library(dplyr)

##############################################
##############################################
##############################################

withinIndividualDissimilarity <- function(asv.relab, metadata, collection.method=NULL){

    if (is.null(collection.method)){
        metadata <- metadata}
    else{
        metadata <- metadata[metadata$Sample.collection.method == collection.method,]
	}
    sample.ids <- rownames(metadata)
    asvs <- asv.relab[,sample.ids]

    # get the patient ids to iterate over
    patient.ids <- unique(metadata$Patient.ID)

    # iterate over patients and calculate
    # the bray.curtis distance

    result.p <- list()
    for (i in 1:length(patient.ids)){
        metadata.p <- metadata[metadata$Patient.ID == patient.ids[i],]

        # don't compute on single samples
        if (nrow(metadata.p) == 1){
	    next}

        asvs.p <- data.frame(asvs[,rownames(metadata.p)])
        rownames(asvs.p) <- rownames(asvs)
        colnames(asvs.p) <- rownames(metadata.p)
        asvs.p <- phyloseq::otu_table(asvs.p, taxa_are_rows=TRUE)
        metadata.p <- sample_data(metadata.p)
        phyob.p <- phyloseq::merge_phyloseq(asvs.p, metadata.p)
	diss <- phyloseq::distance(phyob.p, method="bray")
        diss <- as.data.frame(as.matrix(diss))
        res <- data.frame(dissimilarity = diss[2:nrow(diss),1])
        res$sample <- patient.ids[i]
        result.p[[i]] <- res
    }
    res.all <- bind_rows(result.p)
    return(res.all)
}

##############################################
##############################################
##############################################


betweenIndividualDissimilarity <- function(asv.relab, metadata, collection.method="Brush", tissue.location="Caecum"){

    metadata <- metadata[metadata$Sample.collection.method == collection.method
                         & metadata$Tissue.location == tissue.location,]
                         
    sample.ids <- rownames(metadata)
    asvs <- asv.relab[,sample.ids]

    asvs <- phyloseq::otu_table(asvs, taxa_are_rows=TRUE)
    metadata <- phyloseq::sample_data(metadata)
    phyob <- phyloseq::merge_phyloseq(asvs, metadata)

    diss <- phyloseq::distance(phyob, method="bray")
    res <- data.frame(dissimilarity = as.vector(diss))
    res$sample <- paste0(collection.method, ":", tissue.location)
    return(res)
    }