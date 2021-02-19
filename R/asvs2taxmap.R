############################################
############################################
############################################
# Take a list of annotated ASVs and convert
# to taxmap table for input into phyloseq
############################################
############################################
############################################

asvs2taxmap <- function(asvs){

    taxmap <- matrix(nrow=length(asvs), ncol=6)
    for (i in 1:length(asvs)){
        taxonomy <- unlist(strsplit(asvs[i], ";"))
	taxonomy[1] <- gsub("ASV[0-9]*:", "", taxonomy[1])
	taxmap[i,] <- taxonomy
    }
    	
    taxmap <- as.data.frame(taxmap)
    colnames(taxmap) <- c("Phylum",
                          "Class",
			  "Order",
			  "Family",
			  "Genus",
			  "Species")
    rownames(taxmap) <- asvs
    taxmap[] <- lapply(taxmap, as.character)
    return(taxmap)
}

