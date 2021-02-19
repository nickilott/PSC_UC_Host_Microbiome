library(plyr)

buildEnrichedSets <- function(enrichment.list, modules=NULL, fdr=0.05){

    # build a single dataframe with each module and
    # significant pathways from a list of enrichment
    # results output by goseq

    enriched.sets <- list()
    for (i in 1:length(enrichment.list)){
        module <- names(enrichment.list[i])
        enrichment.res <- enrichment.list[[i]]
	
        sig <- enrichment.res[enrichment.res$over_represented_fdr < fdr,]
              if (nrow(sig) == 0){
	                 sig <- data.frame(category=NA,
	                 over_represented_pvalue=NA,
			             under_represented_pvalue=NA,
			             numDEInCat=NA,
			             numInCat=NA,
			             term=NA,
			             ontology=NA,
			             over_represented_fdr=NA)
              }
        sig$module <- module
        sig$tissue <- gsub("\\.[0-9]*", "", module)
        enriched.sets[[i]] <- sig
     }
    enriched.sets <- ldply(enriched.sets, data.frame)
    if(!(is.null(modules))){
      enriched.sets <- enriched.sets[enriched.sets$module %in%  modules,]
    }
    
    # for plotting purposes need to transform pathway names into
    # numbers so that they are just placed omn an arbitrary axis
    pathways <- unique(enriched.sets$term)
    pathways <- seq(1, length(pathways), 1)
    names(pathways) <- unique(enriched.sets$term)
    enriched.sets$yaxis <- unlist(pathways[enriched.sets$term])

    return(list(enriched.sets=na.omit(enriched.sets), labels=pathways[!(is.na(names(pathways)))]))
    
}

#' Plot module enrichments
#'
#' Plot the modules and enriched gene sets
#' 
#' @param enriched.sets full data frame of modules and their sigificantly enriched pathways
#' @return ggplot object 
#' @import plyr ggplot2
#' @examples
#' plotModuleEnrichments(enriched.sets)
#' @export

plotModuleEnrichments <- function(enriched.sets, labels){

    p1 <- ggplot(enriched.sets, aes(x=module, y=yaxis, size=-log10(over_represented_fdr), color=tissue)) + geom_point()
    p2 <- p1 + theme_bw() + scale_y_continuous(breaks=unlist(labels, use.names=FALSE), labels=names(labels))
    p2 <- p2 + theme(axis.text.x=element_text(angle=90)) + scale_color_manual(values=c("red3", "blue3"))
return(p2)
    
}
