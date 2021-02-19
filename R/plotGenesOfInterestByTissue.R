# require plotGeneOfInterest from MIGTranscriptomeExplorer

# Plot genes of interest - returns a list of grobs
plotGenesOfInterestByTissue <- function(gene.names,
                                        results.table,
                                        expression.matrix,
                                        metadata,
                                        colours,
                                        export.dir="export.dir"){
  grobs <- list()
  for (i in 1:length(gene.names)){
    gene <- gene.names[[i]]
    gene_id <- results.table$gene_id[results.table$gene_name == gene]
    
    p <- plotGeneOfInterest(gene, expression.matrix[gene_id,],
                            metadata,
                            variable="group") + scale_color_manual(values=colours) + geom_segment(aes(x=0, y=max(value)+1, xend=3.5, yend=max(value)+1), color=blues9[3], size=7) + geom_segment(aes(x=3.5, y=max(value)+1, xend=6.5, yend=max(value)+1), color=blues9[6], size=7) + geom_segment(aes(x=6.5, y=max(value)+1, xend=10, yend=max(value)+1), color=blues9[9], size=7) + geom_vline(xintercept = c(3.5,6.5), linetype="dashed") + theme(text=element_text(size=12)) + theme(legend.position="none")
    
    ggsave(paste0(export.dir, "/", gene, ".pdf"), height=3, width=4)
    # ggsave(paste0(export.dir, gene, ".png"), height=4, width=6)
    grobs[[i]] <- p
  }
  return(grobs)
  }