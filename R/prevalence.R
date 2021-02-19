###############
###############
###############

# prevalence
pprev <- function(counts){

    nsamples <- ncol(counts)
    prev <- (rowSums(counts > 0)/nsamples)*100
    prev.df <- data.frame(Taxon=rownames(counts), Prevalence=prev)
    return(prev.df)
    }

###############
###############
###############

# plotting
plotPrev <- function(prev.df, as.is=TRUE){

    if (as.is==TRUE){
        prev.df <- prev.df[order(prev.df$Prevalence, decreasing=TRUE),]
        p1 <- ggplot(prev.df, aes(x=Taxon, y=Prevalence))
        p2 <- p1 + geom_bar(stat="identity")
        p3 <- p2 + theme_bw()
        p4 <- p3 + theme(axis.text.x=element_text(angle=90))
    }else{
        p1 <- ggplot(prev.df, aes(x=Prevalence))
        p2 <- p1 + geom_histogram()
        p3 <- p2 + theme_bw()
        p4 <- p3
    }
    return(p4)
}
