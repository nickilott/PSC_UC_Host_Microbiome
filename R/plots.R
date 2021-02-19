#################################################################################
#################################################################################
# Plotting functions
#################################################################################
#################################################################################

library(ggplot2)
library(reshape)

################
################
################

plotBar <- function(x, colname="ASV"){
    p1 <- ggplot(x, aes_string(x="variable", y="value", fill=colname))
    p2 <- p1 + geom_bar(stat="identity")
    p3 <- p2 + scale_y_continuous(labels = scales::percent)
    p4 <- p2 + theme_bw()
    p5 <- p4 + theme(axis.text.x=element_text(angle=90)) + ylab("% of reads")
    return(p5)
}

################
################
################

plotNumberOfASVs <- function(dat){

    nasvs <- data.frame(colSums(dat > 0))
    nasvs$sample <- rownames(nasvs)
    nasvs.m <- melt(nasvs)

    # number of ASVs
    p1 <- ggplot(nasvs.m, aes(x=sample, y=value))
    p2 <- p1 + geom_bar(stat="identity")
    p3 <- p2 + theme_bw()
    p4 <- p3 + theme(axis.text.x=element_text(angle=90))
    p5 <- p4 + ylab("number of ASVs")
    return(p5)
}

################
################
################

plotSampleComponents <- function(ordinate.output){

    p1 <- ggplot(ordinate.output, aes(x=variable, y=value, colour=covariate, group=sample))
    p2 <- p1 + geom_point() + geom_line()
    p3 <- p2 + facet_wrap(~type, scale="free")
    return(p3)
}

################
################
################

plotPrev <- function(prev.df, as.is=TRUE){

    if (as.is==TRUE){
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

