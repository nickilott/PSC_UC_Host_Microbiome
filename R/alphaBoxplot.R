

alphaBoxplot <- function(diversity.df){
    p1 <- ggplot(diversity.df, aes(x=collection.method, y=alpha.diversity, colour=library.size))
    p2 <- p1 + geom_boxplot(outlier.alpha=0)
    p3 <- p2 + geom_jitter(width=0.2)
    p4 <- p3 + facet_wrap(~tissue.location)
    p5 <- p4 + theme_bw()
    return(p5)
    }