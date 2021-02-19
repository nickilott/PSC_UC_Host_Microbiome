################################################################
################################################################
################################################################
# plot library size vs. number of ASVs in brushes vs biopsies
################################################################
################################################################
################################################################

librarySizeVsAlphaDiversity <- function(dat){

    cor.brush <- cor.test(dat$library.size[dat$collection.method=="Brush"],
                          dat$alpha.diversity[dat$collection.method=="Brush"])
    cor.brush.p <- cor.brush$p.value
    cor.brush.r2 <- cor.brush$estimate^2

    cor.biopsy <- cor.test(dat$library.size[dat$collection.method=="Biopsy"],
                          dat$alpha.diversity[dat$collection.method=="Biopsy"])
    cor.biopsy.p <- cor.biopsy$p.value
    cor.biopsy.r2 <- cor.biopsy$estimate^2

    label.brush <- paste0("Brush: r2 = ", round(cor.brush.r2, 2), ", ", "p = ", round(cor.brush.p, 2))
    label.biopsy <- paste0("Biopsy: r2 = ", round(cor.biopsy.r2, 2), ", ", "p = ", round(cor.biopsy.p, 2))

    p1 <- ggplot(dat, aes(x=library.size, y=alpha.diversity, colour=collection.method))
    p2 <- p1 + geom_point()
    p3 <- p2 + theme_bw()
    p4 <- p3 + scale_colour_manual(values=c("grey", "purple"))
    p5 <- p4 + geom_smooth(method="lm")
    p6 <- p5 + annotate("text", x=10, y=800, label=label.brush)
    p7 <- p6 + annotate("text", x=10, y=750, label=label.biopsy)
    return(p7)
    }