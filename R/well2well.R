#########################################################
#########################################################
#########################################################
# Assess the well to well contamination based on position
# of positive control samples on a plate
#########################################################
#########################################################
#########################################################

library(ggplot2)

getWellXaxis <- function(well.position){

    xaxis <- unlist(strsplit(well.position, "[A-Z]"))[seq(2, length(well.position)*2, 2)]
    return(as.numeric(xaxis))
}

#########################################################
#########################################################
#########################################################

getWellYaxis <- function(well.position){

    yaxis <- unlist(strsplit(well.position, "[0-9]*"))[seq(1, length(well.position)*2, 2)]
    return(yaxis)
}

#########################################################
#########################################################
#########################################################

letter2number <- function(){

    yaxis.map <- seq(1, 8, 1)
    names(yaxis.map) <- c("A", "B", "C", "D", "E", "F", "G", "H")
    return(yaxis.map)
}

#########################################################
#########################################################
#########################################################

plotWell2Well <- function(layout){

    # columns in the data frame are:
    # 1. well.position
    # 2. relative.abundance
    # 3. group
    # 4. plate

    layout$xaxis <- getWellXaxis(layout$well.position)
    layout$yaxis <- getWellYaxis(layout$well.position)

    layout$yaxis <- factor(layout$yaxis, levels=c("H", "G", "F", "E", "D", "C", "B", "A"))

    # plot
    p1 <- ggplot(layout, aes(x=yaxis, y=xaxis, colour=group, size=relative.abundance))
    p2 <- p1 + geom_point()
    p3 <- p2 + coord_flip() + theme_bw()
    p4 <- p3 + xlab("") + ylab("")
    p5 <- p4 + scale_y_continuous(breaks = seq(1,12,1), labels = seq(1,12,1))
    p6 <- p5 + facet_wrap(~plate)
    return(p6)

}