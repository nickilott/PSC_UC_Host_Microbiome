
threewayvenn <- function(df){

    list1 <- df$test.id[df[,2]  < 0.05]
    list2 <- df$test.id[df[,3]  < 0.05]
    list3 <- df$test.id[df[,4]  < 0.05]

    names.list <- unlist(strsplit(colnames(df), "*_padj.*"))
    names.list <- names.list[c(2,3,4)]

    name1 <- names.list[1]
    name2 <- names.list[2]
    name3 <- names.list[3]

    tovenn <- list(list1, list2, list3)
    names(tovenn) <- c(name1, name2, name3)

    v <- venn.diagram(tovenn, filename=NULL,
                      fill=c("red4", "blue4", "green4"),
                      alpha=c(0.5,0.5,0.5),
                      fontfamily=rep("sans", 7),
                      cat.fontfamily=rep("sans", 3),
                      cat.pos=c(0, 180, 180),
                      cat.just=list(c(1,1), c(0,1), c(1,1)),
                      margin=c(0.2,0.2,0.2,0.2))
    grid.draw(v)
}
																						    