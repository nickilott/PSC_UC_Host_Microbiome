###############################################
###############################################
###############################################
# Functions for running anovas for differential
# abundance analysis
###############################################
###############################################
###############################################

library(dplyr)
library(car)

changeFormat <- function(x){
    return(as.numeric(as.character(x)))
}

aovStandard <- function(row, metadata, independent.variable=NULL, covariate=NULL){

    asv <- rownames(row)
    iv <- unlist(metadata[,independent.variable])
    aov.result <- aov(unlist(row) ~ as.factor(iv))
    if (!(is.null(covariate))){
        covariate <- metadata[,covariate]
        aov.result <- aov(unlist(row) ~ as.factor(iv) + as.factor(covariate))
        aov.result <- car::Anova(aov.result, type="III")
        print(aov.result)
    }
    
    tukey.res <- as.data.frame(TukeyHSD(aov.result)[[1]])
    tukey.res <- data.frame(t(tukey.res))[4,]
    rownames(tukey.res) <- asv
    
    aov.result.tab <- data.frame(summary(aov.result)[[1]])
    aov.result.tab <- aov.result.tab[1,]
    aov.result.tab$test.id <- asv
    rownames(aov.result.tab) <- asv
    
    aov.result.tab <- cbind(aov.result.tab, tukey.res)
    return(aov.result.tab)
}

###############################################
###############################################
###############################################

runAovStandard <- function(mat, metadata, independent.variable=NULL, covariate=NULL){

    #aov.results <- matrix(nrow=nrow(mat), ncol=6)
    aov.results <- list()
    for (i in 1:nrow(mat)){
        result <- aovStandard(mat[i,], metadata, independent.variable, covariate=covariate)
        aov.results[[i]] <- as.data.frame(result)
    }
    
    #aov.results <- as.data.frame(aov.results)
    aov.results <- bind_rows(aov.results)
    
    colnames(aov.results)[1:6] <- c("Df", "Sum.Sq", "Mean.Sq", "F.value", "pvalue", "test.id")

    aov.results$padj <- p.adjust(aov.results$pvalue, method="fdr")
    return(aov.results)
}


###############################################
###############################################
###############################################

aovInteraction <- function(row, metadata, independent.variable=NULL, covariate=NULL){

    asv <- rownames(row)
    iv <- unlist(metadata[,independent.variable])
    cv <- unlist(metadata[,covariate])
    
    aov.result <- aov(unlist(row) ~ as.factor(iv)*as.factor(cv))
    aov.result.tab <- data.frame(summary(aov.result)[[1]])

    main.iv.f <- paste0("main.", independent.variable, ".F")
    main.iv.p <- paste0("main.", independent.variable, ".Pval")
    main.cv.f <- paste0("main.", covariate, ".F")
    main.cv.p <- paste0("main.", covariate, ".Pval")
    interaction.f <- "interaction.F"
    interaction.p <- "interaction.Pval"

    aov.result.tab <- data.frame(main.iv.f=aov.result.tab[1,4],
                                 main.iv.p=aov.result.tab[1,5],
				 main.cv.f=aov.result.tab[2,4],
				 main.cv.p=aov.result.tab[2,5],
				 interaction.f=aov.result.tab[3,4],
				 interaction.p=aov.result.tab[3,5])

    aov.result.tab$test.id <- asv
    return(aov.result.tab)
}

###############################################
###############################################
###############################################

runAovInteraction <- function(mat, metadata, independent.variable=NULL, covariate=NULL){

    aov.results <- matrix(nrow=nrow(mat), ncol=7)
    for (i in 1:nrow(mat)){
        result <- aovInteraction(mat[i,], metadata, independent.variable, covariate)
        aov.results[i,1] <- result[,1]
        aov.results[i,2] <- result[,2]
        aov.results[i,3] <- result[,3]
        aov.results[i,4] <- result[,4]
        aov.results[i,5] <- result[,5]
        aov.results[i,6] <- result[,6]
    }
    
    aov.results <- as.data.frame(aov.results)
    main.iv.f <- paste0("main.", independent.variable, ".F")
    main.iv.p <- paste0("main.", independent.variable, ".Pval")
    main.cv.f <- paste0("main.", covariate, ".F")
    main.cv.p <- paste0("main.", covariate, ".Pval")
    interaction.f <- "interaction.F"
    interaction.p <- "interaction.Pval"

    colnames(aov.results) <- c(main.iv.f, main.iv.p, main.cv.f, main.cv.p, interaction.f, interaction.p)

    main.iv.padj <- paste0("main.", independent.variable, ".padj")
    main.cv.padj <- paste0("main.", covariate, ".padj")
    interaction.padj <- "interaction.padj"

    aov.results[,main.iv.padj] <- p.adjust(aov.results[,main.iv.p], method="fdr")
    aov.results[,main.cv.padj] <- p.adjust(aov.results[,main.cv.p], method="fdr")
    aov.results[,interaction.padj] <- p.adjust(aov.results[,interaction.p], method="fdr")

    rownames(aov.results) <- rownames(mat)
    return(aov.results)
}

