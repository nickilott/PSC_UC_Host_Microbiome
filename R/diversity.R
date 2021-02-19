#######################################################
#######################################################
#######################################################
# A set of functions for using with diversity metrics
#######################################################
#######################################################
#######################################################

library(vegan)

multiFactorKruskalTest <- function(data){

    #####################################################################
    # this function takes a data frame as the argument. The data frame
    # must have a column that is called "alpha.diversity". This is just so
    # that any of the diversity metrics can be used from phyloseq
    ######################################################################

    results <- matrix(nrow=ncol(data), ncol=3)
    colnames(results) <- c("~factor", "chi.squared", "p.value")
    for (i in 1:ncol(data)){
        if (colnames(data)[i] == "alpha.diversity"){next}
	res <- kruskal.test(alpha.diversity ~ as.factor(get(colnames(data)[i])), data=data)
	chi.squared <- res$statistic[[1]]
	pvalue <- res$p.value

        # add the values to the
	results[i,1] <- colnames(data)[i]
	results[i,2] <- chi.squared
	results[i,3] <- pvalue

    }
    return(na.omit(as.data.frame(results)))
}

#######################################################
#######################################################
#######################################################

pairwiseAdonis <- function(dat, groups=NULL, method="bray", nperm=1000){

    # NB variable to test differences must be called "group"
    # NB groups must be ordered correctly with input data frame
    # there is no test for this

    if (is.null(groups)){
        stop("must provide groups vector metadata")
    }

    # maker sure groups are characters to
    # begin with at least
    groups$group <- as.character(groups$group)

    # get the pairwise comparisons
    groups.to.test <- as.character(unique(groups$group))
    groups.to.test <- data.frame(t(combn(groups.to.test, 2)))

    # result container
    results <- matrix(nrow = nrow(groups.to.test), ncol=3)
    labels <- c()
    for (i in 1:nrow(groups.to.test)){
        to.take <- which(groups$group == groups.to.test[i,1] | groups$group == groups.to.test[i,2])

        # labels for output
        label <- paste0(groups.to.test[i,1], "-", groups.to.test[i,2])
        labels <- append(labels, label)

        dat.subset <- dat[,to.take]
        group.subset <- data.frame(group=groups[to.take,])

        adonis.result <- adonis(t(dat.subset) ~ as.factor(group), data=group.subset, method=method, permutations=nperm)
	adonis.result <- data.frame(adonis.result$aov.tab)
        results[i, 1] <- adonis.result$R2[1]
        results[i, 2] <- adonis.result$F.Model[1]
        results[i, 3] <- adonis.result$Pr..F.[1]
    }
    colnames(results) <- c("R2", "F", "P.value")
    rownames(results) <- labels
    return(as.data.frame(results))
    }