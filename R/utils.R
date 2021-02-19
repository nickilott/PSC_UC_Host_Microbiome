#######################################
# utlity functions
#######################################

getNumberOfASVs <- function(dat){

    # return the number of ASVs per
    # sample
    nasvs <- data.frame(colSums(dat > 0))
    colnames(nasvs) <- "count"
    return(nasvs)
}

getLibrarySize <- function(dat){

    # return read count for each sample
    library.size <- data.frame(colSums(dat))
    colnames(library.size) <- "nreads"
    return(library.size)
}

getShortNames <- function(longnames, type="ASV", level="phylum"){

    choices <- c("phylum", "class", "order", "family", "genus", "species")
    if (!(level %in% choices)){
        stop("level must be one of phylum, class, order, family, genus, species")
    }


    shortnames <- unlist(strsplit(longnames, ";"))
    end <- length(shortnames)

    if (type=="ASV"){
       by = 6
       }
    else {
       by = 5
    }

    # return vector of shortened names
    if (level == "phylum"){
        start <- 1
        pattern <- "*p__"
    }
    else if (level == "class"){
        start <- 2
        pattern <- "c__"
    }
    else if (level == "order"){
        start <- 3
        pattern <- "o__"
    }
    else if (level == "family"){
        start <- 4
        pattern <- "f__"
    }
    else if (level == "genus"){
        start <- 5
        pattern <- "g__"
    }
    else if (level == "species"){
        start <- 6
        pattern <- "s__"
    }
    
    shortnames <- shortnames[seq(start, end, by)]
    shortnames <- gsub(pattern, "", shortnames)

    # hackaround some names having "genus cluster" nomenclature
    # for clostridium genera
    shortnames <- gsub("Clostridium .*", "Clostridium", shortnames)
    return(shortnames)
}

