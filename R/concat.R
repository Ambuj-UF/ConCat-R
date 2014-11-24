################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                       {Mrinal Mishra, VTT Technical Research Center, Finland}                                #
# "read.nex" and write.nex are Simplified NEXUS data parser by Johan Nylander, nylander @ scs.fsu.edu          #
# "read.nex" is a Simplified NEXUS data parser by Johan Nylander, nylander @ scs.fsu.edu                       #
# WARNING: This is parser reads a restricted nexus format, see the link below for details                      #
# http://svitsrv25.epfl.ch/R-doc/library/ape/html/read.nexus.data.html &                                       #
# http://www.inside-r.org/packages/cran/ape/docs/write.nexus.data                                              #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################



library(seqinr)


read.nexus.data <- function (file)
{
    find.ntax <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bntax", x[i], ignore.case = TRUE))) {
                ntax <- as.numeric(sub("(.+?)(ntax\\s*\\=\\s*)(\\d+)(.+)",
                "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        ntax
    }
    
    find.nchar <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bnchar", x[i], ignore.case = TRUE))) {
                nchar <- as.numeric(sub("(.+?)(nchar\\s*\\=\\s*)(\\d+)(.+)",
                "\\3", x[i], perl = TRUE, ignore.case = TRUE))
                break
            }
        }
        nchar
    }
    
    find.matrix.line <- function (x)
    {
        for (i in 1:NROW(x)) {
            if(any(f <- grep("\\bmatrix\\b", x[i], ignore.case = TRUE))) {
                matrix.line <- as.numeric(i)
                break
            }
        }
        matrix.line
    }
    
    trim.whitespace <- function (x)
    {
        gsub("\\s+", "", x)
    }
    
    trim.semicolon <- function (x)
    {
        gsub(";", "", x)
    }
    
    X <- scan(file = file, what = character(), sep = "\n",
    quiet = TRUE, comment.char = "[", strip.white = TRUE)
    ntax <- find.ntax(X)
    nchar <- find.nchar(X)
    matrix.line <- find.matrix.line(X)
    start.reading <- matrix.line + 1
    Obj <- list()
    length(Obj) <- ntax
    i <- 1
    pos <- 0
    tot.nchar <- 0
    tot.ntax <- 0
    
    for (j in start.reading:NROW(X)) {
        Xj <- trim.semicolon(X[j])
        if(Xj == "") {
            break
        }
        if(any(jtmp <- grep("\\bend\\b", X[j], perl = TRUE, ignore.case = TRUE))) {
            break
        }
        ts <- unlist(strsplit(Xj, "(?<=\\S)(\\s+)(?=\\S)", perl = TRUE))
        if (length(ts) > 2) {
            stop("nexus parser does not handle spaces in sequences or taxon names (ts>2)")
        }
        if (length(ts) !=2) {
            stop("nexus parser failed to read the sequences (ts!=2)")
        }
        Seq <- trim.whitespace(ts[2])
        Name <- trim.whitespace(ts[1])
        nAME <- paste(c("\\b", Name, "\\b"), collapse = "")
        if (any(l <- grep(nAME, names(Obj)))) {
            tsp <- strsplit(Seq, NULL)[[1]]
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[l]][p] <- tsp[k]
                chars.done <- k
            }
        }
        else {
            names(Obj)[i] <- Name
            tsp <- strsplit(Seq, NULL)[[1]]
            for (k in 1:length(tsp)) {
                p <- k + pos
                Obj[[i]][p] <- tsp[k]
                chars.done <- k
            }
        }
        
        tot.ntax <- tot.ntax + 1
        if (tot.ntax == ntax) {
            i <- 1
            tot.ntax <- 0
            tot.nchar <- tot.nchar + chars.done
            if (tot.nchar == nchar*ntax) {
                print("ntot was more than nchar*ntax")
                break
            }
            pos <- tot.nchar
        }
        else {
            i <- i + 1
        }
    }
    
    if (tot.ntax != 0) {
        cat("ntax:",ntax,"differ from actual number of taxa in file?\n")
        stop("nexus parser did not read names correctly (tot.ntax!=0)")
    }
    
    for (i in 1:length(Obj)) {
        if (length(Obj[[i]]) != nchar) {
            cat(names(Obj[i]),"has",length(Obj[[i]]),"characters\n")
            stop("nchar differ from sequence length (length(Obj[[i]])!=nchar)")
        }
    }
    Obj <- lapply(Obj, tolower)
    Obj
}

# Nexus data parser.
#
# Version: 09/13/2006 09:06:33 AM CEST
#
# By:      Johan Nylander, nylander @ scs.fsu.edu
#
# TODO:    Standard data, mixed data, nice indent
#------------------------------------------------------------------

"write.nex" <- function (x, file, format = "dna", datablock = TRUE,
interleaved = TRUE, charsperline = NULL,
gap = NULL, missing = NULL)
{
    
    indent          <- "  "  # Two blanks
    maxtax          <- 5     # Max nr of taxon names to be printed on a line
    defcharsperline <- 80    # Default nr of characters per line if interleaved
    defgap          <- "-"   # Default gap character
    defmissing      <- "?"   # Default missing data character
    
    ntax <- length(x)
    nchars <- length(x[[1]])
    zz <- file(file, "w")
    
    if (is.null(names(x))) {
        names(x) <- as.character(1:ntax)
    }
    
    "fcat" <- function (..., file = zz)
    {
        cat(..., file = file, sep = "", append = TRUE)
    }
    
    "find.max.length" <- function (x)
    {
        max <- 0
        for (i in 1:length(x)) {
            val <- length((strsplit(x[i], split = NULL))[[1]])
            if (val > max) {
                max <- val
            }
        }
        max
    }
    
    "print.matrix" <- function(x, dindent = "    ")
    {
        Names <- names(x)
        printlength <- find.max.length(Names) + 2
        if (interleaved == FALSE) {
            for (i in 1:length(x)) {
                sequence <- paste(x[[i]], collapse = "")
                taxon <- Names[i]
                thestring <- sprintf("%-*s%s%s", printlength, taxon, dindent, sequence)
                fcat(indent, indent, thestring, "\n")
            }
        }
        else {
            ntimes <- ceiling(nchars/charsperline)
            start <- 1
            end <- charsperline
            for (j in 1:ntimes) {
                for (i in 1:length(x)) {
                    sequence <- paste(x[[i]][start:end], collapse = "")
                    taxon <- Names[i]
                    thestring <- sprintf("%-*s%s%s", printlength, taxon, dindent, sequence)
                    fcat(indent, indent, thestring, "\n")
                }
                if (j < ntimes) {
                    fcat("\n")
                }
                start <- start + charsperline
                end <- end + charsperline
                if (end > nchars) {
                    end <- nchars
                }
            }
        }
    }
    
    fcat("#NEXUS\n[Data written by write.nexus.data.R,", " ", date(),"]\n")
    
    NCHAR <- paste("NCHAR=", nchars, sep = "")
    NTAX <- paste("NTAX=", ntax, sep = "")
    
    if (format == "dna") {
        DATATYPE <- "DATATYPE=DNA"
    }
    if (format == "protein") {
        DATATYPE <- "DATATYPE=PROTEIN"
    }
    
    if (is.null(charsperline)) {
        if (nchars < defcharsperline) {
            charsperline <- nchars
            interleaved <- FALSE
        }
        else {
            if (nchars > defcharsperline) {
                charsperline <- defcharsperline
            }
        }
    }
    
    if (is.null(missing)) {
        MISSING <- paste("MISSING=", defmissing, sep = "")
    }
    else {
        MISSING <- paste("MISSING=", missing, sep = "")
    }
    
    if (is.null(gap)) {
        GAP <- paste("GAP=", defgap, sep = "")
    }
    else {
        GAP <- paste("GAP=", gap, sep = "")
    }
    
    if (interleaved == TRUE) {
        INTERLEAVE <- "INTERLEAVE=YES"
    }
    if (interleaved == FALSE) {
        INTERLEAVE <- "INTERLEAVE=NO"
    }
    
    if (datablock == TRUE) {
        fcat("BEGIN DATA;\n")
        fcat(indent,"DIMENSIONS", " ", NTAX, " ", NCHAR, ";\n")
        if (format %in% c("dna", "protein")) {
            fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", GAP, " ", INTERLEAVE, ";\n") # from FranÁois Michonneau (2009-10-02)
        }
        fcat(indent,"MATRIX\n")
        print.matrix(x)
        fcat(indent, ";\n")
        fcat("END;\n\n")
    }
    else {
        fcat("BEGIN TAXA;\n")
        fcat(indent, "DIMENSIONS", " ", NTAX, ";\n")
        fcat(indent, "TAXLABELS\n")
        fcat(indent, indent)
        j <- 0
        for (i in 1:ntax) {
            fcat(names(x[i]), " ")
            j <- j + 1
            if (i == ntax) {
                fcat("\n", indent, ";\n")
            }
            else {
                if (j == maxtax) {
                    fcat("\n", indent, indent)
                    j <- 0
                }
            }
        }
        fcat("END;\n\n")
        fcat("BEGIN CHARACTERS;\n")
        fcat(indent, "DIMENSIONS", " ", NCHAR, ";\n")
        if (format %in% c("dna", "protein")) {
            fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", DATATYPE, " ", INTERLEAVE, ";\n")
        }
        fcat(indent,"MATRIX\n")
        print.matrix(x)
        fcat(indent, ";")
        fcat("\nEND;\n\n")
    }
    close(zz)
}

##   Write DNA Sequences in a File

## Copyright 2003-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

write.dna <- function(x, file, format = "interleaved", append = FALSE,
nbcol = 6, colsep = " ", colw = 10, indent = NULL,
blocksep = 1)
{
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    phylip <- if (format %in% c("interleaved", "sequential")) TRUE else FALSE
    if (inherits(x, "DNAbin")) x <- as.character(x)
    aligned <- TRUE
    if (is.matrix(x)) {
        N <- dim(x)
        S <- N[2]
        N <- N[1]
        xx <- vector("list", N)
        for (i in 1:N) xx[[i]] <- x[i, ]
        names(xx) <- rownames(x)
        x <- xx
        rm(xx)
    } else {
        N <- length(x)
        S <- unique(unlist(lapply(x, length)))
        if (length(S) > 1) aligned <- FALSE
    }
    if (is.null(names(x))) names(x) <- as.character(1:N)
    if (is.null(indent))
    indent <- if (phylip) 10 else  0
    if (is.numeric(indent))
    indent <- paste(rep(" ", indent), collapse = "")
    if (format == "interleaved") {
        blocksep <- paste(rep("\n", blocksep), collapse = "")
        if (nbcol < 0) format <- "sequential"
    }
    zz <- if (append) file(file, "a") else file(file, "w")
    on.exit(close(zz))
    if (phylip) {
        if (!aligned)
        stop("sequences must have the same length for
        interleaved or sequential format.")
        cat(N, " ", S, "\n", sep = "", file = zz)
        if (nbcol < 0) {
            nb.block <- 1
            nbcol <- totalcol <- ceiling(S/colw)
        } else {
            nb.block <- ceiling(S/(colw * nbcol))
            totalcol <- ceiling(S/colw)
        }
        ## Prepare the sequences in a matrix whose elements are
        ## strings with `colw' characters.
        SEQ <- matrix("", N, totalcol)
        for (i in 1:N) {
            X <- paste(x[[i]], collapse = "")
            for (j in 1:totalcol)
            SEQ[i, j] <- substr(X, 1 + (j - 1)*colw, colw + (j - 1)*colw)
        }
        ## Prepare the names so that they all have the same nb of chars
        max.nc <- max(nchar(names(x)))
        ## always put a space between the sequences and the taxa names
        fmt <- paste("%-", max.nc + 1, "s", sep = "")
        names(x) <- sprintf(fmt, names(x))
    }
    switch(format, "interleaved" = {
        ## Write the first block with the taxon names
        colsel <- if (nb.block == 1) 1:totalcol else 1:nbcol
        for (i in 1:N) {
            cat(names(x)[i], file = zz)
            cat(SEQ[i, colsel], sep = colsep, file = zz)
            cat("\n", file = zz)
        }
        ## Write eventually the other blocks
        if (nb.block > 1) {
            for (k in 2:nb.block) {
                cat(blocksep, file = zz)
                endcolsel <- if (k == nb.block) totalcol else nbcol + (k - 1)*nbcol
                for (i in 1:N) {
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1)*nbcol):endcolsel], sep = colsep, file = zz)
                    cat("\n", file = zz)
                }
            }
        }
        
    }, "sequential" = {
        if (nb.block == 1) {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, ], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        } else {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
                for (k in 2:nb.block) {
                    endcolsel <- if (k == nb.block) totalcol else nbcol + (k - 1)*nbcol
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1)*nbcol):endcolsel], sep = colsep, file = zz)
                    cat("\n", file = zz)
                }
            }
        }
    }, "fasta" = {
        for (i in 1:N) {
            cat(">", names(x)[i], file = zz, sep = "")
            cat("\n", file = zz)
            X <- paste(x[[i]], collapse = "")
            S <- length(x[[i]])
            totalcol <- ceiling(S/colw)
            if (nbcol < 0) nbcol <- totalcol
            nb.lines <- ceiling(totalcol/nbcol)
            SEQ <- character(totalcol)
            for (j in 1:totalcol)
            SEQ[j] <- substr(X, 1 + (j - 1) * colw, colw + (j - 1) * colw)
            for (k in 1:nb.lines) {
                endsel <-
                if (k == nb.lines) length(SEQ) else nbcol + (k - 1)*nbcol
                cat(indent, file = zz)
                cat(SEQ[(1 + (k - 1)*nbcol):endsel], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        }
    })
}





insertRow <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,setNames(newrow, names(existingDF)))
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)
}


baseConCat <- function(dataFileExtension, fileFormat){
    files = list.files(pattern=dataFileExtension)
    dataObject = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    header_col = c('FileName', 'Species', 'Sequence')
    names(dataObject) = header_col
    
    for (file in files) {
        data = read.alignment(file, fileFormat)
        x = 1
        while (x <= length(data$seq)){
            dataObject = insertRow(dataObject, data.frame("FileName" = file,
            "Species" = data$nam[x], "Sequence" = data$seq[x]), 1)
            x = x + 1
        }
    }
    
    newDataObject = new.env()
    
    for (filename in files){
        newDataObject[[filename]] = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    }
    
    for (filename in files){
        for (i in 1:length(dataObject$FileName)){
            if (dataObject[i,][1] == filename){
                newDataObject[[filename]] = insertRow(newDataObject[[filename]], dataObject[i,], 1)
            }
        }
    }
    
    row_to_find <- data.frame(b="cat")
    
    speciesAll = unique(dataObject$Species)
    
    
    for (filename in files){
        rFlag = FALSE
        tmp <- sapply(newDataObject[[filename]][1,][3], as.character)
        if (substr(tmp, getLength(tmp), getLength(tmp)) == '\r') {
            rFlag = TRUE
        }
        
        for (species in speciesAll){
            row_to_find = data.frame("Species"=species)
            if (isTRUE(nrow(merge(row_to_find,newDataObject[[filename]]))>0) == FALSE){
                
                if (rFlag == TRUE) {
                    missingObject = paste(paste(rep("?", getLength(tmp)-1), collapse=""), "\r", sep="")
                }
                else {
                    missingObject = paste(rep("?", getLength(tmp)), collapse="")
                }
                
                newDataObject[[filename]] = insertRow(newDataObject[[filename]],
                data.frame("FileName" = filename, "Species" = species,
                "Sequence" = missingObject), 1)
            }
        }
    }
    
    concatFrame = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    
    for (sname in speciesAll){
        seqObject = c()
        for (filename in files){
            for (i in 1:length(newDataObject[[filename]]$Sequence)){
                if (newDataObject[[filename]][i,][2] == sname){
                    seqObject = c(seqObject, newDataObject[[filename]][i,][3])
                }
            }
        }
        
        concatFrame = insertRow(concatFrame, data.frame("FileName" = "MasterConcat",
        "Species" = sname, "Sequence" = do.call(paste, c(as.list(seqObject), sep=""))), 1)
    }
    
    concatFrame = concatFrame[-nrow(concatFrame),]
    concatFrame = concatFrame[-1,]
    
    return(concatFrame)
}


nexConCat <- function(dataFileExtension, fileFormat) {
    nexFiles = list.files(pattern=dataFileExtension)
    speciesAll = c()
    dataObject = data.frame("FileName" = 'Test', "Species" = 'Test', "Sequence" = 'Test')
    for (filename in nexFiles) {
        nexData = tryCatch({read.nexus.data(filename)}, finally = {print(filename)})
        speciesAll = c(speciesAll, names(nexData))
        for (i in 1:length(nexData)) {
            dataObject = insertRow(dataObject, data.frame("FileName" = filename, "Species" = names(nexData[i]), "Sequence" = paste(nexData[[names(nexData[i])]], collapse="")), 1)
        }
        
    }
    
    speciesAll = unique(speciesAll)
    newDataObject = new.env()
    
    for (filename in nexFiles){
        newDataObject[[filename]] = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    }
    
    for (filename in nexFiles){
        for (i in 1:length(dataObject$FileName)){
            if (dataObject[i,][1] == filename){
                newDataObject[[filename]] = insertRow(newDataObject[[filename]], dataObject[i,], 1)
            }
        }
    }
    
    
    for (filename in nexFiles){
        rFlag = FALSE
        tmp <- sapply(newDataObject[[filename]][1,][3], as.character)
        
        if (substr(tmp, getLength(tmp), getLength(tmp)) == '\r') {
            rFlag = TRUE
        }
        
        for (species in speciesAll){
            row_to_find = data.frame("Species"=species)
            if (isTRUE(nrow(merge(row_to_find,newDataObject[[filename]]))>0) == FALSE){
                if (rFlag == TRUE) {
                    missingObject = paste(paste(rep("?", getLength(tmp)-1), collapse=""), "\r", sep="")
                }
                else {
                    missingObject = paste(rep("?", getLength(tmp)), collapse="")
                }
                
                newDataObject[[filename]] = insertRow(newDataObject[[filename]], data.frame("FileName" = filename,
                "Species" = species, "Sequence" = missingObject), 1)
            }
        }
    }
    
    concatFrame = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    
    for (sname in speciesAll){
        seqObject = c()
        for (filename in nexFiles){
            for (i in 1:length(newDataObject[[filename]]$Sequence)){
                if (newDataObject[[filename]][i,][2] == sname){
                    seqObject = c(seqObject, newDataObject[[filename]][i,][3])
                }
            }
        }
        
        concatFrame = insertRow(concatFrame, data.frame("FileName" = "MasterConcat", "Species" = sname, "Sequence" = do.call(paste, c(as.list(seqObject), sep=""))), 1)
    }
    
    concatFrame = concatFrame[-nrow(concatFrame),]
    concatFrame = concatFrame[-1,]
    
    return(concatFrame)
}


concat <- function (ext, form, outform,writeData=TRUE) {
    if (form == 'nexus') {
        outData = nexConCat(ext, 'nexus')
    }
    else {
        outData = baseConCat(ext, form)
    }
    if (writeData == TRUE) {
        write.fasta(as.list(outData$Sequence), outData$Species, nbchar = 60, "Output.fas", open = 'w')
    }
    fasta <- read.alignment("Output.fas", format = "fasta")
    seq_list=fasta$seq
    seq_name=fasta$nam
    seq=list()
    seq1=list()
    for (i in 1:length(seq_list)) {
        seq[[i]]<-sapply(c(strsplit(sapply(seq_list[[i]],as.character),"")),as.character)
        names(seq)[i]<-as.character(fasta$nam[i])
        seq1[[i]]<-sapply(c(strsplit(sapply(seq_list[[i]],as.character),"")),as.character)
        names(seq1)[i]<-as.character(paste(unlist(sapply(strsplit(sapply(seq_name[[i]],as.character),""),as.list)[1:10]),collapse=""))
    }
    

    if (outform == 'nexus') {
        outData= seq
        write.nex(outData, file= "Output.nex", interleaved = TRUE, charsperline = 100)
        unlink("Output.fas", recursive = FALSE)
    }
    if (outform=='phylip_relaxed'){
        outData=seq
        write.dna(outData,file="Output_phy_relaxed.phy" , format = "interleaved")
        unlink("Output.fas", recursive = FALSE)
    }
    if (outform=='phylip_interleaved'){
        outData=seq1
        write.dna(outData,file="Output_phy_interleaved.phy" , format = "interleaved")
        unlink("Output.fas", recursive = FALSE)
    }
    if (outform=='phylip_sequential'){
        outData=seq1
        write.dna(outData,file="Output_phy_sequential.phy" , format = "sequential")
        unlink("Output.fas", recursive = FALSE)
    }
    return(outData)
}




#x=concat('.nex', 'nexus','phylip_interleaved', writeData=TRUE)




