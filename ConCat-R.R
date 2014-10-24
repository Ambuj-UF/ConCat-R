################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# "read.nex" and "write.nexus" is a Simplified NEXUS data parser by Johan Nylander, nylander @ scs.fsu.edu     #
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

#Under development phase


library(seqinr)

#bwt <- function(s) {
#Apply Burrows-Wheeler transform to input string

#    stopifnot(grepl('\0', s) == TRUE)
#    s = paste(s, "\0")
#    table = c()
#    for (i in c(1:length(s))){
#        table = c(table, paste(substr(x, i, length(s)), substr(x, 1, length(s))))
#    }
#
#    table = sort(table)
#    last_column = c()
#    for (row in table){
#        last_column = c(last_column, paste(substr(row, length(row), length(row)), substr(row, 1, length(row))))
#    }

#    return(do.call(paste, c(as.list(last_column), sep="")))
#}

#ibwt <- function(r, *args) {
#Inverse Burrows-Wheeler transform. args is the original index if it was not indicated by a null byte
#Working on it

#}


read.nex <- function (file)
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


write.nexus <- function (x, file, format = "dna", datablock = TRUE,
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
    
    fcat <- function (..., file = zz)
    {
        cat(..., file = file, sep = "", append = TRUE)
    }
    
    find.max.length <- function (x)
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
    
    print.matrix <- function(x, dindent = "    ")
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
            fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", GAP, " ", INTERLEAVE, ";\n") # from François Michonneau (2009-10-02)
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


insertRow <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,setNames(newrow, names(existingDF)))
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)
}


ConCat <- function(dataFileExtension, fileFormat){
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
        nexData = read.nex(filename)
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

concat <- function (extension, inpform) {
    if (inpform == 'nexus') {
        outData = nexConCat(extension, 'nexus')
        write.nexus.data(outData, file="FuncOutput.nex", format = "dna", datablock = TRUE, interleaved = TRUE)
    }
    else {
        outData = ConCat(extension, inpform)
        write.fasta(as.list(outData$Sequence), outData$Species, nbchar = 60, "FuncOutput.fas", open = 'w')
        return(outData)
    }
    
}



