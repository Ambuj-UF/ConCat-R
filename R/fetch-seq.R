################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
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

require(RCurl)
require(XML)
require(reutils)
require(seqinr)



trim <- function (x) gsub("^\\s+|\\s+$", "", x)

searchcds <- function(gene, group) {
    if (group != "None") {
        term = paste(paste(gene, "[sym]", sep=""), paste(group, "[orgn]", sep=""), sep=" ")
        e <- esearch(term, "gene")
    }
    
    else {
        term = paste(gene, "[sym]", sep="")
        e <- esearch(term, "gene")
    }
    
    return(e)
}

fetchIDs <- function(url) {
    webpage <- getURL(url)
    webpage <- readLines(tc <- textConnection(webpage)); close(tc)
    pagetree <- htmlTreeParse(webpage, error=function(...){}, useInternalNodes = TRUE)
    x <- xpathSApply(pagetree, "//*", xmlValue)
    x <- unlist(strsplit(x, "\n"))
    x <- gsub("\t","",x)
    x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
    x <- x[!(x %in% c("", "|"))]

    idList = c()
    for (lines in x) {
        if(grepl('→', sapply(lines, as.character)) == TRUE & grepl(':', sapply(lines, as.character)) == FALSE) {
            idList = c(idList, trim(strsplit(lines, '→')[[1]][1]))
        }
        else if (grepl(' - Gene - NCBI', sapply(lines, as.character)) == TRUE) {
            spName = strsplit(strsplit(sapply(lines, as.character), split="\\[")[[1]][2], split="\\]")[[1]][1]
        }
    }
    
    
    return(list("id" = idList, "species" = spName))
}


fetchSeq <- function(ID, type) {
    if (type == 'cds' | type == 'CDS') {
        x <- efetch(ID, "nuccore", rettype = "fasta_cds_na", retmode = "text")
    }
    else if(type == 'amino acid') {
        x <- efetch(ID, "nuccore", rettype = "fasta_cds_aa", retmode = "text")
    }
    else {x <- efetch(ID, "nuccore", rettype = "fasta", retmode = "text")}

    return(content(x))
    
}


longestSeq <- function(seqList) {
    retSeq = seqList[[1]]
    for (seq in seqList) {
        seqLength = length(strsplit(seq, '\n')[[1]])
        seqObj = strsplit(seq, '\n')[[1]][2:seqLength]
        retSeqLength = length(strsplit(retSeq, '\n')[[1]])
        retSeqObj = strsplit(retSeq, '\n')[[1]][2:retSeqLength]
        if (getLength(paste(seqObj, collapse="")) > getLength(paste(retSeqObj, collapse=""))) {
                retSeq = seq
        }
    }

    return(retSeq)
}


fetch <- function(geneList, orgn="None") {
    for (geneName in geneList) {
        cat("Extracting sequence for ", geneName, "")
        ids = searchcds(geneName, orgn)
        
        totalSeqList = c()
        idData = c()
        
        for (x in 1:1000) {
            if (!is.na(ids[x])) {lenID = x}
                else {break}
        }
        
        for (x in 1:lenID) {
            
            seqList=c()
            url = paste("http://www.ncbi.nlm.nih.gov/gene/", ids[x], sep="")
            retObj = fetchIDs(url)
            idList = unique(retObj["id"]$id)
            cat(retObj["species"]$species)
            cat("\n")
            for (inID in idList) {
                seqList = c(seqList, fetchSeq(inID, 'cds'))
            }
            
            lseq = longestSeq(seqList)
            newSeq = strsplit(sapply(lseq, as.character), '\n')
            newSeq[[1]][1] = paste(strsplit(sapply(retObj["species"]$species, as.character), "\\(")[[1]][1], ">", sep="")
            newSeq = unlist(newSeq)
            newSeq = paste(newSeq, sep="\n")
            totalSeqList = c(totalSeqList, newSeq)
            
        }
        
        sink(paste(geneName, ".fas"))
        for (x in 1:length(totalSeqList)) {cat(sapply(totalSeqList[x], as.character))}
        sink()
    }
}

fetch(c('ASPM'), "Primates")

