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

#' require(RCurl)
#' require(XML)


trim <- function (x) gsub("^\\s+|\\s+$", "", x)

searchcds <- function(gene, group=NULL) {
    if (group != NULL) {
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
    }
    
    return(idList)
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
    retSeq = seqList[1]
    for (seq in seqList) {
        if (getLength(paste(strsplit((seq, '\n')[[1]][2:length(strsplit(seq, '\n')[[1]])], collapse="")) >
            getLength(paste(strsplit((retSeq, '\n')[[1]][2:getLength(strsplit(retSeq, '\n')[[1]])]), collapse=""))) {
                retSeq = seq
        }
    }
    
    return(retSeq)
}


for (geneName in geneList) {
    ids = searchcds(geneName, group = orgn)
    totalSeqList = c()
    for (id in ids) {
        seqList=c()
        url = paste("http://www.ncbi.nlm.nih.gov/gene/", id, sep="")
        idList = unique(fetchIDs(url))
        for (inID in idList) {
            seqList = c(seqList, fetchSeq(inID, type))
        }
        
        totalSeqList = c(totalSeqList, longestSeq(seqList))
    }
    
    sink(paste(geneName, ".fas"))
    
    for (obj in totalSeqList) {cat(obj)}
}


