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

searchcds <- function(gene, group=NULL) {
    if (group != NULL) {
        term = paste(paste(gene, "[sym]", sep=""), paste(group, "[orgn]", sep=""), sep=" ")
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
    for (lines %in% x) {
        if ('→' %in% lines and !(':' %in% lines)) {
            idList = c(idList, str_replace_all(lines.strsplit(lines, '→')[1], " ", ""))
        }
    }
    
    return(idList)
}


