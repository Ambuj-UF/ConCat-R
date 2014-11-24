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

# Calculates fast evolving sites

require(seqinr)

<<<<<<< HEAD

#' Class \code{"fast-evol"}
#' @export
#' usage = fevol(file, cutoff)
fevol <- function(file, cutoff=0.9) {
    seqData = read.alignment(file, "fasta")
=======
fevol <- function(filename, cutoff=0.9) {
    seqData = read.alignment(filename, "fasta")
>>>>>>> FETCH_HEAD
    seqList = list()
    for (pos in 1:getLength(seqData$seq[1])) {
        
        seqObj = c()
        for (i in 1:length(seqData$nam)) {
            seqObj = c(seqObj, substring(seqData[i]$seq, pos, pos))
        }
        
        seqList[[sapply(pos, as.character)]] = seqObj[seqObj != "" & seqObj != "?"]
    }
    
    
    OVlist = list()
    retOV = c()
    for (i in 1:length(seqList)) {
        
        counter = 1
        incounter = 0
        inStore = c()
        for (obj in seqList[[names(seqList[i])]]) {
            inStore = c(inStore, obj)
            for (inObj in seqList[[names(seqList[i])]][!(seqList[[names(seqList[i])]] %in% inStore)]) {
                if (obj != inObj) {incounter = incounter + 1}
            }
        }
        
        posVecLen = length(seqList[[names(seqList[i])]])
        
        k = ((posVecLen)^2 - posVecLen)/2
        OVlist[[names(seqList[i])]] = incounter/k
        if (incounter/k >= as.numeric(cutoff)) {
            retOV = c(retOV, i)
        }
    }
    
    return(retOV)
    
}


x = fevol("MCPH1.fas", 0.7)







