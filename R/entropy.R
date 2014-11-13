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

# Calculates alignment entropy

require(seqinr)

shannon <- function(Obj) {
    totEntVect = c()
    entList = list()
    for (i in 1:getLength(Obj$seq[1])) {
        entList[[i]] = list()
        objects = c()
        for (j in 1:length(Obj$nam)) {
            objects = c(objects, substring(Obj[j]$seq, i, i))
        }
        
        objects = objects[!objects %in% c('?')]
        uniqObj = unique(objects)
        
        freqVector = c()
        for (uniqData in uniqObj) {
            counter = 0
            for (inData in objects) {
                if (inData == uniqData) {
                    counter = counter + 1
                }
            }
            freqVector = c(freqVector, counter/length(objects))
        }
        
        colEntropyVector = c()
        for (values in freqVector) {
            colEntropyVector = c(colEntropyVector, values*log(values))
        }
        
        colEntropy = -sum(colEntropyVector)
        totEntVect = c(totEntVect, colEntropy)
    }
    
    entropy = sum(totEntVect)/getLength(Obj$seq[1])
    return(entropy)
}



shannonProt <- function (Obj) {
    totEntVect = c()
    entList = list()
    for (i in 1:getLength(Obj$seq[1])) {
        entList[[i]] = list()
        objects = c()
        for (j in 1:length(Obj$nam)) {
            objects = c(objects, substring(Obj[j]$seq, i, i))
        }
        
        objects = objects[!objects %in% c('?')]
        
        aminoList = list()
        aminoList$a = c('d', 'e')
        aminoList$b = c('r', 'k')
        aminoList$i = c('i', 'v')
        aminoList$l = c('l', 'm')
        aminoList$f = c('f', 'w')
        aminoList$n = c('n', 'q')
        aminoList$s = c('s', 't')
        aminoList$c = c('c')
        aminoList$h = c('h')
        aminoList$u = c('u')
        aminoList$g = c('g')
        aminoList$p = c('p')
        
        for (i in length(objects)) {
             for (j in 1:length(aminoList)) {
                 if (objects[[i]] %in% aminoList[j][[names(aminoList[j])]]) {
                     objects[[i]] = names(aminoList[j])
                 }
            }
        }
        
        uniqObj = unique(objects)
        
        freqVector = c()
        for (uniqData in uniqObj) {
            counter = 0
            for (inData in objects) {
                if (inData == uniqData) {
                    counter = counter + 1
                }
            }
            freqVector = c(freqVector, counter/length(objects))
        }
        
        colEntropyVector = c()
        for (values in freqVector) {
            colEntropyVector = c(colEntropyVector, values*log(values))
        }
        
        colEntropy = -sum(colEntropyVector)
        totEntVect = c(totEntVect, colEntropy)
    }
    
    entropy = sum(totEntVect)/getLength(Obj$seq[1])
    return(entropy)
}


entcal <- function(file, type="DNA") {
    if (type == "DNA" | type == "dna") {
        alignment = read.alignment("MCPH1.fas", "fasta")
        ent = shannon(alignment)
    }
    else {
        alignment = read.alignment("MCPH1.fas", "fasta")
        ent = shannonProt(alignment)
    }
    
    return(ent)
}


