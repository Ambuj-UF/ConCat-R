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

abs <- function (number) {
    if (number >= 0) {
        number = number
    }
    else {
        number = - number
    }
    return(number)
}


fastaToNexus <- function(aliObj) {
    retAli = list()
    for (i in 1:length(aliObj$nam)) {
        retAli[[aliObj$nam[i]]] = strsplit(sapply(aliObj$seq[i][[1]], as.character), "")
    }
    
    return(retAli)
}


rcvdna <- function(filename){
    seqData = read.alignment(filename, 'fasta')
    seqData = fastaToNexus(seqData)
    numA = 0; numC = 0; numG = 0; numT = 0; numGap = 0
    numAL = list(); numCL = list(); numGL = list(); numTL = list(); numGapL = list()
    for (i in 1:length(seqData)){
        sequence = seqData[[names(seqData[i])]]
        numAL[[names(seqData[i])]] = 0; numCL[[names(seqData[i])]] = 0; numGL[[names(seqData[i])]] = 0; numTL[[names(seqData[i])]] = 0; numGapL[[names(seqData[i])]] = 0;
        for (charData in sequence[[1]]) {
            if (charData == 'a') {
                numA = numA + 1
                numAL[[names(seqData[i])]] = numAL[[names(seqData[i])]] + 1
            }
            else if (charData == 'c') {
                numC = numC + 1
                numCL[[names(seqData[i])]] = numCL[[names(seqData[i])]] + 1
            }
            else if (charData == 'g') {
                numG = numG + 1
                numGL[[names(seqData[i])]] = numGL[[names(seqData[i])]] + 1
            }
            else if (charData == 't') {
                numT = numT + 1
                numTL[[names(seqData[i])]] = numTL[[names(seqData[i])]] + 1
            }
            else if (charData == '-') {
                numGap = numGap + 1
                numGapL[[names(seqData[i])]] = numGapL[[names(seqData[i])]] + 1
            }
        }
    }
    
    rcvCal = 0; nTaxa = length(names(seqData))
    
    for (i in 1:length(seqData)) {
        rcvCal = rcvCal + abs(numAL[[names(seqData[i])]] - (numA/nTaxa)) + abs(numGL[[names(seqData[i])]] - (numG/nTaxa)) + abs(numCL[[names(seqData[i])]] - (numC/nTaxa)) + abs(numTL[[names(seqData[i])]] - (numT/nTaxa)) + abs(numGapL[[names(seqData[i])]] - (numGap/nTaxa))
    }

    totalRCV = rcvCal/(length(seqData)*length(unlist(seqData[[names(seqData[1])]])))
    return(totalRCV)
}


rcvprot <- function(filename) {
    seqData = read.alignment(filename)
    seqData = fastaToNexus(seqData)
    numA = 0; numB = 0; numI = 0; numL = 0; numF = 0; numN = 0;
    numS = 0; numC = 0; numH = 0; numU = 0; numG = 0; numP = 0; numGap = 0
    
    numAL = list(); numBL = list(); numIL = list(); numLL = list();
    numFL = list(); numNL = list(); numSL = list(); numCL = list();
    numHL = list(); numUL = list(); numGL = list(); numPL = list(); numGapL = list()
    
    for (i in 1:length(seqData)){
        sequence = seqData[[names(seqData[i])]]

        numAL[[names(seqData[i])]] = 0; numBL[[names(seqData[i])]] = 0;numIL[[names(seqData[i])]] = 0;numLL[[names(seqData[i])]] = 0;
        numFL[[names(seqData[i])]] = 0; numNL[[names(seqData[i])]] = 0;numSL[[names(seqData[i])]] = 0;numCL[[names(seqData[i])]] = 0;
        numHL[[names(seqData[i])]] = 0; numUL[[names(seqData[i])]] = 0;numGL[[names(seqData[i])]] = 0;numPL[[names(seqData[i])]] = 0;numGapL[[names(seqData[i])]] = 0;
        
        for (charData in sequence) {
            if (charData == 'd' | charData == 'e') {
                numA = numA + 1
                numAL[[names(seqData[i])]] = numAL[[names(seqData[i])]] + 1
            }
            else if (charData == 'r' | charData == 'k') {
                numB = numB + 1
                numBL[[names(seqData[i])]] = numBL[[names(seqData[i])]] + 1
            }
            else if (charData == 'i' | charData == 'v') {
                numI = numI + 1
                numIL[[names(seqData[i])]] = numIL[[names(seqData[i])]] + 1

            }
            else if (charData == 'l' | charData == 'm') {
                numL = numL + 1
                numLL[[names(seqData[i])]] = numLL[[names(seqData[i])]] + 1
            }
            else if (charData == 'f' | charData == 'w') {
                numF = numF + 1
                numFL[[names(seqData[i])]] = numFL[[names(seqData[i])]] + 1
            }
            else if (charData == 'n' | charData == 'q') {
                numN = numN + 1
                numNL[[names(seqData[i])]] = numNL[[names(seqData[i])]] + 1
            }
            else if (charData == 's' | charData == 't') {
                numS = numS + 1
                numSL[[names(seqData[i])]] = numSL[[names(seqData[i])]] + 1
            }
            else if (charData == 'c') {
                numC = numC + 1
                numCL[[names(seqData[i])]] = numCL[[names(seqData[i])]] + 1
            }
            else if (charData == 'h') {
                numH = numH + 1
                numHL[[names(seqData[i])]] = numHL[[names(seqData[i])]] + 1
            }
            else if (charData == 'u') {
                numU = numU + 1
                numUL[[names(seqData[i])]] = numUL[[names(seqData[i])]] + 1
            }
            else if (charData == 'g') {
                numG = numG + 1
                numGL[[names(seqData[i])]] = numGL[[names(seqData[i])]] + 1
            }
            else if (charData == 'p') {
                numP = numP + 1
                numPL[[names(seqData[i])]] = numPL[[names(seqData[i])]] + 1
            }
            else if (charData == '-') {
                numGap = numGap + 1
                numGapL[[names(seqData[i])]] = numGapL[[names(seqData[i])]] + 1
            }
        }
    }
    
    rcvCal = 0; nTaxa = length(names(seqData))
    
    for (i in 1:length(seqData)) {
        rcvCal = rcvCal + abs(numAL[[names(seqData[i])]] - (numA/nTaxa)) + abs(numBL[[names(seqData[i])]] - (numB/nTaxa)) + abs(numIL[[names(seqData[i])]] - (numI/nTaxa)) + abs(numLL[[names(seqData[i])]] - (numL/nTaxa))+ abs(numFL[[names(seqData[i])]] - (numF/nTaxa)) + abs(numNL[[names(seqData[i])]] - (numN/nTaxa))+ abs(numSL[[names(seqData[i])]] - (numS/nTaxa))+ abs(numCL[[names(seqData[i])]] - (numC/nTaxa))+ abs(numHL[[names(seqData[i])]] - (numH/nTaxa))+ abs(numUL[[names(seqData[i])]] - (numU/nTaxa))+ abs(numGL[[names(seqData[i])]] - (numG/nTaxa))+ abs(numPL[[names(seqData[i])]] - (numP/nTaxa)) + abs(numGapL[[names(seqData[i])]] - (numGap/nTaxa))
    }
    
    totalRCV = rcvCal/(length(seqData)*length(unlist(seqData[[names(seqData[1])]])))
    
    return(totalRCV)
    
}



rcv <- function(file, type="DNA") {
    if (type == "DNA" | type == "dna") {
        rcvValue = rcvdna(file)
    }
    
    else {rcvValue = rcvprot(file)}
}
    









