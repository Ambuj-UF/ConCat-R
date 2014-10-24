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


abs <- function (number) {
    if (number >= 0) {
        pass
    }
    else {
        number = - number
    }
    return(number)
}


rcvdna <- functions (filename) {
    seqData = read.nexus(filename)
    numA = 0; numC = 0; numG = 0; numT = 0; numGap = 0
    numAL = list(); numCL = list(); numGL = list(); numTL = list(); numGapL = list();
    for (seq in seqData){
        sequence = seqData[[names(seq)]]
        numAL[[names(seq)]] = 0; numCL[[names(seq)]] = 0; numGL[[names(seq)]] = 0; numTL[[names(seq)]] = 0; numGapL[[names(seq)]] = 0;
        for (charData in sequence) {
            if (charData == 'a') {
                numA = numA + 1
                numAL[[names(seq)]] = numAL[[names(seq)]] + 1
            }
            else if (charData == 'c') {
                numC = numC + 1
                numCL[[names(seq)]] = numCL[[names(seq)]] + 1
            }
            else if (charData == 'g') {
                numG = numG + 1
                numGL[[names(seq)]] = numGL[[names(seq)]] + 1
            }
            else if (charData == 't') {
                numT = numT + 1
                numTL[[names(seq)]] = numTL[[names(seq)]] + 1
            }
            else if (charData == '-') {
                numGap = numGap + 1
                numGapL[[names(seq)]] = numGapL[[names(seq)]] + 1
            }
        }
    }
    
    rcvCal = 0; nTaxa = length(names(seqData))
    
    for (species in names(seqData) {
        rcvCal = rcvCal + abs(numAL[[species]] - (numA/nTaxa)) + abs(numGL[[species]] - (numG/nTaxa)) + abs(numCL[[species]] - (numC/nTaxa)) + abs(numTL[[species]] - (numT/nTaxa)) + abs(numGapL[[species]] - (numGap/nTaxa))
    }
    
    totalRCV = rcvCal/(length(nTaxa)*getLength(seqData[[names(seqData[1])]]))
    
    return(totalRCV)
    
}


rcvprot <- functions (filename) {
    seqData = read.nexus(filename)
    numA = 0; numB = 0; numI = 0; numL = 0; numF = 0; numN = 0;
    numS = 0; numC = 0; numH = 0; numU = 0; numG = 0; numP = 0; numGap = 0
    
    numAL = list(); numBL = list(); numIL = list(); numLL = list();
    numFL = list(); numNL = list(); numSL = list(); numCL = list();
    numHL = list(); numUL = list(); numGL = list(); numPL = list(); numGapL = list()
    
    for (seq in seqData){
        sequence = seqData[[names(seq)]]

        numAL[[names(seq)]] = 0; numBL[[names(seq)]] = 0;numIL[[names(seq)]] = 0;numLL[[names(seq)]] = 0;
        numFL[[names(seq)]] = 0;numNL[[names(seq)]] = 0;numSL[[names(seq)]] = 0;numCL[[names(seq)]] = 0;
        numHL[[names(seq)]] = 0;numUL[[names(seq)]] = 0;numGL[[names(seq)]] = 0;numPL[[names(seq)]] = 0;numGapL[[names(seq)]] = 0;
        
        for (charData in sequence) {
            if (charData == 'd' or charData == 'e') {
                numA = numA + 1
                numAL[[names(seq)]] = numAL[[names(seq)]] + 1
            }
            else if (charData == 'r' or charData == 'k') {
                numB = numB + 1
                numBL[[names(seq)]] = numBL[[names(seq)]] + 1
            }
            else if (charData == 'i' or charData == 'v') {
                numI = numI + 1
                numIL[[names(seq)]] = numIL[[names(seq)]] + 1
            }
            else if (charData == 'l' or charData == 'm') {
                numL = numL + 1
                numLL[[names(seq)]] = numLL[[names(seq)]] + 1
            }
            else if (charData == 'f' or charData == 'w') {
                numF = numF + 1
                numFL[[names(seq)]] = numFL[[names(seq)]] + 1
            }
            else if (charData == 'n' or charData == 'q') {
                numN = numN + 1
                numNL[[names(seq)]] = numNL[[names(seq)]] + 1
            }
            else if (charData == 's' or charData == 't') {
                numS = numS + 1
                numSL[[names(seq)]] = numSL[[names(seq)]] + 1
            }
            else if (charData == 'c') {
                numC = numC + 1
                numCL[[names(seq)]] = numCL[[names(seq)]] + 1
            }
            else if (charData == 'h') {
                numH = numH + 1
                numHL[[names(seq)]] = numHL[[names(seq)]] + 1
            }
            else if (charData == 'u') {
                numU = numU + 1
                numUL[[names(seq)]] = numUL[[names(seq)]] + 1
            }
            else if (charData == 'g') {
                numG = numG + 1
                numGL[[names(seq)]] = numGL[[names(seq)]] + 1
            }
            else if (charData == 'p') {
                numP = numP + 1
                numPL[[names(seq)]] = numPL[[names(seq)]] + 1
            }
            else if (charData == '-') {
                numGap = numGap + 1
                numGapL[[names(seq)]] = numGapL[[names(seq)]] + 1
            }
        }
    }
    
    rcvCal = 0; nTaxa = length(names(seqData))
    
    for (species in names(seqData) {
        rcvCal = rcvCal + abs(numAL[[species]] - (numA/nTaxa)) + abs(numBL[[species]] - (numB/nTaxa)) + abs(numIL[[species]] - (numI/nTaxa)) + abs(numLL[[species]] - (numL/nTaxa))+ abs(numFL[[species]] - (numF/nTaxa)) + abs(numNL[[species]] - (numN/nTaxa))+ abs(numSL[[species]] - (numS/nTaxa))+ abs(numCL[[species]] - (numC/nTaxa))+ abs(numHL[[species]] - (numH/nTaxa))+ abs(numUL[[species]] - (numU/nTaxa))+ abs(numGL[[species]] - (numG/nTaxa))+ abs(numPL[[species]] - (numP/nTaxa)) + abs(numGapL[[species]] - (numGap/nTaxa))
    }
    
    totalRCV = rcvCal/(length(nTaxa)*getLength(seqData[[names(seqData[1])]]))
    
    return(totalRCV)
    
}






