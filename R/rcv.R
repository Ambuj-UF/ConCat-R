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
        numAL[[names(seq)]] = 0; numCL[[names(seq)]] = 0; numGL[[names(seq)]] = 0; numGL[[names(seq)]] = 0; numAL[[names(seq)]] = 0;
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
    
    valData$A = numAL; valData$C = numCL; valData$G = numGL; valData$T = numTL; valData$Gap = numGapL
    
    for (species in names(seqData) {
        rcvCal = rcvCal + abs(numAL[[species]] - (numA/nTaxa)) + abs(numGL[[species]] - (numG/nTaxa)) + abs(numCL[[species]] - (numC/nTaxa)) + abs(numTL[[species]] - (numT/nTaxa)) + abs(numGapL[[species]] - (numGap/nTaxa))
    }
    
    totalRCV = rcvCal/(length(nTaxa)*getLength(seqData[[names(seqData[1])]]))
    
    return(totalRCV)
    
}





