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


fevol <- function (filename, cutoff) {
    seqData = read.nexus(filename)
    seqList = list()
    for (pos in 1:getLength(seqData[[names(seqData)[1]]])) {
        seqObj = c()
        for (i in 1:length(seqData)) {
            seqObj = c(seqObj, seqData[[names(seqData)[i]]][pos])
        }
        seqList[[pos]] = seqObj
    }
    
    OVlist = list()
    retOV = c()
    
    for (x in seqList) {
        counter = 1
        incounter = 0
        for (obj in x[[names(x)]]) {
            for (inObj in x[[names(x)]]) {
                if (obj != inObj) {incounter = incounter + 1}
            }
        }
        k = (length(x[[names(x)]])^2 - length(x[[names(x)]]))/2
        OVlist[[names(x)]] = incounter/k
        if (incounter/k >= cutoff) {
            retOV = c(retOV, names(x))
        }
    }
    
    return(retOV)
    
}







