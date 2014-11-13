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

# RY coding
# Usage rycode(alignment-Object, type=3)


rycode <- function (seqData, type=3) {
    if (type == 3) {
        for (i in 1:length(seqData)) {
            counter = 1
            for (char in unlist(seqData[[names(seqData[i])]])) {
                if (counter%%3 == 0) {
                    if (seqData[[names(seqData[i])]][[1]][counter] == 'a' | seqData[[names(seqData[i])]][[1]][counter] == 'g') {
                        seqData[[names(seqData[i])]][[1]][counter] = 'r'
                    }
                    else if(seqData[[names(seqData[i])]][[1]][counter] == 't' | seqData[[names(seqData[i])]][[1]][counter] == 'c') {
                        seqData[[names(seqData[i])]][[1]][counter] = 'y'
                    }
                    
                }
                
                counter = counter + 1
            }
        }
    }
    
    else if (type == 'all') {
        for (i in 1:length(seqData)) {
            for (char in unlist(seqData[[names(seqData[i])]])) {
                if (seqData[[names(seqData[i])]][[1]][counter] == 'a' | seqData[[names(seqData[i])]][[1]][counter] == 'g') {
                    seqData[[names(seqData[i])]][[1]][counter] = 'r'
                }
                else if(seqData[[names(seqData[i])]][[1]][counter] == 't' | seqData[[names(seqData[i])]][[1]][counter] == 'c') {
                    seqData[[names(seqData[i])]][[1]][counter] = 'y'
                }
            }
        }
    }
    
    sink("RYcoded.fas")
    
    for (i in 1:length(seqData)) {
        cat(paste(">",sapply(names(seqData[i]), as.character), sep=""),"\n")
        cat(paste(unlist(seqData[[names(seqData[i])]])), sep="")
        cat("\n\n")
    }
    
    sink()
}














