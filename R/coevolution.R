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

# Calculates coevlving amino acid sites


Poisson_dist <- function (seq1, seq2) {
    dist = 0
    gap = 0
    for (i in 1:length(seq1)){
        if (seq1[[i]] == seq2[[i]] && seq1[[i]] != '-' && seq2[[i]] != '-') {
            dist = dist + 1
        }
        if (seq1[[i]] == '-' || seq2[[i]] == '-') {
            gap = gap + 1
        }
    }
    return(Poisson(dist, length(seq1) - gap))
}

Poisson <- function (distance, long) {
    pdist = 0
    if ((1 - (distance/long)) != 0) {
        pdist = -log(1 - (distance/long))
    }
    else { pdist = 0}
    
    return(pdist)
}


optimize <- function (seqObj) {
    distance = c()
    seqNameVec = c()
    for (i in 1:length(names(seqObj))) {
        seqNameVec = c(seqNameVec, names(seqObj[i]))
        inSeqObj = seqObj[!names(seqObj[i]) %in% seqNameVec]
        seq1 = seqObj[i][names(seqObj[i])]
        for (j in 1:length(names(seqNameVec))) {
            seq2 = seqObj[j][names(seqObj[j])]
            distance = c(distance, Poisson_dist(seq1, seq2))
        }
    }
    
    sorted_dist = sort(distance)
}








