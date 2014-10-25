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

# taxa authority functions


excTaxaNex <- functions (alignmentObject, authorityFile) {
    taxExc = read.table(authorityFile,header=FALSE,sep="\n")
    for (taxa in taxExc) {
        alignmentObject[[taxa]] = NULL
    }
    
    return(alignmentObject)
}



incTaxaNex <- functions (alignmentObject, authorityFile) {
    taxInc = read.table(authorityFile,header=FALSE,sep="\n")
    newTaxaList = list()
    for (taxa in taxInc) {
        newTaxaList[[taxa]] = alignmentObject[[taxa]]
    }
    
    return(newTaxaList)
}


excTaxaOd <- functions (alignmentObject, authorityFile) {
    taxExc = read.table(authorityFile,header=FALSE,sep="\n")
    nameString = sapply(alignmentObject$Species[1], as.character)
    newTaxaExc = c()
    if (substr(nameString, getLength(nameString), getLength(nameString)) == '\r') {
        for (i in 1:length(taxInc$V!)) {
            newTaxaInc = c(newTaxaInc, taxInc[i,])
            newTaxaInc = c(newTaxaInc, paste(taxInc[i,], '\r', sep=''))
        }
    }
    dataObject = alignmentObject[!(alignmentObject$Species) %in% taxExc,]
    return(alignmentObject)
}


incTaxaOd <- functions (alignmentObject, authorityFile) {
    taxInc = read.table(authorityFile,header=FALSE,sep="\n")
    nameString = sapply(alignmentObject$Species[1], as.character)
    newTaxaInc = c()
    if (substr(nameString, getLength(nameString), getLength(nameString)) == '\r') {
        for (i in 1:length(taxExc$V!)) {
            newTaxaExc = c(newTaxaExc, taxExc[i,])
            newTaxaExc = c(newTaxaExc, paste(taxExc[i,], '\r', sep=''))
        }
    }
    dataObject = alignmentObject[(alignmentObject$Species) %in% newTaxaExc,]
    return(alignmentObject)
}






