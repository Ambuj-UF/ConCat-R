################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
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

"""Under development phase"""


library(seqinr)

files = list.files(pattern = '*.fas')

insertRow <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,setNames(newrow, names(existingDF)))
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)
}

dataObject = data.frame('Test','Test','Test')
header_col = c('FileName', 'Species', 'Sequence')
names(dataObject) = header_col

for (file in files) {
    data = read.alignment(file, "fasta")
    x = 1
    while (x <= length(data$seq)){
        dataObject = insertRow(dataObject, data.frame("FileName" = file, "Species" = data$nam[x], "Sequence" = data$seq[x]), 1)
        x = x + 1
    }
}

newDataObject = new.env()

for (filename in files){
    newDataObject[[filename]] = data.frame('Test','Test','Test')
}

for (filename in files){
    for (data in dataObject){
        if (data$FileName == filename){
            newDataObject[[filename]] = insertRow(newDataObject[[filename]], data, 1)
    }
}
    
row_to_find <- data.frame(b="cat")

speciesAll = unique(dataObject$Species)

for (filename in files){
    for (species in speciesAll){
        row_to_find = data.frame("Species"=species)
        if (isTRUE(nrow(merge(row_to_find,newDataObject[filename]))>0) == FALSE){
            newDataObject[[filename]] = insertRow(newDataObject[[filename]], data.frame("FileName" = filename, "Species" = species, "Sequence" = paste(rep(strsplit("?","")[[1]], each=length(newDataObject[[filename]]$Sequence[1])), collapse="")), 1)
        }
    }
}

seqObject = c()
concatFrame = data.frame('Test','Test','Test')
    
for (sname in speciesAll){
    for (filename in files){
        for (data in newDataObject[[filename]]){
            if (data[2] == sname){
                seqObject = c(seqObject, data[3])
            }
        }
    }
    concatFrame = insertRow(concatFrame, data.frame("Filename" = "MasterConcat", "Species" = sname, "Sequence" = do.call(paste, c(as.list(sdata), sep=""))), 1)
}


fileConn<-file("output.fasta")
for (data in concatFrame){
    writeLines(c(">",data$Species), fileConn)
    writeLines(c(data$Sequence), fileConn)
    writeLines(c("\n\n"), fileConn)
}

close(fileConn)









