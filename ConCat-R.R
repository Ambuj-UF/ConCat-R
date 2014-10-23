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

#Under development phase


library(seqinr)

#bwt <- function(s) {
    #Apply Burrows-Wheeler transform to input string
    
#    stopifnot(grepl('\0', s) == TRUE)
#    s = paste(s, "\0")
#    table = c()
#    for (i in c(1:length(s))){
#        table = c(table, paste(substr(x, i, length(s)), substr(x, 1, length(s))))
#    }
#
#    table = sort(table)
#    last_column = c()
#    for (row in table){
#        last_column = c(last_column, paste(substr(row, length(row), length(row)), substr(row, 1, length(row))))
#    }
    
#    return(do.call(paste, c(as.list(last_column), sep="")))
#}

#ibwt <- function(r, *args) {
    #Inverse Burrows-Wheeler transform. args is the original index if it was not indicated by a null byte
    #Working on it
    
#}

insertRow <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,setNames(newrow, names(existingDF)))
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)
}


ConCat <- function(dataFileExtension){
    files = list.files(pattern=dataFileExtension)
    dataObject = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    header_col = c('FileName', 'Species', 'Sequence')
    names(dataObject) = header_col

    for (file in files) {
        data = read.alignment(file, "fasta")
        x = 1
        while (x <= length(data$seq)){
            dataObject = insertRow(dataObject, data.frame("FileName" = file,
                "Species" = data$nam[x], "Sequence" = data$seq[x]), 1)
            x = x + 1
        }
    }

    newDataObject = new.env()

    for (filename in files){
        newDataObject[[filename]] = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    }
    
    for (filename in files){
        for (i in 1:length(dataObject$FileName)){
            if (dataObject[i,][1] == filename){
                newDataObject[[filename]] = insertRow(newDataObject[[filename]], dataObject[i,], 1)
            }
        }
    }
    
    row_to_find <- data.frame(b="cat")

    speciesAll = unique(dataObject$Species)
    

    for (filename in files){
        for (species in speciesAll){
            row_to_find = data.frame("Species"=species)
            if (isTRUE(nrow(merge(row_to_find,newDataObject[[filename]]))>0) == FALSE){
                tmp <- sapply(newDataObject[[filename]][1,][3], as.character)
                newDataObject[[filename]] = insertRow(newDataObject[[filename]],
                    data.frame("FileName" = filename, "Species" = species,
                        "Sequence" = paste(rep(strsplit("?","")[[1]],
                            each=getLength(tmp)), collapse="")), 1)
            }
        }
    }

    concatFrame = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    
    for (sname in speciesAll){
        seqObject = c()
        for (filename in files){
            for (i in 1:length(newDataObject[[filename]]$Sequence)){
                if (newDataObject[[filename]][i,][2] == sname){
                    seqObject = c(seqObject, newDataObject[[filename]][i,][3])
                }
            }
        }

        concatFrame = insertRow(concatFrame, data.frame("FileName" = "MasterConcat",
            "Species" = sname, "Sequence" = do.call(paste, c(as.list(seqObject), sep=""))), 1)
    }

    concatFrame = concatFrame[-nrow(concatFrame),]
    concatFrame = concatFrame[-1,]

    return(concatFrame)
}



z = ConCat('*.fas')

sink("outfile.txt")


for (i in 1:length(z$Sequence)){
  cat(paste(paste(">",sapply(z[i,][2], as.character)),sapply(z[i,][3], as.character),sep='\n'))
   cat("\n\n")
}

sink()
