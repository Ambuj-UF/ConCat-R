################################################################################################################
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                       {Mrinal Mishra, VTT Technical Research Center, Finland}                                #
# "read.nex" and write.nex are Simplified NEXUS data parser by Johan Nylander, nylander @ scs.fsu.edu          #
# "read.nex" is a Simplified NEXUS data parser by Johan Nylander, nylander @ scs.fsu.edu                       #
# WARNING: This is parser reads a restricted nexus format, see the link below for details                      #
# http://svitsrv25.epfl.ch/R-doc/library/ape/html/read.nexus.data.html &                                       #
# http://www.inside-r.org/packages/cran/ape/docs/write.nexus.data                                              #
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



library(seqinr)
library(ape)


insertRow <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,setNames(newrow, names(existingDF)))
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)
}

################################################################################################
# Non-nexus concatenation.
# Annotations handling - None

baseConCat <- function(dataFileExtension, fileFormat){
    files = list.files(pattern=dataFileExtension)
    dataObject = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    header_col = c('FileName', 'Species', 'Sequence')
    names(dataObject) = header_col
    
    for (file in files) {
        cat(file)
        cat("\n")
        data = read.alignment(file, fileFormat)
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
        rFlag = FALSE
        tmp <- sapply(newDataObject[[filename]][1,][3], as.character)
        if (substr(tmp, getLength(tmp), getLength(tmp)) == '\r') {
            rFlag = TRUE
        }
        
        for (species in speciesAll){
            row_to_find = data.frame("Species"=species)
            if (isTRUE(nrow(merge(row_to_find,newDataObject[[filename]]))>0) == FALSE){
                
                if (rFlag == TRUE) {
                    missingObject = paste(paste(rep("?", getLength(tmp)-1), collapse=""), "\r", sep="")
                }
                else {
                    missingObject = paste(rep("?", getLength(tmp)), collapse="")
                }
                
                newDataObject[[filename]] = insertRow(newDataObject[[filename]],
                data.frame("FileName" = filename, "Species" = species,
                "Sequence" = missingObject), 1)
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


################################################################################################
# Nexus data concatenation.
# Annotations handling - Yes

nexConCat <- function(dataFileExtension, fileFormat) {
    nexFiles = list.files(pattern=dataFileExtension)
    speciesAll = c()
    dataObject = data.frame("FileName" = 'Test', "Species" = 'Test', "Sequence" = 'Test')
    for (filename in nexFiles) {
        cat(filename)
        cat("\n")
        nexData = read.nexus.data(filename)
        speciesAll = c(speciesAll, names(nexData))
        for (i in 1:length(nexData)) {
            dataObject = insertRow(dataObject, data.frame("FileName" = filename, "Species" = names(nexData[i]), "Sequence" = paste(nexData[[names(nexData[i])]], collapse="")), 1)
        }
        
    }
    
    speciesAll = unique(speciesAll)
    newDataObject = new.env()
    
    for (filename in nexFiles){
        newDataObject[[filename]] = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    }
    
    for (filename in nexFiles){
        for (i in 1:length(dataObject$FileName)){
            if (dataObject[i,][1] == filename){
                newDataObject[[filename]] = insertRow(newDataObject[[filename]], dataObject[i,], 1)
            }
        }
    }
    
    
    for (filename in nexFiles){
        rFlag = FALSE
        tmp <- sapply(newDataObject[[filename]][1,][3], as.character)
        
        if (substr(tmp, getLength(tmp), getLength(tmp)) == '\r') {
            rFlag = TRUE
        }
        
        for (species in speciesAll){
            row_to_find = data.frame("Species"=species)
            if (isTRUE(nrow(merge(row_to_find,newDataObject[[filename]]))>0) == FALSE){
                if (rFlag == TRUE) {
                    missingObject = paste(paste(rep("?", getLength(tmp)-1), collapse=""), "\r", sep="")
                }
                else {
                    missingObject = paste(rep("?", getLength(tmp)), collapse="")
                }
                
                newDataObject[[filename]] = insertRow(newDataObject[[filename]], data.frame("FileName" = filename,
                "Species" = species, "Sequence" = missingObject), 1)
            }
        }
    }
    
    concatFrame = data.frame("FileName"='Test',"Species"='Test',"Sequence"='Test')
    
    for (sname in speciesAll){
        seqObject = c()
        for (filename in nexFiles){
            for (i in 1:length(newDataObject[[filename]]$Sequence)){
                if (newDataObject[[filename]][i,][2] == sname){
                    seqObject = c(seqObject, newDataObject[[filename]][i,][3])
                }
            }
        }
        
        concatFrame = insertRow(concatFrame, data.frame("FileName" = "MasterConcat", "Species" = sname, "Sequence" = do.call(paste, c(as.list(seqObject), sep=""))), 1)
    }
    
    concatFrame = concatFrame[-nrow(concatFrame),]
    concatFrame = concatFrame[-1,]
    
    return(concatFrame)
}


concat <- function (ext, form, outform='fasta',writeData=TRUE) {
    if (form == 'nexus') {
        outData = nexConCat(ext, 'nexus')
    }
    else {
        outData = baseConCat(ext, form)
    }
    if (writeData == TRUE) {
        write.fasta(as.list(outData$Sequence), outData$Species, nbchar = 60, "Output.fas", open = 'w')
    }
    
    if (writeData == TRUE & outform != 'fasta') {
        fasta <- read.alignment("Output.fas", format = "fasta")
    
        if (outform != 'fasta') {
            seq_list=fasta$seq
            seq_name=fasta$nam
            seq=list()
            for (i in 1:length(seq_list)) {
                if (outform == 'nexus' | outform=='phylip_relaxed') {
                    seq[[i]]<-sapply(c(strsplit(sapply(seq_list[[i]],as.character),"")),as.character)
                    names(seq)[i]<-as.character(fasta$nam[i])
                }
                else if (outform=='phylip_interleaved' | outform=='phylip_sequential') {
                    seq[[i]]<-sapply(c(strsplit(sapply(seq_list[[i]],as.character),"")),as.character)
                    names(seq)[i]<-as.character(paste(unlist(sapply(strsplit(sapply(seq_name[[i]],as.character),""),as.list)[1:10]),collapse=""))
                }
            }
        }
    
    
        if (outform == 'nexus') {
            outData= seq
            write.nex.data(outData, file= "Output.nex", interleaved = TRUE, charsperline = 100)
            unlink("Output.fas", recursive = FALSE)
        }
        if (outform=='phylip_relaxed'){
            outData=seq
            write.dna(outData,file="Output_phy_relaxed.phy" , format = "interleaved")
            unlink("Output.fas", recursive = FALSE)
        }
        if (outform=='phylip_interleaved'){
            outData=seq
            write.dna(outData,file="Output_phy_interleaved.phy" , format = "interleaved")
            unlink("Output.fas", recursive = FALSE)
        }
        if (outform=='phylip_sequential'){
            outData=seq
            write.dna(outData,file="Output_phy_sequential.phy" , format = "sequential")
            unlink("Output.fas", recursive = FALSE)
        }
    }
    
    return(outData)
}




#x=concat('.nex', 'nexus','phylip_interleaved', writeData=TRUE)




