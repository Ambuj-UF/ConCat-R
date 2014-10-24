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


