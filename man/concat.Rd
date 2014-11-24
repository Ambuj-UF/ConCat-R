\name{concat}
\alias{concat}
\title{Alignment concatenation function}
\usage{
concat(ext, form, outform, writeData = TRUE)
}
\arguments{
    \item{ext}{Input file extension. example: ".nex"}
    \item{form}{Input sequence alignment format. It can be nexus, phylip, fasta and clustal}
    \item{outform}{Output sequence alignment format. It can be nexus, phylip-relaxed,
    phylip-interleived, phylip-sequential, fasta and clustal}
    \item{writeData}{If TRUE (default = TRUE) the program creates a concatenated alignment file. 
    Otherwise it return the concatenated dataframe element}
}

\description{
Performs sequence alignment concatenation. It adds "?" for missing taxon sequence data.
}

\author{
Ambuj Kumar \email{ambuj@ufl.edu}
Mrinal Mishra \email{mrnlmishra@gmail.com}
}