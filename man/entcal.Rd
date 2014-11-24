\name{entcal}
\alias{entcal}
\title{Class \code{"entropy"}}
\usage{
entcal(file, type = "DNA")
}
\arguments{
    \item{file}{Sequence alignment file name specified by either a variable of mode character,
    or a double-quoted string; if \code{file = ""} (the default) then the
    program terminates with message}
    \item{type}{Alignment type. Can be DNA or Protein}
}

\description{
Calculates entropy value for the given input alignment file.
}

\author{Ambuj Kumar \email{ambuj@ufl.edu}}