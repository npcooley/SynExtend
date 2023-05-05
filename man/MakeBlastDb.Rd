\name{MakeBlastDb}
\alias{MakeBlastDb}
\title{
Create a BLAST Database from R
}
\description{
Wrapper to create \href{https://blast.ncbi.nlm.nih.gov/Blast.cgi}{BLAST} databases for subsequent queries using the commandline BLAST tool directly from R.
Can operate on an \code{\link[Biostrings]{XStringSet}} or a \code{FASTA} file. 

This function requires the BLAST+ commandline tools, which can be downloaded 
\href{https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download}{here}.
}
\usage{
MakeBlastDb(seqs, dbtype=c('prot', 'nucl'), 
          dbname=NULL, dbpath=NULL,
          extraArgs='', createDirectory=FALSE, 
          verbose=TRUE)
}
\arguments{
  \item{seqs}{
    Sequence(s) to create a BLAST database from. This can be either an 
    \code{\link[Biostrings]{XStringSet}} or a path to a FASTA file.
  }
  \item{dbtype}{
    Either \code{'prot'} for amino acid input, or \code{'nucl'} for nucleotide input.
  }
  \item{dbname}{
    Name of the resulting database. If not provided, defaults to
    a random string prefixed by \code{blastdb}.
  }
  \item{dbpath}{
    Path where database should be created. If not provided, 
    defaults to \code{\link{TMPDIR}}.
  }
  \item{extraArgs}{
    Additional arguments to be passed to the query
    executed on the command line. This should be a single character string.
  }
  \item{createDirectory}{
    Should a directory be created for the database if it doesn't exist? If \code{FALSE}, the function will throw an error instead of creating a directory.
  }
  \item{verbose}{
    Should output be displayed?
  }
}
\details{
This offers a quick way to create BLAST databases from R. This function essentially wraps the \code{makeblastdb} commandline function. All arguments supported by \code{makeblastdb} are supported in the \code{extraArgs} argument.
}
\value{
Returns a length 2 named character vector specifying the name of the BLAST database and the path to it.
}
\seealso{
\code{\link{BlastSeqs}}
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}
\examples{
#
}
