\name{BlastSeqs}
\alias{BlastSeqs}
\title{
Run BLAST queries from R
}
\description{
Wrapper to run \href{https://blast.ncbi.nlm.nih.gov/Blast.cgi}{BLAST} queries using the commandline BLAST tool directly from R.
Can operate on an \code{\link[Biostrings]{XStringSet}} or a \code{FASTA} file. 

This function requires the BLAST+ commandline tools, which can be downloaded 
\href{https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download}{here}.
}
\usage{
BlastSeqs(seqs, BlastDB, 
              blastType=c('blastn', 'blastp', 'tblastn', 'blastx', 'tblastx'),
              extraArgs='', verbose=TRUE)
}
\arguments{
  \item{seqs}{
    Sequence(s) to run BLAST query on. This can be either an 
    \code{\link[Biostrings]{XStringSet}} or a path to a FASTA file.
  }
  \item{BlastDB}{
    Path to FASTA file in a pre-built BLAST Database. These can be built using either \code{\link{MakeBlastDb}} from R or
    the commandline \code{makeblastdb} function from BLAST+. For more information
    on building BLAST DBs, see \href{https://www.ncbi.nlm.nih.gov/books/NBK569841/}{here}.
  }
  \item{blastType}{
    Type of BLAST query to run. See 'Details' for more information on available types.
  }
  \item{extraArgs}{
    Additional arguments to be passed to the BLAST query
    executed on the command line. This should be a single character string.
  }
  \item{verbose}{
    Should output be displayed?
  }
}
\details{
BLAST implements multiple types of search. Available types are the following:
  \itemize{
    \item{
      \code{blastn}: Nucleotide sequences against database of nucleotide sequences
    }
    \item{
      \code{blastp}: Protein sequences against database of protein sequences
    }
    \item{
      \code{tblastn}: Protein sequences against translated database of nucleotide sequences
    }
    \item{
      \code{blastx}: Translated nucleotide sequences against database of protein sequences
    }
    \item{
      \code{tblastx}: Translated nucleotide sequences against translated database of nucleotide sequences
    }
}
Different BLAST queries require different inputs. The function will throw an error
if the input data does not match expected input for the requested query type.

Input sequences for \code{blastn},
\code{blastx}, and \code{tblastx} should be nucleotide data.

Input sequences for
\code{blastp} and \code{tblastn} should be amino acid data.

Database for \code{blastn}, \code{tblastn}, \code{tblastx} should be nucleotide data.

Database for \code{blastp} and \code{blastx} should be amino acid data.
}
\value{
Returns a data frame (\link[base]{data.frame}) of results of the BLAST query.
}
\seealso{
\code{\link{MakeBlastDb}}
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}
\examples{
#
}
