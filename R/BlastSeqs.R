##### -- Function to run BLAST on query sequences ------------------------------
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

BlastSeqs <- function(seqs, BlastDB, 
                      blastType=c('blastn', 'blastp', 'tblastn', 'blastx', 'tblastx'),
                      extraArgs='', verbose=TRUE){
  if (length(blastType) > 1) blastType <- blastType[1]
  stopifnot("'seqs' must be XStringSet or path to a .fasta file"=
              any(is(seqs, 'DNAStringSet'),  
                is(seqs, 'AAStringSet'),
                is(seqs, 'XStringSet'),
                is(seqs, 'character')))
  
  stopifnot("'BlastDB' file does not exist"=file.exists(BlastDB))
  
  stopifnot('blastType must be a single string'=all(is(blastType, 'character'), length(blastType)==1))
  stopifnot('Undefined blast type'=(blastType %in% c('blastn', 'blastp', 'tblastn', 'blastx', 'tblastx')))
  bcheck <- pipe(paste0('which ', blastType))
  stopifnot('Requested blast not found. Check if the BLAST commandline tool is installed!'=
              length(readLines(bcheck)) > 0)
  close(bcheck)
  if (verbose) cat(blastType, 'installation verified...\n\n')
         
  if (!is(seqs, 'character')){
    if (is(seqs, 'AAStringSet') && blastType %in% c('blastn', 'blastx', 'tblastx'))
      stop("'", blastType, "'expects nucleotide data as input (received AAStringSet!)")
    if (is(seqs, 'DNAStringSet') && blastType %in% c('blastp', 'tblastn'))
      stop("'", blastType, "' expects protein data as input (received DNAStringSet!)")
    if (verbose) cat('Input is XStringSet, creating temporary .fasta file...\n\n')
    f <- tempfile()
    writeXStringSet(seqs, filepath=f, format='FASTA')
  } else {
    stopifnot('File does not exist'=file.exists(seqs))
    stopifnot(readBStringSet(seqs, format='fasta', nrec=1)) #error message never shows but works correctly
    if (verbose) cat('Input file located...\n\n')
    f <- seqs
  }
  
  query <- paste0(blastType, ' -query ', f, ' -db ', BlastDB, ' -outfmt 6 ', extraArgs)
  if (verbose) {
    cat('Running query:\n', query, '\n\n', sep='')
    t0 <- Sys.time()
  }
  blastres <- pipe(query)
  blastResults <- tryCatch({
      b <- read.table( blastres )
      colnames( b ) <- c( "QueryID",  "SubjectID", "Perc.Ident",
                                   "Alignment.Length", "Mismatches", "Gap.Openings", "Q.start", "Q.end",
                                   "S.start", "S.end", "E", "Bits" )
      
      if (verbose) {
        cat('Success!\nOutput head:\n\n') 
        print(head(b))
        print(difftime(Sys.time(), t0))
      }
      return(b)
    },
    error= function(e) NULL,
    finally= function(){ close(blastres) })
  
  if (!is.null(blastResults)) return(blastResults)
  else invisible(blastResults)
}
