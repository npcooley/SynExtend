###### -- fastq dump wrapper --------------------------------------------------
# a very barebones wrapper for fastq-dump and has no overhead checking
# author: nicholas cooley
# contact: npc19@pitt.edu / npcooley@gmail.com

FastQFromSRR <- function(SRR,
                         ARGS = list("--gzip" = NULL,
                                     "--skip-technical" = NULL,
                                     "--readids" = NULL,
                                     "--read-filter" = "pass",
                                     "--dumpbase" = NULL,
                                     "--split-3" = NULL,
                                     "--clip" = NULL),
                         KEEPFILES = FALSE) {
  
  
  DIR <- tempdir()
  w1 <- unname(sapply(X = ARGS,
                      FUN = function(x) {
                        is.null(x)
                      },
                      USE.NAMES = FALSE))
  # w2 <- sapply(X = ARGS[!w1],
  #              FUN = function(x) {
  #                paste(x,
  #                      collapse = " ")
  #              })
  
  # construct query:
  FASTQDUMP <- paste("fastq-dump",
                     paste(names(ARGS)[w1],
                           collapse = " "),
                     paste(names(ARGS)[!w1],
                           unlist(ARGS[!w1]),
                           collapse = " "),
                     paste("--outdir", DIR),
                     paste("--accession", SRR))
  print(paste("Running fastq-dump command:\n",
              FASTQDUMP))
  system(command = FASTQDUMP)
  
  SeqsFound <- list.files(path = DIR,
                          pattern = SRR)
  
  if (length(SeqsFound) > 0) {
    ans <- lapply(X = SeqsFound,
                  FUN = function(x) {
                    readQualityScaledDNAStringSet(filepath = paste0(DIR,
                                                                    "/",
                                                                    x))
                  })
    if (KEEPFILES) {
      
      for (m1 in seq_along(SeqsFound)) {
        system(command = paste("mv",
                               paste(paste0(DIR,
                                            "/",
                                            SeqsFound[m1])),
                               SeqsFound[m1]))
      } # end m1 loop
    } else {
      system(command = paste("rm",
                             paste(paste0(DIR,
                                          "/",
                                          SeqsFound))))
    } # end logical check
    
    return(ans)
  } else {
    print("No seqs returned.")
  }
  
}

