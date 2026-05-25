###### -- Reframe genecalls into subfeature bounds ----------------------------
# reframe a genecall DataFrame into an explicitly square representation of the 
# feature calls with a key to reference back to the source DataFrame,
# while placing 

###### -- NOTES ---------------------------------------------------------------
# assume object is already appropriately sorted
# convert to data.frame to make behavior expectations simpler internally

###### -- FUNCTION ------------------------------------------------------------

FrameDownward <- function(genecalls) {
  
  TempSubFeatures <- unlist(genecalls$Range)
  TempFeatureBlocks <- lengths(genecalls$Range)
  SubKey <- unlist(lapply(X = TempFeatureBlocks,
                          FUN = function(x) {
                            seq(x)
                          }))
  SubFeatureBlocks <- data.frame("Index" = rep(x = genecalls$Index,
                                               times = TempFeatureBlocks),
                                 "Strand" = rep(x = genecalls$Strand,
                                                times = TempFeatureBlocks),
                                 "Left" = start(TempSubFeatures),
                                 "Right" = end(TempSubFeatures),
                                 "Key" = rep(x = seq_along(TempFeatureBlocks),
                                             times = TempFeatureBlocks),
                                 "SubKey" = SubKey,
                                 "Phase" = unlist(genecalls$Phase))
  SubFeatureBlocks <- SubFeatureBlocks[order(SubFeatureBlocks$Index,
                                             SubFeatureBlocks$Left), ]
  rownames(SubFeatureBlocks) <- NULL
  
  res <- SubFeatureBlocks
  attr(x = res,
       which = "origin") <- genecalls
  
  return(res)
}
