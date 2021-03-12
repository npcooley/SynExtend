# SynExtend 1.3.5

* `PairSummaries` now translates sequences based on `transl_table` attributes provided by gene calls
* `PairSummaries` now uses a generic model for predicting PID
* `gffToDataFrame` now parses out the `transl_table` attribute

# SynExtend 1.3.2

* Minor changes to `NucleotideOverlap`
* Major changes to `PairSummaries` - can now take in objects of class `Genes` build by the DECIPHER function `FindGenes()`

# SynExtend 0.99.30

* Vignette and help files edited for clarity

# SynExtend 0.99.1

* `SynExtend` submitted to Bioconductor
* Added function `gffToDataFrame`
* Added function `NucleotideOverlap`
* Added function `PairSummaries`