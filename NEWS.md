# SynExtend 1.10.1
* Memory leak bugfix

# SynExtend 1.9.21
* Adds new `RFDist` function to calculate Robinson-Foulds Distance
* Adds normalization for `GeneralizedRF` to make the distance between 0 and 1
* Minor bugfixes
* Documentation for new functions

# SynExtend 1.9.20
* Adds new `GeneralizedRF` function to calculated information-theoretic Generalized Robinson-Foulds distance between two dendrograms.
* Documentation for new function 
* New ProtWeaver predictor based off of `GeneralizedRF` metric
* New internal C source code for `GeneralizedRF`

# SynExtend 1.9.19
* Adds new `DPhyloStatistic` function to calculate the D-statistic for a binary state against a phylogeny following Fritz and Purvis (2009).
* Documentation for new function
* new internal C source code for `DPhyloStatistic`
* new internal C source code for random utility functions, currently only has functions to generate random numbers

# SynExtend 1.9.18
* Various internal improvements to presence/absence profile methods

# SynExtend 1.9.17
* Adds new prediction algorithm `GainLoss`
* Adds new internal C implementation of dendrograms, significantly faster than R dendrograms
* `ProtWeaver` methods `Behdenna` and `GainLoss` can now infer a species tree when possible
* Updates `Jaccard` and `Hamming` methods to use C implementations for distance calculation
* Adds `HammingGL` method to calculate Hamming distance of gain/loss events
* Minor bugfixes to `ProtWeaver` methods relating to subsetting
* Updates to various `man` pages

# SynExtend 1.9.16
* Removes `flatdendrapply`, function was already included in SynExtend
* minor bugfixes to `ProtWeaver`

# SynExtend 1.9.15
* Edits to `SelectByK`, function can work as intended, but is still too conservative at false positive removal.

# SynExtend 1.9.14
* Adds new function `flatdendrapply` for more options to apply functions to dendrograms. Function is used in `SuperTree`.
* Adds new function `SuperTree` to construct a species tree from a set of gene trees.
* Adds new dataset `SuperTreeEx` for `SuperTree` and `flatdendrapply` examples.

# SynExtend 1.9.13
* `SelectByK` function argument `ClusterSelect` switched to `ClusterScalar`. Cluster number selection now performed by fitting sum of total within cluster sum of squares to a right hyperbola and taking the ceiling of the half-max. Scalar allows a user to pick different tolerances to select more, or less clusters. Plotting behavior updated.

# SynExtend 1.9.12
* `simMat` class now supports empty indexing (`s[]`)
* `simMat` class now supports logical accession (`s[c(T,F,T),]`)

# SynExtend 1.9.11
* Added the function `SelectByK` that allows for quick removal of false positive predicted pairs based on a relatively simple k-means approach. Function is currently designed for use on the single genome-to-genome pairwise comparison, and not on an all-vs-all many genomes scale, though it may provide acceptable results on that scale.

# SynExtend 1.9.10
* New `simMat` class for `dist`-like similarity matrices that can be manipulated like base matrices
* Major update to `ProtWeaver` internals
  * All internal calls use `simMat` objects whenever possible to decrease memory footprint
  * Note `ContextTree` and `ProfDCA` require matrices internally
* `ProtWeb` objects now inherit from `simMat`
  * `ProtWeb.show` and `ProtWeb.print` now display predictions in a more natural way
  * `GetProtWebData()` deprecated; `ProtWeb` now inherits `as.matrix.simMat` and `as.data.frame.simMat` 
* New documentation pages for `simMat` class
* `GetProtWebData` documentation page reworked into `ProtWeb` documentation file.
* Fixes new bug in `Method='Hamming'` introduced in SynExtend 1.9.9

# SynExtend 1.9.9
* Fixes minor bug in `Method='Hamming'`
* Moves some code around

# SynExtend 1.9.8
* Major refactor to file structure of `ProtWeaver` to make individual files more manageable
* Adds new documentation files for individual prediction streams of `predict.ProtWeaver`

# SynExtend 1.9.7
* `BlockReconciliation` now returns a an object of class `PairSummaries`.

# SynExtend 1.9.6
* Fixes an error where warnings were mistakenly output to the user

# SynExtend 1.9.5
* Moves platform-specific files in `src/` (originally added by mistake)

# SynExtend 1.9.4
* Lots of bugfixes to `ResidueMI.ProtWeaver`
* `predict.ProtWeaver` now correctly labels rows/columns with gene names, not numbers
* `predict.ProtWeaver` now correctly handles `Subset` arguments
* `predict.ProtWeaver(..., Subset=3)` will correctly predict for all pairs involving gene `3` (or for any gene `x`, as long as `Subset` is a length 1 character or integer vector).

# SynExtend 1.9.2
* Adds residue MI method to `ProtWeaver`
* Various bugfixes for `ProtWeaver`

# SynExtend 1.7.14
* Various improvements for `GenRearrScen`, improves consistency and output formatting
* Major bugfix for `ProtWeaver` methods using dendrogram objects
* `ProtWeaver` now correctly guards against non-bifurcating dendrograms in methods that expect it

# SynExtend 1.7.13
* Introduces new `ProtWeaver` class to predict functional association of genes from COGs or gene trees. This implements many algorithms commonly used in the literature, such as MirrorTree and Inverse Potts Models.
* `predict(ProtWeaverObject)` returns a `ProtWeb` class with information on predicted associations. 
* Adds `BlastSeqs` to run BLAST queries on sequences stored as an `XStringSet` or `FASTA` file.

# SynExtend 1.7.12

* Updates to `ExtractBy` function. Methods and inputs simplified and adjusted, and significant improvements to speed.

# SynExtend 1.7.11

* Updated `NucleotideOverlaps` to now correctly registers hits in genes with a large degree of overlap with the immediately preceding gene.
* Fixed aberrant behavior in `BlockExpansion` where contigs with zero features could cause an error in expansion attempts.

# SynExtend 1.7.10

* `BlockReconciliation` now allows for setting either block size or mean PID for reconciliation precedence.

# SynExtend 1.7.9

* Added retention thresholds to `BlockReconciliation`.

# SynExtend 1.7.8

* `BlockExpansion` cases corrected for zero added rows.

# SynExtend 1.7.7

* Improvements to `BlockExpansion` and `BlockReconciliation` functions.

# SynExtend 1.7.5

* Began integration of `DECIPHER`'s `ScoreAlignment` function.

# SynExtend 1.7.4

* Fixed a bug in `PairSummaries` function.

# SynExtend 1.7.3

* Added `BlockExpansion` function.

# SynExtend 1.7.2

* Adjustment in how `PairSummaries` handles default translation tables and GFF derived gene calls.

# SynExtend 1.7.1

* Large changes under the hood to `PairSummaries`.
* Failure to accurately assign neighbors in some cases should now be fixed.
* Extraction of genomic features is now faster.
* `OffSetsAllowed` argument now defaults to `FALSE`. This argument may be dropped in the future in favor of a more complex function post-summary.
* Small edits to `SequenceSimilarity`

# SynExtend 1.5.4

* Added the function `SubSetPairs` that allows for easy trimming of predicted pairs based on conflicting predictions and / or prediction statistics.
* Added the function `EstimageGenomeRearrangements` that generates rearrangement scenarios of large scale genomic events using the double cut and join model.

# SynExtend 1.5.3

* Added the function `SequenceSimilarity` and made improvements to runtime in `DisjointSet`.

# SynExtend 1.4.1

* Fixed a small bug in consensus scores in `PairSummaries` where features facing on different strands had their score computed incorrectly.

# SynExtend 1.3.15

* Changes to concensus score in `PairSummaries`.

# SynExtend 1.3.14

* Major changes to the `PairSummaries` function and minor changes to `NucleotideOverlaps`, `ExtractBy`, and `FindSets`. Adjustments to the model that `PairSummaries` calls on to predict PIDs.

# SynExtend 1.3.13

* `ExtractBy` function has been added. Allows extraction of feature sequences into `XStingSet`s organized by the a `PairSummaries` object or the single linkage clusters implied by pairings within the `PairSummaries` objects.
* `DisjointSet` function added to extract single linkage clusters from a `PairSummaries` object.

# SynExtend 1.3.12

* `PairSummaries` now computes 4-mer distance between predicted pairs.

# SynExtend 1.3.11

* `PairSummaries` now returns a column titled Adjacent that provides the number of directly adjacent neighbor pairs to a predicted pair. Gap filling code adjusted.
* The function `FindSets` has been added and performs single linkage clustering on a pairs list as represented by vectors of integers using the Union-Find algorithm. Long term this function will have a larger wrapper function for user ease of access but will remain exposed.

# SynExtend 1.3.10

* `NucleotideOverlap` now passes it's GeneCalls object forward, allowing `PairSummaries` to forego inclusion of that object as an argument.

# SynExtend 1.3.9

* Minor vignette and suggested package changes.

# SynExtend 1.3.8

* `PairSummaries` now allows users to fill in specific matching gaps in blocks of predicted pairs with the arguments `AllowGaps` and `OffSetsAllowed`.

# SynExtend 1.3.7

* Adjustments to progress bars in both `PairSummaries` and `NucleotideOverlap`.
* PID prediction models in `PairSummaries` adjusted.

# SynExtend 1.3.6

* Contig name matching has been implemented. Scheme expects users to follow NCBI contig naming and gff formats, accepting contig names from gffs directly, and removing the first whitespace and everything thereafter from FASTA headers. Contig name matching can be disabled if users wish, using the argument `AcceptContigNames`, but ensuring that the correct contigs in GeneCalls objects are matched to the appropriate contigs in Synteny objects are then the user's responsibility.

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
