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
