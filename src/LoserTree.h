#ifndef LOSERTREE_H
#define LOSERTREE_H

/*
 * This is intended to be a type-unaware LoserTree (Tournament Tree)
 * implementation that supports dynamic refilling of bins for external sorting.
 *
 * The expected usage is the following:
 *  - Initialize tree with LT_alloc()
 *  - for each block to merge:
 *      | read some number of elements into a buffer b
 *      | call LT_fillBin() to assign the bin to b
 *  - initialize tree values with LT_initGame()
 *  - while elements remain in blocks:
 *      | call LT_popOutput(tree)
 *      | if output buffer is full, call LT_dumpOutput
 *      | if tree->empty_bin != -1, refill the empty bin
 *          (e.g., read new elements into buffer, then LT_fillBin() again)
 *      | call LT_updateTree(tree)
 *
 * If you don't care about dynamic bin refilling this could be greatly
 * streamlined, but because I'm assuming block size >> bin size, this
 * added complexity is necessary.
 *
 * For processing with output streamed to a file, use LT_runFileGame(). This
 * will keep running games until a bin is emptied, at which point the index
 * of the bin is returned. The calling function should then call LT_refillBin().
 * If the bin does not need to be refilled, the caller still needs to call
 * LT_refillBin(tree, i, 0, NULL).
 *
 * NOTE: the responsibility of keeping track of the buffers is placed on the
 *  calling function, NOT the LoserTree struct. LoserTree->bins is simply an
 *  array of pointers that will iterate along allocated void* memory. The memory
 *  pointed to by each LoserTree->bins[i] will NOT be alloc'd or free'd, only
 *  the pointers themselves.
 */

typedef struct LoserTree {
  int nbins;
  int full_bins;
  int empty_bin;
  int output_size;
  int cur_output_i;
  size_t e_size; // element size, in bytes
  int *binsize;
  void **bins;
  void *output;
  int *values;
  long nwritten;
  int (*compare)(const void *a, const void *b);
} LoserTree;

LoserTree* LT_alloc(int nbins, int output_size, size_t element_size,
                    int (*compare)(const void *a, const void *b));
void LT_fillBin(LoserTree *tree, int bin, int nelem, void *input);
void LT_initGame(LoserTree *tree);
void LT_popOutput(LoserTree *tree);
void LT_updateTree(LoserTree *tree);
void LT_refillBin(LoserTree *tree, int bin, int nelem, void *input);
size_t LT_dumpOutput(LoserTree *tree, void *output_buffer);
size_t LT_fdumpOutput(LoserTree *tree, FILE *f);
size_t LT_fdumpOutputInplace(LoserTree *tree, size_t block_end,
                            FILE *f, long int *remaining, long int **offsets);
int LT_runFileGame(LoserTree *tree, FILE *f);
int LT_runInplaceFileGame(LoserTree *tree, size_t block_end,
                          FILE *f, long int *remaining, long int **offsets);
void LT_free(LoserTree *tree);

#endif
