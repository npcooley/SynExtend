#include "SEutils.h"
#include "LoserTree.h"

/*
 * See comments on expected usage in LoserTree.h
 */

LoserTree* LT_alloc(int nbins, int output_size, size_t element_size,
										int (*compare)(const void *a, const void *b)){
	// going to fill bins later, for now just leave them

	LoserTree *tree = malloc(sizeof(LoserTree));

	// we need an even power of two for the number of bins
	int actual_bins = 1;
	while(actual_bins < nbins) actual_bins <<= 1;

	// allocating the bins
	tree->nbins = actual_bins;
	tree->full_bins = 0;

	// values will always hold indices, so we know it'll be ints
	int *values = malloc(sizeof(int) * (actual_bins) * 2);
	int *binsize = malloc(sizeof(int) * actual_bins);
	void **bins = malloc(sizeof(void *) * actual_bins);

	for(int i=0; i<actual_bins; i++){
		binsize[i] = 0;
		bins[i] = NULL;
		values[i] = -1;
		values[i+actual_bins] = i; // fill second half of array with indices
	}
	tree->binsize = binsize;
	tree->bins = bins;
	tree->values = values;

	tree->empty_bin = -1;
	tree->output_size = output_size;
	tree->output = malloc(output_size*element_size);
	tree->cur_output_i = 0;
	tree->e_size = element_size;

	tree->compare = compare;

	tree->nwritten = 0;

	return tree;
}

static void LT_playgame(LoserTree *tree, int *a, int *b){
	// helper function
	// assigns the index of smaller value (loser) to a, and the other to b
	if(!tree->binsize[*b]) return;
	if(!tree->binsize[*a] || tree->compare(tree->bins[*a], tree->bins[*b]) > 0){
		// need to swap the values
		int tmp = *a;
		*a = *b;
		*b = tmp;
	}
	return;
}

static int LT_playRecursiveGameAtNodeI(LoserTree *tree, int i){
	if(i >= tree->nbins) return (i - tree->nbins);

	int left = LT_playRecursiveGameAtNodeI(tree, 2*i); // left
	int right = LT_playRecursiveGameAtNodeI(tree, 2*i+1); // right

	int smaller=left, larger=right;
	if(!tree->binsize[right]){
		// if right doesn't exist, it counts as infinity (larger)
		larger = right;
		smaller = left;
	} else if(!tree->binsize[left]){
		// if left doesn't exist, it counts as infinity (larger)
		larger = left;
		smaller = right;
	} else {
		// otherwise, larger is left if cmp(a,b) > 0
		larger = tree->compare(tree->bins[left], tree->bins[right]) < 0 ? right : left;
		smaller = larger == left ? right : left;
	}

	// now "left" is always the index of the smaller element
	// and "right" is always the index of the larger element

	// in loser trees, the smaller element goes on, larger stays
	tree->values[i] = larger;
	return smaller;
}

void LT_fillBin(LoserTree *tree, int bin, int nelem, void *input){
	// fills a bin of the LT with some amount of data
	// input data should be preallocated, no memory copying will be done
	if(bin > tree->nbins)
		error("Attempted to fill out-of-bounds bin in LoserTree!");
	if(tree->binsize[bin] == 0 && nelem) tree->full_bins++;
	tree->binsize[bin] = nelem;
	tree->bins[bin] = nelem ? input : NULL;
	if(nelem && tree->empty_bin == bin)
		tree->empty_bin = -1;
	return;
}

void LT_refillBin(LoserTree *tree, int bin, int nelem, void *input){
	// should be called when the last popOutput() emptied a bin
	// then the top element will still be the (new) empty bin
	if(nelem) LT_fillBin(tree, bin, nelem, input);
	LT_updateTree(tree);
	return;
}

void LT_initGame(LoserTree *tree){
	tree->values[0] = LT_playRecursiveGameAtNodeI(tree, 1);
	return;
}

void LT_popOutput(LoserTree *tree){
	if(tree->output_size <= tree->cur_output_i)
		error("Tried to pop output from LoserTree but buffer is full!");
	size_t size = tree->e_size;
	int cur_min = tree->values[0];
	if(!tree->binsize[cur_min])
		error("Tried to pop LoserTree output from an empty bin!");

	// void_deref(v, i, size) is equivalent to v[i] for a void*
	void *to_write = void_deref(tree->output, (tree->cur_output_i)++, size);

	// copy the top value into the output array
	memcpy(to_write, tree->bins[cur_min], size);


	if(--(tree->binsize[cur_min])){
		// if there's still elements in the bin, advance the pointer by one
		tree->bins[cur_min] = void_deref(tree->bins[cur_min], 1, size);
		tree->empty_bin = -1;
	} else {
		// otherwise the bin is empty, mark the signal value so we can populate it
		// in the calling function
		tree->empty_bin = cur_min;
		tree->full_bins--;
		tree->bins[cur_min] = NULL;
	}

	// note that I'm not going to update the tree here
	// we'll do that in LT_updateTree
	return;
}

void LT_updateTree(LoserTree *tree){
	// this function will update the tree with the next value in the bin
	// assume we've already called LT_popOutput
	int last_popped = tree->values[0];
	int cur_node = last_popped + tree->nbins;
	int v1 = last_popped;
	int v2;

	while(cur_node) {
		v2 = tree->values[cur_node];
		LT_playgame(tree, &v1, &v2);
		tree->values[cur_node] = v2;
		cur_node /= 2; // this moves to parent node
	};

	tree->values[cur_node] = v1;

	return;
}

int LT_runFileGame(LoserTree *tree, FILE *f){
	/*
	 * this function will run games until a bin is emptied
	 * when this happens, control returns to the caller to
	 * allow the bin to be refilled
	 */
	int retval = -1;
	while(tree->full_bins){
		LT_popOutput(tree);
		if(tree->cur_output_i == tree->output_size)
			LT_fdumpOutput(tree, f);
		if(tree->empty_bin != -1){
			retval = tree->empty_bin;
			tree->empty_bin = -1;
			return retval;
		}
		LT_updateTree(tree);
	}
	return -1;
}

int LT_runInplaceFileGame(LoserTree *tree, size_t block_end,
													FILE *f, long int *remaining, long int **offsets){
	/*
	 * Same as LT_runFileGame(), but does it in-place
	 */
	int retval = -1;
	while(tree->full_bins){
		LT_popOutput(tree);
		if(tree->cur_output_i == tree->output_size)
			LT_fdumpOutputInplace(tree, block_end, f, remaining, offsets);
		if(tree->empty_bin != -1){
			retval = tree->empty_bin;
			tree->empty_bin = -1;
			return retval;
		}
		LT_updateTree(tree);
	}
	return -1;
}

size_t LT_dumpOutput(LoserTree *tree, void *output_buffer){
	size_t nbytes = tree->e_size * tree->cur_output_i;
	memcpy(output_buffer, tree->output, nbytes);
	tree->cur_output_i = 0;
	return nbytes;
}

size_t LT_fdumpOutput(LoserTree *tree, FILE *f){
	// assume f is a valid file
	//printf("\n\tWriting %d values\n", tree->cur_output_i);
	size_t to_write = tree->cur_output_i;
	if (!to_write) return 0;
	size_t nwrote = fwrite(tree->output, tree->e_size, to_write, f);
	if(nwrote != to_write)
		error("Failed to write to file! (tried to write %zu elements, wrote %zu elements)",
			to_write, nwrote);
	tree->cur_output_i = 0;
	tree->nwritten += nwrote;
	return nwrote;
}

void reorganize_blocks(LoserTree *tree, size_t block_end, FILE *f,
 												long int *remaining, long int **offsets){
	size_t size = tree->e_size;
	int output_size = tree->cur_output_i;
	long int *offs = *offsets;
	// use as much space as possible, minimize r/w calls
	void *scratch_buf = malloc(size*tree->output_size);
	int nbins = tree->nbins;
	long int write_start, write_end, read_start, read_end, to_read;

	int last_bin = nbins-1;
	while(!remaining[last_bin]) last_bin--;
	write_end = block_end;
	for(int i=last_bin; i>=0; i--){
		// last bin is always in the right place
		// have to move all the rest of the bins down
		if(remaining[i]){
			read_start = offs[i];
			read_end = offs[i] + remaining[i];
			//printf("Moving bin %d from %ld-%ld to end at %ld\n", i, read_start, read_end, write_end);
			while(read_end != offs[i]){
				R_CheckUserInterrupt();
				to_read = output_size;
				if(read_end < offs[i] + to_read) to_read = read_end - offs[i];
				read_start = read_end - to_read;
				write_start = write_end - to_read;
				fseek(f, read_start*size, SEEK_SET);
				fread(scratch_buf, size, to_read, f);
				fseek(f, write_start*size, SEEK_SET);
				fwrite(scratch_buf, size, to_read, f);
				read_end = read_start;
				write_end = write_start;
			}
			offs[i] = write_start;
		}
	}

	free(scratch_buf);
	return;
}


size_t LT_fdumpOutputInplace(LoserTree *tree, size_t block_end,
														FILE *f, long int *remaining, long int **offsets){
	/*
	 * Function to dump output in-place
	 * Requires a bunch of extra values so we can keep track of stuff
	 * Input Variables:
	 *	-       tree: LoserTree structure
	 *	-	 f_reading: file pointer for reading values
	 *	-  f_writing: file pointer for writing values
	 *	-      start: starting line for the set of all blocks in current iteration
	 *	- remaining: pointer to int* containing # of elements remaining per block
	 *	-   offsets: pointer to int* with start position of each block
	 *
	 * The goal is essentially to insert the sorted block at the top of the area.
	 * We have k blocks, and we only pull values from the top of each block.
	 * In other words, if the block is size n and has m remaining, we only need
	 * to save the remaining n-m values.
	 *
	 * To do this, I'm going to allocate a second output buffer equal in size to
	 * the first output buffer, which will hold S elements. We copy off any
	 * remaining unseen elements to the second buffer from the first block of S
	 * elements in the file. Then we write the output buffer over the first S
	 * elements of the file. Now the second output buffer becomes our main output
	 * buffer, and we repeat the process with the second S elements. Once the sum
	 * of the sizes of the two output buffers is at most S, we write both to the
	 * same block and return.
	 *
	 * The crucial part here is making sure we update the `offsets` value
	 * appropriately. `remaining` is not going to change in one copy, since we
	 * just move values around. We also have to copy the final block so it's
	 * aligned with the very end of the block, so that we don't have issues with
	 * gaps in the middle of a block (in case of a partial copy during move).
	 */

	size_t size = tree->e_size;
	size_t start = tree->nwritten;
	int output_size = tree->cur_output_i;
	int nbins = tree->nbins;
	long int *offs = *offsets;

	if(!output_size) return start;
	int first_bin = 0;
	while(first_bin < nbins && !remaining[first_bin]) first_bin++;

	if(first_bin < nbins && offs[first_bin] < (start+output_size))
		reorganize_blocks(tree, block_end, f, remaining, offsets);

	// reset the writing pointer
	fseek(f, (tree->nwritten)*size, SEEK_SET);
	LT_fdumpOutput(tree, f);

	return 0;
}

void LT_free(LoserTree *tree){
	free(tree->bins);
	free(tree->binsize);
	free(tree->output);
	free(tree->values);
	free(tree);
	return;
}

// helper functions
void LT_print(LoserTree *tree){
	int nbins = tree->nbins;
	int nnodes = nbins*2;

	printf("%2d ", tree->values[0]);
	for(int i=1; i<nbins*2; i++){
	    printf("%2d ", tree->values[i]);
	}
	printf("\n\n");

  int width = 1;
  int i=1;
	printf("%2d ", tree->values[0]);
	while(i < nnodes-1){
		for(int j=0; j<width; j++){
			printf(" %2d ", tree->values[i+j]);
		}
		i += width;
		width *= 2;
		printf("\n");
	}
	printf("Output: %d\n", tree->values[0]);
}