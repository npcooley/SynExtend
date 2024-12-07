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

size_t LT_dumpOutput(LoserTree *tree, void *output_buffer){
	size_t nbytes = tree->e_size * tree->cur_output_i;
	memcpy(output_buffer, tree->output, nbytes);
	tree->cur_output_i = 0;
	return nbytes;
}

size_t LT_fdumpOutput(LoserTree *tree, FILE *f){
	// assume f is a valid file
	//printf("\n\tWriting %d values\n", tree->cur_output_i);
	size_t nbytes = tree->e_size * tree->cur_output_i;
	if (!nbytes) return 0;
	size_t nwrote = fwrite(tree->output, 1, nbytes, f);
	if(nwrote != nbytes)
		error("Failed to write to file! (tried to write %zu bytes, wrote %zu bytes)",
			nbytes, nwrote);
	tree->cur_output_i = 0;
	tree->nwritten += nwrote / tree->e_size;
	return nwrote;
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