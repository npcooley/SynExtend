/*
 * Out of memory clustering with fast label propagation
 * Tentative Name: ShoalFinder
 * Author: Aidan Lakshman
 *
 * This set of functions creates 5 files and a directory (using dummy names):
 *	-   OOMhashes: (directory) mapping between hash value and vertex name, along with index assignment
 * 	-  counts.bin: binary file containing edge counts for each vertex, and later cluster assignments
 * 	-  queue1.bin: out of memory queue for traversing nodes
 *  -  queue2.bin: second queue to improve r/w over jumping around queue1
 *  -     csr.bin: csr-compressed graph structure. Contains n+1 uint64_t values corresponding to vertex
 *		 						 indices, where the k'th value denotes the start position in the file for vertex k.
 *		 						 After the indices follow edge information of the form `d1 w1 d2 w2 d3 w3...`, where
 *		 						 each `d` corresponds to the destination of that edge, and each `w` the weight of
 *		 						 that edge. Indices index into this, e.g., if the first two values are 0 100 then the
 *		 						 outgoing edges from the first vertex are the first 100 edge entries (0-99).
 *  - outfile.tsv: .tsv file returned to R, contains two tab-separated columns (vertex name, cluster)
 *
 * Additional Notes:
 *  - sizeof(char) is guaranteed to be 1 (see https://en.wikipedia.org/wiki/Sizeof). 1 is used instead to
 *		eliminate the extra function call and simplify code somewhat.
 *
 * TODOs:
 *	- At some point it's probably worth refactoring all file accesses into some kind of struct w/ accessors
 *  	-> something like a virtual array object that's actually r/w to disk, could be useful in future
 *		-> this implementation should use mmap (and Windows equivalent when necessary) to improve random r/w
 *		-> this is a very far off wishlist item, not necessary
 *  - switch back to uint16_t for node name length
 *	- check performance of enabling use_limited_nodes
 */

#include "SEutils.h"
#include "SynExtend.h"

#define h_uint uint64_t
#define uint uint32_t
#define l_uint uint64_t

/*
 * common limits are defined in limits.h
 * 	PAGE_SIZE: size in bytes of a page
 *  PATH_MAX: max size in bytes of a path name
 */

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif

#define MAX_NODE_NAME_SIZE 254 // max size of a vertex name (char array will have 2 extra spaces for terminator and flag)
#define NODE_NAME_CACHE_SIZE 4096
#define FILE_READ_CACHE_SIZE 4096 // used for mergesorting files

const int L_SIZE = sizeof(l_uint);
const int LEN_SIZE = sizeof(h_uint);
const int MAX_READ_RETRIES = 10;
const char HASH_FNAME[] = "hashfile";
const char HASH_INAME[] = "hashindex";
const char CONS_TMPNAME[] = "tmpgraph";
const char CONSENSUS_CSRCOPY1[] = "tmpcsr1";
const char CONSENSUS_CSRCOPY2[] = "tmpcsr2";
const char CONSENSUS_CLUSTER[] = "tmpclust";

// set this to 1 if we should sample edges rather than use all of them
// MAX_EDGES_EXACT is a soft cap -- if above this, we sample edges probabalistically
const int use_limited_nodes = 0;
const l_uint MAX_EDGES_EXACT = 20000;
const int PRINT_COUNTER_MOD = 811;
const int PROGRESS_COUNTER_MOD = 3083;


/**********************/
/* Struct Definitions */
/**********************/
typedef struct {
	l_uint ctr1;
	l_uint ctr2;
} double_lu;

typedef struct {
	h_uint strlength;
	char s[MAX_NODE_NAME_SIZE];
	h_uint hash;
	l_uint count;
} msort_vertex_line;

typedef struct ll {
	l_uint id;
	double w;
	struct ll* next;
} ll;

typedef struct {
	h_uint len;
	h_uint hash;
	l_uint index;
} iline;

/**************************/
/* Core Utility Functions */
/**************************/

ll* insert_ll(ll* head, l_uint id, double w){
	ll *tmp = head;
	if(!tmp){
		tmp = malloc(sizeof(ll));
		tmp->id = id;
		tmp->w = w;
		tmp->next=NULL;
		return(tmp);
	}
	while(tmp->id != id && tmp->next) tmp = tmp->next;
	if(tmp->id!=id){
		tmp->next = malloc(sizeof(ll));
		tmp = tmp->next;
		tmp->id = id;
		tmp->w = w;
		tmp->next = NULL;
	} else {
		tmp->w += w;
	}

	return head;
}

static void safe_fread(void *buffer, size_t size, size_t count, FILE *stream){
	size_t found_values = fread(buffer, size, count, stream);
	if(found_values != count){
		// two scenarios:

		// 1. read past the end of the file (throw error and return)
		if(feof(stream))
			error("%s", "Internal error: fread tried to read past end of file.\n");

		// 2. some undefined reading error (retry a few times and then return)
		for(int i=0; i<MAX_READ_RETRIES; i++){
			// if we read a partial value, reset the counter back
			if(found_values) fseek(stream, -1*((int)found_values), SEEK_CUR);

			// try to read again
			found_values = fread(buffer, size, count, stream);
			if(found_values == count) return;
		}

		// otherwise throw an error
		error("Internal error: fread read %zu values (expected %zu).\n", found_values, count);
	}
	return;
};

static void safe_filepath_cat(const char* dir, const char* f, char *fname, size_t fnamesize){
	// fname should be preallocated
	char directory_separator;
#ifdef WIN32
	directory_separator = '\\';
#else
	directory_separator = '/';
#endif
	memset(fname, 0, fnamesize);
	// length is len(dir) + len(f) + 1 (separator) + 1 (terminator)
	snprintf(fname, strlen(dir)+strlen(f)+2, "%s%c%s", dir, directory_separator, f);
	return;
}

void errorclose_file(FILE *f1, FILE *f2, const char* message){
	fclose(f1);
	if(f2) fclose(f2);
	error("%s", message);
}

inline double sigmoid_transform(double w, const double slope){
	// should probably expose these at some point
	const double scale = 0.5;
	const double cutoff = 0.2;
	double r = 1 / (1+exp(-1*slope*(w-scale)));
	return r > cutoff ? r : 0;
}

h_uint hash_string_fnv32(const char *str){
	/*
	 * this is a Fowler-Noll-Vo hash function, it's fast and simple -- see wikipedia for constants
	 * https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
	 *
	 * hashes to a 32-bit uint that we can truncate
	*/
	const uint fnv_prime = 0x01000193;
	uint32_t hash = 0x811c9dc5;
	const char *s = str;
	char tmp;

	while(*s){
		tmp = *s++;
		hash ^= tmp;
		hash *= fnv_prime;
	}

	return (h_uint)hash;
}

h_uint hash_string_fnv(const char *str){
	/*
	 * Same as above, just hashes to a 64-bit uint
	 * only use this if h_uint = uint64_t, otherwise just stick to the 32 bit version (faster)
	*/
	const uint64_t fnv_prime = 0xcbf29ce484222325;
	uint64_t hash = 0x100000001b3;
	const char *s = str;
	char tmp;

	while(*s){
		tmp = *s++;
		hash ^= tmp;
		hash *= fnv_prime;
	}

	return hash;
}

/************************/
/* Comparison Functions */
/************************/


int l_uint_compar(const void* a, const void* b){
	double_lu aa = **(double_lu **)(a);
	double_lu bb = **(double_lu **)(b);
	if(aa.ctr2 - bb.ctr2)
		return aa.ctr2 - bb.ctr2;
	return aa.ctr1 - bb.ctr1;
}

int vertex_name_compar(const void* a, const void* b){
	msort_vertex_line *aa = *(msort_vertex_line **)(a);
	msort_vertex_line *bb = *(msort_vertex_line **)(b);
	if(aa->strlength != bb->strlength)
		return aa->strlength - bb->strlength;
	return strcmp(aa->s, bb->s);
}

int vertex_name_hash_compar(const void* a, const void* b){
	// names are sorted as length > hash > strcmp
	msort_vertex_line *aa = *(msort_vertex_line **)(a);
	msort_vertex_line *bb = *(msort_vertex_line **)(b);
	if(aa->strlength != bb->strlength)
		return aa->strlength - bb->strlength;
	if(aa->hash != bb->hash)
		return aa->hash < bb->hash ? -1 : 1;
	return strcmp(aa->s, bb->s);
}

int nohash_name_cmpfunc(const void *a, const void *b){
	// same as above, but skip the hash comparison
	const char *aa = *(const char **)a;
	const char *bb = *(const char **)b;

	// sort first by string length
	int v1 = strlen(aa), v2 = strlen(bb);
	if (v1 != v2) return v1 - v2;

	return strcmp(aa, bb);
}

/**************************/
/* File Mergesort Helpers */
/**************************/

void precopy_dlu1(const char* f1, const char* f2){
	// write and add index
	double_lu dlu = {1,0};
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");
	while(fread(&dlu.ctr2, L_SIZE, 1, orig)){
		fwrite(&dlu, sizeof(double_lu), 1, copy);
		dlu.ctr1++;
	}
	fclose(orig);
	fclose(copy);
	return;
}

void precopy_dlu2(const char* f1, const char* f2){
	// write flipped version (index, clust => clust, index)
	double_lu dlu = {0,0};
	l_uint prev_ind=0;
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");
	while(fread(&dlu, sizeof(double_lu), 1, orig)){
			prev_ind = dlu.ctr1;
			dlu.ctr1 = dlu.ctr2;
			dlu.ctr2 = prev_ind;
			fwrite(&dlu, sizeof(double_lu), 1, copy);
	}
	fclose(orig);
	fclose(copy);
	return;
}

void postcopy_dlu1(const char* f1, const char* f2){
	double_lu dlu = {0,0};
	l_uint prev_ind=0, ctr=0;
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");
	while(fread(&dlu, sizeof(double_lu), 1, orig)){
			if(prev_ind != dlu.ctr2){
				prev_ind = dlu.ctr2;
				dlu.ctr2 = ++ctr;
			} else {
				dlu.ctr2 = ctr;
			}
			fwrite(&dlu, sizeof(double_lu), 1, copy);
	}
	fclose(orig);
	fclose(copy);
	return;
}

void postcopy_dlu2(const char* f1, const char* f2){
	// write only the cluster into file
	// also can reindex such that the first cluster listed is cluster 1
	// (this can be removed)
	double_lu dlu = {0,0};
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");

	// uncomment these lines to make the first vertex have cluster 1
	/*
	l_uint max_found = 0, offset;
	while(fread(&dlu, sizeof(double_lu), 1, orig))
			if(dlu.ctr1 > max_found) max_found = dlu.ctr1;
	rewind(orig);
	fread(&dlu, sizeof(double_lu), 1, orig);
	offset = max_found - dlu.ctr1;
	rewind(orig);
	*/

	while(fread(&dlu, sizeof(double_lu), 1, orig)){
			// dlu.ctr1 = ((dlu.ctr1 + offset) % max_found) + 1;
			fwrite(&dlu.ctr1, L_SIZE, 1, copy);
	}
	fclose(orig);
	fclose(copy);
	return;
}


void precopy_vertexname(const char* f1, const char* f2){
	msort_vertex_line *aa = malloc(sizeof(msort_vertex_line));
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");
	while(fread(&(aa->strlength), LEN_SIZE, 1, orig)){
		memset(&(aa->s), 0, MAX_NODE_NAME_SIZE);
		safe_fread(&(aa->s), 1, aa->strlength, orig);
		safe_fread(&(aa->count), L_SIZE, 1, orig);
		aa->hash = hash_string_fnv(aa->s);
		fwrite(aa, sizeof(msort_vertex_line), 1, copy);
	}

	fclose(orig);
	fclose(copy);
	free(aa);
}

void postcopy_vertexname(const char* f1, const char* f2){
	msort_vertex_line *cur_elem, *tmp_elem;
	cur_elem = calloc(1, sizeof(msort_vertex_line));
	tmp_elem = calloc(1, sizeof(msort_vertex_line));
	FILE *orig = fopen(f1, "rb");
	FILE *copy = fopen(f2, "wb");

	fread(cur_elem, sizeof(msort_vertex_line), 1, orig);
	while(fread(tmp_elem, sizeof(msort_vertex_line), 1, orig)){
		if(tmp_elem->strlength == cur_elem->strlength && !strcmp(tmp_elem->s, cur_elem->s)){
			cur_elem->count += tmp_elem->count;
		} else {
			fwrite(&(cur_elem->strlength), LEN_SIZE, 1, copy);
			fwrite(&(cur_elem->s), cur_elem->strlength, 1, copy);
			fwrite(&(cur_elem->count), L_SIZE, 1, copy);
			memcpy(cur_elem, tmp_elem, sizeof(msort_vertex_line));
		}
	}

	// dont forget to copy the final element
	fwrite(&(cur_elem->strlength), LEN_SIZE, 1, copy);
	fwrite(&(cur_elem->s), cur_elem->strlength, 1, copy);
	fwrite(&(cur_elem->count), L_SIZE, 1, copy);

	free(cur_elem);
	free(tmp_elem);
	fclose(orig);
	fclose(copy);
	return;
}

void mergesort_clust_file(const char* f, const char* dir, size_t element_size,
															int (*compar)(const void *, const void *),
															void (*precopy)(const char*, const char*),
															void (*postcopy)(const char*, const char*)){
	/*
	 * general file mergesort function
	 * arguments:
	 *	-            f: file to sort
	 *	-          dir: directory to store junk files
	 *  - element_size: size of each element to read/write
	 *  -      *compar: function pointer used in qsort / mergesort. ensure proper casting.
	 *  -     *precopy: function to copy f into the first junk file
	 *  -    *postcopy: function to write final values back into f
	 *  notes:
	 *  - *precopy should open the file (assume it does not exist)
	 *  - *compar will provide void** values, make sure to double dereference
	 */

	// two read pointers, one write pointer
	FILE *f1_r1, *f1_r2, *f2_w;
	char *file1 = malloc(PATH_MAX);
	char *file2 = malloc(PATH_MAX);
	char *finalfile;

	//size_t dlu_size = sizeof(double_lu);

	// create the junk files we'll use
	safe_filepath_cat(dir, "tmp_ms1", file1, PATH_MAX);
	safe_filepath_cat(dir, "tmp_ms2", file2, PATH_MAX);

	// first, we'll use the cache to read in preprocessed sorted blocks of size `element_size`
	l_uint block_size = FILE_READ_CACHE_SIZE;
	l_uint total_lines = 0;

	// allocate space for data in a void*
	//void *read_cache[FILE_READ_CACHE_SIZE];
	void **read_cache = malloc(sizeof(void*) * FILE_READ_CACHE_SIZE);
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) read_cache[i] = malloc(element_size);

	// copy the original file into file 1
	precopy(f, file1);

	// open file, read in chunks, sort locally, write to file
	uint cachectr = 0;
	f1_r1 = fopen(file1, "rb");
	if(!f1_r1) error("%s", "Error opening file obtained from mergesort precopy");
	f2_w = fopen(file2, "wb");
	if(!f2_w) error("%s", "Error opening temporary mergesort file for writing");
	while(fread(read_cache[cachectr], element_size, 1, f1_r1)){
		cachectr++;
		total_lines++;
		if(cachectr == block_size){
			qsort(read_cache, cachectr, sizeof(void*), compar);
			for(int i=0; i<cachectr; i++)
				fwrite(read_cache[i], element_size, 1, f2_w);
			cachectr=0;
		}
	}
	if(cachectr){
		qsort(read_cache, cachectr, sizeof(void*), compar);
		for(int i=0; i<cachectr; i++)
			fwrite(read_cache[i], element_size, 1, f2_w);
	}

	fclose(f1_r1);
	fclose(f2_w);
	finalfile = file2;

	l_uint cur_lines = 0;
	int iter1, iter2, previt1, previt2;
	int flip = 0, cmp;
	void *tmp1 = malloc(element_size);
	void *tmp2 = malloc(element_size);
	char *f1, *f2;
	while(block_size < total_lines){
		// f1 is always the reading file, f2 the writing file
		f1 = flip ? file1 : file2;
		f2 = flip ? file2 : file1;
		flip = !flip;

		f1_r1 = fopen(f1, "rb");
		f1_r2 = fopen(f1, "rb");
		f2_w = fopen(f2, "wb");
		if(!(f1_r1 && f1_r2 && f2_w))
			error("%s", "Error opening temporary files in mergesort");
		// move second pointer forward to second block
		fseek(f1_r2, element_size*block_size, SEEK_CUR);

		// sort file 1 into file 2
		while(cur_lines < total_lines){
			// sort one block from file 1 into file 2 -- iter stores lines remaining
			iter1=total_lines - cur_lines;
			if(iter1 > block_size) iter1 = block_size;
			cur_lines += iter1;

			iter2=total_lines - cur_lines;
			if(iter2 > block_size) iter2 = block_size;
			cur_lines += iter2;

			previt1=iter1+1;
			previt2=iter2+1;
			while(iter1 || iter2){
				if(iter1 && iter1 != previt1){
					safe_fread(tmp1, element_size, 1, f1_r1);
					previt1 = iter1;
				}
				if(iter2 && iter2 != previt2){
					safe_fread(tmp2, element_size, 1, f1_r2);
					previt2 = iter2;
				}

				cmp = compar(&tmp1, &tmp2);
				if(iter1 && (!iter2 || cmp <= 0 )){
					fwrite(tmp1, element_size, 1, f2_w);
					iter1--;
				} else {
					fwrite(tmp2, element_size, 1, f2_w);
					iter2--;
				}
			}
			// advance pointers one block past where we just read:
			// if we move too far it doesn't really matter, we'll catch it on the next part
			fseek(f1_r1, element_size*block_size, SEEK_CUR);
			fseek(f1_r2, element_size*block_size, SEEK_CUR);
		}

		fclose(f1_r1);
		fclose(f1_r2);
		fclose(f2_w);
		cur_lines = 0;
		block_size *= 2;
		finalfile = f2;
	}
	// copy result back into f
	postcopy(finalfile, f);

	// free memory allocations
	free(tmp1);
	free(tmp2);
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) free(read_cache[i]);
	free(read_cache);
	free(file1);
	free(file2);
	return;
}

void copy_csrfile_sig(const char* dest, const char* src, l_uint num_v, const double w){
	l_uint *vbuf = malloc(FILE_READ_CACHE_SIZE*L_SIZE);
	double *wbuf = malloc(FILE_READ_CACHE_SIZE*sizeof(double));
	FILE *fd = fopen(dest, "wb");
	FILE *fs = fopen(src, "rb");

	int to_read=0, remaining=num_v+1;

	// first copy over all the vertex offsets
	while(remaining > 0){
		to_read = remaining > FILE_READ_CACHE_SIZE ? FILE_READ_CACHE_SIZE : remaining;
		remaining -= fread(vbuf, L_SIZE, to_read, fs);
		fwrite(vbuf, L_SIZE, to_read, fd);
	}

	// next copy over vertices and weights
	int cachectr = 0;
	while(fread(&vbuf[cachectr], L_SIZE, 1, fs)){
		fread(&wbuf[cachectr], sizeof(double), 1, fs);
		cachectr++;
		if(cachectr == FILE_READ_CACHE_SIZE){
			for(int i=0; i<cachectr; i++){
				wbuf[i] = w < 0 ? 0 : sigmoid_transform(wbuf[i], w);
				fwrite(&vbuf[i], L_SIZE, 1, fd);
				fwrite(&wbuf[i], sizeof(double), 1, fd);
			}
			cachectr = 0;
		}
	}
	if(cachectr){
		for(int i=0; i<cachectr; i++){
			wbuf[i] = w < 0 ? 0 : sigmoid_transform(wbuf[i], w);
			fwrite(&vbuf[i], L_SIZE, 1, fd);
			fwrite(&wbuf[i], sizeof(double), 1, fd);
		}
	}

	free(vbuf);
	free(wbuf);
	fclose(fd);
	fclose(fs);
	return;
}

void unique_strings_with_sideeffects(char **names, int num_to_sort, int *InsertPoint, uint *counts, int useCounts){
	/*
	 * This code is duplicated a lot, so just putting it here for consistency
	 * This uniques the set of strings **names and stores additional information
	 *  -     InsertPoint: number of unique strings
	 *  -          counts: number of each unique string (ignored if !useCounts)
	 */
	int insert_point = 0;
	uint cur_len;

	// first sort the array
	qsort(names, num_to_sort, sizeof(char*), nohash_name_cmpfunc);

	// next, unique the values
	if(useCounts) counts[0] = !!names[0][MAX_NODE_NAME_SIZE];
	cur_len = strlen(names[0]);
	for(int i=1; i<num_to_sort; i++){
		if(cur_len != strlen(names[i]) || (strcmp(names[i], names[insert_point]) != 0)){
			// if the string is different, save it
			insert_point++;
			if(useCounts) counts[insert_point] = 1;
			if(insert_point != i){
				memcpy(names[insert_point], names[i], MAX_NODE_NAME_SIZE+1);
				names[insert_point][MAX_NODE_NAME_SIZE] = 0;
			}
			cur_len = strlen(names[insert_point]);
		} else if(useCounts) {
			// else it's the same, so increment the corresponding count (last byte stores if it's an edge that should be counted)
			counts[insert_point] += !!names[i][MAX_NODE_NAME_SIZE];
		}
	}
	insert_point++;

	// side effect writes
	*InsertPoint = insert_point;
	return;
}

void batch_write_nodes(char **names, int num_to_sort, FILE *f){
	// just going to write everything to one big file, no checking for uniqueness
	// assume file opened in advance as ab
	// TODO: make strings one byte longer, use this as a flag to determine in vs out edge

	// assuming all the nodes are in a const char* array

	// filename
	char *tmp_s;

	// entries per unique string
	uint *stringcounts = malloc(sizeof(uint) * num_to_sort);
	l_uint tmpcount;

	// strlens
	h_uint cur_len;

	// these are sometimes used for temp storage
	int insert_point;

	unique_strings_with_sideeffects(names, num_to_sort, &insert_point, stringcounts, 1);

	// write all the lines to the end of the file
	for(int j=0; j<insert_point; j++){
		tmp_s = names[j];
		tmpcount = stringcounts[j];
		cur_len = strlen(tmp_s);
		fwrite(&cur_len, LEN_SIZE, 1, f);
		fwrite(tmp_s, 1, cur_len, f);
		fwrite(&tmpcount, L_SIZE, 1, f);
	}


	free(stringcounts);
	return;
}

l_uint node_vertex_file_cleanup(const char* dir, const char* fname, const char* countsfname, const char* indexfname, int v){
	// file is formatted as entries of type "%u%s%lu" (strlen, name, count)
	const size_t line_fixed_size = L_SIZE + LEN_SIZE; // size of line (not including vertex name)

	// sort file by string length and name
	if(v) Rprintf("\tSorting vertex names...");
	mergesort_clust_file(fname, dir, sizeof(msort_vertex_line), vertex_name_hash_compar, precopy_vertexname, postcopy_vertexname);
	if(v) Rprintf("done.\n\tRe-indexing vertices...\n");
	/*
	 * now we populate three files:
	 * - fname will hold the indices for where to find each variable (length, first 4 characters (0-padded), start index)
	 * - countsfname will hold the counts for each variable
	 */

	FILE *orig, *countsfile, *indexfile;
	l_uint tmpcount, curloc = 0, num_v=0;

	char *vertname = malloc(MAX_NODE_NAME_SIZE);
	iline *index_line = calloc(1, sizeof(iline));

	orig = fopen(fname, "rb");
	countsfile = fopen(countsfname, "wb");
	indexfile = fopen(indexfname, "wb");
	if(!orig || !countsfile || !indexfile) error("error opening files while cleaning up vertex hashmap");

	while(fread(&index_line->len, LEN_SIZE, 1, orig)){
		memset(vertname, 0, MAX_NODE_NAME_SIZE);
		safe_fread(vertname, 1, index_line->len, orig);
		safe_fread(&tmpcount, L_SIZE, 1, orig);

		index_line->index = curloc;
		index_line->hash = hash_string_fnv(vertname);

		// write length, hash, start index to index file
		fwrite(index_line, sizeof(iline), 1, indexfile);

		// write count to counts file
		fwrite(&tmpcount, L_SIZE, 1, countsfile);
		curloc += line_fixed_size + index_line->len;
		num_v++;
		if(!(num_v % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%lu vertices found\r", num_v);
			else R_CheckUserInterrupt();
		}
	}

	fclose(orig);
	fclose(countsfile);
	fclose(indexfile);
	free(index_line);
	free(vertname);
	if(v) Rprintf("\t%lu vertices found\n", num_v);
	return num_v;
}

l_uint lookup_node_index(char *name, FILE *findex, FILE *fhash, l_uint num_v){
	// binary search to find value in file

	// findex and hashf should be opened rb
	const size_t line_size = sizeof(iline);
	const h_uint namelen = strlen(name);

	h_uint tmplen;
	h_uint curhash = hash_string_fnv(name);
	int cmpresult;

	l_uint min_line = 0, max_line = num_v;
	int64_t prev_line = 0, cur_line = 0;

	char *vname = malloc(MAX_NODE_NAME_SIZE);
	iline *index_line = malloc(sizeof(iline));

	rewind(findex);
	while(min_line <= max_line){
		// move to center of current area
		cur_line = (max_line + min_line) / 2;
		fseek(findex, (cur_line-prev_line)*line_size, SEEK_CUR);

		safe_fread(index_line, sizeof(iline), 1, findex);
		cmpresult = curhash != index_line->hash; // should be 0 when equal
		if(cmpresult)
			cmpresult = curhash < index_line->hash ? -1 : 1;
		if(index_line->len == namelen && !cmpresult){
			memset(vname, 0, MAX_NODE_NAME_SIZE);
			fseek(fhash, index_line->index, SEEK_SET);
			fread(&tmplen, LEN_SIZE, 1, fhash);
			fread(vname, 1, tmplen, fhash);
			cmpresult = strcmp(name, vname);
			if(!cmpresult){
				free(index_line);
				free(vname);
				return cur_line;
			}
		}

		if(namelen < index_line->len){
			// name comes before current position in array
			max_line = cur_line-1;
		} else if(namelen > index_line->len){
			min_line = cur_line+1;
		} else if(cmpresult < 0){
			max_line = cur_line-1;
		} else {
			min_line = cur_line+1;
		}
		prev_line = cur_line+1;
	}

	// we should never get here if we read in all the vertices correctly
	free(vname);
	free(index_line);
	error("Failed to find vertex %s in hash lookup.", name);
}

void hash_file_vnames_batch(const char* fname, const char* dname, const char *hashfname,
	const char sep, const char line_sep, int v, int is_undirected){
	/*
	 * fname: .tsv list of edges
	 * dname: directory of hash codes
	 * hashfname: file to write vertices to
	 */
	FILE *f = fopen(fname, "rb");
	FILE *hashf = fopen(hashfname, "ab");

	// size + 1 so that there's space for marking if it's an outgoing edge
	char *vname = malloc(MAX_NODE_NAME_SIZE+1);
	char **namecache = malloc(sizeof(char*) * NODE_NAME_CACHE_SIZE);

	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) namecache[i] = malloc(MAX_NODE_NAME_SIZE+1);

	int cur_pos = 0, cachectr=0;
	char c = getc(f);
	l_uint print_counter = 0;

	if(v) Rprintf("\tReading file %s...\n", fname);

	while(!feof(f)){
		// going to assume we're at the beginning of a line
		// lines should be of the form `start end weight` or `start end`
		// separation is by char `sep`
		for(int iter=0; iter<2; iter++){
			memset(vname, 0, MAX_NODE_NAME_SIZE+1);
			cur_pos = 0;
			while(c != sep && c != line_sep){
				vname[cur_pos++] = c;
				c = getc(f);
				if(cur_pos == MAX_NODE_NAME_SIZE-1) // max size has to include the null terminator
					errorclose_file(f, hashf, "Node name is larger than max allowed name size.\n");

				if(feof(f)) errorclose_file(f, hashf, "Unexpected end of file.\n");
			}

			memset(namecache[cachectr], 0, MAX_NODE_NAME_SIZE+1);
			memcpy(namecache[cachectr], vname, strlen(vname));

			// mark if edge is outgoing -- always if undirected, otherwise only if the first node name
			namecache[cachectr][MAX_NODE_NAME_SIZE] = is_undirected ? 1 : !iter;

			cachectr++;
			if(cachectr == NODE_NAME_CACHE_SIZE){
				batch_write_nodes(namecache, cachectr, hashf);
				cachectr = 0;
			}

			// if lines are of the form `start end`, we need to leave c on the terminator
			if(c == sep)
				c = getc(f);
		}

		while(c != line_sep && !feof(f)) c = getc(f);
		if(c == line_sep) c=getc(f);
		print_counter++;
		if(!(print_counter % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%lu lines read\r", print_counter);
			else R_CheckUserInterrupt();
		}
	}
	if(cachectr) batch_write_nodes(namecache, cachectr, hashf);

	if(v) Rprintf("\t%lu lines read\n", print_counter);
	fclose(hashf);
	fclose(f);

	free(vname);
	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) free(namecache[i]);
	free(namecache);
	return;
}

void reformat_counts(const char* curcounts, const char* mastertable, l_uint n_vert, int add_self_loops){
	/*
	 * Creates a new table with cumulative counts
	 * leaves the old table unchanged, this will act as a temporary counts file later
	 */
	const uint l_size = L_SIZE;
	l_uint cumul_total = 0, curcount;
	FILE *tmptab = fopen(curcounts, "rb");
	FILE *mtab = fopen(mastertable, "wb+");
	int self_loop = add_self_loops ? 1 : 0;

	for(l_uint i=0; i<n_vert; i++){
		fwrite(&cumul_total, l_size, 1, mtab);
		safe_fread(&curcount, l_size, 1, tmptab);
		cumul_total += curcount + self_loop; // add an extra count for each node if we add self loops
	}

	// ending position of file
	fwrite(&cumul_total, l_size, 1, mtab);
	fclose(tmptab);
	fclose(mtab);
	return;
}

void add_self_loops_to_csrfile(const char *ftable, l_uint num_v, const double self_weight){
	// If self loops are included, the first entry for each node remains empty
	// here we'll fill it in with the node itself
	// thus, we can just write to whatever the first index of the value is and set the value to 0
	const uint entry_size = L_SIZE + sizeof(double);

	l_uint tmp_pos;

	FILE *mastertab = fopen(ftable, "rb+");
	if(!mastertab) error("%s", "error opening CSR file.\n");

	for(l_uint i=0; i<num_v; i++){
		fseek(mastertab, i*L_SIZE, SEEK_SET);
		safe_fread(&tmp_pos, L_SIZE, 1, mastertab);
		// now we're at position i+1, need to go to num_v+2
		// so we just have to move forward (num_v-i+1) positions to get to the end
		fseek(mastertab, (num_v-i)*L_SIZE, SEEK_CUR);

		// then move forward another [entry] amounts and write the current index to get a self loop
		fseek(mastertab, tmp_pos*entry_size, SEEK_CUR);
		fwrite(&i, L_SIZE, 1, mastertab);
		fwrite(&self_weight, sizeof(double), 1, mastertab);
	}

	fclose(mastertab);
}

void normalize_csr_edgecounts(const char* ftable, l_uint num_v){
	const int entry_size = L_SIZE + sizeof(double);
	double tmp_val, normalizer;
	l_uint start, end;
	FILE *mastertab = fopen(ftable, "rb+");
	if(!mastertab) error("%s", "error opening CSR file.\n");

	start = 0;
	for(l_uint i=0; i<num_v-1; i++){
		normalizer = 0;
		fseek(mastertab, (i+1)*L_SIZE, SEEK_SET);
		safe_fread(&end, L_SIZE, 1, mastertab);
		// pointer is now at position i+2, need to go to num_v+1
		// num_v+1-(i+2) = num_v-i-1
		fseek(mastertab, (num_v-i-1)*L_SIZE, SEEK_CUR);
		fseek(mastertab, start*entry_size, SEEK_CUR);
		for(l_uint j=0; j<(end-start); j++){
			fseek(mastertab, L_SIZE, SEEK_CUR);
			safe_fread(&tmp_val, sizeof(double), 1, mastertab);
			normalizer += tmp_val;
		}

		// now we're at the entry at (end), need to move back to start
		// that means moving back (end+start) spaces
		fseek(mastertab, -1*entry_size*(end-start), SEEK_CUR);

		// guard case where all weights sum to 0
		if(!normalizer) normalizer = 1;

		// finally we overwrite each of the values
		for(l_uint j=0; j<(end-start); j++){
			fseek(mastertab, L_SIZE, SEEK_CUR);
			safe_fread(&tmp_val, sizeof(double), 1, mastertab);
			tmp_val /= normalizer;
			fseek(mastertab, -1*sizeof(double), SEEK_CUR);
			fwrite(&tmp_val, sizeof(double), 1, mastertab);
		}
		start = end;
	}

	fclose(mastertab);
	return;
}

void csr_compress_edgelist_batch(const char* edgefile, const char* indexfname, const char* hashfname,
																	const char* curcountfile, const char* ftable,
																	const char sep, const char linesep, l_uint num_v, int v,
																	const int is_undirected, int has_self_loops, const int ignore_weights){
	/*
	 * This should be called after we've already read in all our files
	 * critically, ensure we're rewritten our ftable file such that it is cumulative counts and not vertex counts
	 *
	 * Error checking can be reduced because we would have caught it earlier
	 */
	const size_t entry_size = L_SIZE + sizeof(double);
	const int self_loop_inc = has_self_loops ? 1 : 0;
	l_uint inds[2], loc, offset;
	double weight;

	char *vname = malloc(MAX_NODE_NAME_SIZE);
	int stringctr=0, itermax = is_undirected ? 2 : 1;
	l_uint print_counter = 0;

	FILE *mastertable, *countstable, *edgelist, *hashfile, *indexfile;
	mastertable = fopen(ftable, "rb+");
	if(!mastertable){
		mastertable = fopen(ftable, "ab+");
		fclose(mastertable);
		mastertable = fopen(ftable, "rb+");
		if(!mastertable) error("%s", "error opening temporary counts file.\n");
	}

	countstable = fopen(curcountfile, "rb+");
	if(!countstable) errorclose_file(mastertable, NULL, "error opening master table file.\n");

	edgelist = fopen(edgefile, "rb");
	if(!edgelist) errorclose_file(countstable, mastertable, "error opening edgelist file.\n");

	hashfile = fopen(hashfname, "rb");
	indexfile = fopen(indexfname, "rb");
	if(v) Rprintf("\tReading file %s...\n", edgefile);

	char c = getc(edgelist);
	while(!feof(edgelist)){
		// read in the two vertex names
		for(int i=0; i<2; i++){
			stringctr = 0;
			memset(vname, 0, MAX_NODE_NAME_SIZE);
			while(c != sep && c != linesep){
				vname[stringctr++] = c;
				c = getc(edgelist);
			}
			// short circuit to skip a length-0 name
			if(stringctr == 0){
				i--;
				c = getc(edgelist);
				continue;
			}
			inds[i] = lookup_node_index(vname, indexfile, hashfile, num_v);

			// advance one past the separator if it isn't linesep
			// it would equal linesep if we don't have weights included
			if(c == sep)
				c = getc(edgelist);
		}
		if(!ignore_weights){
			// read in the weight
			stringctr = 0;
			memset(vname, 0, MAX_NODE_NAME_SIZE);
			c = getc(edgelist);
			while(c != linesep){
				vname[stringctr++] = c;
				c = getc(edgelist);
			}
			weight = atof(vname);
		} else {
			weight = 1.0;
			while(c != linesep)
				c = getc(edgelist);
		}

		// advance one past the separator
		c = getc(edgelist);

		// write the edge
		for(int j=0; j<itermax; j++){
			// get offset for location we'll write to in the counts file
			fseek(countstable, (inds[j])*L_SIZE, SEEK_SET);
			safe_fread(&offset, L_SIZE, 1, countstable);

			// get start location in table file we'll write to
			fseek(mastertable, (inds[j])*L_SIZE, SEEK_SET);
			safe_fread(&loc, L_SIZE, 1, mastertable);

			// write the edge
			fseek(mastertable, (num_v+1)*L_SIZE, SEEK_SET);
			fseek(mastertable, (loc+offset+self_loop_inc-1)*entry_size, SEEK_CUR);
			fwrite(&inds[(j+1)%2], L_SIZE, 1, mastertable);
			fwrite(&weight, sizeof(double), 1, mastertable);

			// decrement the counts file
			fseek(countstable, inds[j]*L_SIZE, SEEK_SET);
			offset--;
			fwrite(&offset, L_SIZE, 1, countstable);
		}
		print_counter++;
		if(!(print_counter % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%lu edges read\r", print_counter);
			else R_CheckUserInterrupt();
		}
	}

	if(v) Rprintf("\t%lu edges read\n", print_counter);
	free(vname);
	fclose(mastertable);
	fclose(countstable);
	fclose(edgelist);
	fclose(indexfile);
	fclose(hashfile);
	return;
}


l_uint update_node_cluster(l_uint ind, l_uint offset, FILE *mastertab, FILE *clusterings){
	/*
	 * Determine number of edges using the table file (next - cur)
	 * If number of edges too large, use some sort of hash to bin edges, then rerun with less
	 * Inputs:
	 * 	- 	      ind: 0-indexed vertex id
	 * 	-      offset: location where edges begin in mastertab
	 *	-   mastertab: file pointer to CSR-compressed graph
	 *	- clusterings: file of current cluster numbers (0=unassigned)
	 */
	R_CheckUserInterrupt();
	const size_t entry_size = L_SIZE + sizeof(double);
	l_uint start, end, num_edges, tmp_cl, tmp_id, zeromaxid=ind+1;
	double tmp_w, acceptance_prob, zeromax=0;

	// move to information for the vertex and read in number of edges
	fseek(mastertab, L_SIZE*ind, SEEK_SET);
	safe_fread(&start, L_SIZE, 1, mastertab);
	safe_fread(&end, L_SIZE, 1, mastertab);

	// if we're above the max edges, subsample the edges we read
	num_edges = end - start;
	// if it has no edges we can't do anything
	if(!num_edges) return ind+1;

	acceptance_prob = fmin((double)MAX_EDGES_EXACT / num_edges, 1.0);
	ll *searchpath = NULL;
	GetRNGstate();
	fseek(mastertab, L_SIZE*offset, SEEK_SET);
	fseek(mastertab, entry_size*start, SEEK_CUR);

	for(int i=0; i<num_edges; i++){
		// skip with probability
		if(use_limited_nodes && unif_rand() > acceptance_prob){
			fseek(mastertab, entry_size, SEEK_CUR);
			continue;
		}

		// read in neighbor and weight
		safe_fread(&tmp_id, L_SIZE, 1, mastertab);
		safe_fread(&tmp_w, sizeof(double), 1, mastertab);

		// get which cluster it belongs to
		fseek(clusterings, L_SIZE*tmp_id, SEEK_SET);
		safe_fread(&tmp_cl, L_SIZE, 1, clusterings);

		/*
		 * this solves a special edge case where uninitialized nodes with a self-loop
		 * can be counted incorrectly if their neighbors have already been assigned to
		 * their node. Thus, if it's a self loop and the cluster is uninitialized,
		 * set it to its own cluster (what it would be initialized to)
		 */
		if(tmp_id == ind && !tmp_cl)
			tmp_cl = ind+1;

		// store value
		if(tmp_cl){
			searchpath = insert_ll(searchpath, tmp_cl, tmp_w);
		} else if(tmp_w > zeromax) {
			zeromax = tmp_w;
			zeromaxid = tmp_id+1; // +1 because vertices 0-indexed and clusters 1-indexed
		}
	}


	// now we've read in all the edges, find the new community
	// (and free along the way)
	// max weight and cluster will be stored in tmp_w and tmp_cl (resp.)
	ll *tmp_ll = searchpath;
	tmp_w = -1;
	while(tmp_ll){
		searchpath = tmp_ll;
		tmp_ll = tmp_ll->next;

		if(searchpath->w > tmp_w || (searchpath->w == tmp_w && unif_rand() < 0.5) ){
			// ties are broken randomly
			tmp_w = searchpath->w;
			tmp_cl = searchpath->id;
		}
		free(searchpath);
	}
	PutRNGstate();
	if(zeromax > tmp_w) tmp_cl = zeromaxid;

	/*
	 * now we have the cluster to choose stored in tmp_cl
	 * Edge cases:
	 * - crazy chance results in every edge being skipped
	 *		-> since zeromax=0 and zeromaxid=id+1 at start, we just cluster it with itself and continue
	 * - all clusters are initialized, searchpath ends up NULL
	 * 		-> we won't free a NULL, and since zeromax > -1 we'll set it to the strongest uninitialized connection
	 */

	// last step is to write it to the cluster file
	fseek(clusterings, L_SIZE*ind, SEEK_SET);
	fwrite(&tmp_cl, L_SIZE, 1, clusterings);

	return tmp_cl;
}

void add_to_queue(l_uint clust, l_uint ind, l_uint n_node, FILE *clust_f, FILE *master_f, FILE *q_f, FILE *ctrq_f){
	l_uint start, end, tmp_ind, tmp_cl, nedge;
	l_uint *buf = malloc(L_SIZE*MAX_EDGES_EXACT);
	double dummy;
	int ctr = 0, found;

	fseek(master_f, ind*L_SIZE, SEEK_SET);
	safe_fread(&start, L_SIZE, 1, master_f);
	safe_fread(&end, L_SIZE, 1, master_f);
	fseek(master_f, (n_node+1)*L_SIZE, SEEK_SET);
	fseek(master_f, start*(L_SIZE+sizeof(double)), SEEK_CUR);

	nedge = end - start;
	for(l_uint i=0; i<nedge; i++){
		safe_fread(&tmp_ind, L_SIZE, 1, master_f);
		safe_fread(&dummy, sizeof(double), 1, master_f);
		fseek(clust_f, L_SIZE*tmp_ind, SEEK_SET);
		safe_fread(&tmp_cl, L_SIZE, 1, clust_f);

		if(tmp_cl && tmp_cl == clust) continue;
		tmp_ind++;
		found = 0;
		for(int j=0; j<ctr; j++){
			if(buf[j] == tmp_ind){
				found = 1;
				break;
			}
		}
		if(!found){
			buf[ctr++] = tmp_ind;
			if(ctr == MAX_EDGES_EXACT) break;
		}
	}

	// iterate over queue file, adding numbers if not already there
	rewind(q_f);
	for(int j=0; j<ctr; j++){
		fseek(ctrq_f, buf[j], SEEK_SET);
		found = getc(ctrq_f);
		if(!found){
			fseek(ctrq_f, -1, SEEK_CUR);
			putc(1, ctrq_f);
		} else {
			buf[j] = 0;
		}
	}

	// this is just in case, it adds a little runtime but it's safer to guard fread errors
	fseek(q_f, 0, SEEK_END);
	for(int j=0; j<ctr; j++){
		if(buf[j]){
			buf[j]--;
			fwrite(&buf[j], L_SIZE, 1, q_f);
		}
	}

	free(buf);
	return;
}

l_uint get_qsize(FILE *q){
	l_uint scratch, ctr=0;
	while(fread(&scratch, L_SIZE, 1, q)) ctr++;
	rewind(q);
	return ctr;
}

void initialize_queue(FILE *q, l_uint maxv, FILE *ctr_file){
	GetRNGstate();
	l_uint j, tmp;
	for(l_uint i=0; i<maxv; i++){
		putc(1, ctr_file);
		j = (l_uint) trunc((i+1) * (unif_rand()));
		if(j < i){
			// guarding edge case where unif_rand() returns 1.0

			// tmp = arr[j]
			fseek(q, L_SIZE*j, SEEK_SET);
			safe_fread(&tmp, L_SIZE, 1, q);

			// arr[j] = i
			fseek(q, -1*L_SIZE, SEEK_CUR);
			fwrite(&i, L_SIZE, 1, q);
		} else {
			tmp = i;
		}

		// arr[i] = tmp
		fseek(q, L_SIZE*i, SEEK_SET);
		fwrite(&tmp, L_SIZE, 1, q);
	}
	PutRNGstate();

	return;
}

void cluster_file(const char* mastertab_fname, const char* clust_fname,
									const char *qfile_f1, const char *qfile_f2, const char *qfile_log,
									l_uint num_v, int max_iterations, int v){
	// main runner function to cluster nodes
	FILE *masterfile = fopen(mastertab_fname, "rb");
	FILE *clusterfile = fopen(clust_fname, "rb+");
	FILE *cur_q, *next_q, *ctr_q;
	const char* queues[] = {qfile_f1, qfile_f2};
	//const char* progress = "\\|/-\\|/-";
	const char* progress[] = {"|o-----|", "|-o----|", "|--o---|", "|---o--|", "|----o-|", "|-----o|",
														"|-----o|", "|----o-|", "|---o--|", "|--o---|", "|-o----|", "|o-----|"};
	const int progbarlen = 12;
	const int progbarcharlen = 8;
	char progpreprint[progbarcharlen+1];
	memset(progpreprint, '\b', progbarcharlen);

	ctr_q = fopen(qfile_log, "wb+");
	l_uint cluster_res, qsize, tmp_ind;
	int print_counter=0, i=0;
	uint statusctr=0;

	// randomly initialize queue and ctr file
	if(v) Rprintf("\tInitializing queues...");
	cur_q = fopen(queues[0], "wb+");
	initialize_queue(cur_q, num_v, ctr_q);
	fclose(cur_q);
	if(v) Rprintf("done.\n\tClustering network:\n");

	if(v) Rprintf("\t0%% complete %s", progress[++statusctr%progbarlen]);
	while(i+1 != max_iterations){
		cur_q = fopen(queues[i%2], "rb+");
		next_q = fopen(queues[(i+1)%2], "wb+");
		qsize = get_qsize(cur_q);
		if(!qsize) break;

		while(fread(&tmp_ind, L_SIZE, 1, cur_q)){
			fseek(ctr_q, tmp_ind, SEEK_SET);
			putc(0, ctr_q);
			cluster_res = update_node_cluster(tmp_ind, num_v+1, masterfile, clusterfile);
			add_to_queue(cluster_res, tmp_ind, num_v, clusterfile, masterfile, next_q, ctr_q);
			print_counter++;
			if(!(print_counter % PROGRESS_COUNTER_MOD)){
				if(v){
					if(max_iterations > 0)
						Rprintf("\r\t%.0f%% complete %s", ((double)(i+1) / max_iterations)*100, progress[++statusctr%progbarlen]);
					else
						Rprintf("\r\t%d iterations %s", i+1, progress[++statusctr%progbarlen]);
				}
				else R_CheckUserInterrupt();
			}
		}

		fclose(cur_q);
		fclose(next_q);
		if(v){
			if(max_iterations > 0)
				Rprintf("\r\t%.0f%% complete %s", ((double)(i+1) / max_iterations)*100, progress[++statusctr%progbarlen]);
			else
				Rprintf("\r\t%d iterations %s", i+1, progress[++statusctr%progbarlen]);
		}
		i++;
	}
	if(v){
		if(max_iterations > 0)
			Rprintf("\r\t100%% complete!                \n");
		else
			Rprintf("\r\tComplete! (%d total iterations)     \n", i+1);
	}

	fclose(ctr_q);
	fclose(masterfile);
	fclose(clusterfile);

	remove(queues[0]);
	remove(queues[1]);
	remove(qfile_log);
	return;
}

void resolve_cluster_consensus(FILE *csr, const char* clustername, l_uint num_v, const double nclust){
	// remember that the csr file has num_v+1 entries, meaning they are at locations (0 -> num_v)
	const size_t entry_size = L_SIZE + sizeof(double);
	const double to_add = 1/nclust;
	FILE *clustf = fopen(clustername, "rb");

	l_uint start=0, end=0, pair, cur_clust, tmp_clust;
	double w;

	// iterate over all nodes
	for(l_uint i=0; i<num_v-1; i++){
		// rewind files to known location, read node data (end index, cluster number)
		fseek(csr, (i+1)*L_SIZE, SEEK_SET);
		safe_fread(&end, L_SIZE, 1, csr);

		fseek(clustf, i*L_SIZE, SEEK_SET);
		safe_fread(&cur_clust, L_SIZE, 1, clustf);

		// advance csr to the beginning of the edge block (see normalize_csr_edgecounts for this math)
		fseek(csr, (num_v-i-1)*L_SIZE, SEEK_CUR);
		fseek(csr, start*entry_size, SEEK_CUR);
		for(l_uint j=0; j<(end-start); j++){
			safe_fread(&pair, L_SIZE, 1, csr);
			safe_fread(&w, sizeof(double), 1, csr);
			fseek(clustf, pair*L_SIZE, SEEK_SET);
			safe_fread(&tmp_clust, L_SIZE, 1, clustf);
			if(tmp_clust == cur_clust){
				// if clusters are the same, add the increment to the value
				w += to_add;
				fseek(csr, -1*sizeof(double), SEEK_CUR);
				fwrite(&w, sizeof(double), 1, csr);
			}
		}
		start = end;
	}

	fclose(clustf);
}

void cluster_oom_single(const char* tabfile, const char* clusteroutfile, const char* dir,
												const char* qfile1, const char* qfile2, const char* qfile3,
												l_uint num_v, int num_iter, int verbose, int is_consensus){
	// runner function to cluster for a single file
	// will be called multiple times for consensus clustering
	if(verbose){
		if(is_consensus) Rprintf("\tClustering...\n");
		else Rprintf("Clustering...\n");
	}
 	cluster_file(tabfile, clusteroutfile, qfile1, qfile2, qfile3, num_v, num_iter, verbose);

 	if(!is_consensus){
	 	// reindex the clusters from 1 to n
	 	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu1, postcopy_dlu1);
	 	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu2, postcopy_dlu2);
 	}
	return;
}


void consensus_cluster_oom(const char* csrfile, const char* clusteroutfile, const char* dir,
													 const char* qfile1, const char* qfile2, const char* qfile3,
													 l_uint num_v, int num_iter, int v,
 													 const double* consensus_weights, const int consensus_len){
	char* tmpcsrfilename1 = malloc(PATH_MAX);
	char* tmpcsrfilename2 = malloc(PATH_MAX);
	char* tmpclusterfile = malloc(PATH_MAX);
	safe_filepath_cat(dir, CONSENSUS_CSRCOPY1, tmpcsrfilename1, PATH_MAX);
	safe_filepath_cat(dir, CONSENSUS_CSRCOPY2, tmpcsrfilename2, PATH_MAX);
	safe_filepath_cat(dir, CONSENSUS_CLUSTER, tmpclusterfile, PATH_MAX);

	FILE *dummyclust, *consensuscsr;
	l_uint *zeroclust = calloc(FILE_READ_CACHE_SIZE, num_v);
	int ntw;

	// tmpcsrfilename2 is going to store the final consensus weights
	copy_csrfile_sig(tmpcsrfilename2, csrfile, num_v, -1);

	consensuscsr = fopen(tmpcsrfilename2, "rb+");
	// now we run clustering over consensus_len times
	for(int i=0; i<consensus_len; i++){
		if(v) Rprintf("Iteration %d of %d:\n", i+1, consensus_len);

		// create csr copy with transformed weights
		if(v) Rprintf("\tTransforming edge weights...\n");
		copy_csrfile_sig(tmpcsrfilename1, csrfile, num_v, consensus_weights[i]);

		// set up a dummy cluster file
		dummyclust = fopen(tmpclusterfile, "wb");
		ntw = num_v;
		while(ntw > 0)
			ntw -= fwrite(zeroclust, L_SIZE, ntw > FILE_READ_CACHE_SIZE ? FILE_READ_CACHE_SIZE : ntw, dummyclust);
		fclose(dummyclust);

		// cluster into dummyclust
		cluster_oom_single(tmpcsrfilename1, tmpclusterfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, v, 1);

		if(v) Rprintf("\tRecording results...\n");
		// lastly, add edge to consensus csr file if they're the same cluster
		resolve_cluster_consensus(consensuscsr, tmpclusterfile, num_v, consensus_len);
	}
	fclose(consensuscsr);

	if(v) Rprintf("Clustering on consensus data...\n");
	cluster_oom_single(tmpcsrfilename2, clusteroutfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, v, 1);

	// reindex clusters from 1 to n
	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu1, postcopy_dlu1);
	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu2, postcopy_dlu2);

	free(tmpcsrfilename1);
	free(tmpcsrfilename2);
	free(tmpclusterfile);
	free(zeroclust);
	return;
}


SEXP R_LPOOM_cluster(SEXP FILENAME, SEXP NUM_EFILES, SEXP TABNAME, SEXP TEMPTABNAME, SEXP QFILES, SEXP OUTDIR, // files
										SEXP SEPS, SEXP CTR, SEXP ITER, SEXP VERBOSE, // control flow
										SEXP IS_UNDIRECTED, SEXP ADD_SELF_LOOPS, // optional adjustments
										SEXP IGNORE_WEIGHTS, SEXP NORMALIZE_WEIGHTS,
										SEXP CONSENSUS_WEIGHTS){
	/*
	 * I always forget how to handle R strings so I'm going to record it here
	 * R character vectors are STRSXPs, which is the same as a list (VECSXP)
	 * Each entry in a STRSXP is a CHARSXP, accessed with STRING_ELT(str, ind)
	 * You can get a `const char*` from a CHARSXP with CHAR()
	 * You can also re-encode with `const char* Rf_translateCharUTF8()`
	 */

	/*
	 * Input explanation:
	 *		 FILENAME: file of edges in format `v1 v2 w`
	 *      TABNAME: stores CSR compression of graph structure
	 *	TEMPTABNAME: first used to count edges, then used to store clusters
	 * 			 QFILES: two files, both used for queues
	 * 		   OUTDIR: directory to store hashed strings
	 *
	 * R_hashedgelist(tsv, csr, clusters, queues, hashdir, seps, 1, iter, verbose)
	 */

	// main files
	const char* dir = CHAR(STRING_ELT(OUTDIR, 0));
	const char* edgefile;
	const char* tabfile = CHAR(STRING_ELT(TABNAME, 0));
	const char* temptabfile = CHAR(STRING_ELT(TEMPTABNAME, 0));

	// queue files
	const char* qfile1 = CHAR(STRING_ELT(QFILES, 0));
	const char* qfile2 = CHAR(STRING_ELT(QFILES, 1));
	const char* qfile3 = CHAR(STRING_ELT(QFILES, 2));

	// required parameters
	const char* seps = CHAR(STRING_ELT(SEPS, 0));
	const int num_edgefiles = INTEGER(NUM_EFILES)[0];
	const int num_iter = INTEGER(ITER)[0];
	const int verbose = LOGICAL(VERBOSE)[0];
	l_uint num_v = (l_uint)(REAL(CTR)[0]);

	// optional parameters
	const int is_undirected = LOGICAL(IS_UNDIRECTED)[0];
	const double self_loop_weight = REAL(ADD_SELF_LOOPS)[0];
	const int add_self_loops = self_loop_weight > 0;
	const int should_normalize = LOGICAL(NORMALIZE_WEIGHTS)[0];
	const int ignore_weights = LOGICAL(IGNORE_WEIGHTS)[0];

	// consensus stuff
	const int consensus_len = length(CONSENSUS_WEIGHTS);
	const double* consensus_w = REAL(CONSENSUS_WEIGHTS);

	// eventually let's just make these tempfiles too called by R
	char *hashfile = malloc(PATH_MAX);
	safe_filepath_cat(dir, HASH_FNAME, hashfile, PATH_MAX);

	char *hashindex = malloc(PATH_MAX);
	safe_filepath_cat(dir, HASH_INAME, hashindex, PATH_MAX);

	// first hash all vertex names
	if(verbose) Rprintf("Building hash table for vertex names...\n");
	for(int i=0; i<num_edgefiles; i++){
		edgefile = CHAR(STRING_ELT(FILENAME, i));
		hash_file_vnames_batch(edgefile, dir, hashfile, seps[0], seps[1], verbose, is_undirected);
	}


	// next, reformat the file to get final counts for each node
	if(verbose) Rprintf("Tidying up internal tables...\n");
	num_v = node_vertex_file_cleanup(dir, hashfile, temptabfile, hashindex, verbose);

 	// next, create an indexed table file of where data for each vertex will be located
 	if(verbose) Rprintf("Reformatting vertex degree file...\n");
 	reformat_counts(temptabfile, tabfile, num_v, add_self_loops);

 	// then, we'll create the CSR compression of all our edges
 	if(verbose) Rprintf("Reading in edges...\n");
 	for(int i=0; i<num_edgefiles; i++){
 		edgefile = CHAR(STRING_ELT(FILENAME, i));
 		csr_compress_edgelist_batch(edgefile, hashindex, hashfile, temptabfile, tabfile, seps[0], seps[1],
 																num_v, verbose, is_undirected, add_self_loops, ignore_weights);
 	}

 	if(add_self_loops && verbose) Rprintf("Adding self loops...\n");
 	if(add_self_loops) add_self_loops_to_csrfile(tabfile, num_v, self_loop_weight);

 	if(verbose && should_normalize) Rprintf("Normalizing edge weights...\n");
	if(should_normalize) normalize_csr_edgecounts(tabfile, num_v);

 	if(consensus_len){
 		consensus_cluster_oom(tabfile, temptabfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, verbose,
 													consensus_w, consensus_len);

 	} else {
 		cluster_oom_single(tabfile, temptabfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, verbose, 0);
 	}

 	free(hashfile);
 	free(hashindex);
	SEXP RETVAL = PROTECT(allocVector(REALSXP, 1));
	REAL(RETVAL)[0] = (double) num_v;
	UNPROTECT(1);
	return RETVAL;
}


SEXP R_LP_write_output(SEXP CLUSTERFILE, SEXP HASHEDDIR, SEXP OUTFILE, SEXP SEPS, SEXP VERBOSE){
	/*
	 * Inputs:
	 * 	CLUSTERFILE: file of clusters (cluster_counts)
	 *  HASHEDFILES: character vector of all hash files
	 *    NUM_FILES: number of hash files
	 *      OUTFILE: file to write to
	 *         SEPS: %c%c, entry_sep, line_sep
	 *
	 * R_write_output_clusters(cluster_counts, list.files(hashdir), length(...), out_tsvpath, seps)
	 */

	const int v = LOGICAL(VERBOSE)[0];
	const char* cfile = CHAR(STRING_ELT(CLUSTERFILE, 0));
	const char* hashdir = CHAR(STRING_ELT(HASHEDDIR, 0));
	const char* outfile = CHAR(STRING_ELT(OUTFILE, 0));
	const char* seps = CHAR(STRING_ELT(SEPS, 0));
	char *hashfname = malloc(PATH_MAX);
	safe_filepath_cat(hashdir, HASH_FNAME, hashfname, PATH_MAX);

	h_uint name_len;
	char buf[MAX_NODE_NAME_SIZE];
	char write_buf[PATH_MAX];
	l_uint clust, num_written=0, junk;
	FILE *outf = fopen(outfile, "w");
	if(v) Rprintf("Writing node clusters to output file...\n");

	FILE *fclusters = fopen(cfile, "rb");
	FILE *hashfile = fopen(hashfname, "rb");
	// clusters and names are aligned, so we just need to iterate both files as once
	while(fread(&name_len, LEN_SIZE, 1, hashfile)){
		memset(buf, 0, MAX_NODE_NAME_SIZE);
		memset(write_buf, 0, PATH_MAX);
		safe_fread(buf, 1, name_len, hashfile);
		safe_fread(&junk, L_SIZE, 1, hashfile);

		// get the cluster for the node name
		safe_fread(&clust, L_SIZE, 1, fclusters);

		// prepare data for output
		snprintf(write_buf, (name_len+3)+L_SIZE, "%s%c%llu%c", buf, seps[0], clust, seps[1]);
		fwrite(write_buf, 1, strlen(write_buf), outf);
		num_written++;

		if(!(num_written % PRINT_COUNTER_MOD)){
			if(v) Rprintf("%lu nodes written.\r", num_written);
			else R_CheckUserInterrupt();
		}
	}
	if(v) Rprintf("%lu nodes written.\n", num_written);
	fclose(fclusters);
	fclose(outf);
	fclose(hashfile);
	free(hashfname);

	return R_NilValue;
}