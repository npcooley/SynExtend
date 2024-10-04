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
 *
 * Optimization Ideas:
 *  - Can we speed up reading edges more?
 * 		- add a cache of "recently used vertices" as [c1,c2,...,\0,ptr] to speed up edge reading
 *		- is it possible to order the nodes by what is most often referenced? seems like more trouble than its worth
 *	- Can we optimize RAM consumption?
 *		- might be worthwhile to compress the trie structure *after* its read in
 *	- Can we speed up clustering?
 *		- is it possible to multithread it? would have to be really careful about r/w
 *			- could write my own write lock
 * 			- fwrite is thread safe but gives no performance benefit (https://stackoverflow.com/questions/26565498/multiple-threads-writing-on-same-file)
 */

#include "SEutils.h"
#include "SynExtend.h"
#include "PrefixTrie.h"
#include <time.h>

#define h_uint uint64_t
#define uint uint_fast32_t
#define strlen_uint uint_least16_t
#define l_uint uint_fast64_t
#define lu_fprint PRIuFAST64
#define w_float float

/*
 * common limits are defined in limits.h
 * 	PAGE_SIZE: size in bytes of a page
 *  PATH_MAX: max size in bytes of a path name
 */

/*
 * TODOs:
 *	- add safeguards around safe_malloc() / calloc() to break if allocation fails
 *	- add graceful failure to cleanup files on abort / interrupt
 *		note Simon's suggestion for this:
 *			static void chkIntFn(void *dummy) {
 *			  R_CheckUserInterrupt();
 *			}
 *
 *			// this will call the above in a top-level context so it won't longjmp-out of your context
 *			bool checkInterrupt() {
 *			  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
 *			}
 *
 *			// your code somewhere ...
 *			if (checkInterrupt()) { // user interrupted ... }
 */

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#ifndef PAGE_SIZE
#define PAGE_SIZE 4096
#endif

#define MAX_NODE_NAME_SIZE 254 // max size of a vertex name (char array will have 2 extra spaces for terminator and flag)
#define NODE_NAME_CACHE_SIZE 40960 // holds pointers char* of size MAX_NODE_NAME_SIZE, 4096 is 1MB
#define FILE_READ_CACHE_SIZE 40960 // used for mergesorting files
// number of entries, so total consumption often multiplied by 8
#define CLUSTER_MIN_WEIGHT 0.01

const int L_SIZE = sizeof(l_uint);
const int W_SIZE = sizeof(w_float);
const int LEN_SIZE = sizeof(strlen_uint);
const int MAX_READ_RETRIES = 10;
const int MAX_WRITE_RETRIES = 10;
const char HASH_FNAME[] = "hashfile";
const char HASH_INAME[] = "hashindex";
const char CONS_TMPNAME[] = "tmpgraph";
const char CONSENSUS_CSRCOPY1[] = "tmpcsr1";
const char CONSENSUS_CSRCOPY2[] = "tmpcsr2";
const char CONSENSUS_CLUSTER[] = "tmpclust";

// set this to 1 if we should sample edges rather than use all of them
// MAX_EDGES_EXACT is a soft cap -- if above this, we sample edges probabalistically
const int use_limited_nodes = 0;
const l_uint MAX_EDGES_EXACT = 50000;
const int PRINT_COUNTER_MOD = 811*13;
const int PROGRESS_COUNTER_MOD = 3083;

l_uint *global_ptr;
char **allocated_filenames;
int num_alloced_filenames;

// fast, non-threadsafe getc() for better performance if opening for reading only
#ifdef WIN32
	inline int getc_unsafe(FILE *f) { return _getc_nolock(f); }
#else
	inline int getc_unsafe(FILE *f) { return getc_unlocked(f); }
#endif

static clock_t clock1, clock2, clock3, clock4;
static l_uint GLOBAL_verts_remaining;


/**********************/
/* Struct Definitions */
/**********************/
typedef struct {
	l_uint ctr1;
	l_uint ctr2;
} double_lu;

typedef struct {
	strlen_uint strlength;
	char s[MAX_NODE_NAME_SIZE];
	h_uint hash;
	l_uint count;
} msort_vertex_line;

typedef struct ll {
	l_uint id;
	float w;
	struct ll* next;
} ll;

typedef struct ll2 {
	float w;
	struct ll2* next;
} ll2;

typedef struct {
	strlen_uint len;
	h_uint hash;
	l_uint index;
} iline;

/**************************/
/* Core Utility Functions */
/**************************/

static inline void debug_print(const char* s){
	Rprintf("&D\t%s\n", s);
}

ll* insert_ll(ll* head, l_uint id, float w){
	ll *tmp = head;
	if(!tmp){
		tmp = safe_malloc(sizeof(ll));
		tmp->id = id;
		tmp->w = w;
		tmp->next=NULL;
		return(tmp);
	}
	while(tmp->id != id && tmp->next) tmp = tmp->next;
	if(tmp->id!=id){
		tmp->next = safe_malloc(sizeof(ll));
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

static size_t safe_fwrite(void *buffer, size_t size, size_t count, FILE *stream){
	size_t written_values = fwrite(buffer, size, count, stream);
	if(written_values != count){
		// same scenarios as in safe_fread
		if(feof(stream))
			error("%s", "Internal error: fread tried to read past end of file.\n");

		for(int i=0; i<MAX_WRITE_RETRIES; i++){
			// if we read a partial value, reset the counter back
			if(written_values) fseek(stream, -1*((int)written_values), SEEK_CUR);

			// try to read again
			written_values = fwrite(buffer, size, count, stream);
			if(written_values == count) return count;
		}

		// otherwise throw an error
		error("Internal error: fwrite wrote %zu values (expected %zu).\n", written_values, count);
	}
	return written_values;
}

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

static char* create_filename(const char* dir, const char* to_append){
	size_t nchar1 = strlen(to_append);
	size_t nchar2 = strlen(dir) + nchar1 + 3;
	char *newstr = safe_malloc(nchar2);
	safe_filepath_cat(dir, to_append, newstr, nchar1);
	// save the result somewhere for later
	allocated_filenames[num_alloced_filenames++] = newstr;
	return newstr;
}


void close_all_files(FILE **f, int num_files){
	for(int i=0; i<num_files; i++) fclose(f[i]);
	if(num_files) free(f);
}

void alloced_cleanup(void){
	for(int i=0; i<num_alloced_filenames; i++){
		Rprintf("%d\n", i);
		remove(allocated_filenames[i]);
		free(allocated_filenames[i]);
	}
	Rprintf("Final remove\n");
	if(allocated_filenames) free(allocated_filenames);
	Rprintf("Done.");
}

void errorclose_file(FILE **f, int num_files, const char* message){
	close_all_files(f, num_files);
	alloced_cleanup();
	error("%s", message);
}

FILE *safe_fopen(const char *fname, const char* mode, FILE **opened_files, int *num_open_files, const char *errmsg){
	FILE *f = fopen(fname, mode);
	if(!f) errorclose_file(opened_files, *num_open_files, errmsg);
	opened_files[(*num_open_files)++] = f;
	return f;
}
void graceful_error(const char *errmsg){
	alloced_cleanup();
	error("%s", errmsg);
}

inline float sigmoid_transform(float w, const double slope){
	// should probably expose these at some point
	const float scale = 0.5;
	const float cutoff = 0.2;
	float r = 1 / (1+exp(-1*slope*(w-scale)));
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

int index_compar(const void *a, const void *b){
	return global_ptr[*(int*)a] - global_ptr[*(int*)b];
}

int l_uint_index_compar(const void *a, const void *b){
	l_uint aa = global_ptr[*(l_uint*)a];
	l_uint bb = global_ptr[*(l_uint*)b];
	return (aa > bb) - (aa < bb);
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
		safe_fwrite(&dlu, sizeof(double_lu), 1, copy);
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
			safe_fwrite(&dlu, sizeof(double_lu), 1, copy);
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
			safe_fwrite(&dlu, sizeof(double_lu), 1, copy);
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
			safe_fwrite(&dlu.ctr1, L_SIZE, 1, copy);
	}
	fclose(orig);
	fclose(copy);
	return;
}

void mergesort_clust_file(const char* f, const char* dir, size_t element_size,
															int (*compar)(const void *, const void *),
															void (*precopy)(const char*, const char*),
															void (*postcopy)(const char*, const char*),
															const int verbose){
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
	 *  - *compar will provide void** values, make sure to float dereference
	 */

	// two read pointers, one write pointer
	FILE *f1_r1, *f1_r2, *f2_w;
	char *file1 = safe_malloc(PATH_MAX);
	char *file2 = safe_malloc(PATH_MAX);
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
	void **read_cache = safe_malloc(sizeof(void*) * FILE_READ_CACHE_SIZE);
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) read_cache[i] = safe_malloc(element_size);

	// copy the original file into file 1
	if(verbose) Rprintf("\t\tCopying source file...");
	precopy(f, file1);
	if(verbose) Rprintf("done.\n");

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
				safe_fwrite(read_cache[i], element_size, 1, f2_w);
			cachectr=0;
		}
		if(verbose && total_lines % PROGRESS_COUNTER_MOD == 0)
			Rprintf("\r\t\tLocal sort: %" lu_fprint " lines processed", total_lines);
	}
	if(cachectr){
		qsort(read_cache, cachectr, sizeof(void*), compar);
		for(int i=0; i<cachectr; i++)
			safe_fwrite(read_cache[i], element_size, 1, f2_w);
	}
	if(verbose) Rprintf("\r\t\tLocal sort: %" lu_fprint " lines processed\n", total_lines);

	fclose(f1_r1);
	fclose(f2_w);
	finalfile = file2;

	l_uint cur_lines = 0;
	int iter1, iter2, previt1, previt2;
	int flip = 0, cmp;
	void *tmp1 = safe_malloc(element_size);
	void *tmp2 = safe_malloc(element_size);
	char *f1, *f2;
	l_uint num_iter = 1, total_iter=0;
	if(verbose){
		l_uint tmp_bs = block_size;
		while(tmp_bs < total_lines){
			num_iter++;
			tmp_bs <<= 1;
		}
		total_iter = num_iter-1;
	}
	num_iter = 1;

	while(block_size < total_lines){
		// need an interrupt here or we can brick on larger graphs
		if(verbose){
			Rprintf("\r\t\tIteration %" lu_fprint "/%" lu_fprint ": 0.00%% complete   ", num_iter, total_iter);
		} else {
			R_CheckUserInterrupt();
		}

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
					safe_fwrite(tmp1, element_size, 1, f2_w);
					iter1--;
				} else {
					safe_fwrite(tmp2, element_size, 1, f2_w);
					iter2--;
				}
			}
			// advance pointers one block past where we just read:
			// if we move too far it doesn't really matter, we'll catch it on the next part
			fseek(f1_r1, element_size*block_size, SEEK_CUR);
			fseek(f1_r2, element_size*block_size, SEEK_CUR);
			if(verbose){
				Rprintf("\r\t\tIteration %" lu_fprint "/%" lu_fprint ": %04.2f%% complete  ",
					num_iter, total_iter, 100*(double)cur_lines / total_lines);
			} else {
				R_CheckUserInterrupt();
			}
		}

		fclose(f1_r1);
		fclose(f1_r2);
		fclose(f2_w);
		cur_lines = 0;
		block_size *= 2;
		finalfile = f2;
		num_iter++;
	}
	if(verbose && total_iter) Rprintf("\r\t\tIteration %" lu_fprint "/%" lu_fprint ": 100.00%% complete  \n",
					total_iter, total_iter);
	if(verbose) Rprintf("\t\tCopying sorted results...\n");
	// copy result back into f
	postcopy(finalfile, f);

	// free memory allocations
	free(tmp1);
	free(tmp2);
	for(int i=0; i<FILE_READ_CACHE_SIZE; i++) free(read_cache[i]);
	free(read_cache);
	remove(file1);
	remove(file2);
	free(file1);
	free(file2);
	return;
}

void copy_csrfile_sig(const char* dest, const char* src, l_uint num_v, const double w){
	l_uint *vbuf = safe_malloc(FILE_READ_CACHE_SIZE*L_SIZE);
	float *wbuf = safe_malloc(FILE_READ_CACHE_SIZE*sizeof(float));
	FILE *fd = fopen(dest, "wb");
	FILE *fs = fopen(src, "rb");

	int to_read=0, remaining=num_v+1;

	// first copy over all the vertex offsets
	while(remaining > 0){
		to_read = remaining > FILE_READ_CACHE_SIZE ? FILE_READ_CACHE_SIZE : remaining;
		remaining -= fread(vbuf, L_SIZE, to_read, fs);
		safe_fwrite(vbuf, L_SIZE, to_read, fd);
	}

	// next copy over vertices and weights
	int cachectr = 0;
	while(fread(&vbuf[cachectr], L_SIZE, 1, fs)){
		safe_fread(&wbuf[cachectr], sizeof(float), 1, fs);
		cachectr++;
		if(cachectr == FILE_READ_CACHE_SIZE){
			for(int i=0; i<cachectr; i++){
				wbuf[i] = w < 0 ? 0 : sigmoid_transform(wbuf[i], w);
				safe_fwrite(&vbuf[i], L_SIZE, 1, fd);
				safe_fwrite(&wbuf[i], sizeof(float), 1, fd);
			}
			cachectr = 0;
		}
	}
	if(cachectr){
		for(int i=0; i<cachectr; i++){
			wbuf[i] = w < 0 ? 0 : sigmoid_transform(wbuf[i], w);
			safe_fwrite(&vbuf[i], L_SIZE, 1, fd);
			safe_fwrite(&wbuf[i], sizeof(float), 1, fd);
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
			// else it's the same, so increment the corresponding count
			// (last byte stores if it's an edge that should be counted)
			counts[insert_point] += !!names[i][MAX_NODE_NAME_SIZE];
		}
	}
	insert_point++;

	// side effect writes
	*InsertPoint = insert_point;
	return;
}

h_uint hash_file_vnames_trie(const char* fname, prefix *trie, h_uint next_index,
	const char sep, const char line_sep, int v, int is_undirected){
	/*
	 * fname: .tsv list of edges
	 * dname: directory of hash codes
	 * hashfname: file to write vertices to
	 */
	FILE *f = fopen(fname, "rb");
	FILE **opened_files = malloc(sizeof(FILE*));
	opened_files[0] = f;

	// size + 1 so that there's space for marking if it's an outgoing edge
	//char *vname = safe_malloc(MAX_NODE_NAME_SIZE+1);
	char *vname;
	char **namecache = safe_malloc(sizeof(char*) * NODE_NAME_CACHE_SIZE);
	uint *str_counts = safe_malloc(sizeof(uint) * NODE_NAME_CACHE_SIZE);
	int num_unique;

	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) namecache[i] = safe_malloc(MAX_NODE_NAME_SIZE+1);

	int cur_pos = 0, cachectr=0;
	char c = getc_unsafe(f);
	l_uint print_counter = 0;

	if(v) Rprintf("\tReading file %s...\n", fname);

	while(!feof(f)){
		// going to assume we're at the beginning of a line
		// lines should be of the form `start end weight` or `start end`
		// separation is by char `sep`
		for(int iter=0; iter<2; iter++){
			vname = namecache[cachectr];
			memset(vname, 0, MAX_NODE_NAME_SIZE+1);
			cur_pos = 0;
			while(c != sep && c != line_sep){
				vname[cur_pos++] = c;
				c = getc_unsafe(f);
				if(cur_pos == MAX_NODE_NAME_SIZE-1) // max size has to include the null terminator
					errorclose_file(opened_files, 1, "Node name is larger than max allowed name size.\n");

				if(feof(f)) errorclose_file(opened_files, 1, "Unexpected end of file.\n");
			}

			// mark if edge is outgoing -- always if undirected, otherwise only if the first node name
			vname[MAX_NODE_NAME_SIZE] = is_undirected ? 1 : !iter;
			str_counts[cachectr] = vname[MAX_NODE_NAME_SIZE];

			cachectr++;
			if(cachectr == NODE_NAME_CACHE_SIZE){
				unique_strings_with_sideeffects(namecache, cachectr, &num_unique, str_counts, TRUE);
				for(int i=0; i<num_unique; i++)
					next_index = insert_into_trie(namecache[i], trie, next_index, str_counts[i]);
				cachectr = 0;
			}

			// if lines are of the form `start end`, we need to leave c on the terminator
			if(c == sep)
				c = getc_unsafe(f);
		}

		while(c != line_sep && !feof(f)) c = getc_unsafe(f);
		if(c == line_sep) c=getc_unsafe(f);
		print_counter++;
		if(!(print_counter % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%" lu_fprint " lines read\r", print_counter);
			else R_CheckUserInterrupt();
		}
	}

	if(cachectr){
		unique_strings_with_sideeffects(namecache, cachectr, &num_unique, str_counts, TRUE);
		for(int i=0; i<num_unique; i++)
				next_index = insert_into_trie(namecache[i], trie, next_index, str_counts[i]);
	}

	if(v) Rprintf("\t%" lu_fprint " lines read\n", print_counter);
	close_all_files(opened_files, 1);

	for(int i=0; i<NODE_NAME_CACHE_SIZE; i++) free(namecache[i]);
	free(namecache);
	return next_index;
}

l_uint reindex_trie_and_write_counts(prefix *trie, FILE* csrfile, l_uint cur_index){
	// Vertices are 0-indexed
	// assume that we've already opened the file, since we'll call this recursively
	if(!trie) return cur_index;
	uint8_t bits_remaining = trie->count1 + trie->count2;
	uint8_t ctr = 0;
	if(trie->bmap1 & 1){
		leaf *l = (leaf*)(trie->child_nodes[ctr++]);
		l->index = cur_index++;
		safe_fwrite(&(l->count), L_SIZE, 1, csrfile);
	}
	while(ctr < bits_remaining)
		cur_index = reindex_trie_and_write_counts(trie->child_nodes[ctr++], csrfile, cur_index);

	return cur_index;
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
		safe_fwrite(&cumul_total, l_size, 1, mtab);
		safe_fread(&curcount, l_size, 1, tmptab);
		cumul_total += curcount + self_loop; // add an extra count for each node if we add self loops
	}

	// ending position of file
	safe_fwrite(&cumul_total, l_size, 1, mtab);
	fclose(tmptab);
	fclose(mtab);
	return;
}

void add_self_loops_to_csrfile(const char *ftable, const char* fweights, const char* fneighbors, l_uint num_v, float self_weight){
	// If self loops are included, the first entry for each node remains empty
	// here we'll fill it in with the node itself
	// thus, we can just write to whatever the first index of the value is and set the value to 0
	l_uint tmp_pos, prev_pos, pos_change;

	int num_open_files = 0;
	FILE **opened_files = safe_malloc(sizeof(FILE*)*3);
	FILE *offsets = fopen(ftable, "rb+");
	if(!offsets) graceful_error("error opening CSR file.\n");
	opened_files[num_open_files++] = offsets;

	FILE *weights = fopen(fweights, "rb+");
	if(!weights) errorclose_file(opened_files, num_open_files, "error opening weights file.\n");
	opened_files[num_open_files++] = weights;

	FILE *neighbors = fopen(fneighbors, "rb+");
	if(!neighbors) errorclose_file(opened_files, num_open_files, "error opening neighbors file.\n");
	opened_files[num_open_files++] = neighbors;

	prev_pos = 0;
	for(l_uint i=0; i<num_v; i++){
		// read start of next position from offsets
		safe_fread(&tmp_pos, L_SIZE, 1, offsets);

		// get the change in position needed
		pos_change = tmp_pos - prev_pos;

		// move the other pointers
		if(pos_change){
			fseek(weights, W_SIZE*pos_change, SEEK_CUR);
			fseek(neighbors, L_SIZE*pos_change, SEEK_CUR);
		}

		// write the value
		safe_fwrite(&i, L_SIZE, 1, neighbors);
		safe_fwrite(&self_weight, W_SIZE, 1, weights);

		// store previous position (+1 because we wrote one)
		prev_pos = tmp_pos+1;
	}

	close_all_files(opened_files, num_open_files);
}

void normalize_csr_edgecounts_batch(const char* ftable, const char* fweights, l_uint num_v, int verbose){
	float normalizer;
	l_uint start, end, num_entries;
	int num_open_files = 0;
	FILE **opened_files = safe_malloc(sizeof(FILE*)*2);
	FILE *offsets = safe_fopen(ftable, "rb", opened_files, &num_open_files, "error opening offsets file.\n");
	FILE *weights = safe_fopen(fweights, "rb+", opened_files, &num_open_files, "error opening weights file.\n");

	if(verbose) Rprintf("\tNodes remaining: %" lu_fprint "", num_v);

	// switching to just reading all edges in
	w_float *weights_arr = NULL;

	safe_fread(&end, L_SIZE, 1, offsets);
	for(l_uint i=0; i<num_v; i++){
		normalizer = 0;
		start = end;
		safe_fread(&end, L_SIZE, 1, offsets);
		num_entries = end - start;
		if(!num_entries) continue;

		// read in all the weights
		// note realloc on NULL is the same as malloc
		weights_arr = safe_realloc(weights_arr, W_SIZE*num_entries);
		safe_fread(weights_arr, W_SIZE, num_entries, weights);

		fseek(weights, (-1*num_entries) * W_SIZE, SEEK_CUR);

		// normalize the weights
		for(l_uint j=0; j<num_entries; j++)
			normalizer += fabsf(weights_arr[j]);
		normalizer = normalizer == 0 ? 1 : normalizer;
		for(l_uint j=0; j<num_entries; j++)
			weights_arr[j] /= normalizer;
		safe_fwrite(weights_arr, W_SIZE, num_entries, weights);

		if(i % PROGRESS_COUNTER_MOD == 0){
			if(verbose) Rprintf("\r\tNodes remaining: %" lu_fprint "                   ", num_v-i-1);
			else R_CheckUserInterrupt();
		}
	}

	if(verbose) Rprintf("\r\tNodes remaining: Done!               \n");

	if(weights_arr) free(weights_arr);
	close_all_files(opened_files, num_open_files);
	return;
}

void inflate_csr_edgecounts(FILE *offsets, FILE *weights, FILE *q, l_uint num_v, double exponent){
	GetRNGstate();
	float tmp_val, normalizer;
	float random_nudge = CLUSTER_MIN_WEIGHT / 10;
	l_uint start, end, i, num_edges;
	rewind(q);
	float *random_nudge_arr = NULL;
	w_float *weights_arr = NULL;

	start = 0;
	while(fread(&i, L_SIZE, 1, q)){
		//Rprintf("\tnode %llu\n", i);
		normalizer = 0;
		fseek(offsets, i*L_SIZE, SEEK_SET);
		safe_fread(&start, L_SIZE, 1, offsets);
		safe_fread(&end, L_SIZE, 1, offsets);
		// going to add a small nudge to the weights to bump it out of steady states
		num_edges = end - start;
		if(!num_edges) continue;
		random_nudge_arr = safe_realloc(random_nudge_arr, sizeof(float) * num_edges);
		weights_arr = safe_realloc(weights_arr, W_SIZE * num_edges);
		// random number in [-random_nudge, random_nudge]
		for(l_uint j=0; j<(end-start); j++)
			random_nudge_arr[j] = unif_rand() * (random_nudge*2) - random_nudge;

		// read in the weights
		fseek(weights, start*W_SIZE, SEEK_SET);
		safe_fread(weights_arr, W_SIZE, num_edges, weights);

		// raise to power and get sum
		normalizer = 0;
		for(l_uint j=0; j<num_edges; j++){
			tmp_val = weights_arr[j] + random_nudge_arr[j];
			if(exponent != 1 && tmp_val > 0)
				tmp_val = pow(tmp_val, exponent);
			tmp_val = tmp_val < CLUSTER_MIN_WEIGHT ? 0 : tmp_val;
			normalizer += fabsf(tmp_val);
			weights_arr[j] = tmp_val;
		}

		// guard case where normalizer == 0
		normalizer = normalizer ? normalizer : 1;

		// normalize
		for(l_uint j=0; j<num_edges; j++)
			weights_arr[j] /= normalizer;

		fseek(weights, -1*(num_edges*W_SIZE), SEEK_CUR);
		safe_fwrite(weights_arr, W_SIZE, num_edges, weights);
	}
	PutRNGstate();

	free(weights_arr);
	free(random_nudge_arr);
	return;
}

void csr_compress_edgelist_trie_batch(const char* edgefile, prefix *trie, const char* curcountfile,
																	const char* ftable, const char* fweight, const char* fneighbor,
																	const char sep, const char linesep, l_uint num_v, int v,
																	const int is_undirected, int has_self_loops, const int ignore_weights){
	/*
	 * This should be called after we've already read in all our files
	 * critically, ensure we're rewritten our ftable file such that it is cumulative counts and not vertex counts
	 *
	 * Error checking can be reduced because we would have caught it earlier
	 *
	 * File inputs:
	 * 	ftable: offsets file
	 *  fweight: file where edge weights will be written
	 *  fneighbors: file where edge end nodes will be written
	 *  curcountfile: counts file, becomes clustering file later
	 */
	l_uint *node1, *node2, *locations;
	float *weights;
	int cachectr = 0;
	int *indexes;

	const int self_loop_inc = has_self_loops ? 1 : 0;
	l_uint inds[2], loc, offset;
	float weight;

	// additional variables for writing to cache
	l_uint prev_ind, cur_ind;

	// allocate space for all the stuff
	int stringctr=0;
	l_uint print_counter = 0;
	char *vname = safe_malloc(MAX_NODE_NAME_SIZE);
	node1 = safe_malloc(FILE_READ_CACHE_SIZE * L_SIZE);
	node2 = safe_malloc(FILE_READ_CACHE_SIZE * L_SIZE);
	locations = safe_malloc(FILE_READ_CACHE_SIZE * L_SIZE);
	weights = safe_malloc(FILE_READ_CACHE_SIZE * sizeof(float));
	indexes = safe_malloc(FILE_READ_CACHE_SIZE * sizeof(int));
	global_ptr = node1;

	FILE **opened_files = safe_malloc(sizeof(FILE*)*5);
	int num_open_files = 0;
	FILE *offsetstable, *countstable, *edgelist, *weightstable, *neighbortable;
	offsetstable = fopen(ftable, "rb+");
	if(!offsetstable){
		offsetstable = fopen(ftable, "ab+");
		fclose(offsetstable);
		offsetstable = safe_fopen(ftable, "rb+", opened_files, &num_open_files, "error opening offsets file.\n");
	}

	countstable = safe_fopen(curcountfile, "rb+", opened_files, &num_open_files, "error opening remaining counts file.\n");
	edgelist = safe_fopen(edgefile, "rb", opened_files, &num_open_files, "error opening edges file.\n");

	// initialize neighbors and weights table
	weightstable = safe_fopen(fweight, "wb", opened_files, &num_open_files, "error initializing weights file.\n");
	neighbortable = safe_fopen(fneighbor, "wb", opened_files, &num_open_files, "error initializing neighbors file\n");

	// note: weights/neighbors don't need to be initialized
	// see https://stackoverflow.com/questions/31642389/fseek-a-newly-created-file

	if(v) Rprintf("\tReading file %s...\n", edgefile);

	char c = getc_unsafe(edgelist);
	while(!feof(edgelist)){
		// read in the two vertex names
		for(int i=0; i<2; i++){
			stringctr = 0;
			//memset(vname, 0, MAX_NODE_NAME_SIZE);
			while(c != sep && c != linesep){
				vname[stringctr++] = c;
				c = getc_unsafe(edgelist);
			}
			// short circuit to skip a length-0 name
			if(stringctr == 0){
				i--;
				c = getc_unsafe(edgelist);
				continue;
			}
			vname[stringctr] = 0;
			inds[i] = find_index_for_prefix(vname, trie);

			// advance one past the separator if it isn't linesep
			// it would equal linesep if we don't have weights included
			if(c == sep)
				c = getc_unsafe(edgelist);
		}
		if(!ignore_weights){
			// read in the weight
			stringctr = 0;
			memset(vname, 0, MAX_NODE_NAME_SIZE);
			c = getc_unsafe(edgelist);
			while(c != linesep){
				vname[stringctr++] = c;
				c = getc_unsafe(edgelist);
			}
			weight = atof(vname);
		} else {
			weight = 1.0;
			while(c != linesep)
				c = getc_unsafe(edgelist);
		}

		// add to cache
		node1[cachectr] = inds[0];
		node2[cachectr] = inds[1];
		weights[cachectr] = weight;
		cachectr++;
		if(is_undirected){
			node1[cachectr] = inds[1];
			node2[cachectr] = inds[0];
			weights[cachectr] = weight;
			cachectr++;
		}

		// advance one past the separator
		c = getc_unsafe(edgelist);

		if((cachectr+is_undirected) >= FILE_READ_CACHE_SIZE || (feof(edgelist) && cachectr)){
			// get the order of the edges to write
			for(int i=0; i<cachectr; i++)
				indexes[i] = i;
			qsort(indexes, cachectr, sizeof(int), index_compar);

			// get the offset and start location for where we write to in the counts file
			prev_ind = 0;
			offset = 0;
			rewind(offsetstable);
			rewind(countstable);
			for(int i=0; i<=cachectr; i++){
				if(i < cachectr) cur_ind = node1[indexes[i]];
				if(i==0 || i==cachectr || prev_ind != cur_ind){
					// have to be sure to write the final entry as well
					// (this happens when i==cachectr)

					// get new offset and start
					if(i){
						// if not on the first entry, update counts
						// note counts should already be at the right spot
						safe_fwrite(&offset, L_SIZE, 1, countstable);
						fseek(countstable, -1*L_SIZE, SEEK_CUR);
						if(i==cachectr) break; // exit if writing final entry
					}

					// at this point the entry should be at position prev_ind (or 0)
					// this should always be positive since we already sorted the array
					fseek(countstable, (cur_ind-prev_ind)*L_SIZE, SEEK_CUR);
					safe_fread(&offset, L_SIZE, 1, countstable);
					fseek(countstable, -1*L_SIZE, SEEK_CUR);
					fseek(offsetstable, (cur_ind-prev_ind)*L_SIZE, SEEK_CUR);
					safe_fread(&loc, L_SIZE, 1, offsetstable);
					fseek(offsetstable, -1*L_SIZE, SEEK_CUR);
					prev_ind = cur_ind;
				}
				offset--; // decrement first because offsets are 1 larger than needed
				locations[i] = loc+offset+self_loop_inc;
			}

			// now we have all the offsets to write to in [locations]
			// i think now it's better to increment back to front for
			// locality and read direction on the hard drive

			// start by putting the pointer at the beginning of the edge fields
			offset = 0;
			loc = 0;
			rewind(weightstable);
			rewind(neighbortable);
			for(int i=cachectr-1; i >= 0; i--){
				// (un)signed casting issues? seems to work ok on onlinegdb.com
				if((locations[i] - loc) != 0){
					fseek(neighbortable, (locations[i]-loc)*L_SIZE, SEEK_CUR);
					fseek(weightstable, (locations[i]-loc)*W_SIZE, SEEK_CUR);
				}
				safe_fwrite(&(node2[indexes[i]]), L_SIZE, 1, neighbortable);
				safe_fwrite(&(weights[indexes[i]]), sizeof(float), 1, weightstable);
				loc = locations[i] + 1; // new location
			}
			cachectr = 0;
		}

		print_counter++;
		if(!(print_counter % PRINT_COUNTER_MOD)){
			if(v) Rprintf("\t%" lu_fprint " edges read\r", print_counter);
			else R_CheckUserInterrupt();
		}
	}

	if(v) Rprintf("\t%" lu_fprint " edges read\n", print_counter);
	free(node1);
	free(node2);
	free(locations);
	free(weights);
	free(indexes);
	free(vname);
	close_all_files(opened_files, num_open_files);
	return;
}


l_uint update_node_cluster(l_uint ind,
													FILE *offsets, FILE *clusterings,
													FILE *weightsfile, FILE *neighborfile){
	/*
	 * Determine number of edges using the table file (next - cur)
	 * If number of edges too large, use some sort of hash to bin edges, then rerun with less
	 * Inputs:
	 * 	- 	       ind: 0-indexed vertex id
	 * 	-       offset: location where edges begin in mastertab
	 *	-    	 offsets: file pointer to offsets in CSR format
	 *	-  clusterings: file of current cluster numbers (0=unassigned)
	 *	-  weightsfile: file pointer to weights for each edge
	 *  - neighborfile: file pointer to end nodes for each edge
	 */
	R_CheckUserInterrupt();
	l_uint start, end, num_edges, tmp_clust;
	w_float *weights_arr;
	l_uint *clusts;
	l_uint *indices;

	// move to information for the vertex and read in number of edges
	fseek(offsets, L_SIZE*ind, SEEK_SET);
	safe_fread(&start, L_SIZE, 1, offsets);
	safe_fread(&end, L_SIZE, 1, offsets);

	num_edges = end - start;
	// if it has no edges we can't do anything
	if(!num_edges) return ind+1;
	weights_arr = safe_malloc(W_SIZE*num_edges);
	clusts = safe_malloc(L_SIZE*num_edges);

	// read in the edges
	fseek(weightsfile, start*W_SIZE, SEEK_SET);
	fseek(neighborfile, start*L_SIZE, SEEK_SET);
	safe_fread(weights_arr, W_SIZE, num_edges, weightsfile);
	safe_fread(clusts, L_SIZE, num_edges, neighborfile);

	// read in the clusters
	rewind(clusterings);
	start = 0;
	for(l_uint i=0; i<num_edges; i++){
		end = clusts[i];
		if(end - start)
			fseek(clusterings, (end-start)*L_SIZE, SEEK_CUR);
		safe_fread(&tmp_clust, L_SIZE, 1, clusterings);
		tmp_clust = tmp_clust ? tmp_clust : clusts[i]+1;
		clusts[i] = tmp_clust;
		start = end+1;
	}

	indices = safe_malloc(L_SIZE*num_edges);
	for(l_uint i=0; i<num_edges; i++) indices[i] = i;
	global_ptr = clusts;
	qsort(indices, num_edges, L_SIZE, l_uint_index_compar);

	float max_weight=0, cur_weight=0;
	l_uint max_clust=0, cur_clust=ind+1;
	// TODO: store also the nodes that need to be updated and just update here
	// l_uint max_offset, cur_offset = 0;
	for(l_uint i=0; i<num_edges; i++){
		if(clusts[indices[i]] != cur_clust){
			if(max_weight < cur_weight){
				max_weight = cur_weight;
				max_clust = cur_clust;
				//max_offset = cur_offset;
			}
			//cur_offset = i;
			cur_clust = clusts[indices[i]];
			cur_weight = 0;
		}
		cur_weight += weights_arr[indices[i]];
	}
	if(max_weight < cur_weight){
		max_weight = cur_weight;
		max_clust = cur_clust;
	}

	// have to actually write the new cluster
	fseek(clusterings, ind*L_SIZE, SEEK_SET);
	safe_fwrite(&max_clust, L_SIZE, 1, clusterings);

	free(weights_arr);
	free(clusts);
	free(indices);
	rewind(offsets);
	rewind(clusterings);
	rewind(neighborfile);
	rewind(weightsfile);
	return max_clust;
}

void add_to_queue(l_uint clust, l_uint ind, l_uint n_node, FILE *clust_f, FILE *offset_f,
									FILE *weight_f, FILE *neighbor_f, FILE *q_f, FILE *ctrq_f){
	l_uint start, end, tmp_ind, tmp_cl, nedge;
	// TODO: make buf a minheap or something instead
	// LL would also work better for dynamic sizing
	l_uint *buf = safe_malloc(L_SIZE*MAX_EDGES_EXACT);
	w_float dummy;
	int ctr = 0, found;

	fseek(offset_f, ind*L_SIZE, SEEK_SET);
	safe_fread(&start, L_SIZE, 1, offset_f);
	safe_fread(&end, L_SIZE, 1, offset_f);

	fseek(neighbor_f, start*L_SIZE, SEEK_SET);
	fseek(weight_f, start*W_SIZE, SEEK_SET);

	nedge = end - start;
	for(l_uint i=0; i<nedge; i++){
		safe_fread(&tmp_ind, L_SIZE, 1, neighbor_f);
		safe_fread(&dummy, W_SIZE, 1, weight_f);
		fseek(clust_f, L_SIZE*tmp_ind, SEEK_SET);
		safe_fread(&tmp_cl, L_SIZE, 1, clust_f);
		if((tmp_cl && tmp_cl == clust) || dummy < CLUSTER_MIN_WEIGHT) continue;
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
		// buf[j]-1 added in solution for debug_1.R
		fseek(ctrq_f, buf[j]-1, SEEK_SET);
		found = getc(ctrq_f);
		if(!found){
			fseek(ctrq_f, -1, SEEK_CUR);
			putc(1, ctrq_f);
			GLOBAL_verts_remaining++;
		} else {
			buf[j] = 0;
		}
	}

	fseek(q_f, 0, SEEK_END);
	//safe_fwrite(buf, L_SIZE, ctr, q_f);
	for(int j=0; j<ctr; j++){
		if(buf[j]){
			buf[j]--;
			safe_fwrite(&(buf[j]), L_SIZE, 1, q_f);
		}
	}

	free(buf);
	rewind(q_f);
	return;
}

l_uint get_qsize(FILE *q){
	l_uint scratch, ctr=0;
	while(fread(&scratch, L_SIZE, 1, q)) ctr++;
	rewind(q);
	return ctr;
}


/*
l_uint get_qsize(FILE *q){
	fseek(q, 0, SEEK_END);
	return ftell(q) / L_SIZE;
}
*/

l_uint get_qsize_v(FILE *q){
	l_uint scratch, ctr=0;
	while(fread(&scratch, L_SIZE, 1, q)){
		Rprintf("%" lu_fprint " ", scratch);
	}
	rewind(q);
	Rprintf("\n");
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
			safe_fwrite(&i, L_SIZE, 1, q);
		} else {
			tmp = i;
		}

		// arr[i] = tmp
		fseek(q, L_SIZE*i, SEEK_SET);
		safe_fwrite(&tmp, L_SIZE, 1, q);
	}
	PutRNGstate();

	return;
}

void shuffle_queue(FILE *q){
	GetRNGstate();
	//l_uint maxv = get_qsize(q);
	// when we get here, queue1 should be exhausted
	// so all vertices remaining should be in queue2
	l_uint maxv = GLOBAL_verts_remaining;
	l_uint j, tmp1, tmp2;
	for(l_uint i=0; i<maxv; i++){
		j = (l_uint) trunc((i+1) * (unif_rand()));
		if(j < i){
			// guarding edge case where unif_rand() returns 1.0

			// tmp1 = arr[j]
			fseek(q, L_SIZE*j, SEEK_SET);
			safe_fread(&tmp1, L_SIZE, 1, q);

			// tmp2 = arr[i]
			fseek(q, L_SIZE*i, SEEK_SET);
			safe_fread(&tmp2, L_SIZE, 1, q);

			// arr[i] = tmp1
			fseek(q, -1*L_SIZE, SEEK_CUR);
			safe_fwrite(&tmp1, L_SIZE, 1, q);

			// arr[j] = tmp2
			fseek(q, L_SIZE*j, SEEK_SET);
			safe_fwrite(&tmp2, L_SIZE, 1, q);
		}
	}
	rewind(q);
	PutRNGstate();
	return;
}

void cluster_file(const char* offsets_fname, const char* clust_fname,
									const char* weights_fname, const char* neighbor_fname,
									const char *qfile_f1, const char *qfile_f2, const char *qfile_log,
									l_uint num_v, int max_iterations, int v, double inflation, const int should_shuffle){
	GLOBAL_verts_remaining = num_v;
	int num_open_files = 0;
	FILE **opened_files = safe_malloc(sizeof(FILE *)*7);


	// main runner function to cluster nodes
	FILE *offsetsfile = safe_fopen(offsets_fname, "rb", opened_files, &num_open_files, "error opening offsets file");
	FILE *clusterfile = safe_fopen(clust_fname, "rb+", opened_files, &num_open_files, "error opening cluster file");
	FILE *weightsfile = safe_fopen(weights_fname, "rb+", opened_files, &num_open_files, "error opening weights file");
	FILE *neighborfile = safe_fopen(neighbor_fname, "rb", opened_files, &num_open_files, "error opening neighbors file");
	FILE *cur_q, *next_q, *ctr_q;
	const char* queues[] = {qfile_f1, qfile_f2};
	//const char* progress = "\\|/-\\|/-";
	const char* progress[] = {"|o-----|", "|-o----|", "|--o---|", "|---o--|", "|----o-|", "|-----o|",
														"|-----o|", "|----o-|", "|---o--|", "|--o---|", "|-o----|", "|o-----|"};
	const int progbarlen = 12;
	const int progbarcharlen = 8;
	const float MIN_UPDATE_PCT = 0.001;
	char progpreprint[progbarcharlen+1];
	memset(progpreprint, '\b', progbarcharlen);

	ctr_q = safe_fopen(qfile_log, "wb+", opened_files, &num_open_files, "error opening queue counter file.");
	l_uint cluster_res, tmp_ind;
	int print_counter=0, i=0;
	uint statusctr=0;
	float pct_complete = 0, prev_pct=0;

	// randomly initialize queue and ctr file
	if(v) Rprintf("\tInitializing queues...");
	cur_q = safe_fopen(queues[0], "wb+", opened_files, &num_open_files, "error opening queue.");
	initialize_queue(cur_q, num_v, ctr_q);
	fclose(cur_q);
	opened_files[--num_open_files] = NULL;
	if(v) Rprintf("done.\n\tClustering network:\n");

	if(v) Rprintf("\t0%% complete %s", progress[++statusctr%progbarlen]);
	while(i+1 != max_iterations){
		cur_q = safe_fopen(queues[i%2], "rb+", opened_files, &num_open_files, "error opening queue file(s).");
		next_q = safe_fopen(queues[(i+1)%2], "wb+", opened_files, &num_open_files, "error opening queue file(s).");
		if(v){
			pct_complete = max_iterations ? (i+1) / max_iterations : ((float)(num_v - GLOBAL_verts_remaining)) / num_v;
			if(pct_complete != 1 && fabs(pct_complete - prev_pct) < MIN_UPDATE_PCT){
				pct_complete = prev_pct;
			} else {
				prev_pct = pct_complete;
			}
			Rprintf("\r\t%0.1f%% complete %s   ", (pct_complete)*100, progress[++statusctr%progbarlen]);
		}
		//if(!qsize) break;
		// if no vertices are left, break
		if(!GLOBAL_verts_remaining) break;
		if(should_shuffle) shuffle_queue(cur_q);

		/*
		if(pct_complete > 0.99){
			Rprintf("\n**DEBUG: %" lu_fprint "\n", qsize);
			get_qsize_v(cur_q);
		}
		*/
		while(fread(&tmp_ind, L_SIZE, 1, cur_q)){
			GLOBAL_verts_remaining--;
			fseek(ctr_q, tmp_ind, SEEK_SET);
			putc(0, ctr_q);
			cluster_res = update_node_cluster(tmp_ind, offsetsfile, clusterfile, weightsfile, neighborfile);
			add_to_queue(cluster_res, tmp_ind, num_v, clusterfile, offsetsfile, weightsfile, neighborfile, next_q, ctr_q);
			print_counter++;
			if(!(print_counter % PROGRESS_COUNTER_MOD)){
				if(v){
					pct_complete = max_iterations ? (i+1) / max_iterations : ((float)(num_v - GLOBAL_verts_remaining)) / num_v;
					if(pct_complete != 1 && fabs(pct_complete - prev_pct) < MIN_UPDATE_PCT){
						pct_complete = prev_pct;
					} else {
						prev_pct = pct_complete;
					}
				 	Rprintf("\r\t%0.1f%% complete %s", (pct_complete)*100, progress[++statusctr%progbarlen]);
				} else{
				 	R_CheckUserInterrupt();
				}
			}
		}

		fclose(cur_q);
		opened_files[--num_open_files] = NULL;
		i++;
		// every 2 iterations, apply inflation operator
		// case where inflation == 1.0 is handled in this function
		if(inflation > 1.0)// && i%2 == 0)
		inflate_csr_edgecounts(offsetsfile, weightsfile, next_q, num_v, pow(inflation, 1+(i/10)));
		fclose(next_q);
		opened_files[--num_open_files] = NULL;
	}
	if(v){
		if(max_iterations > 0)
			Rprintf("\r\t100%% complete!                \n");
		else
			Rprintf("\r\tComplete! (%d total iterations)     \n", i+1);
	}

	close_all_files(opened_files, num_open_files);

	remove(queues[0]);
	remove(queues[1]);
	remove(qfile_log);
	return;
}

void resolve_cluster_consensus(FILE *csr, const char* clustername, l_uint num_v, const float nclust){
	// remember that the csr file has num_v+1 entries, meaning they are at locations (0 -> num_v)
	const size_t entry_size = L_SIZE + sizeof(float);
	const float to_add = 1/nclust;
	FILE *clustf = fopen(clustername, "rb");

	l_uint start=0, end=0, pair, cur_clust, tmp_clust;
	float w;

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
			safe_fread(&w, sizeof(float), 1, csr);
			fseek(clustf, pair*L_SIZE, SEEK_SET);
			safe_fread(&tmp_clust, L_SIZE, 1, clustf);
			if(tmp_clust == cur_clust){
				// if clusters are the same, add the increment to the value
				w += to_add;
				fseek(csr, -1*sizeof(float), SEEK_CUR);
				safe_fwrite(&w, sizeof(float), 1, csr);
			}
		}
		start = end;
	}

	fclose(clustf);
}

void cluster_oom_single(const char* tabfile, const char* clusteroutfile,
												const char* weightsfile, const char* neighborfile, const char* dir,
												const char* qfile1, const char* qfile2, const char* qfile3,
												l_uint num_v, int num_iter, int verbose, int is_consensus, double inflation, const int should_shuffle){
	// runner function to cluster for a single file
	// will be called multiple times for consensus clustering
	if(verbose){
		if(is_consensus) Rprintf("\t");
		Rprintf("Clustering...\n");
	}
 	cluster_file(tabfile, clusteroutfile, weightsfile, neighborfile, qfile1, qfile2, qfile3, num_v, num_iter, verbose, inflation, should_shuffle);

 	if(!is_consensus){
 		if(verbose) Rprintf("\tReindexing clusters...\n");
	 	// reindex the clusters from 1 to n
	 	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu1, postcopy_dlu1, verbose);
	 	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu2, postcopy_dlu2, verbose);
 	}
	return;
}


void consensus_cluster_oom(const char* csrfile, const char* clusteroutfile,
													 const char* weightsfile, const char* neighborfile, const char* dir,
													 const char* qfile1, const char* qfile2, const char* qfile3,
													 l_uint num_v, int num_iter, int v, double inflation,
 													 const double* consensus_weights, const int consensus_len, const int should_shuffle){
	char* tmpcsrfilename1 = safe_malloc(PATH_MAX);
	char* tmpcsrfilename2 = safe_malloc(PATH_MAX);
	char* tmpclusterfile = safe_malloc(PATH_MAX);
	safe_filepath_cat(dir, CONSENSUS_CSRCOPY1, tmpcsrfilename1, PATH_MAX);
	safe_filepath_cat(dir, CONSENSUS_CSRCOPY2, tmpcsrfilename2, PATH_MAX);
	safe_filepath_cat(dir, CONSENSUS_CLUSTER, tmpclusterfile, PATH_MAX);

	FILE *dummyclust, *consensuscsr;
	l_uint *zeroclust = safe_calloc(FILE_READ_CACHE_SIZE, num_v);
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
			ntw -= safe_fwrite(zeroclust, L_SIZE, ntw > FILE_READ_CACHE_SIZE ? FILE_READ_CACHE_SIZE : ntw, dummyclust);
		fclose(dummyclust);

		// cluster into dummyclust
		//cluster_oom_single(tmpcsrfilename1, tmpclusterfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, v, 1, inflation, should_shuffle);

		if(v) Rprintf("\tRecording results...\n");
		// lastly, add edge to consensus csr file if they're the same cluster
		resolve_cluster_consensus(consensuscsr, tmpclusterfile, num_v, consensus_len);
	}
	fclose(consensuscsr);

	if(v) Rprintf("Clustering on consensus data...\n");
	//cluster_oom_single(tmpcsrfilename2, clusteroutfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, v, 1, inflation, should_shuffle);

	if(v) Rprintf("Reindexing clusters...\n");
	// reindex clusters from 1 to n
	if(v) Rprintf("\tSorting Iteration 1/2:\n");
	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu1, postcopy_dlu1, v);
	if(v) Rprintf("\tSorting Iteration 2/2:\n");
	mergesort_clust_file(clusteroutfile, dir, sizeof(double_lu), l_uint_compar, precopy_dlu2, postcopy_dlu2, v);

	free(tmpcsrfilename1);
	free(tmpcsrfilename2);
	free(tmpclusterfile);
	free(zeroclust);
	return;
}

l_uint write_output_clusters_trie(FILE *outfile, FILE *clusterfile, prefix *trie,
																char *s, int cur_pos, char *write_buf, const size_t num_bytes,
																const char *seps, l_uint num_v, int verbose){
	// we re-indexed the trie according to a DFS
	// thus if we just traverse the same way, we'll get the names in order
	// the hardest part here is ensuring we keep track of the character string
	if(!num_v) return 0;
	uint8_t bits_remaining = trie->count1 + trie->count2;
	uint8_t ctr = 0;
	uint8_t current_bit = 0;
	uint64_t bitmap = trie->bmap1;
	while(ctr < trie->count1){
		// iterate over first bitmap
		if(bitmap & 1){
			if(current_bit == 0){
				// read in the cluster
				l_uint tmpval = ((leaf *)(trie->child_nodes[0]))->index;
				fseek(clusterfile, tmpval*L_SIZE, SEEK_SET);
				safe_fread(&tmpval, L_SIZE, 1, clusterfile);

				// prepare data for output
				// using tmpval+1 so that clusters start at 1 and not 0
				s[cur_pos] = 0;
				snprintf(write_buf, num_bytes, "%s%c%" lu_fprint "%c", s, seps[0], tmpval+1, seps[1]);
				safe_fwrite(write_buf, 1, strlen(write_buf), outfile);
				num_v--;
				if(num_v % PROGRESS_COUNTER_MOD == 0){
					if(verbose)
						Rprintf("\r\tNodes remaining: %" lu_fprint "                    ", num_v);
					else
						R_CheckUserInterrupt();
				}
			} else {
				// set character and recur
				s[cur_pos] = current_bit + CHAR_OFFSET;
				num_v = write_output_clusters_trie(outfile, clusterfile, trie->child_nodes[ctr],
																		s, cur_pos+1, write_buf, num_bytes, seps, num_v, verbose);
			}
			ctr++;
		}
		bitmap >>= 1;
		current_bit++;
	}

	bitmap = trie->bmap2;
	current_bit = 0;
	while(ctr < bits_remaining){
		// iterate over the second bitmap
		if(bitmap&1){
			// don't have to safecheck the bit==0 case
			s[cur_pos] = current_bit + MIN_VALUE_BMAP2;
			num_v = write_output_clusters_trie(outfile, clusterfile, trie->child_nodes[ctr],
																		s, cur_pos+1, write_buf, num_bytes, seps, num_v, verbose);
			ctr++;
		}
		bitmap >>= 1;
		current_bit++;
	}

	return num_v;
}

void test_getcounts(const char* countsfile, l_uint num_v){
	FILE *f = fopen(countsfile, "rb");
	l_uint tmp, ctr=0;
	for(l_uint i=0; i<num_v; i++){
		if(ctr++ == 20){
			ctr = 0;
			Rprintf("\n");
		}
		safe_fread(&tmp, L_SIZE, 1, f);
		Rprintf("%llu ", tmp);
	}
	Rprintf("\n");
	fclose(f);
	return;
}

void test_writtenedges(const char* offsets, const char* weights, const char *neighbors, l_uint num_v){
	FILE *f1, *f2, *f3;
	f1 = fopen(offsets, "rb");
	f2 = fopen(weights, "rb");
	f3 = fopen(neighbors, "rb");
	l_uint start, end, node;
	float weight_v;

	safe_fread(&end, L_SIZE, 1, f1);
	for(l_uint i=0; i<num_v; i++){
		start = end;
		Rprintf("\t");
		safe_fread(&end, L_SIZE, 1, f1);
		for(l_uint i=start; i<end; i++){
			safe_fread(&weight_v, sizeof(float), 1, f2);
			safe_fread(&node, L_SIZE, 1, f3);
			Rprintf("%llu (%0.3f); ", node, weight_v);
		}
		Rprintf("\n");
	}

	fclose(f1);
	fclose(f2);
	fclose(f3);
	return;
}

SEXP R_LPOOM_cluster(SEXP FILENAME, SEXP NUM_EFILES, SEXP TABNAME, // files
										SEXP TEMPTABNAME, SEXP QFILES, SEXP OUTDIR, SEXP OUTFILE,	// more files
										SEXP SEPS, SEXP CTR, SEXP ITER, SEXP VERBOSE, // control flow
										SEXP IS_UNDIRECTED, SEXP ADD_SELF_LOOPS, // optional adjustments
										SEXP IGNORE_WEIGHTS, SEXP NORMALIZE_WEIGHTS,
										SEXP CONSENSUS_WEIGHTS, SEXP INFLATION_POW, SEXP SHUFFLE_QUEUES){
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

	// TODOs:
	//	- remove the following args:
	//			QFILES, TABNAME, TEMPTABNAME

	// main files
	num_alloced_filenames = 0;
	allocated_filenames = safe_malloc(sizeof(char*) * 7);
	const char* dir = CHAR(STRING_ELT(OUTDIR, 0));
	const char* tabfile = create_filename(dir, "tabfile.bin");
	const char* temptabfile = create_filename(dir, "temptabfile.bin");
	const char* weightsfile = create_filename(dir, "weights.bin");
	const char* neighborfile = create_filename(dir, "neighbors.bin");
	const char* qfile1 = create_filename(dir, "queue1.bin");
	const char* qfile2 = create_filename(dir, "queue2.bin");
	const char* qfile3 = create_filename(dir, "queue3.bin");
	const char* outfile = CHAR(STRING_ELT(OUTFILE, 0));
	const char* edgefile;

	// required parameters
	const char* seps = CHAR(STRING_ELT(SEPS, 0));
	const int num_edgefiles = INTEGER(NUM_EFILES)[0];
	const int num_iter = INTEGER(ITER)[0];
	const int verbose = LOGICAL(VERBOSE)[0];
	l_uint num_v = (l_uint)(REAL(CTR)[0]);

	// optional parameters
	const int is_undirected = LOGICAL(IS_UNDIRECTED)[0];
	const float self_loop_weight = REAL(ADD_SELF_LOOPS)[0];
	const double inflation = REAL(INFLATION_POW)[0];
	const int add_self_loops = self_loop_weight > 0;
	const int should_normalize = LOGICAL(NORMALIZE_WEIGHTS)[0];
	const int ignore_weights = LOGICAL(IGNORE_WEIGHTS)[0];
	const int should_shuffle_queue = LOGICAL(SHUFFLE_QUEUES)[0];

	// consensus stuff
	const int consensus_len = length(CONSENSUS_WEIGHTS);
	const double* consensus_w = REAL(CONSENSUS_WEIGHTS);

	// first hash all vertex names
	if(verbose) Rprintf("Building trie for vertex names...\n");
	prefix *trie = initialize_trie();
	for(int i=0; i<num_edgefiles; i++){
		edgefile = CHAR(STRING_ELT(FILENAME, i));
		num_v = hash_file_vnames_trie(edgefile, trie, num_v, seps[0], seps[1], verbose, is_undirected);
	}

	// next, reformat the file to get final counts for each node
	if(verbose) Rprintf("Tidying up internal tables...\n");
	//num_v = node_vertex_file_cleanup(dir, hashfile, temptabfile, hashindex, verbose);
	FILE *tempcounts = fopen(temptabfile, "wb");
	if(!tempcounts) error("could not open file %s\n", temptabfile);
	num_v = reindex_trie_and_write_counts(trie, tempcounts, 0);
	fclose(tempcounts);
	if(verbose) Rprintf("\tFound %" lu_fprint " unique vertices!\n", num_v);

 	// next, create an indexed table file of where data for each vertex will be located
 	if(verbose) Rprintf("Reformatting vertex degree file...\n");
 	reformat_counts(temptabfile, tabfile, num_v, add_self_loops);

 	// then, we'll create the CSR compression of all our edges
 	if(verbose) Rprintf("Reading in edges...\n");
 	for(int i=0; i<num_edgefiles; i++){
 		edgefile = CHAR(STRING_ELT(FILENAME, i));
 		csr_compress_edgelist_trie_batch(edgefile, trie,
 																temptabfile, tabfile, weightsfile, neighborfile,
 																seps[0], seps[1], num_v, verbose,
 																is_undirected, add_self_loops, ignore_weights);
 	}

 	if(add_self_loops && verbose) Rprintf("Adding self loops...\n");
 	if(add_self_loops) add_self_loops_to_csrfile(tabfile, weightsfile, neighborfile, num_v, self_loop_weight);

 	if(verbose && should_normalize) Rprintf("Normalizing edge weights...\n");
	if(should_normalize) normalize_csr_edgecounts_batch(tabfile, weightsfile, num_v, verbose);

 	if(consensus_len){
 		error("not implemented");
 		//consensus_cluster_oom(tabfile, temptabfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, verbose,
 		//											inflation, consensus_w, consensus_len, should_shuffle_queue);

 	} else {
 		cluster_oom_single(tabfile, temptabfile, weightsfile, neighborfile, dir, qfile1, qfile2, qfile3, num_v, num_iter, verbose, 0, inflation, should_shuffle_queue);
 	}


 	// have to allocate resources for writing out
 	FILE *results = fopen(outfile, "wb");
 	if(!results) graceful_error("Failed to open output file.");
 	FILE *clusters = fopen(temptabfile, "rb");
 	if(!clusters){
 		fclose(results);
 		graceful_error("Internal error while opening CSR compressed graph.");
 	}
 	char vert_name_holder[MAX_NODE_NAME_SIZE];
 	char write_buffer[PATH_MAX];
 	if(verbose) Rprintf("Writing clusters to file...\n\tVertices remaining: %" lu_fprint "", num_v);
 	write_output_clusters_trie(results, clusters, trie, vert_name_holder, 0, write_buffer, (size_t)PATH_MAX, seps, num_v, verbose);
 	if(verbose) Rprintf("\r\tNodes remaining: Done!               \n");
 	fclose(results);
 	fclose(clusters);

 	alloced_cleanup();
	free_trie(trie);
	return R_NilValue;
}
