#ifndef PREFIXTRIE_H
#define PREFIXTRIE_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <R_ext/Error.h>
#include "SEutils.h"

// testing includes
//#include <stdio.h>

#define trie_uint uint64_t
#define FALSE 0
#define TRUE 1
#define CHAR_OFFSET 31
#define MIN_VALUE_BMAP2 87

typedef struct leaf {
	trie_uint count; //(making this uint32_t can save a lot of space)
	trie_uint index;
} leaf;

typedef struct prefix {
	uint64_t bmap1 : 56; // 0, 32-86
	uint8_t count1 : 8;  // counts in this bitmap
	uint64_t bmap2 : 42; // 87-127
	uint8_t count2 : 8;   // counts in this bitmap
	void **child_nodes;  // 0 will be leaf, else will be prefix
} prefix;

prefix *initialize_trie(void);
leaf *find_node_for_prefix(char *s, prefix *trie);
trie_uint insert_into_trie(char *s, prefix *trie, trie_uint ctr, trie_uint to_add);
trie_uint find_index_for_prefix(char *s, prefix *trie);
void get_trie_child_data(prefix *node, trie_uint *index, trie_uint *count);
void free_trie(prefix *trie);
#endif
