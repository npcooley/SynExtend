#ifndef PREFIXTRIE_H
#define PREFIXTRIE_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define trie_uint uint64_t
#define FALSE 0
#define TRUE 1

typedef struct leaf {
	trie_uint count; //(making this uint32_t can save a lot of space)
	trie_uint index;
} leaf;

typedef struct prefix {
	char plen;
	char *s;
	struct prefix *next;
	struct prefix *prev;
	void *child; // can be either prefix or leaf
} prefix;

prefix *initialize_trie(void);
prefix *find_node_for_prefix(char *s, prefix *trie);
trie_uint insert_into_trie(char *s, prefix *trie, trie_uint ctr);
void free_trie(prefix *trie);

#endif