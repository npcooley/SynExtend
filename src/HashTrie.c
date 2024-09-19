#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define trie_uint uint64_t
#define FALSE 0
#define TRUE 1


// MAKE SURE WE INITIALIZE TREE WITH:
// alloc_new_terminal_node(0, NULL)

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

leaf *alloc_leaf(trie_uint index) {
	leaf *node = malloc(sizeof(leaf));
	node->count = 0;
	node->index = index;
	return node;
}

prefix *alloc_prefix() {
	prefix *p = malloc(sizeof(prefix));
	p->plen = 0;
	p->s = NULL;
	p->next = NULL;
	p->prev = NULL;
	p->child = NULL;
	return p;
}

void assign_prefix(prefix *node, char *s) {
	if(s) {
		char plen = (char) strlen(s);
		node->plen = plen;
		node->s = malloc(plen+1);
		node->s[plen] = 0;
		if(plen) memcpy(node->s, s, plen);
	} else {
		node->plen = 0;
		node->s = calloc(1,1);
	}
	return;
}

prefix *alloc_new_terminal_node(char *s) {
	prefix *node = alloc_prefix();
	prefix *child = alloc_prefix();
	child->child = alloc_leaf(0);
	node->child = child;
	assign_prefix(node, s);
	assign_prefix(child, NULL);
	return node;
}

void insert_before_prefix(prefix *node, char *s) {
	prefix *prev = node->prev;
	prefix *newnode = alloc_new_terminal_node(s);
	node->prev = newnode;
	newnode->next = node;
	newnode->prev = prev;
	if(prev) prev->next = newnode;
	return;
}

int compare_prefix(char *s1, char *s2, int maxlen) {
	int prefix_len = 0;
	while(*s1 && *s2 && *s1==*s2 && prefix_len < maxlen) {
		s1++;
		s2++;
		prefix_len++;
	}
	return prefix_len;
}

void split_prefix(prefix *node, char pos) {
	// keep the first pos chars here, send the rest below
	prefix *newchild = alloc_prefix();
	assign_prefix(newchild, &(node->s[pos]));

	if(pos != node->plen) {
		// reallocate and add null-terminating byte if size changed
		node->s = realloc(node->s, pos);
		node->s[pos] = 0;
		node->plen = pos;

		// re-sort position of prefix in array
		// can only move backwards by losing characters
		prefix *prev = node->prev;
		if(prev && prev->plen && strcmp(prev->s, node->s)) {
			// this means we have to move this node
			prev->next = node->next;
			node->next->prev = prev;
			while(prev->prev && strcmp(prev->prev->s, node->s))
				prev = prev->prev;
			// now either at beginning of array (prev has no previous value)
			// or prev->prev->s < node->s
			// either way, we're in the right spot
			node->next = prev;
			node->prev = prev->prev;
			prev->prev = node;
			if(node->prev) node->prev->next = node;
		}
	}

	// add newchild between node and node->child
	newchild->child = node->child;
	prefix *newstart = alloc_prefix();
	assign_prefix(newstart, NULL);
	newstart->child = alloc_leaf(0);
	newstart->next = newchild;
	newchild->prev = newstart;
	node->child = newstart;

	return;
}

prefix *find_node_for_prefix(char *s, prefix *trie) {
	/*
	 * Check the following cases, in order:
	 * 1. trie is empty. Just insert the string as a prefix, and return.
	 * 2. s and trie->s have length 0. Increment leaf and return.
	 * 3. s comes before trie->s. Insert s prior to current node.
	 * 4. s comes after trie->s.
	 *		a. if prefixes match, recur on suffix at child
	 *		b. if prefixes don't match, continue to next node.
	 */
	char l = (char) strlen(s);
	int cmp = 0;

	// first compare against root
	if(trie->plen == 0 && l == 0)
		return trie;

	prefix *cur = trie;
	prefix *prev = cur->prev;
	while(cur) {
		// check the prefix
		cmp = compare_prefix(s, cur->s, cur->plen);

		// some or all of prefix matches
		if(cmp) {
			/*
			 * three cases:
			 * 1. prefix matches the entire string
			 * 2. prefix matches all of prefix but part of s.
			 * 3. prefix matches all of s but part of prefix.
			 * 1,2 are the same because of how we handle length 0 strings
			 */
			// case 2,3: split the node into two nodes
			if(cmp < cur->plen) split_prefix(cur, cmp);
			return find_node_for_prefix(s+cmp, cur->child);
		}

		// string comes before current smallest prefix
		if(strcmp(s, cur->s) < 0) {
			// assign new values before current node
			insert_before_prefix(cur, s);
			return cur->prev->child;
		}

		prev = cur;
		cur = cur->next;
	}

	// if we're here, we didn't find a prefix in the entire linked list
	cur = alloc_new_terminal_node(s);
	prev->next = cur;
	cur->prev = prev;
	return cur->child;
}

trie_uint insert_into_trie(char *s, prefix *trie, trie_uint ctr) {
	prefix *node = find_node_for_prefix(s, trie);
	leaf *l = (leaf*)node->child;
	if(!l->count)
		l->index = ctr++;
	l->count++;
	return ctr;
}

void free_trie(prefix *trie) {
	if(!trie) return;
	free_trie(trie->next);
	if(trie->plen) {
		free_trie(trie->child);
	} else {
		free(trie->child);
	}
	free(trie);

	return;
}