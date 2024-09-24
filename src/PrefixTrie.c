#include "PrefixTrie.h"

leaf *alloc_leaf() {
	leaf *node = malloc(sizeof(leaf));
	node->count = 0;
	node->index = 0;
	return node;
}

prefix *alloc_prefix() {
	prefix *p = malloc(sizeof(prefix));
	p->bmap1 = 0;
	p->bmap2 = 0;
	p->count1 = 0;
	p->count2 = 0;
	p->child_nodes = NULL;
	return p;
}

prefix *initialize_trie(){
	return alloc_prefix();
}

void *insert_into_node_terminal(prefix *node){
	if(!(node->bmap1 & 1)){
		// create a leaf if this node wasn't previously terminal
		uint8_t total_children = node->count1 + node->count2;
		node->bmap1 |= 1;
		node->count1++;
		leaf *newleaf = alloc_leaf();
		void **new_ptr_holder = malloc(sizeof(void*)*(total_children+1));
		new_ptr_holder[0] = newleaf;
		if(total_children){
			memcpy(&(new_ptr_holder[1]), node->child_nodes, sizeof(void*)*total_children);
			free(node->child_nodes);
		}
		node->child_nodes = new_ptr_holder;
	}
	return node->child_nodes[0];
}

void *insert_into_node_nonterminal(prefix *node, char s){
	uint8_t USE_HIGHER = s >= MIN_VALUE_BMAP2;
	// idx is either s-31 or s-86 depending on if we're using
	// the higher or lower bitfield
	uint8_t idx = s-(!USE_HIGHER*CHAR_OFFSET)-(USE_HIGHER*MIN_VALUE_BMAP2);
	uint8_t current_bit = 0;
	uint8_t insert_point = 0;
	uint64_t bitfield = node->bmap1;

	if(USE_HIGHER){
		insert_point = node->count1;
		bitfield = node->bmap2;
	}


	while(current_bit != idx){
		current_bit++;
		insert_point += bitfield & 1;
		bitfield >>= 1;
	}

	if(bitfield & 1){
		return node->child_nodes[insert_point];
	} else {
		uint8_t total_children = node->count1 + node->count2;
		// got to current bit BUT it doesn't exist

		// first create a new holder with an extra space
		void **new_ptr_holder = malloc(sizeof(void*)*(total_children + 1));

		// then copy first <insert_point> pointers
		if(insert_point)
			memcpy(new_ptr_holder, node->child_nodes, sizeof(void*)*insert_point);

		// then add in the new pointer
		prefix *new_child = alloc_prefix();
		new_ptr_holder[insert_point] = new_child;

		// then copy the remaining pointers
		total_children -= insert_point;
		if(total_children)
			memcpy(&(new_ptr_holder[insert_point+1]),
							&(node->child_nodes[insert_point]),
							sizeof(void*)*total_children);

		// reassign the child pointer array
		if(node->child_nodes)
			free(node->child_nodes);
		node->child_nodes = new_ptr_holder;

		// trying to avoid weirdness from potential auto-cast to 32bit
		bitfield = 1;
		bitfield <<= idx;
		if(USE_HIGHER){
			node->count2++;
			node->bmap2 |= bitfield;
		} else {
			node->count1++;
			node->bmap1 |= bitfield;
		}
		return new_child;
	}
}


leaf *find_node_for_prefix(char *s, prefix *trie){
	prefix *tmp = trie;
	while(*s){
		if(*s < 31){
			free_trie(trie);
			error("Labels must contain ASCII values in the range 32-127 (received %u)", (uint8_t)(*s));
		}
		tmp = (prefix *)insert_into_node_nonterminal(tmp, *s++);
	}
	return (leaf *)insert_into_node_terminal(tmp);
}

trie_uint insert_into_trie(char *s, prefix *trie, trie_uint ctr, trie_uint to_add) {
	leaf *node = find_node_for_prefix(s, trie);
	if(!node->count)
		node->index = ctr++;
	node->count += to_add;
	return ctr;
}

trie_uint find_index_for_prefix(char *s, prefix *trie){
	leaf *l = find_node_for_prefix(s, trie);
	return l->index;
}

void free_trie(prefix *trie) {
	if(!trie) return;

	uint8_t ctr = 0;
	uint8_t bits_remaining = trie->count1 + trie->count2;

	if(trie->bmap1 & 1)
		free((leaf*)(trie->child_nodes[ctr++]));

	while(ctr < bits_remaining)
		free_trie(trie->child_nodes[ctr++]);
	if(trie->child_nodes)
		free(trie->child_nodes);
	free(trie);
	return;
}

/*
int main(int argc, char *argv[]){
	int NCHAR = 255;
	char *buffer = calloc(NCHAR, 1);
	char *trash;
	char sep1 = '\t';
	char sep2 = '\n';

	if(argc<3) perror("No data file provided!");
	FILE *f = fopen(argv[1], "rb");
	int LINES_TO_READ = strtol(argv[2], &trash, 10);

	int mod_print = 3113;
	int ctr = 0;
	int pos;
	char c;
	trie_uint idx = 1;
	trie_uint p_idx = 1;
	int s = 0;
	prefix *trie = alloc_prefix();
	while(ctr != LINES_TO_READ){
		pos = 0;
		memset(buffer, 0, NCHAR);
		c = getc(f);
		while(!feof(f) && c != sep1){
			buffer[pos++] = c;
			c = getc(f);
		}
		if(s){
			while(!feof(f) && c != sep2)
				c = getc(f);
			s = 0;
		} else {
			s = 1;
		}
		if(feof(f)) break;
		p_idx = idx;
		idx = insert_into_trie(buffer, trie, idx, 1);
		ctr++;
		if(ctr % mod_print == 0)
			printf("\rInserted %d nodes", ctr);
	}
	printf("\rInserted %d nodes\nDone inserting! %llu unique nodes.\n", ctr, idx-1);

	free_trie(trie);
	free(buffer);
	fclose(f);
	return 0;
}
*/