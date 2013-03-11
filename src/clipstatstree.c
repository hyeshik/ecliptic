#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "clipstatstree.h"


int
clipstats_compare(struct clipstats *a, struct clipstats *b)
{
    if (a->depth < b->depth)
        return -1;
    if (a->depth > b->depth)
        return 1;

    if (a->basetype < b->basetype)
        return -1;
    if (a->basetype > b->basetype)
        return 1;

    if (a->value < b->value)
        return -1;
    if (a->value > b->value)
        return 1;

    return 0;
}

SPLAY_GENERATE(clipstatstree, clipstats, node, clipstats_compare)


CLIPSTATS_TREES *
clipstatstrees_new(void)
{
    CLIPSTATS_TREES *trees;

    trees = malloc(sizeof(CLIPSTATS_TREES));
    if (trees != NULL) {
        SPLAY_INIT(&trees->mod);
        SPLAY_INIT(&trees->del);
        SPLAY_INIT(&trees->moddel);
        SPLAY_INIT(&trees->entropy);
        trees->bulktail = NULL;
    }

    return trees;
}

void
clipstatstrees_destroy(CLIPSTATS_TREES *trees)
{
    struct node_bulk *blockptr, *prevptr;

    for (blockptr = trees->bulktail; blockptr != NULL; blockptr = prevptr) {
        prevptr = blockptr->prev;
        free(blockptr);
    }

    free(trees);
}

static struct clipstats *
clipstatstrees_alloc_node(CLIPSTATS_TREES *trees)
{
    assert(trees->bulktail == NULL || trees->bulktail->consumed <= BULKNODE_UNIT);

    /* prepare new bulk */
    if (trees->bulktail == NULL || trees->bulktail->consumed >= BULKNODE_UNIT) {
        struct node_bulk *newbulk;
        newbulk = malloc(sizeof(struct node_bulk));
        if (newbulk == NULL)
            return NULL;

        memset(newbulk, 0, sizeof(struct node_bulk));
        newbulk->prev = trees->bulktail;
        trees->bulktail = newbulk;
    }

    return &trees->bulktail->slots[trees->bulktail->consumed++];
}

int
clipstats_add(CLIPSTATS_TREES *bundle,
              struct clipstatstree *tree, uint32_t depth,
              uint8_t basetype, float value)
{
    struct clipstats refnode, *found;

    refnode.depth = depth;
    refnode.basetype = basetype;
    refnode.value = value;

    found = SPLAY_FIND(clipstatstree, tree, &refnode);
    if (found == NULL) {
        struct clipstats *newnode;

        newnode = clipstatstrees_alloc_node(bundle);
        if (newnode == NULL)
            return -1;

        memset(newnode, 0, sizeof(struct clipstats));
        newnode->depth = depth;
        newnode->basetype = basetype;
        newnode->value = value;
        newnode->count = 1;

        SPLAY_INSERT(clipstatstree, tree, newnode);
    }
    else
        found->count++;

    return 0;
}

CLIPSTATS_ARRAY *
clipstatsarray_new(uint32_t initial_size)
{
    CLIPSTATS_ARRAY *arr;

    arr = malloc(CLIPSTATS_ARRAY_SIZE(initial_size));
    if (arr == NULL)
        return NULL;

    arr->allocated = initial_size;
    arr->end = 0;

    return arr;
}

void
clipstatsarray_destroy(CLIPSTATS_ARRAY *arr)
{
    free(arr);
}

CLIPSTATS_ARRAY *
clipstatsarray_mergetree(CLIPSTATS_ARRAY *arr, struct clipstatstree *tree)
{
    CLIPSTATS_ARRAY *merged;
    struct clipstats *node;
    struct clipstats *arrptr, *arrend;

    merged = clipstatsarray_new(arr->allocated);
    if (merged == NULL)
        return NULL;

#define EXPAND_IF_NEEDED(a) do {                                \
    if ((a)->end >= (a)->allocated) {                           \
        CLIPSTATS_ARRAY *expanded;                              \
        expanded = clipstatsarray_new((a)->allocated * 2);      \
        if (expanded == NULL)                                   \
            goto onError;                                       \
        memcpy(expanded->nodes, (a)->nodes,                     \
               sizeof(struct clipstats) * (a)->end);            \
        expanded->end = (a)->end;                               \
        free(a);                                                \
        (a) = expanded;                                         \
    }                                                           \
} while (0)

    arrptr = arr->nodes;
    arrend = &arr->nodes[arr->end];

    SPLAY_FOREACH(node, clipstatstree, tree) {
        while (arrptr < arrend && clipstats_compare(arrptr, node) < 0) {
            /* feed node from the existing array */
            EXPAND_IF_NEEDED(merged);
            memcpy(&merged->nodes[merged->end++], arrptr++,
                   sizeof(struct clipstats));
        }

        if (arrptr < arrend && clipstats_compare(arrptr, node) == 0)
            arrptr->count += node->count;
        else {
            EXPAND_IF_NEEDED(merged);
            memcpy(&merged->nodes[merged->end++], node,
                   sizeof(struct clipstats));
        }
    }

    /* flush from the existing array */
    while (arrptr < arrend) {
        EXPAND_IF_NEEDED(merged);
        memcpy(&merged->nodes[merged->end++], arrptr++,
               sizeof(struct clipstats));
    }

    return merged;

  onError:
    if (merged != NULL)
        free(merged);

    return NULL;
}
