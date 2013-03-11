#ifndef _CLIPSTATSTREE_H_
#define _CLIPSTATSTREE_H_

#include <inttypes.h>
#include "freebsd-tree.h"

#define BULKNODE_UNIT      256*1024 /* 8 MB */

struct clipstats {
    SPLAY_ENTRY(clipstats) node;
    uint64_t count;
    uint32_t depth;
    uint8_t basetype;
    float value;
};

SPLAY_HEAD(clipstatstree, clipstats);

struct node_bulk {
    struct node_bulk *prev;
    struct clipstats slots[BULKNODE_UNIT];
    int consumed;
};

typedef struct {
    struct node_bulk *bulktail;
    struct clipstatstree mod;
    struct clipstatstree del;
    struct clipstatstree moddel;
    struct clipstatstree entropy;
} CLIPSTATS_TREES;

typedef struct {
    uint32_t allocated;
    uint32_t end;
    struct clipstats nodes[1];
} CLIPSTATS_ARRAY;
#define CLIPSTATS_ARRAY_SIZE(n) \
    (sizeof(uint32_t) * 2 + sizeof(struct clipstats) * (n))

extern CLIPSTATS_TREES *clipstatstrees_new(void);
extern void clipstatstrees_destroy(CLIPSTATS_TREES *hdl);
extern int clipstats_add(CLIPSTATS_TREES *bundle, struct clipstatstree *tree,
                         uint32_t depth, uint8_t basetype, float value);
extern int clipstats_compare(struct clipstats *a, struct clipstats *b);
extern CLIPSTATS_ARRAY *clipstatsarray_new(uint32_t initial_size);
extern void clipstatsarray_destroy(CLIPSTATS_ARRAY *arr);
extern CLIPSTATS_ARRAY *clipstatsarray_mergetree(CLIPSTATS_ARRAY *arr, struct clipstatstree *tree);

SPLAY_PROTOTYPE(clipstatstree, clipstats, node, clipstats_compare)

#endif
