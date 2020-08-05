#ifndef FILEIO_INCLUDED
#define FILEIO_INCLUDED

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "vector.c"

#ifndef UINT
#define UINT unsigned int
#endif

UINT **read_edges(FILE *fp, int n) {
    Vector *vo_neighbours = vector_setup(n);
    Vector *vi_neighbours = vector_setup(n);
    UINT edge_count = n * 2;
    int a, b;
    while (fscanf(fp, "%d\t%d", &a, &b) != EOF) {
        assert(a < n);
        assert(b < n);
        vector_add(b, &vo_neighbours[a]);
        vector_add(a, &vi_neighbours[b]);
        edge_count++;
    }
#ifdef DEBUG
    fprintf(stderr, "Read %d edges from input file.\n", edge_count - n * 2);
    fprintf(stderr, "Allocating %ld * 2 bytes for edges\n", edge_count * sizeof(UINT));
#endif
    UINT *out_edges = malloc(edge_count * sizeof(UINT));
#ifdef DEBUG
    assert(out_edges != NULL);
#endif
    UINT *in_edges = malloc(edge_count * sizeof(UINT));
#ifdef DEBUG
    assert(in_edges != NULL);
#endif

    UINT **edges = malloc(2 * sizeof(UINT *));
#ifdef DEBUG
    assert(edges != NULL);
#endif
    edges[0] = out_edges;
    edges[1] = in_edges;
    out_edges[0] = edge_count;
    in_edges[0] = edge_count;
    UINT o_idx = n;
    UINT i_idx = n;
    for (int node = 1; node < n; node++) {
        out_edges[node] = o_idx;
        in_edges[node] = i_idx;
        out_edges[o_idx++] = vo_neighbours[node].size;
        in_edges[i_idx++] = vi_neighbours[node].size;
        for (int e = 0; e < vo_neighbours[node].size; e++) {
            out_edges[o_idx++] = vo_neighbours[node].vec[e];
        }
        for (int e = 0; e < vi_neighbours[node].size; e++) {
            in_edges[i_idx++] = vi_neighbours[node].vec[e];
        }
    }
    
    for (int i = 1; i < n; i++) {
        free(vo_neighbours[i].vec);
        free(vi_neighbours[i].vec);
    }
    free(vi_neighbours);
    free(vo_neighbours);
    return edges;
}

void output_universe(int x, int y, UINT *universe, double n) {
    int pad = 2;
    while (n >= 10) {
        pad++;
        n /= 10;
    }
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            printf("%*d", pad, universe[x * i + j]);
        }
        printf("\n");
    }
}

#endif