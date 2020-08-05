#ifndef PAIRWISE_DISTANCES_INCLUDED
#define PAIRWISE_DISTANCES_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "bitvector.c"
#include "queue.c"

#ifndef UINT
#define UINT unsigned int
#endif

void hsort(UINT *arr, UINT n);
void getk(int idx, int k, int p, int n, UINT *o_neighbours, UINT *i_neighbours, 
          UINT *dists, __uint64_t *tsen, Queue *qa, Queue *qb, 
          UINT *out);

UINT **setup_distances(int n, int k, int start, int end, int p, int id, int ntasks, UINT *o_neighbours, UINT *i_neighbours) {
#ifdef DEBUG
    fprintf(stderr, "Task %d is calculating distances %d - %d\n", id, start, end);
#endif

    int e_count = (end - start) * k * 4;;
    
    UINT *locals = malloc(e_count * sizeof(UINT));
    __uint64_t *tsen = bv_make(n);
    UINT *dists = malloc(2 * n * sizeof(UINT));
    Queue *qa = queue_make();
    Queue *qb = queue_make();

    for (int i = 0; i < end - start; i++) {
        getk(i + start, k, p, n, o_neighbours, i_neighbours, dists, tsen, qa, qb, locals + (i * k * 4));
        if (id == 0)
            fprintf(stderr, "Got distances for approx %.2f%% of edges\r", ((double) i * 100) / end);

    }
    if (id == 0) fprintf(stderr, "\n");

#ifdef DEBUG
    fprintf(stderr, "Task %d, done calculating distances\n", id);
#endif

    free(tsen);
    free(dists);
    free(qa->vec);
    free(qa);
    free(qb->vec);
    free(qb);

#ifdef DEBUG
    fprintf(stderr, "Task %d has %d elements\n", id, e_count);
#endif
    int *elems = malloc(ntasks * sizeof(int));
    MPI_Allgather(&e_count, 1, MPI_INT, elems, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef DEBUG
    if (id == 0) {
        fprintf(stderr, "Elements by id:\n    ");
        for (int i = 0; i < ntasks; i++) {
            fprintf(stderr, "(%d, %d) ", i, elems[i]);
        }
        fprintf(stderr, "\n");
    }
#endif
    int off = 0;
    int *offsets = malloc(ntasks * sizeof(int));
    for (int i = 0; i < ntasks; i++) {
        offsets[i] = off;
        off += elems[i];
    }

    UINT *global = malloc(off * sizeof(UINT));
    MPI_Allgatherv(locals, e_count, MPI_UNSIGNED, global, elems, offsets, MPI_UNSIGNED, MPI_COMM_WORLD);
    free(locals);
    free(elems);
    free(offsets);
    UINT **g_pointers = malloc(n * sizeof(UINT *));
    g_pointers[0] = global;
    for (int i = 0; i < n - 1; i++) {
        g_pointers[i + 1] = global + (4 * k * i);
    }
    return g_pointers;
}

void getk(int idx, int k, int p, int n, UINT *o_neighbours, UINT *i_neighbours, 
          UINT *dists, __uint64_t *seen, Queue *qa, Queue *qb, 
          UINT *out) {
    bv_reset(seen);
    for (int i = 2; i < n * 2; i+=2) {
        dists[i] = i / 2;
        dists[i + 1] = 0;
    }
    queue_reset(qa);
    queue_reset(qb);

    queue_offer(0, qa);
    queue_offer(idx, qa);

    int max_w = 0;

    while (qa->size > 0 || qb->size > 0) {
        UINT n_id;
        UINT n_w;
        if (qb->size == 0) {
            n_w = queue_pop(qa);
            n_id = queue_pop(qa);
        } else if (qa->size == 0) {
            n_w = queue_pop(qb);
            n_id = queue_pop(qb);
        } else if (queue_peek(qa) < queue_peek(qb)) {
            n_w = queue_pop(qa);
            n_id = queue_pop(qa);
        } else {
            n_w = queue_pop(qb);
            n_id = queue_pop(qb);
        }
        if (bv_get(n_id, seen)) continue;
        bv_set(n_id, 1, seen);
        dists[n_id * 2 + 1] = n_w;
        max_w = n_w > max_w ? n_w : max_w;
        int start = o_neighbours[n_id];
        int end = start + o_neighbours[start++];
        while (start < end) {
            UINT ne_id = o_neighbours[start++];
            if (bv_get(ne_id, seen) == 0) {
                queue_offer(n_w + 1, qa);
                queue_offer(ne_id, qa);
            }
        }
        start = i_neighbours[n_id];
        end = start + i_neighbours[start++];
        while (start < end) {
            UINT ne_id = i_neighbours[start++];
            if (bv_get(ne_id, seen) == 0) {
                queue_offer(n_w + p, qb);
                queue_offer(ne_id, qb);
            }
        }
    }
    max_w++;
    for (int i = 3; i < n * 2; i += 2) {
        if (dists[i] == 0 && dists[i-1] != idx) {
            dists[i] = max_w;
        }
    }
    hsort(dists, n-1);
    memcpy(out, dists + 4, k * 2 * sizeof(UINT));
    memcpy(out + (2 * k), dists + ((n - k) * 2), k * 2 * sizeof(UINT));
}

UINT parent(UINT idx) {
    return idx >> 1;
}

UINT l_child(UINT idx) {
    return idx << 1;
}

UINT r_child(UINT idx) {
    return (idx << 1) + 1;
}

void swap(UINT *arr, UINT a, UINT b) {
    UINT i_t = arr[a * 2];
    UINT v_t = arr[a * 2 + 1];
    arr[a * 2] = arr[b * 2];
    arr[a * 2 + 1] = arr[b * 2 + 1];
    arr[b * 2] = i_t;
    arr[b * 2 + 1] = v_t;
}

void heapify(UINT *arr, UINT i, UINT n) {
    while (1) {
        UINT l = l_child(i);
        UINT r = r_child(i);
        UINT greatest = i;
        if (l <= n && arr[l * 2 + 1] > arr[greatest * 2 + 1]) {
            greatest = l;
        }
        if (r <= n && arr[r * 2 + 1] > arr[greatest * 2 + 1]) {
            greatest = r;
        }
        if (greatest != i) {
            swap(arr, i, greatest);
            i = greatest;
        } else {
            break;
        }
    }
}

void hsort(UINT *arr, UINT n) {
    for (UINT i = n / 2; i > 0; i--) {
        heapify(arr, i, n);
    }
    while (n > 1) {
        swap(arr, 1, n);
        n--;
        heapify(arr, 1, n);
    }
}

#endif