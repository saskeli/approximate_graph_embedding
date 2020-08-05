#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <mpi.h>

#include "prog_io.c"
#include "pairwise_distances.c"
#include "parse_args.c"
#include "move_calculations.c"

#ifndef UINT
#define UINT unsigned int
#endif

void seed(int x, int y, int n, UINT *universe, UINT *coords);

int main(int argc, char *argv[]) {
    //basic setup
    int id, rc, ntasks;

    rc = MPI_Init(&argc, &argv);
    rc |= MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    rc |= MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "MPI initialization failed\n");
        exit(1);
    }

    int n, k, x, y, p, interval, iters, report_interval;
    FILE *fd = NULL;
    UINT *o_neighbours;
    UINT *i_neighbours;

    if (id == 0) {
        struct arguments *args = parse_arguments(argc, argv, id);
        n = args->n + 1;
        k = args->k;
        x = args->x;
        y = args->y;
        p = args->p;
        interval = args->interval;
        iters = args->iters;

        if (k == 0) {
            while (n) {
                k++;
                n >>= 1;
            }
            n = args->n + 1;
#ifdef DEBUG
            fprintf(stderr, "Calculated value for k: %d\n", k);
#endif
        }

        if (x == 0 && y == 0) {
            x = 2;
            while (x * x <= n * 2) {
                x += 1;
            }
            y = x;
#ifdef DEBUG
            fprintf(stderr, "Calculated values for (x, y): (%d, %d)\n", x, y);
#endif
        } else if (x == 0) {
            x = 2;
            while (x * y <= n * 2) {
                x += 1;
            }
#ifdef DEBUG
            fprintf(stderr, "Calculated value for x. (x, y): (%d, %d)\n", x, y);
#endif
        } else if (y == 0) {
            y = 2;
            while (x * y <= n * 2) {
                y += 1;
            }
#ifdef DEBUG
            fprintf(stderr, "Calculated value for y. (x, y): (%d, %d)\n", x, y);
#endif
        }

        if (p == 0) p = 3;
    
        fd = fopen(args->file_path, "r");
        if (fd == NULL) {
            fprintf(stderr, "Error reading file\n");
            return 1;
        }
#ifdef DEBUG
        fprintf(stderr, "Opened input file.\n");
#endif
        //reading edges from file
        UINT **edges = read_edges(fd, n);
        o_neighbours = edges[0];
        i_neighbours = edges[1];
        free(edges);
        fclose(fd);
#ifdef DEBUG
        if (o_neighbours[0] < 101) {
            fprintf(stderr, "Out edges:\n");
            for (UINT node = 1; node < n; node++) {
                fprintf(stderr, "node %d out:\n", node);
                UINT nloc = o_neighbours[node];
                for (UINT e = 1; e <= o_neighbours[nloc]; e++) {
                    fprintf(stderr, "%3d", o_neighbours[nloc + e]);
                }
                fprintf(stderr, "\n");
            }
            fprintf(stderr, "In edges:\n");
            for (UINT node = 1; node < n; node++) {
                fprintf(stderr, "node %d in:\n", node);
                UINT nloc = i_neighbours[node];
                for (UINT e = 1; e <= i_neighbours[nloc]; e++) {
                    fprintf(stderr, "%3d", i_neighbours[nloc + e]);
                }
                fprintf(stderr, "\n");
            }
        }
#endif
        free(args);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&y, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    //calculating network partitioning.
    int partition = n / ntasks;
    int overflow = n - partition * ntasks;
    int start = id * partition + (overflow > id ? id : overflow);
    if (start == 0) start++;
    int end = (id + 1) * partition + (overflow > id + 1 ? id + 1 : overflow);

    //calculating the required distances in the graph
    UINT ne_size;
    if (id == 0) {
        ne_size = o_neighbours[0];
    }
#ifdef DEBUG
    if (id == 0) fprintf(stderr, "Broadcasting edges\n");
#endif
    MPI_Bcast(&ne_size, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    if (id != 0) {
        o_neighbours = malloc(ne_size * sizeof(UINT));
        i_neighbours = malloc(ne_size * sizeof(UINT));
    }
    MPI_Bcast(o_neighbours, ne_size, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(i_neighbours, ne_size, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    if (id == 0) fprintf(stderr, "Calculating pairwise distances\n");
    UINT **neighbours = setup_distances(n, k, start, end, p, id, ntasks, o_neighbours, i_neighbours);
    free(o_neighbours);
    free(i_neighbours);
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) {
        fprintf(stderr, "Calculated pairwise distances\n");
        if (n <= 20) {
            for (int i = 1; i < n; i++) {
                fprintf(stderr, "%d + %d neighbours for node %d\n\t", k, k, i);
                for (int j = 0; j < 2 * k; j++) {
                    fprintf(stderr, "(%d, %d) ", neighbours[i][j * 2], neighbours[i][j * 2 + 1]);
                }
                fprintf(stderr, "\n");
            }
        }
    }
#endif

    UINT *universe;
    UINT *coords = malloc(n * 2 * sizeof(UINT));
    if (id == 0) {
        universe = calloc(x * y, sizeof(UINT));
        seed(x, y, n, universe, coords);
    } else {
        universe = malloc(x * y * sizeof(UINT));
    }
    
    UINT *all_moves;
    if (id == 0) all_moves = malloc(n * 2 * sizeof(UINT));
    UINT *task_moves = malloc((end - start) * 2 * sizeof(UINT));
    UINT elems = (end - start) * 2;
    UINT *all_elems;

    if (id == 0) all_elems = malloc(ntasks * sizeof(UINT));
    MPI_Gather(&elems, 1, MPI_UNSIGNED, all_elems, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    int *offsets;
    if (id == 0) {
        UINT off = 2;
        offsets = malloc(ntasks * sizeof(UINT));
        for (int i = 0; i < ntasks; i++) {
            offsets[i] = off;
            off += all_elems[i];
        }
    }

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr, "%d is ready for iterating\n", id);
#endif
    if (id == 0)
        report_interval = interval > 0 ? interval : 60;
    for (int i = 0; i < iters; i++) {
        MPI_Bcast(universe, x * y, MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Bcast(coords, n * 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        calculate_moves(start, end, universe, coords, task_moves, neighbours, k, x, y, i + 1);

        MPI_Gatherv(task_moves, elems, MPI_UNSIGNED, all_moves, all_elems, offsets, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        if (id == 0) {
            execute_moves(all_moves, n, universe, coords, neighbours, k, x, y);
#ifndef DEBUG
            if (i % report_interval == 0) {
#endif
                int cost = 0;
                for (int i = 1; i < n; i++) {
                    cost += calc_cost(i, coords, neighbours, k);
                }
                fprintf(stderr, "iter %d cost: %d.\n", i, cost);
#ifndef DEBUG
            }

            if (interval > 0 && i % interval == 0) {
#endif
                printf("Iter %d:\n", i);
                output_universe(x, y, universe, n);
#ifndef DEBUG
            }
#endif
        }
    }
#ifndef DEBUG
    if (id == 0) {
        if (interval > 0 && (iters - 1) % interval != 0) {
            if (interval > 0) printf("Final universe:\n");
            output_universe(x, y, universe, n);
        }
    }
#endif
#ifdef DEBUG
    fprintf(stderr, "%d is done iterating.\n", id);
#endif
    free(universe);
    free(neighbours[0]);
    free(neighbours);
    free(coords);
    if (id == 0) {
        free(all_elems);
        free(all_moves);
        free(offsets);
    }
    free(task_moves);

    MPI_Finalize();
    return 0;
}

void seed(int x, int y, int n, UINT *universe, UINT *coords) {
#ifdef DEBUG
    fprintf(stderr, "Seeding universe of size %d x %d, with %d nodes.\n", x, y, n-1);
#endif
    int u = x * y;
    UINT *ptr = universe;
    for (UINT i = 1; i < n; i++) {
        universe[i] = i;
    }
    srand(time(NULL));
#ifdef DEBUG
    srand(1337);
#endif
    for (UINT i = 0; i < u; i++) {
        UINT target = i + (rand() % (u - i));
        UINT t = universe[target];
        universe[target] = universe[i];
        universe[i] = t;
        if (t > 0) {
            coords[t * 2] = i / x;
            coords[t * 2 + 1] = i % y;
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Created universe.\n");
    if (u < 100) {
        for (UINT i = 0; i < x; i++) {
            for (UINT j = 0; j < y; j++) {
                fprintf(stderr, "%3d", universe[i * x + j]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "node coordinates:\n");
        for (UINT i = 1; i < n; i++) {
            fprintf(stderr, "%3d: (%2d, %2d)\n", i, coords[i * 2], coords[i * 2 + 1]);
        }
    }
#endif
}