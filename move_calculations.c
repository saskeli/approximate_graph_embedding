#ifndef MOVE_CALCULATIONS_INCLUDED
#define MOVE_CALCULATIONS_INCLUDED

#include <stdlib.h>

#ifndef UINT
#define UINT unsigned int
#endif

void maybe_swap(UINT *coords, UINT *universe, int x, int y, UINT xa, UINT ya, UINT xb, UINT yb);

int calc_cost(int n, UINT *coords, UINT **neighbours, int k) {
    int cost = 0;
    int x = coords[n * 2];
    int y = coords[n * 2 + 1];
    UINT *cn = neighbours[n];
    for (int i = 0; i < 2 * k; i+=2) {
        int n_id = cn[i];
        int delta = abs(coords[n_id * 2] - x) + abs(coords[n_id * 2 + 1] - y);
        delta -= cn[i + 1];
        cost += delta < 0 ? 0 : delta;
    }
    for (int i = k * 2; i < 4 * k; i+=2) {
        int n_id = cn[i];
        int delta = cn[i + 1];// * cn[i + 1];
        delta -= abs(coords[n_id * 2] - x) + abs(coords[n_id * 2 + 1] - y);
        delta = delta < 0 ? 0 : delta;
        cost += delta * delta;
    }
    return cost;
}

void calculate_moves(UINT start, UINT end, UINT *universe, UINT *coords, UINT *moves, UINT **neighbours, int k, int x, int y, UINT limit_movement) {
    int idx = 0;
    limit_movement = x * y / limit_movement;
    limit_movement = limit_movement < 1 ? 1 : limit_movement;
    for (UINT n = start; n < end; n++) {
        UINT cx = coords[n * 2];
        UINT nx = cx;
        UINT cy = coords[n * 2 + 1];
        UINT ny = cy;
        int c_cost = calc_cost(n, coords, neighbours, k);
        int limit = cx - limit_movement;
        limit = limit < 0 ? 0 : limit;
        for (UINT i = limit; i < cx; i++) {
            maybe_swap(coords, universe, x, y, cx, cy, i, cy);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, i, cy, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = i;
            }
        }
        limit = cx + limit_movement + 1;
        limit = limit > x ? x : limit;
        for (UINT i = cx + 1; i < limit; i++) {
            maybe_swap(coords, universe, x, y, cx, cy, i, cy);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, i, cy, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = i;
            }
        }
        limit = cy - limit_movement;
        limit = limit < 0 ? 0 : limit;
        for (UINT j = limit; j < cy; j++) {
            maybe_swap(coords, universe, x, y, cx, cy, cx, j);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, cx, j, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = cx;
                ny = j;
            }
        }
        limit = cy + limit_movement + 1;
        limit = limit > y ? y : limit;
        for (UINT j = cy + 1; j < limit; j++) {
            maybe_swap(coords, universe, x, y, cx, cy, cx, j);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, cx, j, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = cx;
                ny = j;
            }
        }
        UINT i = cx - 1;
        UINT j = cy - 1;
        for (int _i = 0; _i < limit; _i++) {
            if (i >= x || j >= y) break;
            maybe_swap(coords, universe, x, y, cx, cy, i, j);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, i, j, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = i;
                ny = j;
            }
            i--;
            j--;
        }
        i = cx + 1;
        j = cy + 1;
        for (int _i = 0; _i < limit; _i++) {
            if (i >= x || j >= y) break;
            maybe_swap(coords, universe, x, y, cx, cy, i, j);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, i, j, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = i;
                ny = j;
            }
            i++;
            j++;
        }
        i = cx + 1;
        j = cy - 1;
        for (int _i = 0; _i < limit; _i++) {
            if (i >= x || j >= y) break;
            maybe_swap(coords, universe, x, y, cx, cy, i, j);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, i, j, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = i;
                ny = j;
            }
            i++;
            j--;
        }
        i = cx - 1;
        j = cy + 1;
        for (int _i = 0; _i < limit; _i++) {
            if (i >= x || j >= y) break;
            maybe_swap(coords, universe, x, y, cx, cy, i, j);
            int n_cost = calc_cost(n, coords, neighbours, k);
            if (universe[cx * x + cy] > 0) {
                n_cost += calc_cost(universe[cx * x + cy], coords, neighbours, k);
            }
            maybe_swap(coords, universe, x, y, i, j, cx, cy);
            if (n_cost < c_cost) {
                c_cost = n_cost;
                nx = i;
                ny = j;
            }
            i--;
            j++;
        }

        moves[idx * 2] = nx;
        moves[idx * 2 + 1] = ny;
        idx++;
    }
}

void maybe_swap(UINT *coords, UINT *universe, int x, int y, UINT xa, UINT ya, UINT xb, UINT yb) {
    coords[universe[xa * x + ya] * 2] = xb;
    coords[universe[xa * x + ya] * 2 + 1] = yb;
    if (universe[xb * x + yb] > 0) {
        coords[universe[xb * x + yb] * 2] = xa;
        coords[universe[xb * x + yb] * 2 + 1] = ya;
        UINT t = universe[xb * x + yb];
        universe[xb * x + yb] = universe[xa * x + ya];
        universe[xa * x + ya] = t;
    } else {
        universe[xb * x + yb] = universe[xa * x + ya];
        universe[xa * x + ya] = 0;
    }
}

void execute_moves(UINT *moves, int n, UINT *universe, UINT *coords, UINT **neighbours, int k, int x, int y) {
    for (int i = 1; i < n; i++) {
        UINT cx = coords[i * 2];
        UINT cy = coords[i * 2 + 1];
        maybe_swap(coords, universe, x, y, cx, cy, moves[i * 2], moves[i * 2 + 1]);
    }
}

#endif