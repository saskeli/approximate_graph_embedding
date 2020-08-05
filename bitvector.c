#ifndef BITVECTOR_INCLUDED
#define BITVECTOR_INCLUDED

#include <stdlib.h>
#include <stdio.h>
#ifdef DEBUG
#include <assert.h>
#endif
#include <string.h>

__uint64_t *bv_make(int size) {
#ifdef DEBUG
    assert(size >= 0);
#endif
    int vlen = 2 + (size >> 6);
    __uint64_t *vec = calloc(vlen, sizeof(__uint64_t));
    vec[0] = vlen - 1;
    return vec;
}

int bv_get(int index, __uint64_t *vec) {
#ifdef DEBUG
    assert(index >= 0);
#endif
    int v_idx = 1 + (index >> 6);
    int v_off = index & 63;
#ifdef DEBUG
    assert(v_idx == 1 + (index / 64));
    assert(v_off == index % 64);
#endif
    return (vec[v_idx] >> v_off) & 1;
}

int bv_set(int index, int v, __uint64_t *vec) {
    __uint64_t mask = 1;
#ifdef DEBUG
    assert(index >= 0);
#endif
    int v_idx = 1 + (index >> 6);
    mask <<= (index & 63);
#ifdef DEBUG
    assert(v_idx == 1 + (index / 64));
    assert(mask == ((__uint64_t) 1) << (index % 64));
#endif
    if (v == 0) {
        vec[v_idx] &= ~mask;
    } else {
        vec[v_idx] |= mask;
    }
    return index;
}

void bv_reset(__uint64_t *vec) {
    memset(vec + 1, 0, vec[0] * sizeof(__uint64_t));
}

#endif