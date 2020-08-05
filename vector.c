#ifndef VECTOR_INCLUDED
#define VECTOR_INCLUDED

#ifndef UINT
#define UINT unsigned int
#endif

typedef struct {
    UINT size;
    UINT capacity;
    UINT *vec;
} Vector;

Vector *vector_setup(int n) {
    Vector *vecs = malloc(n * sizeof(Vector));
    for (int i = 1; i < n; i++) {
        vecs[i].capacity = 4;
        vecs[i].size = 0;
        vecs[i].vec = malloc(4 * sizeof(UINT));
    }
    return vecs;
}

void vector_add(UINT e, Vector *v) {
    if (v->capacity == v->size) {
        v->vec = realloc(v->vec, v->capacity * 2 * sizeof(UINT));
        v->capacity *= 2;
    }
    v->vec[v->size++] = e;
}

#endif