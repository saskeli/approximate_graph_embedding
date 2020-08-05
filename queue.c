#ifndef QUEUE_INCLUDED
#define QUEUE_INCLUDED

#ifdef DEBUG
#include <assert.h>
#endif

#ifndef UINT
#define UINT unsigned int
#endif

typedef struct {
    int size;
    int start;
    int end;
    int capacity;
    UINT *vec;
} Queue;

Queue *queue_make() {
    Queue *q = malloc(sizeof(Queue));
    q->size = 0;
    q->start = 0;
    q->end = 0;
    q->vec = malloc(16 * sizeof(UINT));
    q->capacity = 16;
    return q;
}

void queue_offer(UINT v, Queue *q) {
    if (q->size == q->capacity) {
        UINT *nq = malloc((q->capacity << 1) * sizeof(UINT));
        int idx = 0;
        for (int i = q->start; i < q->capacity; i++) {
            nq[idx++] = q->vec[i];
        }
        for (int i = 0; i < q->end; i++) {
            nq[idx++] = q->vec[i];
        }
        q->start = 0;
        q->end = idx;
        q->capacity <<= 1;
        free(q->vec);
        q->vec = nq;
    }
    q->vec[q->end++] = v;
    q->end = q->end == q->capacity ? 0 : q->end;
    q->size++;
}

UINT queue_pop(Queue *q) {
#ifdef DEBUG
    assert(q->size > 0);
#endif
    q->size--;
    UINT v = q->vec[q->start++];
    if (q->start == q->capacity) q->start = 0;
    return v;
}

UINT queue_peek(Queue *q) {
    return q->vec[q->start];
}

void queue_reset(Queue *q) {
    q->start = 0;
    q->end = 0;
    q->size = 0;
}

#endif