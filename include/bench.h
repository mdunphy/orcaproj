#ifndef ORCA_TRIPOLAR_BENCH_H
#define ORCA_TRIPOLAR_BENCH_H

#include <time.h>
// this is only used for benchmarking, and requires compiling with gnu99 instead of c99

static inline double elapsed_ns(struct timespec t0, struct timespec t1)
{
    return (double)(t1.tv_sec - t0.tv_sec) * 1e9 +
           (double)(t1.tv_nsec - t0.tv_nsec);
}

static inline double elapsed_us(struct timespec t0, struct timespec t1)
{
    return elapsed_ns(t0,t1)*1e-3;
}

#endif /* ORCA_TRIPOLAR_BENCH_H */
