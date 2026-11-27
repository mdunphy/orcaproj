#ifndef CHEBCHEB_H
#define CHEBCHEB_H

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include "chebcheb_common.h"
#include "chebcheb_data.h"

typedef struct {
    int    nx;
    double scale_x, shift_x;
    double ax[CHEBCHEB_MAX_AX];
} row_cache_entry_t;

typedef struct {
    double j0;
    uint64_t filled;
    row_cache_entry_t cache[CHEBCHEB_MAX_MODELS];
} row_cache_t;

int find_covering_models(double i, double j, int out_k[2]);

void chebcheb_init_rowcache(row_cache_t *R, double j0);

double eval_cheb_tile(int k, double i, double j, row_cache_t *cache);

void eval_tile_grad(int k, double i, double j, double *fx, double *fy);

void cheb1d_eval_f_and_dfdi_and_dfdj(const cheb1d_edge_t *E,
                                     double j,
                                     double *f,
                                     double *dfdi,
                                     double *dfdj);

#endif /* CHEBCHEB_H */
