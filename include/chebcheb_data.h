#ifndef CHEBCHEB_DATA_H
#define CHEBCHEB_DATA_H
#include <stdint.h>
#include "chebcheb_common.h"

typedef struct {
    int Py, Px;          /* degrees in (y,x) */
    double jmin, jmax;   /* physical y-range */
    double imin, imax;   /* physical x-range */
    uint64_t off_C;      /* offset to coefficients */
    uint64_t cnt_C;      /* number of doubles in C */
} chebcheb_model_t;

typedef struct {
    double jmin, jmax;   /* physical y-range */
    int k0, k1;          /* index of first and last patch on this strip */
} chebcheb_strip_t;

typedef struct {
    int tile_index;    /* index into chebcheb models */
    int x0;            /* -1 (eval left) or +1 (eval right) */
    int Py;            /* polynomial degree in y (Py ≥ 0) */
    double jmin, jmax; /* physical y-range of this tile */

    /* Chebyshev coefficients, length Py+1 */
    double f[CHEBCHEB_MAX_AX];     /*    f(x0, y) series, T type along y */
    double dfdi[CHEBCHEB_MAX_AX];  /* df/di(x0,y) coeffs, T type along y */
    double dfdj[CHEBCHEB_MAX_AX];  /* df/dj(x0,y) coeffs, U type along y */
} cheb1d_edge_t;

// From the generated tables
extern const uint64_t chebcheb_total_doubles;
extern const chebcheb_model_t chebcheb_models[];
extern const chebcheb_strip_t chebcheb_strips[];
extern const int chebcheb_num_models;
extern const int chebcheb_num_strips;
extern const int chebcheb_num_edge_tiles_UL;
extern const int chebcheb_num_edge_tiles_UR;
extern cheb1d_edge_t edge_tiles_UL[];
extern cheb1d_edge_t edge_tiles_UR[];

// Pointer to the loaded coefficients
extern const double *chebcheb_coeffs;

int  chebcheb_load_coeffs(const char *path);
void chebcheb_unload_coeffs(void);
void chebcheb_init_vertical_edges(void);

/* Accessor into the linear coefficient array */
static inline const double *get_C_ptr(int k)
{
    return chebcheb_coeffs + chebcheb_models[k].off_C;
}

#endif /* CHEBCHEB_DATA_H */
