/* pchip.h — Interface for tensor-product PCHIP evaluator. */

#ifndef ORCA_TRIPOLAR_PCHIP_H
#define ORCA_TRIPOLAR_PCHIP_H

#include "pchip_data.h"

#define PCHIP_CHEB_BLEND_WIDTH 1.0

typedef struct {
    int    tile;           /* Tile index; -1 if outside all tiles */
    int    in_edge;        /* Non-zero if the point lies in the stitched edge band */
    double i, j;           /* Query coordinates, i in [0,360], j in [1,220.5] */
    double i_loc, j_loc;   /* Logical coordinates after clamping/scaling */
    int    ii, jj;         /* Integer cell indices */
    double ti, tj;         /* Fractional local coordinates in [0,1] */
    const double *M;       /* Pointer to tile data block */
    const pchip_tile_t *T; /* Pointer to tile metadata */
    double j_blend_min;    /* j value where we start blending */
    double j_blend_max;    /* j value where we stop blending */
} pchip_tile_query_t;


pchip_tile_query_t pchip_locate(double i, double j);

int pchip_eval_tile(const pchip_tile_query_t *Q,
                    double *out,
                    double *out_di,
                    double *out_dj,
                    const double *f_edge,
                    const double *dfdi_edge);

#endif /* ORCA_TRIPOLAR_PCHIP_H */
