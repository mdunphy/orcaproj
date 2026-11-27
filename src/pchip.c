/* pchip.c — Tensor-product PCHIP tile evaluator */

#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>

#include "common.h"
#include "pchip.h"
#include "pchip_data.h"
#include "internal_errors.h"

typedef enum {
    PCHIP_ENDPOINT_STANDARD,
    PCHIP_ENDPOINT_LINEAR
} pchip_endpoint_mode_t;

pchip_tile_query_t pchip_locate(double i, double j)
{
    pchip_tile_query_t Q;
    Q.tile = -1;
    Q.in_edge = 0;
    Q.i = i;
    Q.j = j;
    // Determine which tile contains (i,j)
    for (int t = 0; t < pchip_num_tiles; t++) {
        const pchip_tile_t *T = &pchip_tiles[t];

        if ( j < T->jmin) continue;

        if (j >= T->jmin && j <= T->jmax &&
            i >= T->imin && i <= T->imax)
        {
            Q.tile = t;
            Q.T = T;
            Q.M = pchip_coeffs + T->off;

            double nxm1 = (double)(T->nx - 1);
            double nym1 = (double)(T->ny - 1);

            double i_loc = (i - T->imin) * nxm1 / (T->imax - T->imin);
            double j_loc = (j - T->jmin) * nym1 / (T->jmax - T->jmin);

            if (i_loc <  0.0) i_loc = 0.0;
            if (i_loc > nxm1) i_loc = nxm1;
            if (j_loc <  0.0) j_loc = 0.0;
            if (j_loc > nym1) j_loc = nym1;

            int ii = (int)floor(i_loc);
            int jj = (int)floor(j_loc);

            if (ii >= T->nx - 1) { ii = T->nx - 2; i_loc = nxm1; }
            if (jj >= T->ny - 1) { jj = T->ny - 2; j_loc = nym1; }
            //printf("DBG %.22f %.22f %d %f\n",i,i_loc,ii,i_loc-ii);
            Q.i_loc = i_loc;
            Q.j_loc = j_loc;
            Q.ii = ii;
            Q.jj = jj;
            Q.ti = i_loc - ii;
            Q.tj = j_loc - jj;

            if (t == 0) {
                // UL tile stitched at RIGHT edge
                if (ii == T->nx - 2)
                    Q.in_edge = 1;
            } else {
                // UR tile stitched at LEFT edge
                if (ii == 0)
                    Q.in_edge = 1;
            }
            Q.j_blend_min = T->jmin;
            Q.j_blend_max = T->jmin + PCHIP_CHEB_BLEND_WIDTH;

            return Q;
        }
    }

    // Point not in either tile
    return Q;
}

static inline double pchip_slope_pair(double dL, double dR)
{
    // Fritsch–Butland (1984) monotone cubic slope formula
    if (dL * dR <= 0.0) return 0.0;
//    if (dL * dR <= 0.0)
//        return 0.5 * (dL + dR);   // secant average, not zero
    return (2*dL*dR) / (dL + dR);
}

static inline double pchip_endpoint_slope(double d0, double d1)
{
    // uniform-spacing PCHIP endpoint slope (then limit)
    double m = 0.5 * (3.0*d0 - d1);

    if (m * d0 <= 0.0) return 0.0;
    if (fabs(m) > 3.0 * fabs(d0)) return 3.0 * d0;
    return m;
}


// Simple monotone limiter for a single Hermite segment
static inline void limit_monotone(double fL, double fR, double *mL, double *mR)
{
    double df = fR - fL;
    if (df == 0.0) {
        *mL = *mR = 0.0;
        return;
    }
    if ((*mL) * df <= 0.0) *mL = 0.0;
    if ((*mR) * df <= 0.0) *mR = 0.0;

    double sgn   = (df > 0.0) ? 1.0 : -1.0;
    double limit = 3.0 * fabs(df);

    if (fabs(*mL) > limit) *mL = sgn * limit;
    if (fabs(*mR) > limit) *mR = sgn * limit;
}

static inline double hermite_unit_interval(double y0, double y1, double m0, double m1, double t)
{
    double t2 = t * t;
    double t3 = t2 * t;
    double h00 =  2*t3 - 3*t2 + 1;
    double h10 =      t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 =      t3 - t2;
    return h00*y0 + h10*m0 + h01*y1 + h11*m1;
}


static inline double hermite_unit_interval_deriv(double y0, double y1,
                                                 double m0, double m1,
                                                 double t)
{
    double t2 = t * t;

    double dh00 = 6*t2 - 6*t;
    double dh10 = 3*t2 - 4*t + 1;
    double dh01 = -dh00;
    double dh11 = 3*t2 - 2*t;

    return dh00*y0 + dh10*m0 + dh01*y1 + dh11*m1;
}


static inline void pchip_local_eval_1d(
    const double *base,  /* pointer to the start of the 1D data column or row       */
    int stride,          /* memory stride between consecutive samples               */
    int n,               /* total number of samples in the 1D array                 */
    int k,               /* integer cell index: evaluate between base[k]..base[k+1] */
    double t,            /* local coordinate in [0,1] within that interval          */
    pchip_endpoint_mode_t endpoint_mode, /* linear extension or standrad */
    double *out_f,       /* function value                                          */
    double *out_df       /* derivative of function value, optional (nullable)       */
)
{
    // This function is called for evaluation in I or J
    // endpoint mode should be linear for i, standard otherwise

    size_t s = (size_t)stride;
    if (k < 0) k = 0;
    if (k > n - 2) k = n - 2;

    size_t k0 = (size_t)k * s;
    size_t k1 = (size_t)(k + 1) * s;

    double fk  = base[k0];
    double fk1 = base[k1];
    double dk  = fk1 - fk;

    double m_k, m_k1;

    if (k == 0) {
        // Left endpoint
        double d0 = fk1 - fk;    // f1 - f0
        double d1 = (n > 2) ? (base[(size_t)2 * s] - fk1) : d0;

        if (endpoint_mode == PCHIP_ENDPOINT_STANDARD) {
            // Endpoint slope
            m_k = 0.5 * (3.0*d0 - d1);
            if (m_k * d0 <= 0.0) m_k = 0.0;
        }
        else {
            // Linear continuation
            m_k = d0;
        }
        if (fabs(m_k) > 3.0 * fabs(d0)) m_k = 3.0 * d0;

        // Interior slope
        m_k1 = pchip_slope_pair(d0, d1);
    }
    else if (k == n - 2) {
        // Right endpoint
        double fkm1 = base[(size_t)(k - 1) * s];
        double dNm3 = fk  - fkm1;   // f[n-2] - f[n-3]
        double dNm2 = fk1 - fk;     // f[n-1] - f[n-2]

        // Interior slope
        m_k = pchip_slope_pair(dNm3, dNm2);

        // Endpoint slope
        if (endpoint_mode == PCHIP_ENDPOINT_STANDARD) {
            m_k1 = 0.5 * (3.0*dNm2 - dNm3);
            if (m_k1 * dNm2 <= 0.0) m_k1 = 0.0;
        }
        else {
            m_k1 = dNm2;   // linear continuation
        }
        if (fabs(m_k1) > 3.0 * fabs(dNm2)) m_k1 = 3.0 * dNm2;
    }
    else {
        // Interior
        double fkm1 = base[(size_t)(k - 1) * s];
        double fk2  = base[(size_t)(k + 2) * s];

        m_k  = pchip_slope_pair(fk - fkm1, dk);
        m_k1 = pchip_slope_pair(dk, fk2 - fk1);
    }

    if (out_f)
        *out_f  = hermite_unit_interval(fk, fk1, m_k, m_k1, t);
    if (out_df)
        *out_df = hermite_unit_interval_deriv(fk, fk1, m_k, m_k1, t);
}

static inline void pchip2d_contract_j(
    const pchip_tile_query_t *Q, /* locator output                         */
    double vj[4],       /* output: four 1D values sampled along j          */
    double dvj[4]       /* output, optional: derivatives wrt j.            */
)
{
    int nx = Q->T->nx;
    int ny = Q->T->ny;

    int base;
    if (Q->ii <= 0) {
        base = 0;            // use columns 0,1,2,3
    } else if (Q->ii >= nx - 2) {
        base = nx - 4;       // use columns nx-4,nx-3,nx-2,nx-1
    } else {
        base = Q->ii - 1;    // centered stencil
    }

    for (int k = 0; k < 4; k++) {
        int ic = base + k;
        if (ic < 0)      ic = 0;
        if (ic > nx - 1) ic = nx - 1;

        const double *col0 = Q->M + (size_t)ic;
        pchip_local_eval_1d(col0, nx, ny, Q->jj, Q->tj,
                            PCHIP_ENDPOINT_STANDARD,
                            &vj[k],                 // function value
                            dvj ? &dvj[k] : NULL);  // derivative if requested
    }
}


static int pchip2d_eval(
    const pchip_tile_query_t *Q,   /* locator output                         */
    double *out,         /* output value                                     */
    double *out_di,      /* nullable */
    double *out_dj       /* nullable */
)
{
    if (!Q || !Q->T || !Q->M || !out ) return ORCA_INT_EFAULT;
    int nx = Q->T->nx;
    int ny = Q->T->ny;
    if (nx < 2 || ny < 2) return ORCA_INT_BAD_VALUE;
    double vj[4];
    double dvj[4];

    // Evaluate at j, 4 points
    pchip2d_contract_j(Q, vj, out_dj ? dvj : NULL);


    // evaluate 1D PCHIP along i
    double fi, dfi;
    pchip_local_eval_1d(vj, 1, 4, 1,
                        Q->ti,
                        PCHIP_ENDPOINT_LINEAR,
                        &fi,
                        out_di ? &dfi : NULL);
    *out = fi;
    if (out_di) {
        if (Q->ii == 0 && Q->ti == 0.0) {
            *out_di = (vj[1] - vj[0]);
        }
        else if (Q->ii == nx-2 && Q->ti == 1.0) {
            *out_di = (vj[3] - vj[2]);
        }
        else {
            // interior PCHIP
            *out_di = dfi;
        }
    }

    // compute df/dj from dvj[] via Hermite basis in i-direction
    if (out_dj) {
    double dfj;
    pchip_local_eval_1d(dvj, 1, 4, 1,
                        Q->ti,
                        PCHIP_ENDPOINT_STANDARD,
                        &dfj,
                        NULL);
    *out_dj = dfj;
    }

    return ORCA_INT_OK;
}


int pchip_eval_tile(
    const pchip_tile_query_t *Q,
    double *out,               /* output interpolated value                    */
    double *out_di,            /* nullable */
    double *out_dj,            /* nullable */
    const double *f_edge,      /* optional boundary value at i = imin or imax  */
    const double *dfdi_edge)   /* optional boundary df/di at i = imin or imax  */
{
    if (!Q || !Q->T || !Q->M || !out) return ORCA_INT_EFAULT;
    const pchip_tile_t *T = Q->T;
    int nx = T->nx;
    int ny = T->ny;
    double scale_i = (double)(nx - 1) / (T->imax - T->imin);
    double scale_j = (double)(ny - 1) / (T->jmax - T->jmin);

    /* in_edge=1 means we're in the interval that has provided
       slope from chebyshev tile on one end.
       ==> here we use the provided f & df/di                  */

    /* If we are NOT in_edge, or
       if we are in_edge but with no stitching data,
          THEN we use regular normal 2D eval                     */
    if ((!Q->in_edge ) || (!f_edge || !dfdi_edge)) {
        double fi, dfi, dfj;
        int rc = pchip2d_eval(Q, &fi,
                              out_di ? &dfi : NULL,
                              out_dj ? &dfj : NULL);
        if (rc != ORCA_INT_OK) return rc;
        *out = fi;
        if (out_di) *out_di = dfi * scale_i;
        if (out_dj) *out_dj = dfj * scale_j;
        return rc;
    }

    // If we get here, we must be in one of the stitched tile edge intervals
    double vj[4];
    double dvj[4];

    // Contract in j into vj & dvj
    pchip2d_contract_j(Q, vj, dvj);

    int left_edge  = (Q->ii == 0);
    int right_edge = (Q->ii == nx - 2);

    double fB = *f_edge;
    double dB = (*dfdi_edge) / scale_i;

    double fL=-4, fR=-5, mL=-6, mR=-7, s0=-8, s1=-9;

    // Handle the leftmost interval that is stitched to the chebcheb tiles
    // This is the UR tile
    if (left_edge) {
        fL = fB;
        mL = dB;
        fR = vj[1];
        double f2 = vj[2];
        mR = pchip_slope_pair(fR - fL, f2 - fR);
    }
    // Handle the rightmost interval that is stitched to the chebcheb tiles
    // This is the UL tile
    else if (right_edge) {
        fL = vj[2];
        double f0 = vj[1];
        fR = fB;
        mR = dB;
        s0 = fL - f0;
        s1 = fR - fL;
        mL = pchip_slope_pair(s0, s1);
    }

    limit_monotone(fL, fR, &mL, &mR);
    *out = hermite_unit_interval(fL, fR, mL, mR, Q->ti);

    if (out_di) {
        double dfdt = hermite_unit_interval_deriv(fL, fR, mL, mR, Q->ti);
        *out_di = dfdt * scale_i;
    }
    // df/dj: use Hermite basis in i-direction
    if (out_dj) {
        double dfj_loc;
        // dvj[] is df/d(j_index) at 4 i-columns (from pchip2d_contract_j)
        // Interpolate along i with the same 1D PCHIP used elsewhere.
        pchip_local_eval_1d(dvj, 1, 4, 1, Q->ti, PCHIP_ENDPOINT_STANDARD, &dfj_loc, NULL);
        *out_dj = dfj_loc * scale_j;   // convert index-j derivative to physical df/dj
    }

    return ORCA_INT_OK;
}
