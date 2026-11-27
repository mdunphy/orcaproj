/* blended.c — Blended evaluator/solver stitching PCHIP and Chebyshev regions. */
#include <math.h>
#include <stdio.h>
#include <float.h>
#include "blended.h"
#include "chebcheb.h"
#include "pchip.h"
#include "solver.h"
#include "common.h"
#include "functions.h"
#include "internal_errors.h"

typedef struct {
    double j0, theta_target;
    row_cache_t Rcheb;    // row-cache for Chebyshev evals
} blended_root_ctx_t;

#define CUBIC 1
#define QUINTIC (1-CUBIC)
static inline double weight_function(double x)
{
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    /* cubic 3x^2 - 2x^3 on [0,1], clipped */
    if (CUBIC) return x * x * (3.0 - 2.0 * x);
    /* quintic 6x^5 - 15x^4 + 10x^3 on [0,1], clipped */
    else if (QUINTIC) return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

static inline double weight_function_prime(double x)
{
    /* derivative of weight_function on (0,1), zero outside */
    if (x <= 0.0 || x >= 1.0) return 0.0;
    if (CUBIC) return 6.0 * x * (1.0 - x);
    else if (QUINTIC) return 30.0 * x * x * (1.0 - x) * (1.0 - x);
}

static inline int cheb_eval_blended(double i, double j,
                                    double *out,
                                    double *dout_di,
                                    double *dout_dj,
                                    int *chosen_idx,
                                    row_cache_t *rowcache)
{
    int ks[2];
    int n = find_covering_models(i, j, ks);
    if (n == 0) return 0;

    if (n == 1) {
        int k = ks[0];
        double val = eval_cheb_tile(k, i, j, rowcache);
        if (dout_di || dout_dj) {
            double fx, fy;
            eval_tile_grad(k, i, j, &fx, &fy);
            if (dout_di) *dout_di = fx;
            if (dout_dj) *dout_dj = fy;
        }
        if (out) *out = val;
        if (chosen_idx) *chosen_idx = k;
        return 1;
    }

    int a = ks[0], b = ks[1];
    double va = eval_cheb_tile(a, i, j, rowcache);
    double vb = eval_cheb_tile(b, i, j, rowcache);

    const chebcheb_model_t *A = &chebcheb_models[a];
    const chebcheb_model_t *B = &chebcheb_models[b];

    double jmin = fmax(A->jmin, B->jmin);
    double jmax = fmin(A->jmax, B->jmax);
    double t    = (j - jmin) / (jmax - jmin);
    double w = weight_function(t);
    if (out)
        *out = (1.0 - w)*va + w*vb;

    if (dout_di || dout_dj) {
        double fax, fay;
        double fbx, fby;
        eval_tile_grad(a, i, j, &fax, &fay);
        eval_tile_grad(b, i, j, &fbx, &fby);

        double dw_dj = weight_function_prime(t) / (jmax - jmin);

        /* dC/di = (1-w)*dva/di + w*dvb/di  (w depends only on j) */
        if (dout_di)
            *dout_di = (1.0 - w)*fax + w*fbx;

        /* dC/dj = (1-w)*dva/dj + w*dvb/dj + dw/dj * (vb - va) */
        if (dout_dj)
            *dout_dj = (1.0 - w)*fay + w*fby + dw_dj * (vb - va);
    }
    if (chosen_idx) *chosen_idx = a;
    return 1;
}


static inline void cheb_eval_blended_edge_by_index(
        double j,
        int pchip_tile_idx,
        double *f,
        double *dfdi,
        double *dfdj)
{
    if (pchip_tile_idx == 0) {
        // UL corner tile
        cheb1d_eval_f_and_dfdi_and_dfdj(&edge_tiles_UL[0], j, f, dfdi, dfdj);
        return;
    }
    if (pchip_tile_idx == 1) {
        // UR corner tile
        double fa, dfdia, dfdja;
        double fb, dfdib, dfdjb;
        const cheb1d_edge_t *edge1 = &edge_tiles_UR[0];
        const cheb1d_edge_t *edge2 = &edge_tiles_UR[1];
        cheb1d_eval_f_and_dfdi_and_dfdj(edge1, j, &fa, &dfdia, &dfdja);
        cheb1d_eval_f_and_dfdi_and_dfdj(edge2, j, &fb, &dfdib, &dfdjb);

        // j-blending
        double jmin = fmax(edge1->jmin, edge2->jmin);
        double jmax = fmin(edge1->jmax, edge2->jmax);
        if (jmax <= jmin) {
            *f    = fa;
            *dfdi = dfdia;
            *dfdj = dfdja;
            return;
        }
        double t    = (j - jmin) / (jmax - jmin);
        double w = weight_function(t);
        double dw_dj = weight_function_prime(t) / (jmax - jmin);

        *f    = (1.0-w)*fa      + w*fb;
        *dfdi = (1.0-w)*dfdia   + w*dfdib;
        *dfdj = (1.0-w)*dfdja   + w*dfdjb + dw_dj * (fb - fa);

        return;
    }
    return;  // should never get here
}


// blended_eval_ij_core:
//  - If Rcheb == NULL, we evaluate Chebyshev without cache
//  - If Rcheb != NULL, then Rcheb must point to a row row-cache, which we use
static orca_int_error_t blended_eval_ij_core(double i, double j,
                                             double *theta_out,
                                             double *dtheta_di,
                                             double *dtheta_dj,
                                             row_cache_t *Rcheb)
{
    if (!theta_out)
        return ORCA_INT_EFAULT;
    if ((dtheta_di == NULL) != (dtheta_dj == NULL))
        return ORCA_INT_EFAULT;  // we must have both null or both not-null

    int derivs = (dtheta_di && dtheta_dj);  // flag for derivatives
    if (derivs) {
        *dtheta_di = NAN;
        *dtheta_dj = NAN;
    }

    double P = 0.0, C = 0.0;
    double dPdi = 0.0, dPdj = 0.0;
    double dCdi = 0.0, dCdj = 0.0;

    int has_pchip = 0, has_cheb = 0;
    int pchip_idx = -1;
    int cheb_idx  = -1;

    // Determine if current point is in a Chebyshev tile
    double *C_di = derivs ? &dCdi : NULL;
    double *C_dj = derivs ? &dCdj : NULL;
    has_cheb = cheb_eval_blended(i, j, &C, C_di, C_dj, &cheb_idx, Rcheb);
    if ( i == bounds.I_MIN ) C = bounds.THETA_MIN;
    if ( i == bounds.I_MAX ) C = bounds.THETA_MAX;

    // Determine if the current point is in a PCHIP tile
    pchip_tile_query_t Q = pchip_locate(i,j);
    pchip_idx = Q.tile;
    has_pchip = (pchip_idx >= 0);
    if (has_pchip) {
        orca_int_error_t rc;
        double *P_di = derivs ? &dPdi : NULL;
        double *P_dj = derivs ? &dPdj : NULL;
        if (Q.in_edge) {
            double f,dfdi,dfdj;
            cheb_eval_blended_edge_by_index(j, pchip_idx, &f, &dfdi, &dfdj);
            rc = pchip_eval_tile(&Q,
                                 &P, P_di, P_dj,
                                 &f, &dfdi);
        }
        else {
            rc = pchip_eval_tile(&Q,
                                 &P, P_di, P_dj,
                                 NULL, NULL);
            if ( i == bounds.I_MIN) P = bounds.THETA_MIN;
            if ( i == bounds.I_MAX) P = bounds.THETA_MAX;
        }
        if (rc != ORCA_INT_OK || !isfinite(P)) {
            *theta_out = NAN;
            return ORCA_INT_BAD_VALUE;
        }
    }

    // impossible case
    if (!has_pchip && !has_cheb) {
        *theta_out = NAN;
        return ORCA_INT_NO_TILE_COVERAGE;
    }

    // PCHIP only
    if (has_pchip && !has_cheb) {
        *theta_out = P;
        if (derivs) {
            *dtheta_di = dPdi;
            *dtheta_dj = dPdj;
        }
        return ORCA_INT_OK;
    }

    // CHEB only
    if (!has_pchip && has_cheb) {
        *theta_out = C;
        if (derivs) {
            *dtheta_di = dCdi;
            *dtheta_dj = dCdj;
        }
        return ORCA_INT_OK;
    }

    // Blended region
    double wj = 0, dwj_dj = 0;
    if (has_pchip) {
        /* Cheb -> PCHIP */
        double A = Q.j_blend_min;
        double B = Q.j_blend_max;
        double t = (j - A) / (B - A);
        wj = weight_function(t);
        dwj_dj = weight_function_prime(t) / (B - A);
    }
    double theta = (1.0 - wj)*C + wj*P;

    if ( theta < bounds.THETA_MIN) theta=bounds.THETA_MIN;
    if ( theta > bounds.THETA_MAX) theta=bounds.THETA_MAX;

    *theta_out = theta;
    if (derivs) {
        // dtheta/di: w_j depends only on j
        double dti = (1.0 - wj)*dCdi + wj*dPdi;

        // dtheta/dj: product rule
        double dtj = (1.0 - wj)*dCdj + wj*dPdj + dwj_dj*(P - C);

        *dtheta_di = dti;
        *dtheta_dj = dtj;
    }

    return ORCA_INT_OK;
}


// Direct evaluator; does not use a row cache. Called by inverse.c
orca_int_error_t blended_eval_ij(double i, double j, double *theta_out, double *dtheta_di, double *dtheta_dj)
{
    return blended_eval_ij_core(i, j, theta_out, dtheta_di, dtheta_dj, NULL);
}


static double blended_root_f(double i, void *userdata)
{
    // Evaluator for root finding; we do use a row cache here
    blended_root_ctx_t *ctx = (blended_root_ctx_t*)userdata;

    double theta;
    orca_int_error_t rc = blended_eval_ij_core(i, ctx->j0, &theta, NULL, NULL, &ctx->Rcheb);
    if ( (rc != ORCA_INT_OK) || !isfinite(theta))
        return NAN;

    return theta - ctx->theta_target;
}

// Called by forward.c
orca_int_error_t blended_find_i(double j0, double theta, double lb, double ub,
                   double tol, int maxiter, double *i_out, int *niter)
{
    if (!i_out || !niter)
        return ORCA_INT_EFAULT;

    const double a = lb;
    const double b = ub;
    blended_root_ctx_t ctx;
    ctx.j0 = j0;
    ctx.theta_target = theta;
    chebcheb_init_rowcache(&ctx.Rcheb, j0);

    double fa = blended_root_f(a, &ctx);
    double fb = blended_root_f(b, &ctx);

    // endpoint roots
    if (fabs(fa) < tol) {
        *niter = 0;
        *i_out = a;
        return ORCA_INT_OK;
    }
    if (fabs(fb) < tol) {
        *niter = 0;
        *i_out = b;
        return ORCA_INT_OK;
    }

    if (fa * fb > 0) {
        //printf("blended_find_i %.18e %.18e %.18e %.18e, %.18e\n",  fa, fb, j0, a, b);
        *niter = 0;
        *i_out = NAN;
        return ORCA_INT_ROOT_NO_BRACKET;
    }

    // Now call Brent to solve the root
    double i_star = -888;
    int n_iter = 0;
    int rc = brent_solve(&blended_root_f, &ctx,
                         a, b, tol, maxiter,
                         &i_star, &n_iter, &fa, &fb);

    if (rc == 2) { // rc=2, Brent rejects bracket
        *i_out = (fabs(fa) <= fabs(fb)) ? a : b;
        return ORCA_INT_ROOT_NO_BRACKET;
    }

    if (rc == 1) { // rc=1, max iter but still a valid bracket
        *i_out = i_star;
        return ORCA_INT_ROOT_MAXITER;
    }

    // Success case
    *niter = n_iter;
    //*niter = __builtin_popcountll(ctx.Rcheb.filled); // check how many tiles were collapsed for this solve
    *i_out = i_star;
    return ORCA_INT_OK;
}
