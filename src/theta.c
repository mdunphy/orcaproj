#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "functions.h"
#include "blended.h"
#include "inverse.h"
#include "orca.h"
#include "errors.h"
#include "orca_error_map.h"
#include "theta.h"

#define DBG_TIME 0

#if DBG_TIME
#include "bench.h"
#endif

int ij_to_theta_vector(const double *ig_in,
                       const double *ej_in,
                       size_t n,
                       double *theta_out,
                       double *dtheta_di,
                       double *dtheta_dj,
                       int *s_out)
{
    if (!ig_in || !ej_in || !theta_out || !s_out) return ORCA_NULL_PTR;
    if ((dtheta_di == NULL) != (dtheta_dj == NULL)) return ORCA_INVALID_ARG;
    if (n == 0) return ORCA_BAD_RANGE;
    if (orca_init() != ORCA_OK) return ORCA_INIT_FAILED;

    #ifdef _OPENMP
        #pragma omp parallel for default(none) \
        shared(ig_in, ej_in, theta_out, dtheta_di, dtheta_dj, s_out, n) \
        schedule(dynamic, 64)
    #endif
    for (size_t k = 0; k < n; ++k) {
        double *pdi = dtheta_di ? &dtheta_di[k] : NULL;
        double *pdj = dtheta_dj ? &dtheta_dj[k] : NULL;
        orca_error_t st = ij_to_theta(ig_in[k], ej_in[k], &theta_out[k], pdi, pdj);
        s_out[k] = (int)st;
    }
    return ORCA_OK;
}

// Direct fit evaluator, used for Jacobian etc in Python
int ij_to_theta(double ig, double ej, double *theta_out, double *dtheta_di, double *dtheta_dj)
{
    if (orca_init() != ORCA_OK) return ORCA_INIT_FAILED;

#if DBG_TIME
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
#endif
    if (!theta_out) return ORCA_NULL_PTR;
    *theta_out = NAN;

    int derivs = dtheta_di && dtheta_dj;  // flag for derivatives
    if (derivs) *dtheta_di = *dtheta_dj = NAN; // set to NaN

    // we require ej in [EJ_MIN,EJ_MAX] = [1,220.5] for this function
    if (ig < bounds.IG_MIN || ig > bounds.IG_MAX ||
        ej < bounds.EJ_MIN || ej > bounds.EJ_MAX) {
        return ORCA_BAD_RANGE;
    }

    // Fold i across the seam
    const double i = (ig > bounds.I_MAX) ? (bounds.IG_MAX - ig) : ig;

    // Get theta from direct evaluator
    double THETA=-3, DTHETA_DI=-4, DTHETA_DJ=-5;
    orca_error_t rc_theta;
    orca_int_error_t rc_int_theta;
    if (derivs) {
        rc_int_theta = blended_eval_ij(i, ej, &THETA, &DTHETA_DI, &DTHETA_DJ);
        if (ig > bounds.I_MAX) DTHETA_DJ *= -1;
    }
    else
        rc_int_theta = blended_eval_ij(i, ej, &THETA, NULL, NULL);

    if ((rc_int_theta != ORCA_INT_OK) || !isfinite(THETA)) {
        rc_theta = orca_map_internal_error(rc_int_theta);
        return rc_theta;
    }

    *theta_out = THETA;
    if (derivs) {
        *dtheta_di = DTHETA_DI;
        *dtheta_dj = DTHETA_DJ;
    }

#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t1);
    *theta_out = elapsed_us(t0, t1);
#endif
    return ORCA_OK;
}
