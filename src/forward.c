#include <math.h>
#include <float.h>
#include <stdio.h>

#include "common.h"
#include "functions.h"
#include "blended.h"
#include "solver.h"
#include "forward.h"
#include "xy_to_ej.h"
#include "orca.h"
#include "errors.h"
#include "orca_error_map.h"

// TODO: rework this hack for exporting the timing/iteration counts
#define DBG_ITER 0
#define DBG_TIME 0
#if DBG_TIME
#include "bench.h"
#endif

int lonlat_to_ij_vector(const double *lon_deg,
                        const double *lat_deg,
                        size_t n,
                        double *i_out,
                        double *j_out,
                        int    *s_out)
{
    if (!lon_deg || !lat_deg || !i_out || !j_out || !s_out)  return ORCA_NULL_PTR;
    if (n == 0) return ORCA_BAD_RANGE;
    if (orca_init() != ORCA_OK) return ORCA_INIT_FAILED;

    #ifdef _OPENMP
        #pragma omp parallel for default(none) \
        shared(lon_deg, lat_deg, i_out, j_out, s_out, n) \
        schedule(dynamic, 64)
    #endif
    for (size_t k = 0; k < n; ++k) {
        s_out[k] = (int)lonlat_to_ij(lon_deg[k], lat_deg[k], &i_out[k], &j_out[k]);
    }
    return ORCA_OK;
}



// Forward transform (scalar)
int lonlat_to_ij(double lon_deg,
                 double lat_deg,
                 double *i_out,
                 double *j_out)
{
    if (orca_init() != ORCA_OK) return ORCA_INIT_FAILED;
    if (!i_out || !j_out) return ORCA_NULL_PTR;
#if DBG_TIME
    struct timespec t0, t1, t2;
    clock_gettime(CLOCK_MONOTONIC, &t0);
#endif

    // Input validation
    if (isnan(lon_deg) || isnan(lat_deg)) return ORCA_NAN_INPUT;
    if (!isfinite(lon_deg) || !isfinite(lat_deg)) return ORCA_INF_INPUT;
    if ((lat_deg > 90.0) || (lat_deg < -77.0104755487162)) return ORCA_BAD_RANGE;

    *i_out = -444.0;
    *j_out = NAN;

    // Southern (Mercator) branch
    double ejcheck = 999;
    const double lat0 = 20.08; // slightly above the cutoff for Mercator
    if (lat_deg <= lat0) {
        // We might be on the Mercator grid; compute j to check
        const double j = phi2j(lat_deg);
        ejcheck = fsj(j);
        // Only when ej < 1.0 do we return the Mercator grid
        if (ejcheck <= 1.0) {
            const double lonw = wraplon(lon_deg - rlampo) + 180.0;
            const double i = lam2i(lonw) - 1.0; // The converts from fortran indexing
            *i_out = i;
            *j_out = j;
            return ORCA_OK;
        }
    }

    // Northern ORCA branch

    // 1: Project (lon,lat) to (x,y)
    double x = NAN, y = NAN;
    fslamphixy(lon_deg, lat_deg, &x, &y);

    // 2: Solve ej
    double ej = -2.0;
    int niter_ej = -1;
#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t1);
#endif
    orca_int_error_t rc_int_ej = xy_to_ej(x, y, &ej, &niter_ej);
    orca_error_t rc_ej = orca_map_internal_error(rc_int_ej);
#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t2);
    double ns_j = elapsed_us(t1, t2);
#endif

    if (rc_ej != ORCA_OK) {
#if DBG_TIME
        *j_out = ns_j;
#elif DBG_ITER
        *j_out = niter_ej;
#else
        *j_out = fsej(ej);
#endif
        return rc_ej;
    }

    // 3: Convert ej to jg
    double jg = fsej(ej);
    *j_out = jg;


    // 4: Compute theta
    double A,B,Y0;
    aby0(ej, &A, &B, &Y0);
    //printf("FWD DBG %.20f %.20f %.20f | fs  %.20g %.20g | yy %.20f\n",ej, B, fsb(220.5), fsfun(ej,x,y), fsfun(220.5,x,y), y-Y0);

    double theta;
    int flip = x > 0.0;

    if (B < 1e-12) {
        double arg = (y - Y0) / A;
        if (arg >  1.0) arg =  1.0;
        if (arg < -1.0) arg = -1.0;
        theta = M_PI - asin(arg);
        //printf("XY %.20f %.20f ASIN %.20f %.20f\n",x,y,theta-M_PI/2,arg);
    }
    else {
        double t0 = atan2(B * (y - Y0), A * x);
        theta = (x <= 0.0 && (t0 - M_PI_2) < 0.0)
              ? t0 + 2.0*M_PI
              : (x > 0.0 ? M_PI - t0 : t0);
        //printf("XY %.20f %.20f A %.20f Y0 %.20f ATAN %.20f %.20f\n",x,y,A,Y0,theta-M_PI/2,t0);
    }


    // 5: Solve THETA(i,ej) = theta for i
    double i_val = -333.0;
    int    niter_i = 0;

    const double tol = 1e-16;
    const int    maxiter = 50;

    // Reduce the initial root bracket.
    // TODO: improve this initial bracket
    const double L = bounds.THETA_MIN;
    const double R = bounds.THETA_MAX;
    double lb = -55 + (theta - L)*322/(R - L);
    double ub =   8 + (theta - L)*375/(R - L);
    if (lb < bounds.I_MIN) lb = bounds.I_MIN;
    if (ub > bounds.I_MAX) ub = bounds.I_MAX;
    //lb=bounds.I_MIN; ub=bounds.I_MAX;
#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t1);
#endif
    //printf("FWD %.20f %.20f | B %.20f | ub %.20f \n",theta-M_PI/2,ej,B,ub);
    orca_int_error_t rc_int_i = blended_find_i(ej, theta, lb, ub, tol, maxiter, &i_val, &niter_i);
#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t2);
    double ns_i = elapsed_us(t1, t2);
    *i_out=ns_i;
#endif
    // Failure except for maxiter
    if (rc_int_i != ORCA_INT_OK && rc_int_i != ORCA_INT_ROOT_MAXITER) {
        *i_out = i_val;
        //printf("FWD DBG5 %g %g %d\n",ej,i_val,niter_i);
        return orca_map_internal_error(rc_int_i);
    }

    // 6: Apply mirroring
    if (flip)
        i_val = bounds.IG_MAX - i_val;
    *i_out = i_val;

#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double ns_total = elapsed_us(t0, t1);
    *i_out = ns_i;
    *j_out = ns_j;
#elif DBG_ITER
    *j_out = niter_ej;
    *i_out = niter_i;
#endif
    return ORCA_OK;
}
