#include <math.h>

#include "common.h"
#include "functions.h"
#include "blended.h"
#include "theta.h"
#include "inverse.h"
#include "orca.h"
#include "errors.h"

// TODO: as in forward.c, need to improve this
#define DBG_TIME 0
#if DBG_TIME
#include "bench.h"
#endif

int ij_to_lonlat_vector(const double *i_in,
                        const double *j_in,
                        size_t n,
                        double *lon_out,
                        double *lat_out,
                        int *s_out)
{
    if (!i_in || !j_in || !lon_out || !lat_out || !s_out) return ORCA_NULL_PTR;
    if (n == 0) return ORCA_BAD_RANGE;
    if (orca_init() != ORCA_OK) return ORCA_INIT_FAILED;

    #ifdef _OPENMP
        #pragma omp parallel for default(none) \
        shared(i_in, j_in, lon_out, lat_out, s_out, n) \
        schedule(dynamic, 64)
    #endif
    for (size_t k = 0; k < n; ++k) {
        s_out[k] = (int)ij_to_lonlat(i_in[k], j_in[k], &lon_out[k], &lat_out[k]);
    }
    return ORCA_OK;
}

// Inverse transform, scalar
int ij_to_lonlat(double ig, double jg, double *lon_out, double *lat_out)
{
    if (orca_init() != ORCA_OK) return ORCA_INIT_FAILED;
#if DBG_TIME
    struct timespec t0, t1, t2, t3, t4, t5;
    clock_gettime(CLOCK_MONOTONIC, &t0);
#endif
    if (!lon_out || !lat_out) return ORCA_NULL_PTR;
    *lon_out = NAN; *lat_out = NAN;

    if (isnan(ig) || isnan(jg)) return ORCA_NAN_INPUT;
    if (!isfinite(ig) || !isfinite(jg)) return ORCA_INF_INPUT;
    if (!(ig >= bounds.IG_MIN && ig <= bounds.IG_MAX)) return ORCA_BAD_RANGE;
    if (!(jg >= bounds.JG_MIN && jg <= bounds.JG_MAX)) return ORCA_BAD_RANGE;
    // jg is in [JG_MIN, JG_MAX]

    double ej = fsj(jg); // Coordinate shift from grid j [1, 510.5] to ellipsej, ej in [-289, 220.5]

    // Southern (Mercator) branch
    if (ej <= 1.0) {
        const double lat = fsy0p(jg);
        const double lon = wraplon(fslam(ig+1) + rlampo - 180.0);
        *lon_out = lon;
        *lat_out = lat;
        // printf("REV SH: lon %g lat %g ig %g ii %g jg %g j %g\n", lon, lat, ig, ii, jg, j);
        return ORCA_OK;
    }

    // Northern ORCA branch
#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t1);
#endif
    // Get theta by directly evaluating the fit
    double THETA;
    orca_error_t rc_theta = ij_to_theta(ig, ej, &THETA, NULL, NULL);
    if ((rc_theta != ORCA_OK) || !isfinite(THETA)) {
        return rc_theta;
    }
#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t2);
#endif

    // Convert theta to x,y
    double A,B,Y0;
    aby0(ej, &A, &B, &Y0);

    double x =      B * cos(THETA);
    double y = Y0 + A * sin(THETA);
    if (ig > bounds.I_MAX) x *= -1;

    // Inverse stereographic to lon/lat
    double lon, lat;
    fsxylamphi(x, y, &lon, &lat);

    // Ensure longitude is in [-180, 180)
    lon = wraplon(lon);
    *lon_out = lon;
    *lat_out = lat;

#if DBG_TIME
    clock_gettime(CLOCK_MONOTONIC, &t3);
    double ns_total = elapsed_us(t0, t3);
    double ns_eval = elapsed_us(t1, t2);
    *lon_out = ns_total;
    *lat_out = ns_eval;
#endif
    return rc_theta;
}
