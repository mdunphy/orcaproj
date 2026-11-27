#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "orca.h"
#include "bench.h"

// Select grid mode here
#define GRID_MODE 1

int main(void) {
    double di,dj;
#if GRID_MODE == 1
    di = dj = 0.5; // ORCA05
#elif GRID_MODE == 2
    di = dj = 0.25; // ORCA025 (sort of)
#elif GRID_MODE == 3
    di = dj = 0.1; // Fine grid from grid generator (0.5 / 5 = 0.1)
#else
#error "Invalid GRID_MODE"
#endif

    double i_start = 0.0, i_end = 720.0, j_start = 1.0, j_end = 220.5;
    size_t NI = (size_t)llround((i_end - i_start) / di) + 1;
    size_t NJ = (size_t)llround((j_end - j_start) / dj) + 1;
    size_t NPTS = NI * NJ;

    printf("Grid mode: %d (%.2f°x%.2f°) (NI=%zu, NJ=%zu, N=%zu)\n",GRID_MODE, di, dj, NI, NJ, NPTS);

    double *I  = malloc(NPTS * sizeof(double));
    double *J  = malloc(NPTS * sizeof(double));
    double *lon = malloc(NPTS * sizeof(double));
    double *lat = malloc(NPTS * sizeof(double));
    double *I2  = malloc(NPTS * sizeof(double));
    double *J2  = malloc(NPTS * sizeof(double));
    double *lon2= malloc(NPTS * sizeof(double));
    double *lat2= malloc(NPTS * sizeof(double));
    int    *s   = malloc(NPTS * sizeof(int));

    if (!I || !J || !lon || !lat || !I2 || !J2 || !lon2 || !lat2 || !s) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    struct timespec t0, t1;
    double ns_total, ns_per_point;
    int status;

    // Step 1: Fill I,J arrays
    printf("Initializing I,J arrays...\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (size_t j = 0; j < NJ; j++) {
        double j_val = j_start + dj * (double)j + 290;
        for (size_t i = 0; i < NI; i++) {
            double i_val = i_start + di * (double)i;
            size_t idx = j * NI + i;
            I[idx] = i_val;
            J[idx] = j_val;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    printf("Filled I,J arrays in %.3f s\n\n", ns_total / 1e9);

    // Step 2: Inverse projection
    printf("Inverse projection: ij_to_lonlat_vector\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    status = ij_to_lonlat_vector(I, J, NPTS, lon, lat, s);
    if (status != 0)
        fprintf(stderr, "ij_to_lonlat_vector failed: %d\n", status);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    ns_per_point = ns_total / (double)NPTS;
    printf("Elapsed: %.3f s (%.3f µs/pt)\n\n", ns_total / 1e9, ns_per_point / 1e3);

    // Step 3: Forward projection
    printf("Forward projection: lonlat_to_ij_vector\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    status = lonlat_to_ij_vector(lon, lat, NPTS, I2, J2, s);
    if (status != 0)
        fprintf(stderr, "lonlat_to_ij_vector failed: %d\n", status);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    ns_per_point = ns_total / (double)NPTS;
    printf("Elapsed: %.3f s (%.3f µs/pt)\n\n", ns_total / 1e9, ns_per_point / 1e3);

    // Step 4: Inverse projection again
    printf("Inverse projection: ij_to_lonlat_vector\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    status = ij_to_lonlat_vector(I2, J2, NPTS, lon2, lat2, s);
    if (status != 0)
        fprintf(stderr, "ij_to_lonlat_vector failed: %d\n", status);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    ns_per_point = ns_total / (double)NPTS;
    printf("Elapsed: %.3f s (%.3f µs/pt)\n\n", ns_total / 1e9, ns_per_point / 1e3);

    // Step 5: Compute round-trip errors
    printf("Computing round-trip error...\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    double max_err_i = 0.0;
    double max_err_j = 0.0;
    double max_err_lon = 0.0;
    double max_err_lat = 0.0;
    double max_err_i_nocorners = 0.0;
    double max_err_j_nocorners = 0.0;
    double max_err_lon_nocorners = 0.0;
    double max_err_lat_nocorners = 0.0;
    for (size_t j = 0; j < NJ-1; j++) {
        for (size_t i = 0; i < NI-1; i++) {
            size_t idx = j * NI + i;
            double err_i = fabs(I2[idx] - I[idx]); if (err_i > 719.0) err_i -= 720; // north fold
            double err_j = fabs(J2[idx] - J[idx]);
            double err_lon = fabs(lon2[idx] - lon[idx]);
            double err_lat = fabs(lat2[idx] - lat[idx]);
            // Raw errors
            if (err_i > max_err_i) max_err_i = err_i;
            if (err_j > max_err_j) max_err_j = err_j;
            if (err_lon > max_err_lon) max_err_lon = err_lon;
            if (err_lat > max_err_lat) max_err_lat = err_lat;
            // Mask out corners
            if (!(  ((I[idx] <   2) && (J[idx] > 216)) ||
                    ((I[idx] > 328) && (J[idx] > 210))  )) {
                if (err_i > max_err_i_nocorners) max_err_i_nocorners = err_i;
                if (err_j > max_err_j_nocorners) max_err_j_nocorners = err_j;
                if (err_lon > max_err_lon_nocorners) max_err_lon_nocorners = err_lon;
                if (err_lat > max_err_lat_nocorners) max_err_lat_nocorners = err_lat;
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);

    printf("Round-trip error check done in %.3f s\n", ns_total / 1e9);
    printf("Max |ΔI|   = %.18f\n", max_err_i);
    printf("Max |ΔJ|   = %.18f\n", max_err_j);
    printf("Max |Δlon| = %.18f\n", max_err_lon);
    printf("Max |Δlat| = %.18f\n", max_err_lat);
    printf("Max |ΔI|   = %.18f (no corners)\n", max_err_i_nocorners);
    printf("Max |ΔJ|   = %.18f (no corners)\n", max_err_j_nocorners);
    printf("Max |Δlon| = %.18f (no corners)\n", max_err_lon_nocorners);
    printf("Max |Δlat| = %.18f (no corners)\n", max_err_lat_nocorners);

    // cleanup
    free(I);
    free(J);
    free(lon);
    free(lat);
    free(I2);
    free(J2);
    free(lon2);
    free(lat2);
    free(s);

    return 0;
}
