#include <cstdio>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#include "orca.h"
#include "bench.h"

#define GRID_MODE 2

int main() {
    double di, dj;
#if GRID_MODE == 1
    di = dj = 0.5;
#elif GRID_MODE == 2
    di = dj = 0.25;
#elif GRID_MODE == 3
    di = dj = 0.1;
#else
#error "Invalid GRID_MODE"
#endif

    double i_start = 0.0, i_end = 720.0, j_start = 1.0, j_end = 220.5;

    size_t NI = static_cast<size_t>(std::llround((i_end - i_start) / di)) + 1;
    size_t NJ = static_cast<size_t>(std::llround((j_end - j_start) / dj)) + 1;
    size_t NPTS = NI * NJ;

    std::printf("Grid mode: %d (%.2f°x%.2f°) (NI=%zu, NJ=%zu, N=%zu)\n",
                GRID_MODE, di, dj, NI, NJ, NPTS);

    std::vector<double> I(NPTS), J(NPTS);
    std::vector<double> lon(NPTS), lat(NPTS);
    std::vector<double> I2(NPTS), J2(NPTS);
    std::vector<double> lon2(NPTS), lat2(NPTS);
    std::vector<int>    s(NPTS);

    struct timespec t0{}, t1{};
    double ns_total, ns_per_point;
    int status;

    // Step 1
    std::printf("Initializing I,J arrays...\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (size_t j = 0; j < NJ; j++) {
        double j_val = j_start + dj * static_cast<double>(j) + 290;
        for (size_t i = 0; i < NI; i++) {
            double i_val = i_start + di * static_cast<double>(i);
            size_t idx = j * NI + i;
            I[idx] = i_val;
            J[idx] = j_val;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    std::printf("Filled I,J arrays in %.3f s\n\n", ns_total / 1e9);

    // Step 2
    std::printf("Inverse projection: ij_to_lonlat_vector\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    status = ij_to_lonlat_vector(I.data(), J.data(), NPTS,
                                lon.data(), lat.data(), s.data());
    if (status != 0)
        std::fprintf(stderr, "ij_to_lonlat_vector failed: %d\n", status);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    ns_per_point = ns_total / static_cast<double>(NPTS);
    std::printf("Elapsed: %.3f s (%.3f µs/pt)\n\n",
                ns_total / 1e9, ns_per_point / 1e3);

    // Step 3
    std::printf("Forward projection: lonlat_to_ij_vector\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    status = lonlat_to_ij_vector(lon.data(), lat.data(), NPTS,
                                I2.data(), J2.data(), s.data());
    if (status != 0)
        std::fprintf(stderr, "lonlat_to_ij_vector failed: %d\n", status);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    ns_per_point = ns_total / static_cast<double>(NPTS);
    std::printf("Elapsed: %.3f s (%.3f µs/pt)\n\n",
                ns_total / 1e9, ns_per_point / 1e3);

    // Step 4
    std::printf("Inverse projection: ij_to_lonlat_vector\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    status = ij_to_lonlat_vector(I2.data(), J2.data(), NPTS,
                                lon2.data(), lat2.data(), s.data());
    if (status != 0)
        std::fprintf(stderr, "ij_to_lonlat_vector failed: %d\n", status);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);
    ns_per_point = ns_total / static_cast<double>(NPTS);
    std::printf("Elapsed: %.3f s (%.3f µs/pt)\n\n",
                ns_total / 1e9, ns_per_point / 1e3);

    // Step 5
    std::printf("Computing round-trip error...\n");
    clock_gettime(CLOCK_MONOTONIC, &t0);

    double max_err_i = 0.0;
    double max_err_j = 0.0;
    double max_err_lon = 0.0;
    double max_err_lat = 0.0;
    double max_err_i_nocorners = 0.0;
    double max_err_j_nocorners = 0.0;
    double max_err_lon_nocorners = 0.0;
    double max_err_lat_nocorners = 0.0;

    for (size_t j = 0; j < NJ - 1; j++) {
        for (size_t i = 0; i < NI - 1; i++) {
            size_t idx = j * NI + i;

            double err_i = std::fabs(I2[idx] - I[idx]);
            if (err_i > 719.0) err_i -= 720;

            double err_j = std::fabs(J2[idx] - J[idx]);
            double err_lon = std::fabs(lon2[idx] - lon[idx]);
            double err_lat = std::fabs(lat2[idx] - lat[idx]);

            if (err_i > max_err_i) max_err_i = err_i;
            if (err_j > max_err_j) max_err_j = err_j;
            if (err_lon > max_err_lon) max_err_lon = err_lon;
            if (err_lat > max_err_lat) max_err_lat = err_lat;

            if (!(((I[idx] < 2) && (J[idx] > 216)) ||
                  ((I[idx] > 328) && (J[idx] > 210)))) {

                if (err_i > max_err_i_nocorners) max_err_i_nocorners = err_i;
                if (err_j > max_err_j_nocorners) max_err_j_nocorners = err_j;
                if (err_lon > max_err_lon_nocorners) max_err_lon_nocorners = err_lon;
                if (err_lat > max_err_lat_nocorners) max_err_lat_nocorners = err_lat;
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    ns_total = elapsed_ns(t0, t1);

    std::printf("Round-trip error check done in %.3f s\n", ns_total / 1e9);
    std::printf("Max |ΔI|   = %.18f\n", max_err_i);
    std::printf("Max |ΔJ|   = %.18f\n", max_err_j);
    std::printf("Max |Δlon| = %.18f\n", max_err_lon);
    std::printf("Max |Δlat| = %.18f\n", max_err_lat);
    std::printf("Max |ΔI|   = %.18f (no corners)\n", max_err_i_nocorners);
    std::printf("Max |ΔJ|   = %.18f (no corners)\n", max_err_j_nocorners);
    std::printf("Max |Δlon| = %.18f (no corners)\n", max_err_lon_nocorners);
    std::printf("Max |Δlat| = %.18f (no corners)\n", max_err_lat_nocorners);

    return 0;
}