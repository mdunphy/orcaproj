#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "data.h"

int load_coeff_file(const char *path, uint64_t expected_count, double **out_buf, const double **out_const)
{
    if (*out_const) return 0;

    if (!path) {
        fprintf(stderr, "load_coeff_file: path not provided\n");
        return -1;
    }

    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); return -1; }

    if (fseek(f, 0, SEEK_END) != 0) { perror("fseek"); fclose(f); return -1; }
    long sz = ftell(f);
    if (sz < 0) { perror("ftell"); fclose(f); return -1; }
    rewind(f);

    uint64_t nd = (uint64_t)sz / sizeof(double);

    if ((uint64_t)sz != nd * sizeof(double)) {
        fprintf(stderr, "bin not multiple of sizeof(double)\n");
        fclose(f);
        return -1;
    }

    if (nd != expected_count) {
        fprintf(stderr, "bin doubles=%" PRIu64 " expected=%" PRIu64 "\n", nd, expected_count);
        fclose(f);
        return -1;
    }

    double *buf = (double*)malloc((size_t)nd * sizeof(double));
    if (!buf) { perror("malloc"); fclose(f); return -1; }

    size_t got = fread(buf, sizeof(double), (size_t)nd, f);
    fclose(f);

    if (got != (size_t)nd) {
        fprintf(stderr, "Short read\n");
        free(buf);
        return -1;
    }

    *out_buf = buf;
    *out_const = buf;

    return 0;
}


void unload_coeff_file(double **buf, const double **pub) {
    if (*buf) {
        free(*buf);
        *buf = NULL;
    }
    *pub = NULL;
}
